package PDL::EEG::ICA;
use strict;
use warnings;
use PDL;
use PDL::MatrixOps;
use Exporter 'import';

our $VERSION = '0.01';
our @EXPORT_OK = qw(
    whiten fastica ica_decompose
    identify_by_reference apply_ica ica_operator
    sym_orth inv_sqrt_sym
);

# =====================================================================
# Data convention
#   X : piddle of shape (nt, nch)  -- dim0 = time samples, dim1 = channels
#       (same internal layout as PDL::EEG::GED)
#   sources : (nt, k)              -- one column per independent component
#   unmix   : (nch, k) PDL         -- math (k x nch),  S = unmix x Xc
#   mix     : (k, nch) PDL         -- math (nch x k),  columns = scalp maps
# All matrix products use PDL 'x': (A x B) contracts A dim0 with B dim1.
# =====================================================================

# symmetric matrix inverse square root, via eigen-decomposition
sub inv_sqrt_sym {
    my ($M) = @_;
    my ($U, $s) = eigens_sym($M);
    my $is = $s->clip(1e-15, undef) ** -0.5;
    my $D = zeroes($M->dim(0), $M->dim(0));
    $D->diagonal(0,1) .= $is;
    return $U x $D x $U->transpose;
}

# symmetric orthonormalization  B <- (B B^T)^{-1/2} B
sub sym_orth {
    my ($B) = @_;
    return inv_sqrt_sym($B x $B->transpose) x $B;
}

# ---------------------------------------------------------------------
# whiten: center per channel, PCA-whiten, optional rank reduction.
#   returns hashref: Z (nt,k) whitened sources-basis data,
#                    wh (nch,k) whitening (math k x nch),
#                    dewh (k,nch) de-whitening (math nch x k),
#                    mean (nch), evals (k), k
# opts: keep => n  (keep top-n PCs), rank_tol => frac of max eigenvalue
# ---------------------------------------------------------------------
sub whiten {
    my ($X, %o) = @_;
    my $nt  = $X->dim(0);
    my $nch = $X->dim(1);
    my $mean = $X->average;                    # (nch), mean over time
    my $Xc = $X - $mean->dummy(0);
    my $C = ($Xc x $Xc->transpose) / $nt;      # (nch,nch)
    my ($E, $l) = eigens_sym($C);              # E columns = eigenvectors

    # sort eigenvalues descending
    my $ord = $l->qsorti->slice('-1:0');
    $l = $l->index($ord);
    $E = $E->dice_axis(0, $ord);               # reorder eigenvector columns (dim0)

    my $k = $nch;
    if (defined $o{keep}) { $k = $o{keep}; }
    elsif (defined $o{rank_tol}) {
        my $thr = $o{rank_tol} * $l->max;
        $k = ($l > $thr)->sum->sclr;
        $k = 1 if $k < 1;
    }
    my $lk = $l->slice("0:".($k-1))->copy;
    my $Ek = $E->slice("0:".($k-1).",:")->copy;   # (k, nch): top-k eigenvec columns

    my $ih = $lk ** -0.5;
    my $wh   = $Ek->transpose * $ih->dummy(0); # PDL (nch,k) = math (k x nch)
    # de-whitening (math nch x k) = E_k diag(l^1/2); PDL (k,nch), scale dim0(components)
    my $dewh = $Ek * ($lk ** 0.5);
    my $Z = $wh x $Xc;                         # (nt, k)

    return {
        Z => $Z, wh => $wh, dewh => $dewh,
        mean => $mean, evals => $lk, k => $k,
        Ek => $Ek, lk => $lk,
    };
}

# ---------------------------------------------------------------------
# fastica: symmetric FastICA on whitened data.
#   returns orthonormal rotation B (k,k), iters, conv
# opts: nl => 'tanh'|'gauss', maxit, tol, seed
# ---------------------------------------------------------------------
sub fastica {
    my ($Z, %o) = @_;
    my $k  = $Z->dim(1);
    my $N  = $Z->dim(0);
    my $maxit = $o{maxit} // 300;
    my $tol   = $o{tol}   // 1e-9;
    my $nl    = $o{nl}    // 'tanh';
    srand($o{seed}) if defined $o{seed};   # PDL exports srand; do NOT use PDL::srand (not a package sub in all builds)

    my $B = sym_orth(grandom($k, $k));
    my ($it, $conv) = (0, -1);
    for my $i (1 .. $maxit) {
        $it = $i;
        my $S = $B x $Z;                        # (N,k)
        my ($G, $Gp);
        if ($nl eq 'gauss') {
            my $u2 = $S * $S; my $e = exp(-0.5 * $u2);
            $G = $S * $e; $Gp = (1 - $u2) * $e;
        } else {
            $G = tanh($S); $Gp = 1 - $G * $G;
        }
        my $beta = $Gp->average;                # (k), mean over time
        my $Bn = ($G x $Z->transpose) / $N - $B * $beta->dummy(0);
        $Bn = sym_orth($Bn);
        my $dg = ($Bn x $B->transpose)->diagonal(0,1)->abs;
        $conv = 1 - $dg->min->sclr;
        $B = $Bn;
        last if $conv < $tol;
    }
    return ($B, $it, $conv);
}

# ---------------------------------------------------------------------
# ica_decompose: full pipeline. X (nt,nch) -> decomposition hashref
#   unmix (nch,k), mix (k,nch), sources (nt,k), plus whitening info.
# opts passed through to whiten() and fastica().
# ---------------------------------------------------------------------
sub ica_decompose {
    my ($X, %o) = @_;
    my $w = whiten($X, %o);
    my ($B, $it, $conv) = fastica($w->{Z}, %o);
    my $unmix = $B x $w->{wh};                  # (nch,k)  math k x nch
    # mix = pinv(unmix) = dewh x B^T   (math nch x k) -> PDL (k,nch)
    my $mix = $w->{dewh} x $B->transpose;       # (k,nch) = math nch x k
    my $S = $unmix x ($X - $w->{mean}->dummy(0)); # (nt,k)
    return {
        unmix => $unmix, mix => $mix, sources => $S,
        mean => $w->{mean}, k => $w->{k},
        iters => $it, conv => $conv,
    };
}

# ---------------------------------------------------------------------
# identify_by_reference: rank components by |correlation| with a
#   reference signal (e.g. vEOG/hEOG or a known time course).
#   sources (nt,k), ref (nt) -> list of [idx, corr] sorted by |corr| desc
# ---------------------------------------------------------------------
sub identify_by_reference {
    my ($sources, $ref) = @_;
    my $k = $sources->dim(1);
    my $r0 = $ref - $ref->average;
    my $rn = $r0 / sqrt(($r0*$r0)->average);
    my @out;
    for my $j (0 .. $k-1) {
        my $c = $sources->slice(":,($j)");
        my $c0 = $c - $c->average;
        my $cn = $c0 / sqrt(($c0*$c0)->average + 1e-30);
        push @out, [$j, (($cn*$rn)->average)->sclr];
    }
    return sort { abs($b->[1]) <=> abs($a->[1]) } @out;
}

# ---------------------------------------------------------------------
# ica_operator: build the (nch,nch) linear cleaning operator P such that
#   X_clean = (X - mean) P^T + mean  when zeroing components in @bad.
#   Here returned as removal in source space for exact back-projection.
# apply_ica: X (nt,nch), ica hashref, \@bad -> cleaned X (nt,nch)
#   X_clean = X - mix[:,bad] x unmix[bad,:] x Xc   (back-projection)
# ---------------------------------------------------------------------
sub apply_ica {
    my ($X, $ica, $bad) = @_;
    return $X->copy unless @$bad;
    my $idx = pdl(long, @$bad);
    my $Xc = $X - $ica->{mean}->dummy(0);
    my $S  = $ica->{unmix} x $Xc;               # (nt,k)
    my $Sbad   = $S->dice_axis(1, $idx);        # (nt,nbad)
    my $mixbad = $ica->{mix}->dice_axis(0, $idx); # (nbad,nch)
    my $removed = $mixbad x $Sbad;              # (nt,nch)
    return $X - $removed;
}

1;

__END__

=head1 NAME

PDL::EEG::ICA - FastICA blind source separation for EEG artifact removal (pure PDL)

=head1 DESCRIPTION

Symmetric fixed-point FastICA (Hyvarinen) on PCA-whitened multichannel EEG,
with reference-based component identification and back-projection removal.
Uses only PDL::MatrixOps (LAPACK not required). Data layout is (nt, nch),
matching PDL::EEG::GED. Removal is X - mix[:,bad] (unmix[bad,:] Xc), the same
back-projection form used by the GED spatial filter, so the two methods are
directly comparable on identical data and metrics.

=cut
