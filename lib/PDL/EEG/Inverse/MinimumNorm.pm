package PDL::EEG::Inverse::MinimumNorm;
use strict; use warnings;
use PDL;
use PDL::MatrixOps qw(identity eigens_sym stretcher);  # symmetric pinv, no LA dependency
use Exporter 'import';
our @EXPORT_OK = qw(
    forward_project centering
    avg_reference
    inverse_operator apply_inverse source_estimate source_power
);
our $VERSION = '0.01';

# ---------------------------------------------------------------------------
# Convention (explicit, because dim order is the classic silent bug source):
#   Leadfield $K stores the forward model b = K j with math shape (Ne x Ns):
#     Ne electrodes, Ns sources, one scalar per source (surface-normal
#     constrained, as in the New York Head V_fem_normal).
#   In PDL a math (rows x cols) matrix is stored (cols,rows). So $K has
#     PDL dims (Ns, Ne)  == PDL::IO::NYHead->leadfield ((Nsrc,Nchan)? no):
#   NOTE the NYHead reader returns V_fem_normal as PDL dims (231,74382) =
#     (electrode, source). That is the *transpose* of what this module wants.
#   Pass it through leadfield_for_inverse() below, or ->transpose, so that
#   $K here is (Ns, Ne). All entry points assert Ne < Ns and Ne == length(b),
#   which catches an accidentally-transposed leadfield immediately.
# ---------------------------------------------------------------------------

# Moore-Penrose pseudoinverse of a symmetric PSD matrix via eigendecomposition.
# Robust to the rank-(Ne-1) Gram that average-referenced (CAR) leadfields
# produce -- the constant direction is a true null and must be dropped, not
# inverted. rcond is relative to the largest eigenvalue.
sub _pinv {
    my ($C, $rcond) = @_;
    $rcond //= 1e-9;
    my ($ev, $e) = eigens_sym($C);
    my $tol  = $rcond * $e->abs->max;
    my $keep = ($e->abs > $tol);
    my $einv = $keep * (1.0/($e + (1-$keep)));   # 1/e where kept, 0 in null space
    return $ev x stretcher($einv) x $ev->transpose;
}
*_inv = \&_pinv;   # all call sites use the pseudoinverse

# average-reference (centering) operator H = I - 11^T/n  (n x n)
sub centering {
    my ($n) = @_;
    return identity($n) - ones($n,$n)/$n;
}

# average-reference a leadfield over its electrode axis: each source column
# is centred over electrodes (1^T k_s = 0). Use before building a car operator
# on an electrode subset, so K and data share the same reference.
#   $K dims (Ns,Ne) -> returns (Ns,Ne) with per-source electrode mean removed
sub avg_reference {
    my ($K) = @_;
    my $m = $K->xchg(0,1)->avgover;      # (Ns) per-source mean over electrodes
    return $K - $m;                     # broadcasts over dim0=Ns
}

# b = K j    ($K dims (Ns,Ne); $j dims (Ns) or (Ns,Nt)) -> (Ne) or (Ne,Nt)
sub forward_project {
    my ($K, $j) = @_;
    my ($Ns,$Ne) = $K->dims;
    my $jj = $j->getndims < 2 ? $j->dummy(0) : $j->xchg(0,1); # (1|Nt, Ns)
    my $b  = ($K x $jj);            # b = K j   ($K is (Ns,Ne) = math K)                    # (1|Nt, Ne)
    return $b->getndims==2 && $b->dim(0)==1 ? $b->slice('(0),:')->sever : $b->xchg(0,1)->sever;
}

sub _gram {
    my ($K, $alpha, $ref) = @_;
    my ($Ns,$Ne) = $K->dims;
    my $H = ($ref eq 'car') ? centering($Ne) : identity($Ne);
    my $Kt = $K->transpose;                  # (Ne,Ns)
    my $C  = ($K x $Kt) + $alpha*$H;         # (Ne,Ne)  == K K^T + aH
    return ($C, $H, $Kt);
}

sub _alpha {
    my ($K, %o) = @_;
    return $o{alpha} if defined $o{alpha};
    my ($Ns,$Ne) = $K->dims;
    my $tr = ($K*$K)->sum / $Ne;             # trace(K K^T)/Ne = mean eigen scale
    my $frac = defined $o{reg_frac} ? $o{reg_frac} : 0.05;
    return $frac * $tr;
}

# Build a data-independent inverse operator once; reuse across time points.
#   method => 'mne' | 'sloreta' | 'eloreta'
#   ref    => 'car' (default) | 'none'
#   alpha  => absolute regularizer, or reg_frac => fraction of mean eigenscale
#   (eloreta) max_iter (default 100), tol (default 1e-10)
# returns hashref: { method, ref, alpha, T (Ne,Ns), std (Ns) [sloreta],
#                    w (Ns), iters, rel [eloreta] }
sub inverse_operator {
    my ($K, %o) = @_;
    my $method = lc($o{method} // 'eloreta');
    my $ref    = lc($o{ref}    // 'car');
    my $alpha  = _alpha($K, %o);
    my ($Ns,$Ne) = $K->dims;
    die "leadfield looks transposed: expected (Ns,Ne) with Ns>Ne, got ($Ns,$Ne)"
        if $Ns < $Ne;

    if ($method eq 'mne' or $method eq 'sloreta') {
        my ($C,$H,$Kt) = _gram($K, $alpha, $ref);
        my $Ci = _inv($C);
        my $T  = $Kt x $Ci;                  # (Ne,Ns) = math Ns x Ne  (MNE operator)
        my %op = (method=>$method, ref=>$ref, alpha=>$alpha, T=>$T);
        if ($method eq 'sloreta') {
            # standardizer S_i = k_i^T C^{-1} k_i  (resolution diagonal)
            my $UU = $Ci x $K;               # (Ns,Ne): row s = C^{-1} k_s
            $op{std} = ($K*$UU)->xchg(0,1)->sumover;   # (Ns)
        }
        return \%op;
    }
    elsif ($method eq 'eloreta') {
        my $max = $o{max_iter} // 100;
        my $tol = $o{tol} // 1e-10;
        my $H  = ($ref eq 'car') ? centering($Ne) : identity($Ne);
        my $Kt = $K->transpose;              # (Ne,Ns)
        my $w  = ones($Ns);
        my ($it,$rel) = (0,1);
        for (1..$max) {
            my $Winv = 1.0/$w;
            my $M  = _inv( (($K*$Winv) x $Kt) + $alpha*$H );   # (K W^-1 K^T + aH)^-1
            my $UU = $M x $K;                                  # (Ns,Ne) row s = M k_s
            my $wn = sqrt( ($K*$UU)->xchg(0,1)->sumover );     # (Ns)
            $rel = (abs($wn-$w)->max)/$wn->max;
            $w = $wn; $it++;
            last if $rel < $tol;
        }
        my $Winv = 1.0/$w;
        my $M = _inv( (($K*$Winv) x $Kt) + $alpha*$H );
        my $T = ($Kt * $Winv->dummy(0)) x $M;   # (Ne,Ns) = W^-1 K^T M
        return { method=>'eloreta', ref=>$ref, alpha=>$alpha,
                 T=>$T, w=>$w, iters=>$it, rel=>$rel };
    }
    die "unknown method '$method'";
}

# Apply operator to scalp data b (Ne) or (Ne,Nt) -> source estimate (Ns) or (Ns,Nt).
# For sloreta this returns the *standardized* amplitude jhat/sqrt(std).
sub apply_inverse {
    my ($op, $b) = @_;
    my $T = $op->{T};                        # (Ne,Ns)
    my $bb = $b->getndims < 2 ? $b->dummy(0) : $b->xchg(0,1);  # (1|Nt, Ne)
    my $j  = ($T x $bb);                      # (1|Nt, Ns)
    $j = $j->dim(0)==1 ? $j->slice('(0),:')->sever : $j->xchg(0,1)->sever;  # (Ns)|(Ns,Nt)
    if ($op->{method} eq 'sloreta') {
        my $s = sqrt($op->{std});
        $j = $j->getndims<2 ? $j/$s : $j / $s->dummy(1);
    }
    return $j;
}

# convenience: operator + data in one call
sub source_estimate {
    my ($K, $b, %o) = @_;
    return apply_inverse( inverse_operator($K,%o), $b );
}

# quantity to colour the cortex with: power (>=0). For sloreta this is the
# standardized statistic^2; for mne/eloreta the current-density power.
sub source_power {
    my ($op_or_K, $b, @rest) = @_;
    my $j = ref($op_or_K) eq 'HASH' ? apply_inverse($op_or_K,$b)
                                    : source_estimate($op_or_K,$b,@rest);
    return $j*$j;
}
1;
