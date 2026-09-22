package PDL::EEG::GED;

# GED (Generalized Eigenvalue Decomposition) spatial filtering for EEG,
# with a blink-artifact removal workflow built on top.
#
# Core:   solve  S w = lambda R w  for a pair of channel covariances and
#         remove a chosen component subspace by oblique back-projection.
# Blink:  build S from a peak-normalised blink-averaged evoked and R from
#         blink-free background (both from the SAME recording), identify the
#         blink component by its correlation with the vEOG-evoked time course,
#         and remove it.
#
# The generalized eigenproblem is solved with PDL::LinearAlgebra (LAPACK
# sygvd), whose eigenvectors are R-orthonormal (W' R W = I). That makes the
# forward patterns A = R W biorthogonal to the filters (W' A = I), so removing
# a component subset is an exact oblique projection:
#     X_clean = X - A_sel (W_sel' X).
#
# Data layout everywhere: X is (nch, nt).

use strict;
use warnings;
use PDL;
use PDL::LinearAlgebra qw(msymgeigen);
use Exporter 'import';

our @EXPORT_OK = qw(
    ged_cov ged_operator apply_ged
    detect_blinks blink_free blink_evoked fit_blink remove_blinks
);

# ===========================================================================
# Core GED
# ===========================================================================

# ged_cov($X) -> (nch,nch) time-mean-removed sample covariance.  $X : (nch,nt)
sub ged_cov {
    my ($X) = @_;
    my ($nch, $nt) = $X->dims;
    my $Xi   = $X->transpose;                     # (nt, nch)  view
    my $mean = $Xi->avgover;                       # per-channel mean over time
    $Xi      = $Xi - $mean->dummy(0);              # centre
    my $C    = ($Xi x $Xi->transpose) / ($nt - 1); # (nch,nch)  BLAS
    return ($C + $C->transpose) / 2;               # symmetrise (guard fp)
}

# _shrink($R, $gamma) = (1-g) R + g (tr R / nch) I
sub _shrink {
    my ($R, $g) = @_;
    my $nch = $R->dim(0);
    my $mu  = $R->diagonal(0,1)->avg;
    return (1 - $g) * $R + $g * $mu * identity($nch);
}

# ged_operator($S, $R, %opt) -> operator hashref
#   opt: reg => shrinkage gamma on R (default 0.05)
#   returns { evals (ascending), W (filters), A (patterns), R (Rreg), nch }
sub ged_operator {
    my ($S, $R, %opt) = @_;
    my $g   = defined $opt{reg} ? $opt{reg} : 0.05;
    my $nch = $S->dim(0);
    my $Rr  = $g > 0 ? _shrink($R, $g) : $R;
    my ($evals, $W, $info) = msymgeigen($S, $Rr, 0, 1, 1, 'sygvd');
    die "PDL::EEG::GED: msymgeigen failed (info=$info)\n" if $info != 0;
    my $A = $Rr x $W;                              # forward patterns, W' A = I
    return { evals => $evals, W => $W, A => $A, R => $Rr, nch => $nch };
}

# apply_ged($op, $X, %opt) -> cleaned $X (nch,nt)
#   opt: idx => piddle/arrayref of component indices to remove
#        n   => remove the top-n (largest eigenvalue) components if idx absent (default 1)
sub apply_ged {
    my ($op, $X, %opt) = @_;
    my $nch = $op->{nch};
    my $idx;
    if (defined $opt{idx}) {
        $idx = ref($opt{idx}) eq 'ARRAY' ? pdl(long, @{$opt{idx}}) : $opt{idx}->long;
    } else {
        my $n = defined $opt{n} ? $opt{n} : 1;
        $idx  = ($nch - $n) + sequence(long, $n);  # top-n = last n (ascending evals)
    }
    my $mask = zeroes($nch);
    $mask->index($idx) .= 1;

    my $Xi   = $X->transpose;                      # (nt,nch)
    my $Y    = $op->{W}->transpose x $Xi;          # component activations (nt,nch)
    my $Ysel = $Y * $mask->dummy(0);               # keep only removed comps
    my $Xart = $op->{A} x $Ysel;                   # artifact reconstruction (nt,nch)
    return ($Xi - $Xart)->transpose->sever;        # (nch,nt)
}

# ===========================================================================
# Blink workflow
# ===========================================================================

# detect_blinks($veog, %opt) -> peak sample indices (1-D piddle)
#   opt: sfreq=>1000, thresh_k=>2.5 (x robust SD), refractory_ms=>200
sub detect_blinks {
    my ($veog, %opt) = @_;
    my $sf   = $opt{sfreq} // 1000;
    my $k    = $opt{thresh_k} // 2.5;
    my $refr = int(($opt{refractory_ms} // 200) * $sf / 1000);
    my $vd   = $veog - $veog->avg;
    my $base = 1.4826 * median(abs($vd - median($vd)));
    my $thr  = $k * $base;
    my $mag  = abs($vd);
    my $mask = ($mag > $thr);
    my $nt   = $veog->nelem;
    my $rise = which(($mask == 1) & ($mask->rotate(1) == 0));
    my (@pk); my $last = -1e9;
    for my $s ($rise->list) {
        next if $s - $last < $refr;
        my $e = $s; $e++ while $e < $nt - 1 && $mask->at($e + 1);
        push @pk, $s + $mag->slice("$s:$e")->maximum_ind->sclr;
        $last = $pk[-1];
    }
    return pdl(long, @pk);
}

# blink_free($veog, %opt) -> sample indices free of blinks (for R)
#   opt: sfreq=>1000, thresh_k=>4, guard_ms=>200
sub blink_free {
    my ($veog, %opt) = @_;
    my $sf  = $opt{sfreq} // 1000;
    my $k   = $opt{free_thresh_k} // 4;
    my $w   = int(($opt{guard_ms} // 200) * $sf / 1000);
    my $vd  = $veog - $veog->avg;
    my $thr = $k * 1.4826 * median(abs($vd - median($vd)));
    my $is  = (abs($vd) > $thr);
    my $nt  = $veog->nelem;
    my $c   = $is->double->cumusumover;
    my $lo  = (sequence($nt) - $w)->clip(0, $nt - 1);
    my $hi  = (sequence($nt) + $w)->clip(0, $nt - 1);
    return which((($c->index($hi) - $c->index($lo)) > 0) == 0);
}

# blink_evoked($X, $peaks, %opt) -> (nch,L) blink-averaged evoked
#   opt: half_ms=>150, sfreq=>1000, veog=>idx, normalize=>1 (peak-normalise each
#        epoch by its vEOG peak before averaging = shape-based, amplitude-robust)
sub blink_evoked {
    my ($X, $peaks, %opt) = @_;
    my $sf   = $opt{sfreq} // 1000;
    my $half = int(($opt{half_ms} // 150) * $sf / 1000);
    my $norm = defined $opt{normalize} ? $opt{normalize} : 1;
    my $vi   = $opt{veog};
    my ($nch, $nt) = $X->dims;
    my $sum = zeroes($nch, 2 * $half + 1);
    my $n   = 0;
    for my $p ($peaks->list) {
        my $a = $p - $half; my $b = $p + $half;
        next if $a < 0 || $b >= $nt;
        my $ep = $X->slice(":,$a:$b")->sever;
        if ($norm && defined $vi) {
            my $vpk = abs($ep->slice("($vi),"))->max;
            $ep /= ($vpk > 1e-6 ? $vpk : 1);
        }
        $sum += $ep; $n++;
    }
    die "PDL::EEG::GED: no usable blink epochs\n" if $n == 0;
    my $evk = $sum / $n;
    return wantarray ? ($evk, $n) : $evk;
}

# fit_blink($X, %opt) -> blink GED operator
#   $X : (nch,nt) recording to be cleaned; R (background) is taken from $X.
#   opt: veog=>idx (REQUIRED),
#        fit_source=>$Y  data whose blinks define S (default: $X itself; pass a
#                        dedicated voluntary-blink recording with the SAME channel
#                        order to fit cross-session),
#        reg=>0.05, plus detect/evoked/free options.
#   Adds fields: peaks, n_blinks, veog_corr (per component), blink_comp (index).
#   S = covariance of the peak-normalised blink-averaged evoked of fit_source.
#   R = blink-free covariance of $X.
sub fit_blink {
    my ($X, %opt) = @_;
    my $vi = $opt{veog};
    die "PDL::EEG::GED::fit_blink: veog=>index required\n" unless defined $vi;
    my $src = defined $opt{fit_source} ? $opt{fit_source} : $X;

    my $veog_src = $src->slice("($vi),")->sever;
    my $peaks    = detect_blinks($veog_src, %opt);
    die "PDL::EEG::GED::fit_blink: too few blinks (" . $peaks->nelem . ")\n"
        if $peaks->nelem < 8;

    my $evk  = blink_evoked($src, $peaks, %opt, veog => $vi, normalize => 1);  # S
    my $free = blink_free($X->slice("($vi),")->sever, %opt);                   # R from $X
    my $op   = ged_operator(ged_cov($evk), ged_cov($X->dice_axis(1, $free)),
                            reg => (defined $opt{reg} ? $opt{reg} : 0.05));

    # identify blink component: correlation of each filter's output on the
    # (un-normalised) fit_source blink evoked with its vEOG-evoked time course.
    my $evk_raw = blink_evoked($src, $peaks, %opt, veog => $vi, normalize => 0);
    my $vevk = $evk_raw->slice("($vi),")->flat;
    my $vy   = $vevk - $vevk->avg;
    my $nch  = $op->{nch};
    my $cc   = zeroes($nch);
    for my $k (0 .. $nch - 1) {
        my $w    = $op->{W}->slice("($k),")->flat;
        my $comp = ($w->dummy(1) * $evk_raw)->sumover;
        my $cx   = $comp - $comp->avg;
        $cc->set($k, abs(sum($cx * $vy)) / sqrt(sum($cx**2) * sum($vy**2) + 1e-30));
    }
    $op->{peaks}      = $peaks;
    $op->{n_blinks}   = $peaks->nelem;
    $op->{veog_corr}  = $cc;
    $op->{blink_comp} = $cc->maximum_ind->sclr;
    return $op;
}

# remove_blinks($X, %opt) -> ($subtracted, $removed, $op)
#   High-level entry. $X : (nch,nt) recording to clean.
#   opt: veog=>idx (REQUIRED), fit=>$voluntary (optional data for S, same montage),
#        n=>1 (remove n most blink-like components), plus fit options.
#   $subtracted : cleaned $X (nch,nt);  $removed = $X - $subtracted (the artifact).
sub remove_blinks {
    my ($X, %opt) = @_;
    my %fopt = %opt;
    $fopt{fit_source} = delete $fopt{fit} if defined $fopt{fit};
    my $op  = fit_blink($X, %fopt);
    my $n   = defined $opt{n} ? $opt{n} : 1;
    my $idx;
    if ($opt{components}) { $idx = pdl(long, @{$opt{components}}); }
    else { $idx = $op->{veog_corr}->qsorti->slice("-$n:-1")->long; }  # top-n by vEOG corr
    my $sub = apply_ged($op, $X, idx => $idx);
    my $removed = $X - $sub;
    return wantarray ? ($sub, $removed, $op) : $sub;
}

1;

__END__

=head1 NAME

PDL::EEG::GED - GED spatial filtering and eye-blink removal for EEG

=head1 SYNOPSIS

    use PDL::EEG::GED qw(remove_blinks);

    # $task, $voluntary : (nch, nt) piddles, identical channel order.
    # vEOG is channel index $i in that montage.
    my ($subtracted, $removed, $op) =
        remove_blinks($task, veog => $i, fit => $voluntary);

    # $subtracted (nch,nt) = cleaned data;  $removed = $task - $subtracted.

=head1 DESCRIPTION

Generalized Eigenvalue Decomposition (GED) spatial filtering, with an
eye-blink removal workflow on top.

Given two channel covariances S and R it solves the generalized
symmetric-definite eigenproblem C<S w = lambda R w> and removes a chosen
component subspace from the data by oblique back-projection.

For blink removal, S is the covariance of a peak-normalised blink-averaged
evoked and R is the blink-free background covariance; the blink component is
identified by its correlation with the vEOG-evoked time course and removed.

The eigenproblem is solved with L<PDL::LinearAlgebra> (LAPACK C<sygvd>), whose
eigenvectors are R-orthonormal (C<W' R W = I>). The forward patterns C<A = R W>
are then biorthogonal to the filters (C<W' A = I>), so removing a component
subset is an exact oblique projection C<X_clean = X - A_sel (W_sel' X)>. R is
regularised by shrinkage C<(1-g) R + g (tr R / nch) I> (default C<g = 0.05>).

The blink covariance S is built from the blink-averaged evoked with each epoch
peak-normalised by its vEOG peak before averaging, making the estimated blink
pattern independent of individual blink amplitude.

Data layout everywhere is C<(nch, nt)>.

=head1 FUNCTIONS

Core:

    $C   = ged_cov($X);                       # (nch,nt) -> (nch,nch)
    $op  = ged_operator($S, $R, reg => 0.05);
    $out = apply_ged($op, $X, n => 1);        # or idx => [...]

C<$op> is a hashref: C<evals> (ascending), C<W> (filters), C<A> (patterns),
C<R> (regularised), C<nch>.

Blink workflow:

    $peaks = detect_blinks($veog, sfreq => 1000, thresh_k => 2.5, refractory_ms => 200);
    $idx   = blink_free($veog, sfreq => 1000);
    $evk   = blink_evoked($X, $peaks, veog => $i, normalize => 1);
    $op    = fit_blink($X, veog => $i, fit_source => $Y);  # $Y provides S; R from $X
    ($sub, $rem, $op) = remove_blinks($X, veog => $i, fit => $Y, n => 1);

C<remove_blinks> is the high-level entry: C<$X> (nch,nt) is the recording to
clean, C<fit> is an optional dedicated voluntary-blink recording (same channel
order) used to estimate the blink pattern S; the background R is always taken
from C<$X>. C<$sub> is the cleaned data and C<$rem = $X - $sub> is the removed
artifact. The blink component is identified by correlation with the vEOG-evoked
time course (C<blink_comp>, C<veog_corr>).

=head1 REQUIREMENTS

L<PDL>, L<PDL::LinearAlgebra> (LAPACK C<sygvd> via C<msymgeigen>).

=head1 SEE ALSO

L<PDL::EEG::IO::EDF> (read_edf / clean_edf_label),
F<examples/blink_ged_clean.pl>, F<t/18_ged.t>.

=cut
