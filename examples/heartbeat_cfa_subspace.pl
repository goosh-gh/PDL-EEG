#!/usr/bin/env perl
# examples/heartbeat_cfa_subspace.pl
#
#   perl examples/heartbeat_cfa_subspace.pl REC.edf [options]
#
#     --bn1 L --bn2 L    non-cephalic electrodes; R peaks are detected on BN1-BN2
#     --skip S           ignore the first S seconds                      (10)
#     --band LO,HI       zero-phase Butterworth band, order 4 [Hz]       (0.5,100)
#     --reject UV        per-channel peak-to-peak limit over the epoch   (150)
#     --baseline A,B     baseline window [ms], removed beat by beat      (-75,-50)
#     --windows LIST     NAME:A:B[,...] latency windows [ms]
#                        (P:-250:-90,QRS:-45:55,T:150:450,Late:450:650)
#     --ref LIST         rec (as recorded), avg (average of the head channels) (rec,avg)
#     --nsurr N          surrogate (non-heartbeat) averages               (100)
#     --seed N           seed for the surrogate offsets                   (1)
#     --kmax K           largest subspace size reported                   (6)
#     --exclude LIST     head channels to leave out, comma separated     (none)
#     --slow M           also repeat everything with the first M principal
#                        components of the Late window of the average removed
#                        from the head data (the slow, pulse-like wave)    (1)
#     --out PREFIX       output prefix                 (REC without .edf + _cfasub)
#
# How many dimensions does the cardiac field occupy on the head, and do the
# P, QRS and T waves share them?  For the head channels (19 scalp + A1, A2) the
# script computes, in each reference:
#
#  1. averages: for every latency window, the spatial covariance of the
#     R-locked average, and its generalized eigenvalues against the noise
#     covariance of the non-heartbeat averages (every R moved by 0.25..0.60 s,
#     same windows).  lambda_k is the power of the k-th spatial component in
#     units of the noise in that component; lambda near 1 is noise.
#     The fraction of the window's power carried by the first k principal
#     components (no whitening) is printed next to it.
#  2. subspace sharing: the first k principal components of the QRS window of
#     the average span a subspace; for every other window the share of its
#     power inside that subspace is printed for k = 1..kmax, next to the same
#     share for the non-heartbeat averages (what a subspace of that size takes
#     from data unrelated to the heart).  The cosines of the principal angles
#     between the k-dimensional QRS and T subspaces follow.
#     For each of the first three principal components of P, T and Late, the
#     share of its (unit) pattern that lies inside the k-dimensional QRS
#     subspace (1 = inside, 0 = orthogonal).
#  3. single beats: the QRS windows of all beats stacked, against the Late
#     window of all beats stacked (background EEG); generalized eigenvalues.
#     The number of components well above 1 is how many dimensions the QRS
#     field needs beat by beat (movement of the heart, e.g. with breathing,
#     shows up as more dimensions here than in the average).
#
# With --slow M > 0 everything is repeated after the first M principal
# components of the Late window of the average have been projected out of the
# head data (the slow wave that runs through the whole cardiac cycle; it is
# not in the QRS subspace and otherwise dominates the T, Late and P windows).
#
# Output: PREFIX_<ref>_patterns.tsv  (the QRS and T principal components, per
#         head channel, for topographic maps)

use strict;
use warnings;
use PDL;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use PDL::EEG::IO::EDF   qw(read_edf);
use PDL::EEG::Filter    qw(zp_butter);
use PDL::EEG::Heartbeat qw(detect_r beat_epochs surrogate_r);
use PDL::EEG::Linalg    qw(eigh_sym);

my ($bn1, $bn2, $skip, $band, $reject, $bline) = ('BN1', 'BN2', 10, '0.5,100', 150, '-75,-50');
my ($wins, $refs, $nsurr, $seed, $kmax, $slow, $out)
    = ('P:-250:-90,QRS:-45:55,T:150:450,Late:450:650', 'rec,avg', 100, 1, 6, 1);
my $exclude = '';
GetOptions('exclude=s' => \$exclude, 'bn1=s' => \$bn1, 'bn2=s' => \$bn2, 'skip=f' => \$skip, 'band=s' => \$band,
           'reject=f' => \$reject, 'baseline=s' => \$bline, 'windows=s' => \$wins,
           'ref=s' => \$refs, 'nsurr=i' => \$nsurr, 'seed=i' => \$seed,
           'kmax=i' => \$kmax, 'slow=i' => \$slow, 'out=s' => \$out) or die;
my $file = shift @ARGV or die "usage: $0 REC.edf [options]\n";
($out = $file) =~ s/\.edf$//i, $out .= '_cfasub' unless defined $out;
my ($f_lo, $f_hi) = split /,/, $band;
my @bl  = split /,/, $bline;
my @win = map { my ($n, $a0, $a1) = split /:/; { name => $n, ms => [$a0, $a1] } } split /,/, $wins;
my %W   = map { $_->{name} => $_ } @win;
exists $W{QRS} && exists $W{Late} or die "--windows must contain QRS and Late\n";

# ---- read, select, filter --------------------------------------------------
my $e  = read_edf($file);
my $fs = $e->{fs};
my %ix; @ix{ @{ $e->{labels} } } = 0 .. $#{ $e->{labels} };
my %excl = map { $_ => 1 } split /,/, $exclude;
my @head = grep { exists $ix{$_} && !$excl{$_} }
           qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 C3 Cz C4 T4 T5 P3 Pz P4 T6 O1 O2 A1 A2);
exists $ix{$_} or die "channel '$_' not in $file\n" for $bn1, $bn2;
my $nh = @head;
my $X  = zp_butter($e->{data}->dice_axis(0, pdl(long, [ @ix{@head, $bn1, $bn2} ])),
                   fs => $fs, low => $f_lo, high => $f_hi, order => 4);
my $nt = $X->dim(1);
my $H  = $X->slice("0:" . ($nh - 1) . ",:");
my $ecg = $X->slice("($nh),:") - $X->slice("(" . ($nh + 1) . "),:");

# ---- R peaks and epochs ----------------------------------------------------
my ($r) = detect_r($ecg, $fs);
my $tmin = _ms2s(-400); my $tmax = _ms2s(700);
my ($n0, $n1) = (_rnd($tmin * $fs), _rnd($tmax * $fs));
my $lo = _rnd($skip * $fs) - $n0;
$r = $r->where(($r >= $lo) & ($r + $n1 <= $nt));
my $rr = ($r->slice("1:-1") - $r->slice("0:-2"))->double / $fs;
printf "%s: fs %g Hz, %d head channels, %d beats, RR median %.3f s (%.3f..%.3f)\n",
    $file, $fs, $nh, $r->nelem, $rr->median, $rr->pct(0.05), $rr->pct(0.95);
printf "band %s Hz, baseline %s ms, reject %g uV, windows %s\n", $band, $bline, $reject,
    join(' ', map { "$_->{name} $_->{ms}[0]..$_->{ms}[1]" } @win);
my $nep = $n1 - $n0;
$_->{i} = [ _rnd($_->{ms}[0] * 1e-3 * $fs) - $n0, _rnd($_->{ms}[1] * 1e-3 * $fs) - $n0 - 1 ] for @win;
my @bi  = (_rnd($bl[0] * 1e-3 * $fs) - $n0, _rnd($bl[1] * 1e-3 * $fs) - $n0 - 1);

# surrogate points: far enough from the record edges for a whole epoch
my @sur = map { surrogate_r($r, $nt, $fs, seed => $seed + $_, lo => -$n0, hi => $nt - $n1) } 0 .. $nsurr - 1;

for my $ref (split /,/, $refs) {
  my $H0 = $ref eq 'avg' ? $H - $H->average->dummy(0)
         : $ref eq 'rec' ? $H
         : die "--ref: unknown reference '$ref'\n";
  for my $m ($slow > 0 ? (0, $slow) : (0)) {
    my $Hr = $H0;
    print "\n==== reference: ", ($ref eq 'avg' ? 'average of the head channels' : 'as recorded'),
          ($m ? ", first $m Late PC(s) of the average removed" : ''), "\n";
    if ($m) {
        my ($A0) = _average($H0, $r);
        $A0 = _to_rowspace($A0, $H0);       # per-channel rejection breaks the zero sum of 'avg'
        my ($l0, $l1) = @{ $W{Late}{i} };
        my (undef, $QL) = _pca($A0->slice(":,$l0:$l1"));
        $QL = $QL->slice("0:" . ($m - 1) . ",:");                   # nh x m
        $Hr = $H0 - ($H0 x $QL) x $QL->transpose;
        for my $j (0 .. $m - 1) {
            my $q  = $QL->slice("($j),:");
            my $sg = $q->index($q->abs->maximum_ind)->sclr < 0 ? -1 : 1;
            my $tc = ($A0 x $QL->slice("$j:$j,:"))->flat * $sg;
            print "  removed Late PC", $j + 1, " pattern: ",
                join(' ', map { sprintf '%s %.2f', $head[$_], $sg * $q->at($_) } 0 .. $nh - 1), "\n";
            print "  its time course in the average [uV], every 50 ms from ", 1000 * $tmin, " ms: ",
                join(' ', map { sprintf '%.1f', $tc->at($_) } grep { $_ < $nep } map { $_ * _rnd(0.05 * $fs) } 0 .. 30), "\n";
        }
    }

    my ($A, $E, $good) = _average($Hr, $r);
    printf "beats used per channel: %d..%d\n", $good->mv(1,0)->sumover->minmax;
    my @S = map { [ _average($Hr, $_) ]->[0] } @sur;           # (nh,nep) each
    my $Sall = cat(@S);                                         # (nh,nep,nsurr)

    # rank of the head data in this reference (system reference -> C3, C4 tied)
    my $rank  = _rank($Hr);
    printf "rank of the head channels: %d of %d\n", $rank, $nh;
    # averages with per-channel rejection leave the data's row space slightly
    # (with 'avg' they no longer sum to zero over channels): project them back
    $A    = _to_rowspace($A, $Hr);
    $Sall = _to_rowspace($Sall->clump(1, 2), $Hr)->reshape($nh, $nep, $nsurr);

    # ---- 1. averages: GED against the surrogate noise --------------------
    print "\n-- 1. average: generalized eigenvalues against non-heartbeat averages (lambda_1..)\n";
    print "      and cumulative share of the window power in the first k PCs\n";
    my %pc;
    for my $w (@win) {
        my ($i0, $i1) = @{ $w->{i} };
        my $Xw  = $A->slice(":,$i0:$i1");
        my $Cn  = _cov($Sall->slice(":,$i0:$i1,:")->clump(1, 2));
        my $lam = _ged(_cov($Xw), $Cn, $rank);
        my ($val, $vec) = _pca($Xw);
        $pc{ $w->{name} } = $vec;
        my $cum = $val->cumusumover / $val->sum;
        printf "  %-5s lambda %s\n", $w->{name}, join(' ', map { sprintf '%7.1f', $_ } $lam->slice("0:" . ($kmax - 1))->list);
        printf "  %-5s share  %s\n", '', join(' ', map { sprintf '%7.3f', $_ } $cum->slice("0:" . ($kmax - 1))->list);
    }

    # ---- 2. sharing of the QRS subspace ----------------------------------
    print "\n-- 2. share of each window's power inside the first k PCs of the QRS average\n";
    printf "  %-12s %s\n", 'window', join(' ', map { sprintf '   k=%d', $_ } 1 .. $kmax);
    for my $w (grep { $_->{name} ne 'QRS' } @win) {
        my ($i0, $i1) = @{ $w->{i} };
        for my $src (['heart', $A], ['nonheart', $Sall]) {
            my $Xw = $src->[1]->slice(":,$i0:$i1");
            $Xw = $Xw->clump(1, 2) if $Xw->ndims == 3;
            my $tot = ($Xw**2)->sum;
            my @sh  = map {
                my $Q = $pc{QRS}->slice("0:" . ($_ - 1) . ",:");    # nh x k
                my $P = $Xw x $Q;                                   # (ns x nh)(nh x k)
                (($P**2)->sum / $tot)->sclr
            } 1 .. $kmax;
            printf "  %-5s %-6s %s\n", $w->{name}, $src->[0], join(' ', map { sprintf '%6.3f', $_ } @sh);
        }
    }
    for my $k (1 .. 3) {
        next unless exists $pc{T};
        my $M = $pc{QRS}->slice("0:" . ($k - 1) . ",:")->transpose x $pc{T}->slice("0:" . ($k - 1) . ",:");
        my ($c2) = eigh_sym($M x $M->transpose);
        printf "  cos(principal angles) QRS vs T, k=%d: %s\n", $k,
            join(' ', map { sprintf '%.3f', sqrt($_ > 0 ? $_ : 0) } reverse $c2->list);
    }

    print "  share of each principal component's pattern inside the k-dim QRS subspace:\n";
    for my $wn (grep { $_ ne 'QRS' && exists $pc{$_} } map { $_->{name} } @win) {
        for my $j (0 .. 2) {
            my $v = $pc{$wn}->slice("$j:$j,:");                       # nh x 1
            printf "  %-5s PC%d %s\n", $wn, $j + 1, join(' ', map {
                sprintf '%6.3f', ((($v->transpose x $pc{QRS}->slice("0:" . ($_ - 1) . ",:"))**2)->sum)->sclr
            } 1 .. $kmax);
        }
    }

    # ---- 3. single beats: QRS against Late (background) ------------------
    my ($q0, $q1) = @{ $W{QRS}{i} };
    my ($l0, $l1) = @{ $W{Late}{i} };
    my $gb  = $good->andover->which;                               # beats good in every channel
    my $Eg  = $E->dice_axis(2, $gb);
    my $Cq  = _cov($Eg->slice(":,$q0:$q1,:")->clump(1, 2));
    my $Cl  = _cov($Eg->slice(":,$l0:$l1,:")->clump(1, 2));
    my $lam = _ged($Cq, $Cl, $rank);
    printf "\n-- 3. single beats (%d): QRS window against Late window, lambda_1..: %s\n",
        $gb->nelem, join(' ', map { sprintf '%.2f', $_ } $lam->slice("0:" . ($kmax + 2))->list);
    # the same with the average beat removed: what varies from beat to beat
    my $D   = $Eg - $Eg->mv(2, 0)->average->dummy(2);
    my $lamv = _ged(_cov($D->slice(":,$q0:$q1,:")->clump(1, 2)),
                    _cov($D->slice(":,$l0:$l1,:")->clump(1, 2)), $rank);
    printf "   beat-to-beat part (average beat removed): %s\n",
        join(' ', map { sprintf '%.2f', $_ } $lamv->slice("0:" . ($kmax + 2))->list);

    # ---- patterns --------------------------------------------------------
    my $fn = "${out}_${ref}" . ($m ? "_slow$m" : '') . "_patterns.tsv";
    open my $fh, '>', $fn or die "$fn: $!";
    my @cols = map { my $w = $_; map { "${w}_pc$_" } 1 .. 3 } grep { exists $pc{$_} } qw(QRS T);
    print $fh join("\t", 'ch', @cols), "\n";
    for my $c (0 .. $nh - 1) {
        print $fh join("\t", $head[$c], map { my ($w, $k) = /^(\w+?)_pc(\d)$/;
            sprintf '%.5f', $pc{$w}->at($k - 1, $c) } @cols), "\n";
    }
    close $fh;
    print "wrote $fn\n";
  }
}

# ---------------------------------------------------------------------------
# null space of the head data Y (nh, nt) in its first minute: eigenvectors of the
# covariance below 1e-9 of the largest (avg reference: the constant vector; after
# --slow: also the removed patterns). Returns nh x nnull (or undef) and the rank.
sub _null {
    my $Y = shift;
    my ($ev, $V) = eigh_sym(_cov($Y->slice(":,0:" . (_min(60 * $fs, $Y->dim(1)) - 1))));
    my $nn = ($ev < 1e-9 * $ev->max)->sum->sclr;
    return ($nn ? $V->slice("0:" . ($nn - 1) . ",:") : undef, $nh - $nn);
}
sub _rank { return (_null(shift))[1] }
sub _to_rowspace {
    my ($Z, $Y) = @_;
    my ($N) = _null($Y);
    return defined $N ? $Z - ($Z x $N) x $N->transpose : $Z;
}
sub _min { $_[0] < $_[1] ? $_[0] : $_[1] }

# ($avg (nh,nep), $E (nh,nep,nb) baseline-removed epochs, $good (nh,nb))
sub _average {
    my ($Y, $rp) = @_;
    my $E = beat_epochs($Y, $rp, $n0, $n1)->double;
    $E -= $E->slice(":,$bi[0]:$bi[1]")->mv(1, 0)->average->dummy(1);
    my $Em   = $E->mv(1, 0);
    my $good = ($Em->maximum - $Em->minimum) <= $reject;        # (nh,nb)
    my $n    = $good->mv(1, 0)->sumover->double;
    my $avg  = ($E * $good->dummy(1))->mv(2, 0)->sumover / $n->dummy(1)->clip(1, undef);
    return ($avg, $E, $good);
}

# spatial covariance (nh,nh) of Y (nh, nsamp), no centering (windows are baselined),
# per sample (divided by nsamp) so that windows of different lengths compare
sub _cov { my $Y = shift; return ($Y->transpose x $Y) / $Y->dim(1) }

# principal components of Y (nh, nsamp): (values desc (nh), Q (nh_comp, nh))
# Q is a matrix with one component per column, largest first: component j at
# channel c is $Q->at($j, $c); its first k columns are $Q->slice("0:k-1,:")
sub _pca {
    my $Y = shift;
    my ($val, $vec) = eigh_sym(_cov($Y));                   # ascending, columns
    my $ord = sequence(long, $val->nelem)->slice("-1:0");
    return ($val->index($ord)->clip(0, undef), $vec->dice_axis(0, $ord));
}

# generalized eigenvalues of (Ca, Cb), descending, in the top-$rank subspace of Cb
sub _ged {
    my ($Ca, $Cb, $rank) = @_;
    my ($vb, $Vb) = eigh_sym($Cb);
    my $n   = $vb->nelem;
    my $sel = sequence(long, $rank) + ($n - $rank);         # largest $rank
    my $Wm  = $Vb->dice_axis(0, $sel) / sqrt($vb->index($sel))->dummy(1);   # columns scaled
    my ($lam) = eigh_sym($Wm->transpose x $Ca x $Wm);
    return $lam->slice("-1:0");
}

sub _rnd  { my $v = shift; $v >= 0 ? int($v + 0.5) : -int(-$v + 0.5) }
sub _ms2s { $_[0] * 1e-3 }
