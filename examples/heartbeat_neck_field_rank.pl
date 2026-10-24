#!/usr/bin/env perl
# examples/heartbeat_neck_field_rank.pl
#
# How many independent time courses does the cardiac field have on the chest,
# along the back of the neck, and on the head?  And how many chest components
# are needed to explain the cardiac field on the head?
#
#   perl examples/heartbeat_neck_field_rank.pl REC.edf [options]
#
#     --sets LIST     NAME:CH,CH,...;NAME:...   electrode sets
#                     (chest:ClavL,ClavR,BN2,BN1,Th12;neck:Th12,BN1,Cv45,Cv12;head:HEAD)
#                     HEAD = the 19 scalp channels + A1, A2 present in the file
#     --body LIST     sets whose components are used to explain the head       (chest,neck)
#     --bn1 L --bn2 L R peaks are detected on BN1-BN2                           (BN1, BN2)
#     --band LO,HI    zero-phase Butterworth band, order 4 [Hz]                (0.5,100)
#     --baseline A,B  baseline window [ms] (PR segment), removed beat by beat  (-75,-50)
#     --windows LIST  NAME:A:B[,...] [ms]  (all:-250:450,P:-250:-90,QRS:-45:55,T:150:450)
#     --reject-head UV  per-channel peak-to-peak limit over the epoch, head    (150)
#     --reject-body UV  the same for channels that are not in the head set     (3000)
#     --skip S        ignore the first S seconds                               (10)
#     --edge S        leave out beats within S s of a recording break          (2)
#     --nsurr N       surrogate (non-heartbeat) averages                       (50)
#     --seed N        seed for the surrogate offsets                           (1)
#     --kmax K        largest number of components reported                    (5)
#     --out PREFIX    output prefix                         (REC without .edf + _neckrank)
#
# Recording breaks are the "Block N" annotations; a block that starts with a
# "REC START ... CAL" annotation (calibration) is left out entirely.
#
# Every set is re-referenced to its own average (the system reference and any
# signal common to all channels of the set drop out), so a set of n channels
# has at most n-1 dimensions.
#
#  1. For every set and window: the share of the R-locked average's power in
#     the first k principal components (no whitening), and the generalized
#     eigenvalues against the covariance of the non-heartbeat averages
#     (lambda near 1 = noise level).  The number of lambdas well above 1 is the
#     number of independent cardiac time courses the set sees.
#  2. Head explained by the body sets, window by window: the first k principal
#     time courses of a body set are regressed out of the head average; printed
#     is the share of the head's power they explain, next to the share explained
#     by the head's own first k components (the best any k time courses can do)
#     and the share they explain in the non-heartbeat averages (chance).  The
#     best single body channels and pairs follow (against the system reference,
#     and re-referenced to the set average).  The head's first two PC patterns
#     are printed for each window.
#  3. |cos| between the head's PC patterns of different windows.
#
# Output: PREFIX_patterns.tsv       principal components (window 1) of every set,
#         PREFIX_timecourses.tsv    their time courses in the average,
#         PREFIX_head_patterns.tsv  the head's first two PC patterns in every window
#                                   (sign: the largest |value| positive).
# The sign of a principal component is arbitrary: in the printed head patterns and
# in PREFIX_head_patterns.tsv the channel with the largest |value| is made positive,
# so two windows can show the same topography with opposite signs.

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

my %o = (sets => 'chest:ClavL,ClavR,BN2,BN1,Th12;neck:Th12,BN1,Cv45,Cv12;head:HEAD', body => 'chest,neck',
         bn1 => 'BN1', bn2 => 'BN2', band => '0.5,100', baseline => '-75,-50',
         windows => 'all:-250:450,P:-250:-90,QRS:-45:55,T:150:450',
         reject_head => 150, reject_body => 3000, skip => 10, edge => 2, nsurr => 50, seed => 1, kmax => 5);
GetOptions('sets=s' => \$o{sets}, 'body=s' => \$o{body}, 'bn1=s' => \$o{bn1}, 'bn2=s' => \$o{bn2},
           'band=s' => \$o{band}, 'baseline=s' => \$o{baseline}, 'windows=s' => \$o{windows},
           'reject-head=f' => \$o{reject_head}, 'reject-body=f' => \$o{reject_body},
           'skip=f' => \$o{skip}, 'edge=f' => \$o{edge}, 'nsurr=i' => \$o{nsurr}, 'seed=i' => \$o{seed},
           'kmax=i' => \$o{kmax}, 'out=s' => \$o{out}) or die "bad options (see the header of $0)\n";
my $file = shift @ARGV or die "usage: $0 REC.edf [options]\n";
($o{out} = $file) =~ s/\.edf$//i, $o{out} .= '_neckrank' unless defined $o{out};

my @HEAD = qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 C3 Cz C4 T4 T5 P3 Pz P4 T6 O1 O2 A1 A2);

# ---- read --------------------------------------------------------------------
my $e  = read_edf($file);
my $fs = $e->{fs};
my $nt = $e->{data}->dim(1);
my %ix; @ix{ @{ $e->{labels} } } = 0 .. $#{ $e->{labels} };

my (@sets, %set);
for my $s (split /;/, $o{sets}) {
    my ($name, $list) = split /:/, $s, 2;
    my @ch = $list eq 'HEAD' ? (grep { exists $ix{$_} } @HEAD) : split /,/, $list;
    exists $ix{$_} or die "channel '$_' (set $name) not in $file\n" for @ch;
    @ch >= 2 or die "set $name needs at least 2 channels\n";
    push @sets, $name;
    $set{$name} = [@ch];
}
exists $set{head} or die "--sets needs a set named head\n";
my @body = split /,/, $o{body};
exists $set{$_} or die "--body: no set named $_\n" for @body;
exists $ix{$_} or die "channel '$_' not in $file\n" for $o{bn1}, $o{bn2};

# all channels used, once each; the head channels get the head rejection limit
my %headch = map { $_ => 1 } @{ $set{head} };
my (@use, %pos);
for my $c (map { @{ $set{$_} } } @sets) { next if exists $pos{$c}; $pos{$c} = @use; push @use, $c }
my ($f_lo, $f_hi) = split /,/, $o{band};
my $X = zp_butter($e->{data}->dice_axis(0, pdl(long, [ @ix{@use} ])),
                  fs => $fs, low => $f_lo, high => $f_hi, order => 4);
my $ecg = zp_butter($e->{data}->slice("($ix{$o{bn1}}),:") - $e->{data}->slice("($ix{$o{bn2}}),:"),
                    fs => $fs, low => $f_lo, high => $f_hi, order => 4);
my $limit = pdl(map { $headch{$_} ? $o{reject_head} : $o{reject_body} } @use);

# ---- recording breaks and calibration blocks -------------------------------
my @ev = sort { $a->{onset} <=> $b->{onset} } @{ $e->{events} // [] };
my @breaks = map { _rnd($_->{onset} * $fs) } grep { ($_->{label} // '') =~ /^Block \d+$/ && $_->{onset} > 0 } @ev;
my @recst  = grep { ($_->{label} // '') =~ /^REC START/ } @ev;
my @cal;                                                    # [start, end) samples
for my $j (0 .. $#recst) {
    next unless $recst[$j]{label} =~ /\bCAL\b/;
    my $s0 = _rnd($recst[$j]{onset} * $fs);
    my $s1 = $j < $#recst ? _rnd($recst[$j + 1]{onset} * $fs) : $nt;
    push @cal, [$s0, $s1];
}

# ---- R peaks, epochs ---------------------------------------------------------
my @win = map { my ($n, $a0, $a1) = split /:/; { name => $n, ms => [$a0, $a1] } } split /,/, $o{windows};
my ($n0, $n1) = (_rnd(-0.40 * $fs), _rnd(0.70 * $fs));
my $nep = $n1 - $n0;
$_->{i} = [ _rnd($_->{ms}[0] * 1e-3 * $fs) - $n0, _rnd($_->{ms}[1] * 1e-3 * $fs) - $n0 - 1 ] for @win;
$_->{i}[0] >= 0 && $_->{i}[1] < $nep or die "window $_->{name} outside -400..700 ms\n" for @win;
my @bl = split /,/, $o{baseline};
my @bi = (_rnd($bl[0] * 1e-3 * $fs) - $n0, _rnd($bl[1] * 1e-3 * $fs) - $n0 - 1);

my ($r) = detect_r($ecg, $fs);
my $nr_all = $r->nelem;
my $ok = ($r + $n0 >= _rnd($o{skip} * $fs)) & ($r + $n1 <= $nt);
my $edge = _rnd($o{edge} * $fs);
$ok &= (($r + $n1 <= $_ - $edge) | ($r + $n0 >= $_ + $edge)) for @breaks;
$ok &= (($r + $n1 <= $_->[0]) | ($r + $n0 >= $_->[1])) for @cal;
$r = $r->where($ok);
my $rr = ($r->slice("1:-1") - $r->slice("0:-2"))->double / $fs;
printf "%s: fs %g Hz, %.1f min, %d R peaks on %s-%s, %d kept (skip %g s, %d break(s) +-%g s, %d calibration block(s) out)\n",
    $file, $fs, $nt / $fs / 60, $nr_all, $o{bn1}, $o{bn2}, $r->nelem, $o{skip}, scalar @breaks, $o{edge}, scalar @cal;
printf "RR median %.3f s; band %s Hz, baseline %s ms, reject head %g / body %g uV p-p\n",
    $rr->median, $o{band}, $o{baseline}, $o{reject_head}, $o{reject_body};
printf "set %-6s %s\n", $_, join(' ', @{ $set{$_} }) for @sets;

my ($A, $good) = _average($X, $r);
printf "beats used per channel: %s\n", join(' ', map { "$use[$_] " . $good->slice("($_),:")->sum } 0 .. $#use);
my $lo = -$n0;
my @S = map { (_average($X, surrogate_r($r, $nt, $fs, seed => $o{seed} + $_, lo => $lo, hi => $nt - $n1)))[0] } 0 .. $o{nsurr} - 1;
my $Sall = cat(@S);                                         # (nch, nep, nsurr)

my $K = $o{kmax};
my $ms = (sequence($nep) + $n0) * (1000 / $fs);

# ---- 1. dimensions of every set ---------------------------------------------
print "\n-- 1. share of power in the first k PCs (no whitening) | lambda_k against non-heartbeat averages\n";
my (%pcs);
for my $sn (@sets) {
    my @c  = @{ $set{$sn} };
    my $Ax = _center($A->dice_axis(0, pdl(long, [ @pos{@c} ])));
    my $Sx = _center($Sall->dice_axis(0, pdl(long, [ @pos{@c} ])));
    my $k  = _min($K, @c - 1);
    printf "  %s (%d channels, rank %d after re-referencing to their average)\n", $sn, scalar @c, @c - 1;
    for my $w (@win) {
        my ($i0, $i1) = @{ $w->{i} };
        my ($val, $vec) = _pca($Ax->slice(":,$i0:$i1"));
        my $cum = $val->cumusumover / $val->sum;
        my $lam = _ged(_cov($Ax->slice(":,$i0:$i1")), _cov($Sx->slice(":,$i0:$i1,:")->clump(1, 2)), @c - 1);
        printf "    %-4s rms %7.2f uV | share %s | lambda %s\n", $w->{name},
            sqrt($val->sum / @c), join(' ', map { sprintf '%5.3f', $cum->at($_) } 0 .. $k - 1),
            join(' ', map { sprintf '%8.1f', $lam->at($_) } 0 .. $k - 1);
        $pcs{$sn} = { val => $val, vec => $vec, A => $Ax, S => $Sx } if $w == $win[0];
    }
}

# ---- 2. head explained by the body sets ------------------------------------
# PDL matrices are (columns, rows): an average window (nch, ns) is ns x nch, the PCs (ncomp, nch)
# are nch x ncomp, so  window x PCs  gives the time courses (ncomp, ns) = ns x ncomp, one per column.
my $HA = $pcs{head}{A};
my @hc = @{ $set{head} };
print "\n-- 2. head average explained by time courses of the body sets (share of the head's power)\n";
print "     chance = the same regressors on the non-heartbeat averages; best k = the head's own first k PCs\n";
for my $w (@win) {
    my ($i0, $i1) = @{ $w->{i} };
    my $H  = $HA->slice(":,$i0:$i1");                       # (nh, ns)
    my $HS = $pcs{head}{S}->slice(":,$i0:$i1,:");
    my ($hv, $hq) = _pca($H);
    printf "  window %s (%d..%d ms), head rms %.2f uV\n", $w->{name}, @{ $w->{ms} }, sqrt(($H**2)->avg);
    for my $j (0 .. 1) {
        my $q  = $hq->slice("($j),:");
        my $sg = $q->index($q->abs->maximum_ind)->sclr < 0 ? -1 : 1;
        printf "    head PC%d (%.3f): %s\n", $j + 1, ($hv->at($j) / $hv->sum),
            join(' ', map { sprintf '%s %+.2f', $hc[$_], $sg * $q->at($_) } 0 .. $#hc);
    }
    for my $bn (@body) {
        my @bc = @{ $set{$bn} };
        my $Bw = $pcs{$bn}{A}->slice(":,$i0:$i1");          # set-average reference
        my (undef, $bq) = _pca($Bw);
        my $tc = $Bw x $bq;
        printf "    %-26s %7s %7s %7s\n", "$bn PCs (this window)", 'share', 'chance', 'best k';
        for my $k (1 .. _min($K, @bc - 1)) {
            my $U = $tc->slice("0:" . ($k - 1) . ",:");     # (k, ns) = ns x k
            printf "    %-26s %7.3f %7.3f %7.3f\n", "  PC1..$k", _share($H, $U), _share_s($HS, $U),
                ($hv->slice("0:" . ($k - 1))->sum / $hv->sum)->sclr;
        }
        # single channels and pairs; against the system reference the channels also carry
        # what is common to the whole body set (lost in the set-average reference)
        my $Braw = $A->dice_axis(0, pdl(long, [ @pos{@bc} ]))->slice(":,$i0:$i1");
        for my $kind (['vs system reference', $Braw], ['set-average reference', $Bw]) {
            my @res;
            for my $p (0 .. $#bc) {
                push @res, [ $bc[$p], _share($H, $kind->[1]->slice("$p:$p,:")) ];
                push @res, [ "$bc[$p]+$bc[$_]", _share($H, $kind->[1]->dice_axis(0, pdl(long, [$p, $_]))) ]
                    for $p + 1 .. $#bc;
            }
            my @top = (sort { $b->[1] <=> $a->[1] } @res)[0 .. _min(5, scalar @res) - 1];
            printf "    %-26s %s\n", "  $kind->[0]", join('  ', map { sprintf '%s %.3f', @$_ } @top);
        }
    }
}

# head patterns: how far apart are the windows?  |cos| between unit PC patterns
print "\n-- 3. |cos| between the head's PC patterns of different windows (1 = same topography)\n";
my %hq = map { my ($i0, $i1) = @{ $_->{i} }; ($_->{name} => (_pca($HA->slice(":,$i0:$i1")))[1]) } @win;
for my $x (0 .. $#win) {
    for my $y ($x + 1 .. $#win) {
        my ($qa, $qb) = ($hq{ $win[$x]{name} }, $hq{ $win[$y]{name} });
        printf "  %-4s vs %-4s  %s\n", $win[$x]{name}, $win[$y]{name},
            join('  ', map { my ($ja, $jb) = @$_;
                sprintf 'PC%d.PC%d %.3f', $ja + 1, $jb + 1, abs(($qa->slice("($ja),:") * $qb->slice("($jb),:"))->sum) }
                ([0, 0], [0, 1], [1, 0], [1, 1]));
    }
}

# ---- output -------------------------------------------------------------------
open my $fp, '>', "$o{out}_patterns.tsv" or die "$o{out}_patterns.tsv: $!";
print $fp join("\t", qw(set channel), map { "PC$_" } 1 .. $K), "\n";
for my $sn (@sets) {
    my @c = @{ $set{$sn} };
    for my $ci (0 .. $#c) {
        print $fp join("\t", $sn, $c[$ci], map { sprintf '%.4f', $pcs{$sn}{vec}->at($_, $ci) } 0 .. _min($K, scalar @c) - 1), "\n";
    }
}
close $fp;
open my $ft, '>', "$o{out}_timecourses.tsv" or die "$o{out}_timecourses.tsv: $!";
print $ft join("\t", 'ms', map { my $sn = $_; map { "$sn.PC$_" } 1 .. _min($K, scalar @{ $set{$sn} }) } @sets), "\n";
my %T = map { $_ => ($pcs{$_}{A} x $pcs{$_}{vec}) } @sets;                 # (ncomp, nep)
for my $t (0 .. $nep - 1) {
    print $ft join("\t", sprintf('%.1f', $ms->at($t)),
        map { my $sn = $_; map { sprintf '%.3f', $T{$sn}->at($_, $t) } 0 .. _min($K, scalar @{ $set{$sn} }) - 1 } @sets), "\n";
}
close $ft;
# head PC patterns of every window (the input of nyhead_dipole_gof.pl)
open my $fh, '>', "$o{out}_head_patterns.tsv" or die "$o{out}_head_patterns.tsv: $!";
print $fh join("\t", 'channel', map { ("$_->{name}.PC1", "$_->{name}.PC2") } @win), "\n";
my @hp = map { my ($i0, $i1) = @{ $_->{i} }; my ($v, $q) = _pca($HA->slice(":,$i0:$i1"));
               map { my $c = $q->slice("($_),:"); $c->index($c->abs->maximum_ind)->sclr < 0 ? -$c : $c } 0, 1 } @win;
for my $ci (0 .. $#hc) {
    print $fh join("\t", $hc[$ci], map { sprintf '%.4f', $_->at($ci) } @hp), "\n";
}
close $fh;
print "\nwrote $o{out}_patterns.tsv, $o{out}_timecourses.tsv, $o{out}_head_patterns.tsv\n";

# =================================================================================
sub _average {                                  # (avg (nch,nep), good (nch,nb))
    my ($Y, $rp) = @_;
    my $E = beat_epochs($Y, $rp, $n0, $n1)->double;
    $E -= $E->slice(":,$bi[0]:$bi[1]")->mv(1, 0)->average->dummy(1);
    my $Em   = $E->mv(1, 0);
    my $good = ($Em->maximum - $Em->minimum) <= $limit->dummy(1);
    my $n    = $good->mv(1, 0)->sumover->double;
    my $avg  = ($E * $good->dummy(1))->mv(2, 0)->sumover / $n->dummy(1)->clip(1, undef);
    return ($avg, $good);
}
sub _center { my $Y = shift; return $Y - $Y->average->dummy(0) }     # subtract the mean over channels (dim 0)
sub _cov { my $Y = shift; return ($Y->transpose x $Y) / $Y->dim(1) }
sub _pca {                                      # (values desc, Q: column j = component j)
    my $Y = shift;
    my ($val, $vec) = eigh_sym(_cov($Y));
    my $ord = sequence(long, $val->nelem)->slice("-1:0");
    return ($val->index($ord)->clip(0, undef), $vec->dice_axis(0, $ord));
}
sub _ged {                                      # generalized eigenvalues, descending, top-$rank subspace of Cb
    my ($Ca, $Cb, $rank) = @_;
    my ($vb, $Vb) = eigh_sym($Cb);
    my $n   = $vb->nelem;
    my $sel = sequence(long, $rank) + ($n - $rank);
    my $Wm  = $Vb->dice_axis(0, $sel) / sqrt($vb->index($sel))->dummy(1);
    my ($lam) = eigh_sym($Wm->transpose x $Ca x $Wm);
    return $lam->slice("-1:0");
}
# share of Y's power explained by least squares on time courses U.
# Y (nch, ns) is the matrix ns x nch, U (k, ns) is ns x k; Yh = U (U'U)^-1 U'Y (no ns x ns matrix)
sub _share {
    my ($Y, $U) = @_;
    my $Yh = $U x (inv($U->transpose x $U) x ($U->transpose x $Y));
    return (($Yh**2)->sum / ($Y**2)->sum)->sclr;
}
sub _share_s {                                   # mean share over the surrogate averages (nch, ns, nsurr)
    my ($YS, $U) = @_;
    my @s = map { _share($YS->slice(":,:,($_)"), $U) } 0 .. $YS->dim(2) - 1;
    my $m = 0; $m += $_ for @s; return $m / @s;
}
sub _min { my $m = shift; $_ < $m and $m = $_ for @_; $m }
sub _rnd { my $v = shift; $v >= 0 ? int($v + 0.5) : -int(-$v + 0.5) }
