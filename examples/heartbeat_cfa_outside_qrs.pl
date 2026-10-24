#!/usr/bin/env perl
# examples/heartbeat_cfa_outside_qrs.pl
#
#   perl examples/heartbeat_cfa_outside_qrs.pl REC.edf [options]
#
#     --bn1 L --bn2 L   non-cephalic electrodes; R peaks are detected on BN1-BN2
#     --skip S          ignore the first S seconds                       (10)
#     --band LO,HI      zero-phase Butterworth band, order 4 [Hz]        (0.5,40)
#     --kqrs K          size of the QRS subspace that is removed         (3)
#     --slow M          Late-window principal components removed first  (1)
#     --nsurr N         non-heartbeat averages for the band              (200)
#     --seed N          seed for the surrogate offsets                   (1)
#     --exclude LIST     head channels to leave out, comma separated     (none)
#
# What is left of the R-locked average after the cardiac QRS field is taken out?
# Head channels (19 scalp + A1, A2), average reference, every beat (no rejection,
# so the average stays in the data's row space), baseline -75..-50 ms.
#
#  1. the first M principal components of the Late window (450..650 ms) of the
#     average are projected out of the data (the slow wave that runs through the
#     whole cardiac cycle);
#  2. the first K principal components of the QRS window (-45..55 ms) of what
#     is left are projected out (the QRS field);
#  3. the leading principal component u of the T window (150..450 ms) of the
#     remainder: its pattern, its time course u'R(t) every 25 ms, the pointwise
#     95th percentile of |u'S(t)| over the non-heartbeat averages S (every R
#     moved by 0.25..0.60 s, same projections), and the BN1-BN2 average (the
#     ECG T wave) next to it; the correlation of the two over 50..450 ms; the
#     largest |u'R| over 50..700 ms against the 95th percentile of the same
#     maximum over the non-heartbeat averages.
#
# u is chosen on the heartbeat average itself, so its time course is biased
# toward the T window; the non-heartbeat band does not include that selection.

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

my ($bn1, $bn2, $skip, $band, $kq, $slow, $nsurr, $seed) = ('BN1', 'BN2', 10, '0.5,40', 3, 1, 200, 1);
my $exclude = '';
GetOptions('exclude=s' => \$exclude, 'bn1=s' => \$bn1, 'bn2=s' => \$bn2, 'skip=f' => \$skip, 'band=s' => \$band,
           'kqrs=i' => \$kq, 'slow=i' => \$slow, 'nsurr=i' => \$nsurr, 'seed=i' => \$seed) or die;
my $file = shift @ARGV or die "usage: $0 REC.edf [options]\n";
my ($f_lo, $f_hi) = split /,/, $band;

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
my $H  = $X->slice("0:" . ($nh - 1) . ",:");
$H     = $H - $H->average->dummy(0);
my $nt = $H->dim(1);
my $ecg = $X->slice("($nh),:") - $X->slice("(" . ($nh + 1) . "),:");

my ($n0, $n1) = (_s(-400), _s(700));
my ($r) = detect_r($ecg, $fs);
$r = $r->where(($r >= _s(1000 * $skip) - $n0) & ($r + $n1 <= $nt));
my %w  = (base => [-75, -50], QRS => [-45, 55], T => [150, 450], Late => [450, 650]);
my %wi = map { $_ => [ _s($w{$_}[0]) - $n0, _s($w{$_}[1]) - $n0 - 1 ] } keys %w;

my $A = _mean($H, $r);
my $Rm = $A;
my $Hs = $H;
my @proj;                                                   # projectors applied in order
if ($slow) {
    my (undef, $V) = _pcs(_win($A, 'Late'));
    push @proj, $V->slice("0:" . ($slow - 1) . ",:");
    $Rm = _po($Rm, $proj[-1]);
}
my (undef, $VQ) = _pcs(_win($Rm, 'QRS'));
push @proj, $VQ->slice("0:" . ($kq - 1) . ",:");
my $R1 = $Rm;
$Rm = _po($Rm, $proj[-1]);
my (undef, $VT) = _pcs(_win($Rm, 'T'));
my $u  = $VT->slice("0:0,:");
$u    *= -1 if $u->flat->index($u->flat->abs->maximum_ind)->sclr < 0;
my $tc = ($Rm x $u)->flat;

my $stc = cat(map {
    my $S = _mean($H, surrogate_r($r, $nt, $fs, seed => $seed + $_, lo => -$n0, hi => $nt - $n1));
    $S = _po($S, $_) for @proj;
    ($S x $u)->flat
} 0 .. $nsurr - 1);                                         # (nep, nsurr)
my $bandp = $stc->abs->mv(1, 0)->pctover(0.95);            # (nep)
my $ecgA  = _mean($ecg->dummy(0), $r)->flat;

printf "%s: band %s Hz, %d beats, Late PCs removed %d, QRS subspace k=%d removed\n",
    $file, $band, $r->nelem, $slow, $kq;
print "pattern: ", join(' ', map { sprintf '%s %.2f', $head[$_], $u->flat->at($_) } 0 .. $nh - 1), "\n";
printf "share of the T-window power (after the Late PCs) outside the QRS subspace: %.3f\n",
    ((_win($Rm, 'T')**2)->sum / (_win($R1, 'T')**2)->sum)->sclr;
print "lag(ms)  comp(uV)  |nonheart|95  ECG BN1-BN2(uV)\n";
for (my $ms = -150; $ms <= 675; $ms += 25) {
    my $i = _s($ms) - $n0;
    printf "%6d  %8.2f  %10.2f  %12.1f%s\n", $ms, $tc->at($i), $bandp->at($i), $ecgA->at($i),
        abs($tc->at($i)) > $bandp->at($i) ? '  *' : '';
}
my ($c0, $c1) = (_s(50) - $n0, _s(450) - $n0);
my $ca = $tc->slice("$c0:$c1");   $ca = $ca - $ca->avg;
my $cb = $ecgA->slice("$c0:$c1"); $cb = $cb - $cb->avg;
printf "corr(comp, ECG) over 50..450 ms: %.3f; peak |comp| at %d ms, peak |ECG| (100..450 ms) at %d ms\n",
    (($ca * $cb)->sum / sqrt(($ca**2)->sum * ($cb**2)->sum))->sclr,
    1000 * ($tc->slice("$c0:$c1")->abs->maximum_ind->sclr + $c0 + $n0) / $fs,
    1000 * ($ecgA->slice((_s(100) - $n0) . ":$c1")->abs->maximum_ind->sclr + _s(100)) / $fs;
my ($m0, $m1) = (_s(50) - $n0, _s(700) - $n0 - 1);
printf "max |comp| over 50..700 ms: %.2f uV; 95th percentile of that maximum over the non-heartbeat averages: %.2f uV\n",
    $tc->slice("$m0:$m1")->abs->max->sclr,
    $stc->slice("$m0:$m1")->abs->maximum->pct(0.95)->sclr;

# ---------------------------------------------------------------------------
sub _s    { my $ms = shift; my $v = $ms * 1e-3 * $fs; $v >= 0 ? int($v + 0.5) : -int(-$v + 0.5) }
sub _win  { my ($Z, $n) = @_; $Z->slice(":,$wi{$n}[0]:$wi{$n}[1]") }
# average over every beat (nh, nep), each beat baselined
sub _mean {
    my ($Y, $rp) = @_;
    my $E = beat_epochs($Y, $rp, $n0, $n1)->double;
    $E -= $E->slice(":,$wi{base}[0]:$wi{base}[1]")->mv(1, 0)->average->dummy(1);
    return $E->mv(2, 0)->average;
}
# principal components of Z (nh, ns): values desc, vectors as columns (nh x nh)
sub _pcs {
    my $Z = shift;
    my ($v, $V) = eigh_sym($Z->transpose x $Z);
    my $ord = sequence(long, $v->nelem)->slice("-1:0");
    return ($v->index($ord), $V->dice_axis(0, $ord));
}
sub _po { my ($Z, $Q) = @_; return $Z - ($Z x $Q) x $Q->transpose }
