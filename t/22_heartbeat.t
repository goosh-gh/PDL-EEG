use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::Heartbeat qw(:all);

# Deterministic synthetic ECG (no RNG): 64 beats with a slowly varying RR and R
# amplitude, a T wave 300 ms after each R, 0.25 Hz drift, 50 Hz and 7.3 Hz sines.
my $PI = 4 * atan2(1, 1);
my ($fs, $nt) = (500, 30000);
my $t = sequence(double, $nt) / $fs;
my $x = 0.3 * sin(2 * $PI * 0.25 * $t) + 0.05 * sin(2 * $PI * 50 * $t)
      + 0.02 * sin(2 * $PI * 7.3 * $t);
my @true;
for my $k (0 .. 63) {
    my $c = 1.0 + 0.9 * $k + 0.05 * sin(0.7 * $k);
    my $a = 1 + 0.1 * sin(1.3 * $k);
    $x += $a * exp(-0.5 * (($t - $c) / 0.012) ** 2)
        + 0.25 * exp(-0.5 * (($t - $c - 0.3) / 0.06) ** 2);
    push @true, int($c * $fs + 0.5);
}

# Reference values for this signal from an independent implementation:
# scipy.signal.butter + filtfilt, and detect_r() of hrp_first100.py (scipy 1.18.1).
my @IDX = (5000, 9999, 15000, 20001, 25000);
my %REF = (
    bp => [-0.241612764, -0.176433326, -0.024888856,  0.004606886, -0.012554735],
    hp => [-0.050495358, -0.112319427,  0.025290246,  0.132529752, -0.004125916],
    lp => [ 0.022371938, -0.006903421,  0.102499107,  0.182370058,  0.071445484],
);
my @R_REF = qw(501 965 1424 1872 2310 2742 3179 3624 4083 4551 5015 5474 5922 6359
    6792 7229 7675 8134 8602 9066 9524 9972 10408 10841 11278 11725 12184 12652 13116
    13574 14022 14457 14891 15328 15775 16234 16702 17167 17624 18071 18506 18941
    19377 19825 20285 20752 21218 21674 22121 22555 22990 23427 23875 24335 24803
    25269 25724 26171 26605 27040 27476 27925 28385 28853);

# ---- zp_butter -------------------------------------------------------------
my %F = (bp => [low => 5, high => 20, order => 3],
         hp => [low => 0.5, order => 2],
         lp => [high => 40, order => 4]);
for my $k (sort keys %F) {
    my $y = zp_butter($x, fs => $fs, @{ $F{$k} });
    is($y->nelem, $nt, "zp_butter $k: length kept");
    my $d = ($y->index(pdl(long, \@IDX)) - pdl($REF{$k}))->abs->max;
    cmp_ok($d, '<', 1e-6, "zp_butter $k: interior samples match filtfilt (max diff $d)");
}
{
    my $X2 = cat($x, -2 * $x)->xchg(0, 1);                    # (2, nt)
    my $y2 = zp_butter($X2, fs => $fs, low => 5, high => 20, order => 3);
    is_deeply([$y2->dims], [2, $nt], 'zp_butter: (nch,nt) in, (nch,nt) out');
    my $y1 = zp_butter($x, fs => $fs, low => 5, high => 20, order => 3);
    cmp_ok(($y2->slice("(0)") - $y1)->abs->max, '<', 1e-12, 'zp_butter: row 0 = 1-D result');
    cmp_ok(($y2->slice("(1)") + 2 * $y1)->abs->max, '<', 1e-12, 'zp_butter: linear in the input');
    # a pure 10 Hz sine keeps the analytic Butterworth gain, a 1 Hz sine is removed
    my $g = sub {
        my $f = shift;
        my $y = zp_butter(sin(2 * $PI * $f * $t), fs => $fs, low => 5, high => 20, order => 3);
        sqrt(2 * ($y->slice("5000:24999") ** 2)->avg);
    };
    my ($w, $w1, $w2) = map { sin($PI * $_ / $fs) / cos($PI * $_ / $fs) } 10, 5, 20;
    my $q = ($w * $w - $w1 * $w2) / (($w2 - $w1) * $w);
    cmp_ok(abs($g->(10) - 1 / (1 + $q ** 6)), '<', 1e-6, 'zp_butter: 10 Hz gain = 1/(1+q^6)');
    cmp_ok($g->(1), '<', 1e-3, 'zp_butter: 1 Hz removed by the 5-20 Hz band');
}

# ---- detect_r --------------------------------------------------------------
{
    my ($r, $sign) = detect_r($x, $fs);
    is($sign, 1, 'detect_r: positive R');
    is_deeply([$r->list], \@R_REF, 'detect_r: same R samples as the reference implementation');
    cmp_ok(($r - pdl(long, \@true))->abs->max, '<=', 2, 'detect_r: within 2 samples of the true R times');
    my ($rn, $sn) = detect_r(-$x, $fs);
    is($sn, -1, 'detect_r: negative R recognised');
    is_deeply([$rn->list], \@R_REF, 'detect_r: inverted input gives the same samples');
    my ($rf) = detect_r($x, $fs, sign => 1);
    is_deeply([$rf->list], \@R_REF, 'detect_r: sign => 1');
}

# ---- beat_epochs / beat_average ---------------------------------------------
{
    my $r = pdl(long, \@R_REF);
    my $X = cat($x, 2 * $x + 1, sequence(double, $nt))->xchg(0, 1);      # (3, nt)
    my $E = beat_epochs($X, $r, -100, 200);
    is_deeply([$E->dims], [3, 300, 64], 'beat_epochs: (nch, nep, nbeat)');
    is($E->at(2, 0, 5), $r->at(5) - 100, 'beat_epochs: first sample is r+n0');
    is($E->at(2, 299, 63), $r->at(63) + 199, 'beat_epochs: last sample is r+n1-1');
    is($E->at(0, 100, 7), $x->at($r->at(7)), 'beat_epochs: lag 0 holds the R sample');
    eval { beat_epochs($X, pdl(long, [50]), -100, 200) };
    like($@, qr/runs off the record/, 'beat_epochs: refuses an epoch outside the record');

    my $A = beat_average($X, $r, $fs, tmin => -0.2, tmax => 0.4);
    is_deeply([$A->{avg}->dims], [3, 300], 'beat_average: avg is (nch, nep)');
    is($A->{t_ms}->at(0), -200, 'beat_average: t_ms starts at tmin');
    is($A->{t_ms}->at(100), 0, 'beat_average: t_ms 0 at the R sample');
    my $man = zeroes(double, 300);
    $man += $x->slice(($_ - 100) . ":" . ($_ + 199)) for @R_REF;
    $man /= 64;
    cmp_ok(($A->{avg}->slice("(0)") - $man)->abs->max, '<', 1e-12, 'beat_average: mean = explicit sum / n');
    cmp_ok(($A->{avg}->slice("(1)") - (2 * $man + 1))->abs->max, '<', 1e-12, 'beat_average: per channel');
    is($A->{n}->at(0), 64, 'beat_average: all beats counted');
    my $v = 0; $v += ($x->at($_) - $man->at(100)) ** 2 for @R_REF;
    cmp_ok(abs($A->{sd}->at(0, 100) - sqrt($v / 63)), '<', 1e-12, 'beat_average: sd over beats (n-1)');

    my $B = beat_average($X, $r, $fs, tmin => -0.2, tmax => 0.4, baseline => [-0.2, -0.1]);
    cmp_ok($B->{avg}->slice(":,0:49")->mv(1, 0)->average->abs->max, '<', 1e-9,
           'beat_average: baseline window averages to 0');
    cmp_ok(($B->{avg}->slice("(1)") - 2 * $B->{avg}->slice("(0)"))->abs->max, '<', 1e-9,
           'beat_average: baseline removes a channel offset');

    # one beat made huge in channel 0 only
    my $Y = $X->copy;
    $Y->slice("(0)," . ($r->at(10) + 20)) += 100;
    my $C = beat_average($Y, $r, $fs, tmin => -0.2, tmax => 0.4, reject => pdl(50, 1e30, 1e30));
    is_deeply([$C->{n}->list], [63, 64, 64], 'beat_average: per-channel rejection');
    my $man63 = ($man * 64 - $x->slice(($r->at(10) - 100) . ":" . ($r->at(10) + 199))) / 63;
    cmp_ok(($C->{avg}->slice("(0)") - $man63)->abs->max, '<', 1e-12,
           'beat_average: rejected beat left out of the mean');
    my $M = beat_average($Y, $r, $fs, tmin => -0.2, tmax => 0.4, stat => 'median');
    cmp_ok(abs($M->{avg}->at(0, 120) - $A->{avg}->at(0, 120)), '<', 0.05,
           'beat_average: median ignores one outlying beat');

    my $D = beat_average($X, pdl(long, [50, @R_REF, 29990]), $fs, tmin => -0.2, tmax => 0.4);
    is($D->{r}->nelem, 64, 'beat_average: beats that do not fit are dropped');
}

# ---- surrogate_r -----------------------------------------------------------
{
    my $r  = pdl(long, \@R_REF);
    my $s1 = surrogate_r($r, $nt, $fs, seed => 7);
    my $s2 = surrogate_r($r, $nt, $fs, seed => 7);
    my $s3 = surrogate_r($r, $nt, $fs, seed => 8);
    is_deeply([$s1->list], [$s2->list], 'surrogate_r: same seed, same points');
    isnt(join(',', $s1->list), join(',', $s3->list), 'surrogate_r: another seed differs');
    is($s1->nelem, 64, 'surrogate_r: one point per beat');
    ok($s1->min >= 0 && $s1->max < $nt, 'surrogate_r: inside the record');
    # every surrogate point is 0.25..0.60 s from some R
    my $d = ($s1->dummy(0) - $r->dummy(1))->abs->minimum->double / $fs;
    ok($d->max <= 0.6 + 1e-9, 'surrogate_r: no point farther than max from an R');
    my $s4 = surrogate_r($r, $nt, $fs, seed => 7, min => 0.25, max => 0.30);
    my $d4 = ($s4->dummy(0) - $r->dummy(1))->abs->minimum->double / $fs;
    ok($d4->min >= 0.25 - 1e-9 && $d4->max <= 0.30 + 1e-9, 'surrogate_r: offsets within [min, max]');
    my $s5 = surrogate_r($r, $nt, $fs, seed => 7, lo => 1000, hi => 28000);
    ok($s5->min >= 1000 && $s5->max < 28000, 'surrogate_r: lo/hi respected');
}

# ---- bn_balance ------------------------------------------------------------
{
    my $s  = sin(2 * $PI * 3 * $t->slice("0:999"));
    my $c  = cos(2 * $PI * 5 * $t->slice("0:999"));           # common part
    my ($a, $b) = (0.6 * $s + $c, -1.4 * $s + $c);            # 0.7*a + 0.3*b = c
    my ($w, $rms, $ra, $rb) = bn_balance($a, $b);
    cmp_ok(abs($w - 0.7), '<', 1e-9, 'bn_balance: weight');
    cmp_ok(abs($rms - sqrt(($c ** 2)->avg)), '<', 1e-9, 'bn_balance: residual = common part');
    cmp_ok(abs($ra - sqrt(($a ** 2)->avg)), '<', 1e-12, 'bn_balance: rms of input 1');
    my ($w2) = bn_balance($b, $a);
    cmp_ok(abs($w2 - 0.3), '<', 1e-9, 'bn_balance: swapping inputs gives 1-w');
    eval { bn_balance($a, $a) };
    like($@, qr/identical/, 'bn_balance: identical inputs refused');
}

done_testing;
