use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::Filter qw(:all);

my $PI = 4 * atan2(1, 1);
my ($fs, $nt) = (500, 30000);
my $t = sequence(double, $nt) / $fs;

# ---- zp_butter -------------------------------------------------------------
# Deterministic synthetic ECG (no RNG), the signal of t/22_heartbeat.t: 64 beats
# with a slowly varying RR and R amplitude, a T wave 300 ms after each R,
# 0.25 Hz drift, 50 Hz and 7.3 Hz sines.
my $x = 0.3 * sin(2 * $PI * 0.25 * $t) + 0.05 * sin(2 * $PI * 50 * $t)
      + 0.02 * sin(2 * $PI * 7.3 * $t);
for my $k (0 .. 63) {
    my $c  = 1.0 + 0.9 * $k + 0.05 * sin(0.7 * $k);
    my $am = 1 + 0.1 * sin(1.3 * $k);
    $x += $am * exp(-0.5 * (($t - $c) / 0.012) ** 2)
        + 0.25 * exp(-0.5 * (($t - $c - 0.3) / 0.06) ** 2);
}
# Reference values for this signal from an independent implementation:
# scipy.signal.butter + filtfilt (scipy 1.18.1).
my @IDX = (5000, 9999, 15000, 20001, 25000);
my %REF = (
    bp => [-0.241612764, -0.176433326, -0.024888856,  0.004606886, -0.012554735],
    hp => [-0.050495358, -0.112319427,  0.025290246,  0.132529752, -0.004125916],
    lp => [ 0.022371938, -0.006903421,  0.102499107,  0.182370058,  0.071445484],
);
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
    eval { zp_butter($x, fs => $fs) };
    like($@, qr/give low and\/or high/, 'zp_butter: refuses a call without a corner frequency');
}

# ---- analytic_signal / envelope / inst_phase --------------------------------
# 10 Hz carrier, amplitude 1 + 0.5 sin(2 pi 0.7 t), plus a 23.3 Hz sine.
my $carrier = 2 * $PI * 10 * $t + 0.3;
my $am      = 1 + 0.5 * sin(2 * $PI * 0.7 * $t);
my $z = $am * cos($carrier) + 0.2 * sin(2 * $PI * 23.3 * $t);
{
    # Reference values from scipy.signal.hilbert (scipy 1.18.1), which takes the
    # record as one period; the first two and the last index are at the edges.
    my @J  = (0, 3, 4321, 15017, 20001, 29999);
    my @IM = (0.095520207, 0.507037070, 0.358131828, 0.644149352, 0.223253373, -0.018796449);
    my @EN = (0.960099952, 1.071275789, 1.037667318, 1.198888142, 0.997793543, 0.922975635);
    my ($re, $im) = analytic_signal($z, exact => 1);
    is_deeply([$re->dims], [$nt], 'analytic_signal: (nt) in, (nt) out');
    cmp_ok(($re - $z)->abs->max, '<', 1e-10, 'analytic_signal exact: real part is the input');
    my $j = pdl(long, \@J);
    cmp_ok(($im->index($j) - pdl(\@IM))->abs->max, '<', 1e-8, 'analytic_signal exact: imaginary part = scipy hilbert, edges included');
    cmp_ok((envelope($z, exact => 1)->index($j) - pdl(\@EN))->abs->max, '<', 1e-8, 'envelope exact = |scipy hilbert|');
    # odd length (2999 samples): the other branch of the frequency weights
    my @JO  = (0, 7, 1500, 2998);
    my @IMO = (0.491804255, 1.099669432, 0.220632783, 0.404988116);
    my ($ro, $io) = analytic_signal($z->slice("0:2998"), exact => 1);
    cmp_ok(($io->index(pdl(long, \@JO)) - pdl(\@IMO))->abs->max, '<', 1e-8, 'analytic_signal exact, odd length = scipy hilbert');
    cmp_ok(($ro - $z->slice("0:2998"))->abs->max, '<', 1e-10, 'analytic_signal exact, odd length: real part is the input');
}
{
    # whole cycles of one cosine: the Hilbert transform is the sine
    my $c = cos(2 * $PI * 10 * $t);
    my ($re, $im) = analytic_signal($c, exact => 1);
    cmp_ok(($im - sin(2 * $PI * 10 * $t))->abs->max, '<', 1e-9, 'analytic_signal exact: cos -> sin');
    my ($rc, $ic) = analytic_signal(ones($nt) * 3, exact => 1);
    cmp_ok($ic->abs->max, '<', 1e-9, 'analytic_signal exact: a constant has no imaginary part');
}
{
    # default mode: amplitude and phase of the carrier alone
    my $s   = $am * cos($carrier);
    my $in  = "2500:27499";                                   # 5 s away from both ends
    my $env = envelope($s);
    my $ph  = inst_phase($s);
    cmp_ok(($env - $am)->slice($in)->abs->max, '<', 1e-3, 'envelope: recovers the amplitude modulation (interior)');
    my $d = $ph - $carrier; $d = atan2(sin($d), cos($d));
    cmp_ok($d->slice($in)->abs->max, '<', 1e-3, 'inst_phase: recovers the carrier phase (interior)');
    my $depth = sub { my $e = shift->slice($in); ($e->max - $e->min) / ($e->max + $e->min) };
    cmp_ok(abs($depth->($env) - 0.5), '<', 1e-3, 'envelope: modulation depth 0.5');
    # negative control: no modulation in, none out
    cmp_ok($depth->(envelope(cos($carrier))), '<', 1e-3, 'envelope: an unmodulated carrier gives a flat envelope');
    my ($re) = analytic_signal($s);
    cmp_ok(($re - $s)->abs->max, '<', 1e-10, 'analytic_signal: real part is the input');
    # a channel offset changes neither the envelope nor the imaginary part
    my ($r5, $i5) = analytic_signal($s + 5);
    my ($r0, $i0) = analytic_signal($s);
    cmp_ok(($i5 - $i0)->abs->max, '<', 1e-9, 'analytic_signal: an offset leaves the imaginary part unchanged');
    cmp_ok(($r5 - $s - 5)->abs->max, '<', 1e-9, 'analytic_signal: real part keeps the offset');
}
{
    my $s  = $am * cos($carrier);
    my $X2 = cat($s, -2 * $s)->xchg(0, 1);                    # (2, nt)
    my ($r2, $i2) = analytic_signal($X2);
    my ($r1, $i1) = analytic_signal($s);
    is_deeply([$i2->dims], [2, $nt], 'analytic_signal: (nch,nt) in, (nch,nt) out');
    cmp_ok(($i2->slice("(0)") - $i1)->abs->max, '<', 1e-12, 'analytic_signal: row 0 = 1-D result');
    cmp_ok(($i2->slice("(1)") + 2 * $i1)->abs->max, '<', 1e-12, 'analytic_signal: linear in the input');
    is_deeply([envelope($X2)->dims], [2, $nt], 'envelope: (nch,nt)');
    is_deeply([inst_phase($X2)->dims], [2, $nt], 'inst_phase: (nch,nt)');
    # a length with a large prime factor (29989 is prime) is filled with zeros to 30000
    my ($rp, $ip) = analytic_signal($s->slice("0:29988"));
    is($ip->nelem, 29989, 'analytic_signal: prime length kept');
    cmp_ok(($ip - $i1->slice("0:29988"))->slice("2500:24999")->abs->max, '<', 5e-3, 'analytic_signal: prime length agrees with the full record in the interior');
    eval { analytic_signal(zeroes(2, 3, 4)) };
    like($@, qr/must be \(nt\) or \(nch,nt\)/, 'analytic_signal: refuses 3-D input');
}

done_testing;
