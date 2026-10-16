use strict;
use warnings;
use Test::More;
use PDL;

# PDL::EEG::SynthNoise: reproducible noise for synthetic data.
# The AR(1) reference below is a plain Perl loop written from the definition
#   x[0] = e[0] / sqrt(1 - phi^2),  x[i] = phi * x[i-1] + e[i]
# followed by zero mean and SD rescaling (divisor nt); it shares no code with
# the block-wise implementation in the module.

use_ok('PDL::EEG::SynthNoise', qw(seed_rng gaussian uniform synth_noise));

sub ref_ar1_scaled {
    my ($e, $phi, $sd) = @_;                # $e: Perl array of innovations
    my @x = ($e->[0] / sqrt(1 - $phi * $phi));
    push @x, $phi * $x[-1] + $e->[$_] for 1 .. $#$e;
    my $n = @x;
    my $m = 0; $m += $_ for @x; $m /= $n;
    $_ -= $m for @x;
    my $v = 0; $v += $_ * $_ for @x; $v = sqrt($v / $n);
    return [map { $_ * $sd / $v } @x];
}

sub max_abs_diff {
    my ($x, $ref) = @_;
    my @a1 = $x->list;
    my $m = 0;
    for my $i (0 .. $#a1) { my $d = abs($a1[$i] - $ref->[$i]); $m = $d if $d > $m }
    return $m;
}

# ---- seed_rng: same seed -> same numbers; different seed -> different ----
seed_rng(5);  my $g1 = grandom(1000); my @r1 = map { CORE::rand() } 1 .. 5;
seed_rng(5);  my $g2 = grandom(1000); my @r2 = map { CORE::rand() } 1 .. 5;
seed_rng(6);  my $g3 = grandom(1000);
ok(all($g1 == $g2), 'seed_rng seeds PDL grandom');
is_deeply(\@r1, \@r2, 'seed_rng seeds CORE::rand');
ok(any($g1 != $g3), 'different seeds give different numbers');
eval { seed_rng(-1) };  like($@, qr/non-negative integer/, 'negative seed is rejected');
eval { seed_rng() };    like($@, qr/seed is required/, 'missing seed is rejected');

# ---- gaussian/uniform do not depend on the number of threads ----
{
    my $old = get_autopthread_targ();
    my $n = 2_000_000;                      # above the default autothreading size
    my %got;
    for my $targ (1, 2, 4, 8) {
        set_autopthread_targ($targ);
        seed_rng(9); my $x = gaussian($n);
        seed_rng(9); my $u = uniform($n);
        $got{$targ} = [$x, $u];
    }
    set_autopthread_targ($old);
    ok(!grep({ !all($got{$_}[0] == $got{1}[0]) } 2, 4, 8) &&
       !grep({ !all($got{$_}[1] == $got{1}[1]) } 2, 4, 8),
       'gaussian and uniform: same numbers with 1, 2, 4, 8 threads');
    is(get_autopthread_targ(), $old, 'autothreading setting restored');
    seed_rng(9); my $short = gaussian(1000);
    ok(all($got{1}[0]->slice('0:999') == $short), 'a longer draw starts with the shorter one');
    my $u = $got{1}[1];
    ok($u->min >= 0 && $u->max < 1, 'uniform is in [0, 1)');
}

# ---- white noise: zero mean, exact SD, per channel; dims ----
{
    my $x = synth_noise(5000, nch => 4, sd => 3, seed => 1);
    is_deeply([$x->dims], [5000, 4], 'dims (nt, nch)');
    my $m = $x->average;
    my $s = sqrt(($x * $x)->average);
    ok(all(abs($m) < 1e-12), 'each channel has zero mean');
    ok(all(abs($s - 3) < 1e-12), 'each channel has SD 3 exactly');
    is_deeply([synth_noise(100, seed => 1)->dims], [100], '1-D without nch');
    my $y = synth_noise(5000, nch => 4, sd => 3, seed => 1);
    ok(all($x == $y), 'same seed, same noise');
    my $lag1 = (($x->slice('1:-1') * $x->slice('0:-2'))->average / 9);
    ok(all(abs($lag1) < 0.05), 'white: lag-1 autocorrelation near 0 (got '
       . join(' ', map { sprintf '%.3f', $_ } $lag1->list) . ')');
}

# ---- AR(1): equals the plain loop on every sample ----
for my $case ([0.95, 20_000], [-0.5, 3000], [0.999, 5000], [0.3, 2], [0.9999, 400_000]) {
    my ($phi, $nt) = @$case;
    seed_rng(30); my $e = gaussian($nt);
    my $ref = ref_ar1_scaled([$e->list], $phi, 2.0);
    my $x = synth_noise($nt, kind => 'ar1', ar => $phi, sd => 2.0, seed => 30);
    my $d = max_abs_diff($x, $ref);
    ok($d < 1e-9, "ar1 phi=$phi nt=$nt equals the plain loop on all samples (max diff "
       . sprintf('%.1e', $d) . ')');
}

# ---- AR(1) with several channels: each channel is its own AR(1) ----
{
    my ($nt, $nch, $phi) = (8000, 3, 0.95);
    seed_rng(31); my $e = gaussian($nt, $nch);
    my $x = synth_noise($nt, nch => $nch, kind => 'ar1', ar => $phi, seed => 31);
    my $worst = 0;
    for my $c (0 .. $nch - 1) {
        my $ref = ref_ar1_scaled([$e->slice(":,($c)")->list], $phi, 1);
        my $d = max_abs_diff($x->slice(":,($c)"), $ref);
        $worst = $d if $d > $worst;
    }
    ok($worst < 1e-9, 'ar1 with 3 channels equals the plain loop per channel');
    my $lag1 = ($x->slice('1:-1') * $x->slice('0:-2'))->average;
    ok(all(abs($lag1 - $phi) < 0.03), 'ar1: lag-1 autocorrelation near phi (got '
       . join(' ', map { sprintf '%.3f', $_ } $lag1->list) . ')');
}

# ---- uniform innovations ----
{
    seed_rng(32); my $u = uniform(10_000);
    my $e = [map { ($_ - 0.5) * sqrt(12) } $u->list];
    my $ref = ref_ar1_scaled($e, 0, 1);      # phi = 0: just rescaling (x[0] = e[0])
    my $x = synth_noise(10_000, dist => 'uniform', seed => 32);
    ok(max_abs_diff($x, $ref) < 1e-12, 'white uniform equals rescaled (u - 0.5)');
}

# ---- threshold: uniform white stays below k*1.4826*MAD, Gaussian does not ----
sub crossings {
    my ($v, $k) = @_;
    my $vd = $v - $v->avg;
    my $thr = $k * 1.4826 * median(abs($vd - median($vd)));
    return (abs($vd) > $thr)->sum->sclr;
}
{
    my $n_uni = 0;
    $n_uni += crossings(synth_noise(20_000, dist => 'uniform', seed => $_), 2.5) for 1 .. 10;
    is($n_uni, 0, 'white uniform: no sample above 2.5 x 1.4826 x MAD (10 seeds)');
    my $n_gau = crossings(synth_noise(20_000, seed => 1), 2.5);
    ok($n_gau > 100, "negative control: white Gaussian crosses it ($n_gau samples)");
    my $n_ar = crossings(synth_noise(20_000, kind => 'ar1', dist => 'uniform', seed => 1), 2.5);
    ok($n_ar > 0, "negative control: AR(1) with uniform innovations crosses it ($n_ar samples)");
}

# ---- argument checks ----
for my $bad ([[1], qr/nt must be/], [[100, kind => 'pink'], qr/kind must be/],
             [[100, dist => 'cauchy'], qr/dist must be/], [[100, kind => 'ar1', ar => 1], qr/\|ar\| < 1/],
             [[100, sd => 0], qr/sd must be/], [[100, nch => 0], qr/nch must be/],
             [[100, color => 'white'], qr/unknown option/]) {
    my ($args, $re) = @$bad;
    eval { synth_noise(@$args) };
    like($@, $re, "rejects (@$args)");
}

done_testing();
