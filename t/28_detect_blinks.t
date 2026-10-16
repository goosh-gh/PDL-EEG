use strict;
use warnings;
use Test::More;
use PDL;

eval { require PDL::EEG::GED; PDL::EEG::GED->import('detect_blinks'); 1 }
    or plan skip_all => "PDL::EEG::GED not loadable (needs PDL::LinearAlgebra): $@";
use PDL::EEG::SynthNoise qw(synth_noise);

# Reference written from the description of detect_blinks, in plain Perl arrays:
#   supra-threshold samples: |v - mean(v)| > thresh_k * 1.4826 * MAD
#   a run starts at s when sample s is supra-threshold and the sample before it
#     is not; the sample before 0 is the last sample (the run test wraps)
#   the run continues while the following samples are supra-threshold, up to nt-1
#   peak = first sample of the largest |v - mean| in the run
#   a run starting less than refractory samples after the previous peak is skipped
sub ref_detect {
    my ($v, %opt) = @_;
    my $sf   = $opt{sfreq} // 1000;
    my $k    = $opt{thresh_k} // 2.5;
    my $refr = int(($opt{refractory_ms} // 200) * $sf / 1000);
    my $vd   = $v - $v->avg;
    my $thr  = ($k * 1.4826 * median(abs($vd - median($vd))))->sclr;
    my @mag  = map { abs($_) } $vd->list;
    my $nt   = @mag;
    my @sup  = map { $_ > $thr ? 1 : 0 } @mag;
    my @pk; my $last = -1e9;
    for my $s (0 .. $nt - 1) {
        next unless $sup[$s] && !$sup[($s - 1) % $nt];
        next if $s - $last < $refr;
        my ($best, $bi) = (-1, $s);
        for (my $j = $s; $j < $nt && $sup[$j]; $j++) {
            ($best, $bi) = ($mag[$j], $j) if $mag[$j] > $best;
        }
        push @pk, $bi;
        $last = $bi;
    }
    return @pk;
}

# a quiet trace with bumps of given (start, length, height).
# The noise is white uniform with SD 2 (|v| <= about 3.5): bounded, so on its
# own it never crosses 2.5 x 1.4826 x MAD (see PDL::EEG::SynthNoise). Each trace
# gets its own fixed seed.
my $trace_seed = 280;
sub trace {
    my ($nt, @bumps) = @_;
    my $v = synth_noise($nt, dist => 'uniform', sd => 2, seed => $trace_seed++);
    for my $b3 (@bumps) {
        my ($s, $len, $h) = @$b3;
        my $x = sequence($len);
        my $shape = $h * (1 + 0.5 * sin(3.1 * $x / $len + $s));   # max inside the run
        $v->slice("$s:" . ($s + $len - 1)) .= $shape;
    }
    return $v;
}

my @cases = (
    ['bumps in the middle',      trace(5000, [500, 80, 50], [1500, 120, 40], [3000, 1, 60])],
    ['run at the start',         trace(5000, [0, 90, 50], [2500, 60, 45])],
    ['run at the end',           trace(5000, [1000, 70, 50], [4930, 70, 55])],
    ['runs at both ends (wrap)', trace(5000, [0, 90, 50], [2000, 60, 40], [4950, 50, 55])],
    ['within the refractory',    trace(5000, [1000, 30, 50], [1100, 30, 60], [1400, 30, 45])],
    ['single-sample runs',       trace(5000, map { [300 + 400 * $_, 1, 50 + $_] } 0 .. 10)],
    ['nothing above threshold',  trace(5000)],
);
CORE::srand(28);
for my $i (1 .. 6) {               # random bump trains
    my @b;
    my $s = int(CORE::rand(300));
    while ($s < 19_000) {
        my $len = 1 + int(CORE::rand(150));
        push @b, [$s, $len, 20 + CORE::rand(60)];
        $s += $len + 1 + int(CORE::rand(3000));
    }
    push @cases, ["random train $i", trace(20_000, grep { $_->[0] + $_->[1] <= 20_000 } @b)];
}

for my $c (@cases) {
    my ($name, $v) = @$c;
    next if $name eq 'nothing above threshold';
    my @w = ref_detect($v, refractory_ms => 0);
    ok(@w > 0 && !(grep { abs($v->at($_)) < 10 } @w), "fixture $name: every peak is on a bump, none on the noise");
}

for my $c (@cases) {
    my ($name, $v) = @$c;
    for my $o ([], [refractory_ms => 0], [refractory_ms => 500], [sfreq => 500, thresh_k => 4]) {
        my $tag = "$name" . (@$o ? " (@$o)" : '');
        my $got = detect_blinks($v, @$o);
        my @want = ref_detect($v, @$o);
        is_deeply([$got->list], \@want, "$tag: peaks equal the reference");
        is($got->type->symbol, 'PDL_L', "$tag: long") if @want;
    }
}

is(detect_blinks(trace(5000))->nelem, 0, 'no bump: no peaks (empty, not [0])');

# the reference itself sees the cases it is meant to see
{
    my %c = map { $_->[0] => $_->[1] } @cases;
    my @w = ref_detect($c{'runs at both ends (wrap)'});
    ok(!(grep { $_ < 90 } @w), 'reference: a run at 0 joined to a run at the end does not start a peak');
    ok((grep { $_ >= 4950 } @w), 'reference: the run at the end gives a peak');
    is(scalar(ref_detect($c{'nothing above threshold'})), 0, 'reference: no peaks without a bump');
    is(scalar(ref_detect($c{'within the refractory'}, refractory_ms => 0)), 3, 'reference: 3 bumps, 3 peaks');
    is(scalar(ref_detect($c{'within the refractory'})), 2, 'reference: 200 ms refractory drops one');
}

done_testing;
