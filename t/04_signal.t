use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::Signal qw(detect_square_pulses);

# ---------------------------------------------------------------------------
# Synthetic recording:
#   - 20 EEG channels: modest Gaussian noise (ranges ~ 100-500 uV)
#   - 4 TTL trigger channels: baseline 0, occasional pulses to +3188 uV (rail)
#   - 2 constant-rail markers: stuck at -3133 uV (range ~ 0)
# The detector must pick exactly the 4 trigger channels.
# ---------------------------------------------------------------------------

my $n_samp = 2000;
srandom(42);

my @rows;
# 20 EEG channels, sigma 15..60
for my $k (0 .. 19) {
    my $sigma = 15 + ($k % 10) * 5;
    push @rows, grandom($n_samp) * $sigma;
}
# 4 trigger channels: mostly 0, pulse windows to the rail
my @trig_pos_expected;
for my $k (0 .. 3) {
    my $ch = zeroes(double, $n_samp);
    my $s  = 200 + $k * 300;                     # distinct pulse windows
    $ch->slice("$s:" . ($s + 60)) .= 3188;       # ~60 ms pulse
    $ch->slice(($s + 500) . ":" . ($s + 560)) .= 3188;
    push @rows, $ch;
}
# 2 constant-rail markers
push @rows, (zeroes($n_samp) - 3133) for 1 .. 2;

my $data = cat(@rows)->transpose;                # [n_ch, n_samp]
is($data->dim(0), 26, 'assembled 26 channels');

# trigger channels are positions 20,21,22,23
my %want = map { $_ => 1 } (20, 21, 22, 23);

my $cands = detect_square_pulses($data, fs => 1000, skip_sec => 0, n => 4);
is(scalar @$cands, 4, 'exactly 4 candidates returned');

my %got = map { $_->{pos} => 1 } @$cands;
is_deeply(\%got, \%want, 'the 4 trigger channels are detected (not EEG, not markers)');

# markers must be rejected even without n (range guard)
my $all = detect_square_pulses($data, fs => 1000, skip_sec => 0);
my %all = map { $_->{pos} => 1 } @$all;
ok(!$all{24} && !$all{25}, 'constant-rail markers rejected');
ok(!$all{0} && !$all{5},   'ordinary EEG channels rejected');

done_testing();
