use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::Heartbeat::BeatFeatures qw(beat_features);

# Synthetic lead (no noise): a beat every 0.8 s; QRS a Gaussian (SD 10 ms) at R,
# T a Gaussian (SD 40 ms, 0.3 of the QRS height) at +300 ms. Over a cycle of 5
# beats (phase p = 2 pi k / 5) the beat changes by known amounts:
#   QRS height x (1 + 0.05 cos p), QRS width x (1 + 0.04 sin p),
#   T shift 4 cos p ms, T height x (1 + 0.10 sin p).
# The +-10 s running mean spans whole cycles, so the beat minus its slow mean is
# the change alone.
my $PI = 4 * atan2(1, 1);
for my $fs (500, 1000) {
    my $dur = 120;
    my $nt  = $dur * $fs;
    my $t   = sequence(double, $nt) / $fs;
    my $x   = zeroes(double, $nt);
    my (@r, @p);
    for (my $k = 0, my $tr = 1.0; $tr < $dur - 1; $k++, $tr += 0.8) {
        my $p   = 2 * $PI * $k / 5;
        my $rs  = int($tr * $fs + 0.5);
        my $trs = $rs / $fs;
        my $q   = 0.010 * (1 + 0.04 * sin($p));
        my $tt  = $trs + 0.300 + 0.004 * cos($p);
        my $lo  = $rs - int(0.4 * $fs);
        my $hi  = $rs + int(0.7 * $fs);
        my $s   = $t->slice("$lo:$hi");
        $x->slice("$lo:$hi") += (1 + 0.05 * cos($p)) * exp(-0.5 * (($s - $trs) / $q) ** 2)
                              + 0.3 * (1 + 0.10 * sin($p)) * exp(-0.5 * (($s - $tt) / 0.040) ** 2);
        push @r, $rs;
        push @p, $p;
    }
    my $B = beat_features($x, $fs, r => pdl(long, \@r));
    my $u = $B->{use};
    my $p = pdl(\@p)->index($u);
    my $in = which(($B->{tb}->index($u) > 12) & ($B->{tb}->index($u) < $dur - 12));   # whole cycles in the running mean
    my $fit = sub {                                     # slope of a feature on its known modulation, and the residual RMS
        my ($name, $mod) = @_;
        my $y = $B->{F}{$name}->index($u)->index($in);
        my $m = $mod->index($in);
        my $g = sum($y * $m) / sum($m * $m);
        return ($g, sqrt((($y - $g * $m) ** 2)->avg) / sqrt((($g * $m) ** 2)->avg));
    };
    is_deeply($B->{names}, [qw(R Qsz Qwid Qlat Tsz Tlat Twid RRp RRn)], "fs $fs: feature names");
    is($B->{tb}->nelem, scalar @r, "fs $fs: every beat kept");
    for my $c (['R', 5 * cos($p), 0.02], ['Qsz', 5 * cos($p), 0.05], ['Qwid', 4 * sin($p), 0.05],
               ['Tlat', 4 * cos($p), 0.05], ['Tsz', 10 * sin($p), 0.05]) {
        my ($name, $mod, $tol) = @$c;
        my ($g, $rel) = $fit->($name, $mod);
        note sprintf 'fs %d %-4s gain %.3f, residual %.3f', $fs, $name, $g, $rel;
        cmp_ok(abs($g - 1), '<', $tol, "fs $fs: $name follows its modulation (gain within $tol of 1)");
        cmp_ok($rel, '<', 0.10, "fs $fs: $name residual < 10 %");
    }
    for my $c (['Qlat', 'Qsz'], ['Twid', 'Tsz']) {      # quantities that were not changed stay small
        my ($none, $ref) = @$c;
        my $ratio = $B->{F}{$none}->index($u)->index($in)->abs->max / $B->{F}{$ref}->index($u)->index($in)->abs->max;
        cmp_ok($ratio, '<', 0.2, "fs $fs: $none small (no such change made)");
    }
    cmp_ok(($B->{F}{RRp}->index($u) - 800)->abs->max, '<=', 1000 / $fs + 1e-9, "fs $fs: RRp = 800 ms within a sample");
}
eval { beat_features(zeroes(1000), 500) };
like($@, qr/r => R samples required/, 'beat_features: refuses a call without R samples');

done_testing();
