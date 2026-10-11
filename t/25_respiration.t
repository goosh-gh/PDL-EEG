use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::Heartbeat::Respiration qw(:all);

# Deterministic synthetic recording (no RNG). Breathing phase phi(t) with a rate
# drifting from 0.27 to 0.35 Hz over 200 s. RR interval ending at a beat:
# 0.90 + 0.05 cos(phi) s. R amplitude of two leads: +60 and -120 uV, both times
# 1 + 0.06 cos(phi - 0.5). Each lead also carries a constant offset and a
# 0.05 Hz drift.
my $PI = 4 * atan2(1, 1);
my ($fs, $dur) = (250, 200);
my $phi = sub { my $t = shift; 2 * $PI * (0.27 * $t + 0.0002 * $t * $t) };
my @r = (int(2.0 * $fs));
while (1) {
    my $t0 = $r[-1] / $fs;
    my $t1 = $t0 + 0.90;
    $t1 = $t0 + 0.90 + 0.05 * cos($phi->($t1)) for 1 .. 6;       # RR depends on the phase at the beat it ends at
    last if $t1 > $dur - 2;
    push @r, int($t1 * $fs + 0.5);
}
my $r  = pdl(long, \@r);
my $tb = $r->double / $fs;                                      # beat times on the sample grid
my $nb = $r->nelem;
my $nt = $dur * $fs;
my $t  = sequence(double, $nt) / $fs;
my $gain = 1 + 0.06 * cos($phi->($tb) - 0.5);                   # (nbeat)
my @lead = (7 + 5 * sin(2 * $PI * 0.05 * $t), -3 - 4 * sin(2 * $PI * 0.05 * $t));
for my $k (0 .. $nb - 1) {
    my $g = exp(-0.5 * (($t - $tb->at($k)) / 0.010) ** 2) * $gain->at($k);
    $lead[0] += 60 * $g;
    $lead[1] -= 120 * $g;
}
my $X = cat(@lead)->xchg(0, 1);                                 # (2, nt)
my $deg = sub { $_[0] * 180 / $PI };
my $wrap = sub { my $d = shift; atan2(sin($d), cos($d)) };

# ---- rr_series ---------------------------------------------------------------
my ($t_rr, $rr) = rr_series($r, $fs);
is($rr->nelem, $nb - 1, 'rr_series: nbeat-1 intervals');
is_deeply([map { sprintf '%.6f', $_ } $t_rr->list], [map { sprintf '%.6f', $r[$_] / $fs } 1 .. $nb - 1], 'rr_series: times of beats 1..n-1');
is_deeply([map { sprintf '%.6f', $_ } $rr->list], [map { sprintf '%.6f', ($r[$_] - $r[$_ - 1]) / $fs } 1 .. $nb - 1], 'rr_series: interval ending at each beat');
eval { rr_series(pdl(long, [5]), $fs) };
like($@, qr/at least 2 beats/, 'rr_series: refuses a single beat');

# ---- r_amplitude ---------------------------------------------------------------
my $a = r_amplitude($X, $r, $fs);
is_deeply([$a->dims], [2, $nb], 'r_amplitude: (nlead, nbeat)');
cmp_ok(($a->slice("(0)") - 60 * $gain)->abs->max, '<', 0.2, 'r_amplitude: lead 0 = 60 x gain, offset and drift removed');
cmp_ok(($a->slice("(1)") + 120 * $gain)->abs->max, '<', 0.2, 'r_amplitude: lead 1 = -120 x gain');
my $a1 = r_amplitude($X->slice("(0)"), $r, $fs);
is_deeply([$a1->dims], [$nb], 'r_amplitude: (nt) in, (nbeat) out');
cmp_ok(($a1 - $a->slice("(0)"))->abs->max, '<', 1e-12, 'r_amplitude: 1-D result = row 0');
eval { r_amplitude($X, $r, $fs, pre => [-0.02, 0.01]) };
like($@, qr/pre must lie before R/, 'r_amplitude: refuses a window that reaches R');

# ---- resp_phase / phase_at ---------------------------------------------------
my $in  = which(($t_rr > 20) & ($t_rr < $dur - 20));           # beats away from both ends
my $tin = $t_rr->index($in);
my $p_rr = resp_phase($t_rr, $rr);
is($p_rr->{fs}, 4, 'resp_phase: 4 Hz grid by default');
cmp_ok(abs($p_rr->{t}->at(1) - $p_rr->{t}->at(0) - 0.25), '<', 1e-12, 'resp_phase: grid step 1/fs');
cmp_ok($p_rr->{peak_hz}, '>', 0.26, 'resp_phase RR: spectral peak not below the breathing rates');
cmp_ok($p_rr->{peak_hz}, '<', 0.36, 'resp_phase RR: spectral peak not above the breathing rates');
cmp_ok($p_rr->{band_frac}, '>', 0.8, 'resp_phase RR: most power inside the band');
{
    my $d = $wrap->(phase_at($p_rr, $tin) - $phi->($tin));
    my ($plv, $md) = phase_locking(phase_at($p_rr, $tin), $phi->($tin));
    cmp_ok($plv, '>', 0.98, "RR phase follows the true phase (PLV $plv)");
    cmp_ok(abs($deg->($md)), '<', 3, 'RR phase: mean error below 3 deg');
    cmp_ok($deg->($d->abs->max), '<', 10, 'RR phase: every interior beat within 10 deg');
}
my $mag  = sqrt(($a ** 2)->sumover);                             # (nbeat)
my $p_am = resp_phase($tb, $mag);
{
    my ($plv, $md) = phase_locking(phase_at($p_am, $tin), $phi->($tin) - 0.5);
    cmp_ok($plv, '>', 0.98, "R-size phase follows the true phase minus 0.5 rad (PLV $plv)");
    cmp_ok(abs($deg->($md)), '<', 3, 'R-size phase: mean error below 3 deg');
    my ($plv2, $md2) = phase_locking(phase_at($p_am, $tin), phase_at($p_rr, $tin));
    cmp_ok($plv2, '>', 0.98, 'the two estimates are locked to each other');
    cmp_ok(abs($md2 + 0.5), '<', 0.2, 'R-size phase lags RR phase by the 0.5 rad built in');
}
{
    # negative control: an RR series with no breathing in it (deterministic, k^2 argument)
    my $k   = sequence(double, $nb - 1);
    my $rr0 = 0.90 + 0.05 * sin(0.7 * $k * $k);
    my $p0  = resp_phase($t_rr, $rr0);
    my ($plv) = phase_locking(phase_at($p0, $tin), $phi->($tin));
    cmp_ok($plv, '<', 0.3, "control without breathing: no locking to the true phase (PLV $plv)");
    cmp_ok($p0->{band_frac}, '<', $p_rr->{band_frac}, 'control without breathing: smaller share of power in the band');
}
{
    my $g = $p_rr->{t}->slice("40:49");
    cmp_ok(($wrap->(phase_at($p_rr, $g) - atan2($p_rr->{im}->slice("40:49"), $p_rr->{re}->slice("40:49"))))->abs->max, '<', 1e-12,
           'phase_at: at grid times = angle of the analytic signal');
    my $mx = $p_rr->{series}->slice("100:-100")->maximum_ind->sclr + 100;
    cmp_ok(abs(phase_at($p_rr, $p_rr->{t}->at($mx))->sclr), '<', 0.5, 'phase_at: near 0 at a maximum of the band-passed series');
    my $mn = $p_rr->{series}->slice("100:-100")->minimum_ind->sclr + 100;
    cmp_ok(abs(phase_at($p_rr, $p_rr->{t}->at($mn))->sclr), '>', $PI - 0.5, 'phase_at: near +-pi at a minimum');
    my $T2 = $tin->slice("0:5")->dummy(0, 2);                  # (2, 6)
    is_deeply([phase_at($p_rr, $T2)->dims], [2, 6], 'phase_at: keeps the shape of the times');
    eval { phase_at($p_rr, pdl(0.5)) };
    like($@, qr/times must lie within/, 'phase_at: refuses a time before the grid');
    eval { phase_at($p_rr, pdl($dur + 5)) };
    like($@, qr/times must lie within/, 'phase_at: refuses a time after the grid');
    eval { resp_phase($t_rr, $rr, band => [0.2, 2.5]) };
    like($@, qr/band must satisfy/, 'resp_phase: refuses a band beyond fs/2');
    eval { resp_phase($t_rr->slice("0:9"), $rr->slice("0:9")) };
    like($@, qr/fewer than 3 cycles/, 'resp_phase: refuses a series shorter than 3 cycles of the lower edge');
}

# ---- phase_locking -------------------------------------------------------------
{
    my $p = sequence(double, 360) * 2 * $PI / 360 - $PI;
    my ($plv, $d) = phase_locking($p, $p - 0.4);
    cmp_ok(abs($plv - 1), '<', 1e-12, 'phase_locking: constant difference gives 1');
    cmp_ok(abs($d - 0.4), '<', 1e-12, 'phase_locking: and that difference');
    my ($plv0) = phase_locking($p, zeroes(360));
    cmp_ok($plv0, '<', 1e-12, 'phase_locking: uniform differences give 0');
    my ($plvh) = phase_locking(pdl(0, 0, $PI / 2, $PI / 2), zeroes(4));
    cmp_ok(abs($plvh - sqrt(0.5)), '<', 1e-12, 'phase_locking: two directions 90 deg apart give sqrt(1/2)');
}

done_testing;
