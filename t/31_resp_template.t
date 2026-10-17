use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::Heartbeat::Respiration qw(:all);
use PDL::EEG::SynthNoise qw(seed_rng gaussian);

# Synthetic beats: about 1.1 beats/s, breathing rate drifting from 0.27 to 0.31 Hz
# over 300 s. The per-beat series follows a breath that is not a sinusoid:
#   s(phi) = cos(phi) + 0.45 cos(2 phi - 0.8)
# Its second harmonic (0.54-0.62 Hz) lies above half the beat rate, so after
# resampling and band-passing it folds into the breathing band and distorts
# the Hilbert phase; the template fit works on the beat samples themselves.
my $PI = 4 * atan2(1, 1);
seed_rng(31);
my @t = (1.0);
push @t, $t[-1] + 0.88 + 0.06 * sin(0.05 * $t[-1]) while $t[-1] < 300;
my $tb  = pdl(double, \@t);
my $n   = $tb->nelem;
my $phi_true = 2 * $PI * (0.27 * $tb + 0.04 / 600 * $tb ** 2);
my $shape = sub { my $p = shift; cos($p) + 0.45 * cos(2 * $p - 0.8) };
my $y   = $shape->($phi_true) + 0.15 * gaussian($n);
my $deg = 180 / $PI;
my $err = sub {                                     # circular SD of (estimate - truth) [deg], 20 s trimmed at both ends
    my $p  = shift;
    my $in = which(($tb > 20) & ($tb < $tb->at(-1) - 20));
    return circ_sd(wrap_phase($p->index($in) - $phi_true->index($in))) * $deg;
};

# ---- wrap_phase / unwrap_phase / circ_sd ---------------------------------------
my $x = pdl(-7, -$PI, -1, 0, 1, $PI, 7);
my $w = wrap_phase($x);
ok(all(($w >= -$PI) & ($w < $PI)), 'wrap_phase: within -pi .. pi');
cmp_ok((cos($w) - cos($x))->abs->max + (sin($w) - sin($x))->abs->max, '<', 1e-12, 'wrap_phase: same angle');
my $lin = sequence(50) * 0.9;
cmp_ok((unwrap_phase(wrap_phase($lin)) - $lin)->abs->max, '<', 1e-12, 'unwrap_phase: restores a steadily rising phase');
cmp_ok(circ_sd(pdl(0.3, 0.3, 0.3)), '<', 1e-6, 'circ_sd: constant angle -> 0');
my $alt = pdl(map { $_ % 2 ? 0.4 : -0.4 } 0 .. 99);
cmp_ok(abs(circ_sd($alt) - sqrt(-2 * log(cos(0.4)))), '<', 1e-12, 'circ_sd: +-d alternating -> sqrt(-2 ln cos d)');

# ---- fourier_fit / fourier_eval ------------------------------------------------
my $c = fourier_fit($phi_true, $shape->($phi_true) + 2, 2);
is_deeply([$c->dims], [5, 1], 'fourier_fit: (2K+1, nf)');
my $want = pdl(2, 1, 0, 0.45 * cos(0.8), 0.45 * sin(0.8));
cmp_ok(($c->flat - $want)->abs->max, '<', 1e-10, 'fourier_fit: coefficients of a noise-free series');
cmp_ok((fourier_eval($c, $phi_true)->flat - $shape->($phi_true) - 2)->abs->max, '<', 1e-10, 'fourier_eval: reproduces the series');
my $c2 = fourier_fit($phi_true, cat($shape->($phi_true), -3 * $shape->($phi_true)), 2);
is_deeply([$c2->dims], [5, 2], 'fourier_fit: two series -> two columns');
cmp_ok(($c2->slice(':,(1)') + 3 * $c2->slice(':,(0)'))->abs->max, '<', 1e-10, 'fourier_fit: columns fitted independently');

# ---- phase_template ----------------------------------------------------------------
my $ph = resp_phase($tb, $y, band => [0.15, 0.5]);
my $pH = phase_at($ph, $tb);
my $L  = phase_template($tb, $y, $pH);
my ($eH, $eL) = ($err->($pH), $err->($L->{phi}));
note sprintf 'phase error: Hilbert %.1f deg, template %.1f deg', $eH, $eL;
cmp_ok($eL, '<', 8, 'phase_template: error < 8 deg against the true phase');
cmp_ok($eL, '<', 0.8 * $eH, 'phase_template: below the Hilbert phase error');
my $amp = sqrt($L->{c}->slice('1:-1:2,(0)') ** 2 + $L->{c}->slice('2:-1:2,(0)') ** 2);
my $ratio = $amp->at(1) / $amp->at(0);
note sprintf 'second / first harmonic of the template: %.2f (true 0.45)', $ratio;
cmp_ok(abs($ratio - 0.45), '<', 0.15, 'phase_template: template keeps the non-sinusoidal shape');
my $sine = cos($phi_true) + 0.15 * gaussian($n);
my $Ls = phase_template($tb, $sine, phase_at(resp_phase($tb, $sine, band => [0.15, 0.5]), $tb));
my $as = sqrt($Ls->{c}->slice('1:-1:2,(0)') ** 2 + $Ls->{c}->slice('2:-1:2,(0)') ** 2);
cmp_ok($as->at(1) / $as->at(0), '<', 0.1, 'phase_template: no second harmonic invented for a sinusoidal breath');
eval { phase_template($tb, $y->slice('0:9'), $pH) };
like($@, qr/equal length/, 'phase_template: refuses series of different length');

# ---- phase_template_multi ------------------------------------------------------
my $Y = cat($y, 0.5 * $shape->($phi_true - 0.7) + 0.15 * gaussian($n), gaussian($n));   # (n, 3)
my $M = phase_template_multi($tb, $Y, $pH);
is_deeply([$M->{used}->list], [0, 1], 'phase_template_multi: the pure-noise series is gated out');
my $eM = $err->($M->{phi});
note sprintf 'phase error: multi %.1f deg', $eM;
cmp_ok($eM, '<', 12, 'phase_template_multi: error < 12 deg');
is_deeply([$M->{C}->dims], [5, 2], 'phase_template_multi: one template per used series');

# ---- phase_agreement -------------------------------------------------------------
my $q0 = phase_agreement($tb, $phi_true, wrap_phase($phi_true + 1.3));
cmp_ok($q0->abs->max, '<', 1e-6, 'phase_agreement: a constant difference does not count');
my $alt_n = pdl(map { $_ % 2 ? 0.3 : -0.3 } 0 .. $n - 1);
my $q1 = phase_agreement($tb, $phi_true + $alt_n, $phi_true);
my $mid = which(($tb > 70) & ($tb < $tb->at(-1) - 70));
cmp_ok(($q1->index($mid) - sqrt(-2 * log(cos(0.3))))->abs->max, '<', 1e-3, 'phase_agreement: +-d alternating -> sqrt(-2 ln cos d)');
my $q2 = phase_agreement($tb->slice('0:29'), $pH->slice('0:29'), $L->{phi}->slice('0:29'));
is($q2->isfinite->sum, 0, 'phase_agreement: fewer than min beats -> NaN');
my $q3 = phase_agreement($tb->slice('0:29'), $pH->slice('0:29'), $L->{phi}->slice('0:29'), min => 10);
is($q3->isfinite->sum, 30, 'phase_agreement: min can be lowered');
my $qa = phase_agreement($tb, $pH, $L->{phi}) * $deg;
note sprintf 'q_agree median %.1f deg', $qa->median;
cmp_ok($qa->median, '<', 15, 'phase_agreement: H and template phase of a clean series agree within 15 deg');

done_testing();
