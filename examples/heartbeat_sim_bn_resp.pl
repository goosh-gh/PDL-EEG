#!/usr/bin/env perl
# examples/heartbeat_sim_bn_resp.pl -- synthetic BN1-BN2 ECG with a known breathing curve, written as EDF, to measure how well the
# ECG-derived respiration methods (PDL::EEG::Heartbeat::Respiration) recover the phase and the shape of one breath.
# All random numbers come from PDL::EEG::SynthNoise (seed_rng / gaussian / synth_noise), so a seed gives the same file on
# Linux and on the Mac.
#
# Model
#   breathing   breaths of length T_i = period (1 + cv g_i) (clipped 1.8..6 s); phase u = 0..1 within a breath, linear in time.
#               shape c(u), +1 = end of inspiration, -1 = end of expiration:
#                 sine     -cos(2 pi u)
#                 pause    rise 0-0.40, fall 0.40-0.70, end-expiratory pause 0.70-1 (flat at -1)
#                 plateau  rise 0-0.35, end-inspiratory hold 0.35-0.60 (flat at +1), fall 0.60-1
#               depth d(t) = (1 + am_cv h_i) (1 + am_slow sin(2 pi t / am_period)), h_i per breath (interpolated linearly
#               between breath middles, so the curve has no steps)   -> resp(t) = d(t) c(u(t))
#   beats       RR_k = rr + HRV (AR(1) over beats, SD hrv_sd, coefficient hrv_ar) - rsa * resp(t_k - rsa_lag)
#   each beat   sum of Gaussian waves (time re R [ms], size [uV], width SD [ms]) shaped on the average BN1-BN2 beat of a recording at rest:
#                 P (-195, 22, 15) (-165, 45, 10) (-135, 30, 12), Q (-30, -10, 6), R (0, 185, 10), R' (+30, 38, 7), T (+310, 35, 55)
#               R and R' are scaled by (1 + mod_r * resp(t_k - lag_r)), R' further by (1 + mod_rp * resp(t_k - lag_rp)),
#               T by (1 + mod_t * e_T * resp(t_k - lag_t)) and moved by t_shift * e_T * resp(t_k - lag_tl) ms, where e_T = 1 + t_am_cv h'_i is
#               the T wave's own breath-by-breath coupling (independent of the depth; interpolated like the depth), and each wave group
#               (QRS incl. R', T) by an independent per-beat factor (1 + beat_sd_qrs g), (1 + beat_sd_t g)
#   baseline    base_resp * resp(t - base_lag) [uV] + slow drift (AR(1) sampled at 10 Hz, SD drift_sd uV, interpolated)
#   noise       white (noise_ar 0) or AR(1) with coefficient noise_ar, SD noise_sd uV at fs
# Defaults are matched to BN1-BN2 recordings at rest (beat-by-beat features folded on the phase): R SD ~7 uV, T SD ~6 uV, fundamental of the folded R size ~5 % of the
# mean, residual about the 4-harmonic fit 1.1-1.8 % (mod_r 0.05, beat_sd_qrs 0.008, beat_sd_t 0.15, white noise 7 uV)
# Output
#   <out>.edf      BN1 = +x/2, BN2 = -x/2 (so BN1-BN2 = x), RESP = 100 resp(t) [uV units, i.e. 100 = full breath at depth 1]
#   <out>_beats.tsv  per beat: R sample, time [s], true breath phase u (0..1, 0 = start of inspiration), resp(t_R), depth,
#                  R factor, T factor (the true multipliers, before the per-beat noise and after it), R' extra factor, T shift [ms], e_T
# Presets (--preset, options given after it override): base (sine, nothing else), shape (pause), plateau, rsa, am, baseline, all,
#   multi (several features with their own lags, see below), multi_sine (the same with a sinusoidal breath),
#   multi_drift (multi with the lags of RSA, R', T size and T latency wandering slowly, AR(1) ~50 s, SD lag_drift_sd 0.3 s)
# usage: perl heartbeat_sim_bn_resp.pl --preset all --seed 1 --out sim_all_s1
use strict; use warnings; use FindBin; use lib "$FindBin::Bin/../lib"; use PDL; use Getopt::Long;
use PDL::EEG::SynthNoise qw(seed_rng gaussian synth_noise); use PDL::EEG::IO::EDF qw(write_edf);
my $PI = 3.14159265358979;
my %d = (dur => 300, fs => 1000, period => 3.1, cv => 0.15, shape => 'sine', am_cv => 0, am_slow => 0, am_period => 60,
         rr => 0.98, hrv_sd => 0.02, hrv_ar => 0.5, rsa => 0, rsa_lag => 0, mod_r => 0.05, lag_r => 0, mod_t => 0, lag_t => 0, mod_rp => 0, lag_rp => 0, t_shift => 0, lag_tl => 0, t_am_cv => 0, lag_drift_sd => 0,
         beat_sd_qrs => 0.008, noise_ar => 0, beat_sd_t => 0.15, base_resp => 0, base_lag => 0, drift_sd => 0, noise_sd => 7, seed => 1);
my %preset = (base => {}, shape => { shape => 'pause' }, plateau => { shape => 'plateau' }, rsa => { rsa => 0.03, rsa_lag => 1.0 },
              am => { am_cv => 0.2, am_slow => 0.3 }, baseline => { base_resp => 15, drift_sd => 10 },
              all => { shape => 'pause', rsa => 0.03, rsa_lag => 1.0, am_cv => 0.2, am_slow => 0.3, base_resp => 15, drift_sd => 10, mod_t => 0.05, lag_t => 0.3 },
              # several features with their own phase lags (like BN1-BN2 recordings at rest: R' -32 deg, T size inverted, T latency -40 deg) and a T
              # coupling that varies by itself from breath to breath
              multi => { shape => 'pause', rsa => 0.03, rsa_lag => 1.0, am_cv => 0.2, am_slow => 0.3, drift_sd => 10, mod_rp => 0.3, lag_rp => -0.28,
                         mod_t => -0.12, lag_t => -0.13, t_shift => 8, lag_tl => -0.35, t_am_cv => 0.3 });
$preset{multi_sine} = { %{ $preset{multi} }, shape => 'sine' };
$preset{multi_drift} = { %{ $preset{multi} }, lag_drift_sd => 0.3 };   # the feature lags wander slowly and independently
my ($pre) = grep { $ARGV[$_] eq '--preset' } 0 .. $#ARGV; my %o = %d;
if (defined $pre) { my $p = $ARGV[$pre + 1]; die "unknown preset $p\n" unless $preset{$p}; %o = (%o, %{ $preset{$p} }); splice @ARGV, $pre, 2 }
GetOptions(\%o, (map { "$_=" . ($_ eq 'shape' ? 's' : $_ eq 'seed' ? 'i' : 'f') } keys %d), 'out=s') or die "bad option\n";
die "--out required\n" unless $o{out}; die "unknown shape $o{shape}\n" unless $o{shape} =~ /^(sine|pause|plateau)$/;
seed_rng($o{seed});
my ($fs, $nt) = ($o{fs}, int($o{dur} * $o{fs})); my $t = sequence(double, $nt) / $fs;

# --- breathing (computed on the sample grid; extra breaths before 0 so that lagged look-ups exist) ---
my $nbr = int(($o{dur} + 20) / ($o{period} * (1 - 2 * $o{cv})) ) + 10;
my $T = ($o{period} * (1 + $o{cv} * gaussian(double, $nbr)))->clip(1.8, 6);
my $bs = append(pdl(0), $T->cumusumover) - 10;                                              # breath starts [s], first at -10 s
my $h = gaussian(double, $nbr + 1);
sub shape_c { my ($u, $s) = @_; return -cos(2 * $PI * $u) if $s eq 'sine';
    if ($s eq 'pause') { return ($u < 0.4) * -cos($PI * $u / 0.4) + (($u >= 0.4) & ($u < 0.7)) * cos($PI * ($u - 0.4) / 0.3) - ($u >= 0.7) }
    return ($u < 0.35) * -cos($PI * $u / 0.35) + (($u >= 0.35) & ($u < 0.6)) + ($u >= 0.6) * cos($PI * ($u - 0.6) / 0.4) }   # plateau
sub resp_at { my ($tt) = @_; my ($cyc) = interpol($tt, $bs, sequence(double, $bs->nelem)); my $i = $cyc->floor->long->clip(0, $nbr);
    my $u = $cyc - $cyc->floor; my ($hh) = interpol($cyc - 0.5, sequence(double, $nbr + 1), $h);       # per-breath depth, linear between breath middles (no steps)
    my $dep = (1 + $o{am_cv} * $hh) * (1 + $o{am_slow} * sin(2 * $PI * $tt / $o{am_period}));
    return ($dep * shape_c($u, $o{shape}), $u, $dep) }
my ($resp) = resp_at($t);
# slowly wandering lags (s) of RSA, R', T size, T latency: independent AR(1) at 1 Hz (coefficient 0.98, ~50 s), SD lag_drift_sd; drawn only when used
my @LD = map { zeroes(1) } 0 .. 3; my $ldt = sequence(double, int($o{dur}) + 31) - 15;
if ($o{lag_drift_sd} > 0) { @LD = map { synth_noise($ldt->nelem, kind => 'ar1', ar => 0.98, sd => $o{lag_drift_sd}) } 0 .. 3 }
sub ldrift { my ($j, $tt) = @_; return zeroes($tt->nelem) + 0 if $LD[$j]->nelem == 1; my ($v) = interpol($tt, $ldt, $LD[$j]); return $v }

# --- beats ---
my $nmax = int($o{dur} / ($o{rr} * 0.5)) + 10; my $hrv = $o{hrv_sd} > 0 ? synth_noise($nmax, kind => 'ar1', ar => $o{hrv_ar}, sd => $o{hrv_sd}) : zeroes($nmax);
my @tb; my $tk = 0.5 + 0.5 * gaussian(double, 1)->abs->sclr; my $k = 0;
while ($tk < $o{dur} - 0.6) { push @tb, $tk; my ($rl) = resp_at(pdl($tk) - $o{rsa_lag} - ldrift(0, pdl($tk))); $tk += $o{rr} + $hrv->at($k++) - $o{rsa} * $rl->sclr }
my $tb = pdl(@tb); my $nb = $tb->nelem;
my ($rR) = resp_at($tb - $o{lag_r}); my ($rT) = resp_at($tb - $o{lag_t} - ldrift(2, $tb)); my ($r0, $u0, $dep0) = resp_at($tb);
my $eT = ones($nb);                                                                          # T coupling of its own (drawn only when used, so other presets keep their numbers)
my $fR = 1 + $o{mod_r} * $rR; my $nq = gaussian(double, $nb); my $gR = $fR * (1 + $o{beat_sd_qrs} * $nq);
my ($rP) = resp_at($tb - $o{lag_rp} - ldrift(1, $tb)); my $fRP = 1 + $o{mod_rp} * $rP; my $gRP = $gR * $fRP;
my $nT = gaussian(double, $nb);
if ($o{t_am_cv} > 0) { my $h2 = gaussian(double, $nbr + 1); my ($cyc) = interpol($tb, $bs, sequence(double, $bs->nelem)); my ($hh) = interpol($cyc - 0.5, sequence(double, $nbr + 1), $h2); $eT = 1 + $o{t_am_cv} * $hh }
my $fT = 1 + $o{mod_t} * $eT * $rT; my $gT = $fT * (1 + $o{beat_sd_t} * $nT);
my ($rL) = resp_at($tb - $o{lag_tl} - ldrift(3, $tb)); my $sT = $o{t_shift} * $eT * $rL;                       # T latency change [ms]
my @wave = ([-195, 22, 15, 'P'], [-165, 45, 10, 'P'], [-135, 30, 12, 'P'], [-30, -10, 6, 'QRS'], [0, 185, 10, 'QRS'], [30, 38, 7, 'RP'], [310, 35, 55, 'T']);
my $lag = sequence(double, int(1.0 * $fs)) / $fs * 1000 - 400;                               # -400..+600 ms window
my $x = zeroes(double, $nt);
for my $kb (0 .. $nb - 1) { my $s0 = int($tb->at($kb) * $fs + 0.5); my $frac = $tb->at($kb) * 1000 - $s0 * 1000 / $fs;   # sub-sample R time
    my $i0 = $s0 - 400 * $fs / 1000; my $i1 = $i0 + $lag->nelem - 1; my ($j0, $j1) = ($i0 < 0 ? -$i0 : 0, $i1 >= $nt ? $lag->nelem - 1 - ($i1 - $nt + 1) : $lag->nelem - 1);
    next if $j1 < $j0; my $l = $lag->slice("$j0:$j1") - $frac; my $w = zeroes(double, $l->nelem);
    for my $wv (@wave) { my ($c, $size, $sd, $g) = @$wv; my $f = $g eq 'QRS' ? $gR->at($kb) : $g eq 'RP' ? $gRP->at($kb) : $g eq 'T' ? $gT->at($kb) : 1;
        $c += $sT->at($kb) if $g eq 'T'; $w += $f * $size * exp(-0.5 * (($l - $c) / $sd) ** 2) }
    $x->slice(($i0 + $j0) . ':' . ($i0 + $j1)) += $w }
if ($o{base_resp}) { my ($rb) = resp_at($t - $o{base_lag}); $x += $o{base_resp} * $rb }
if ($o{drift_sd}) { my $n10 = int($o{dur} * 10) + 2; my $dr = synth_noise($n10, kind => 'ar1', ar => 0.995, sd => $o{drift_sd}); my ($dv) = interpol($t, sequence(double, $n10) / 10, $dr); $x += $dv }
$x += ($o{noise_ar} > 0 ? synth_noise($nt, kind => 'ar1', ar => $o{noise_ar}, sd => $o{noise_sd}) : synth_noise($nt, sd => $o{noise_sd})) if $o{noise_sd} > 0;
my $data = cat($x / 2, -$x / 2, 100 * $resp)->transpose;                                     # (3, nt)
write_edf({ data => $data, fs => $fs, labels => [qw(BN1 BN2 RESP)] }, "$o{out}.edf");
open my $fh, '>', "$o{out}_beats.tsv" or die; print $fh "# heartbeat_sim_bn_resp.pl ", join(' ', map { "$_=$o{$_}" } sort grep { $_ ne 'out' } keys %o), "\n";
print $fh join("\t", qw(sample t_s u resp depth fR_true fR_beat fT_true fT_beat fRP_true Tshift_ms eT dlag_rsa dlag_rp dlag_t dlag_tl)), "\n";
my @DL = map { ldrift($_, $tb) } 0 .. 3;
printf $fh "%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.3f\t%.4f\t%.3f\t%.3f\t%.3f\t%.3f\n", int($tb->at($_) * $fs + 0.5), $tb->at($_), $u0->at($_), $r0->at($_), $dep0->at($_), $fR->at($_), $gR->at($_), $fT->at($_), $gT->at($_), $fRP->at($_), $sT->at($_), $eT->at($_), $DL[0]->at($_), $DL[1]->at($_), $DL[2]->at($_), $DL[3]->at($_) for 0 .. $nb - 1;
close $fh; printf "%s: %d beats, %.0f s, breaths %.2f +- %.2f s, shape %s -> %s.edf, %s_beats.tsv\n", $o{out}, $nb, $o{dur}, $T->avg, sqrt((($T - $T->avg) ** 2)->avg), $o{shape}, $o{out}, $o{out};
