#!/usr/bin/env perl
# examples/heartbeat_resp_vs_belt.pl -- how right is the ECG-only respiratory phase (the methods used on BN1-BN2) when checked against breathing belts?
# For one Challenge 2018 PSG record: the ECG lead is treated exactly as BN1-BN2 is (0.5-45 Hz, detect_r, BeatFeatures features, same beat
# selection), the respiratory phase is estimated from the ECG alone with each method below, and compared with the phase of the belt.
# Nothing from the belt goes into the ECG estimates; the belt is used only to score them.
#
# Belt ("true") phase: the belt (--belt ABD, up = inspiration) 0.05-2 Hz; inspiration onsets = belt minima (near the troughs of the Hilbert
#   phase of its 0.1-0.7 Hz band, refined to the minimum within +-0.6 s); linear in time within each breath, 0 = inspiration onset, 360 =
#   next onset; end of inspiration = belt maximum.  Breaths 1.5-8 s with amplitude > 0.3 x median of the 21 around are kept.
# ECG phases (0 = largest R size): H  Hilbert phase of the 0.15-0.5 Hz band of the R size (resp_phase)
#                                  L  phase_template on the R size (window w s, template K 2, i rounds) started from H
#                                  M  phase_template_multi on R size, R' stretch, T latency (w 5, i 3)
# Scores (per method):
#   offset      circular mean of (ECG phase - belt phase): where the largest R size sits in the breath
#   SD          circular SD of the difference after the offset
#   |err|>45/90 share of beats whose error (after the offset) exceeds 45 / 90 deg
#   offset SD   circular SD of the offset across 10-min windows (is one offset enough for a whole recording?)
#   warp        mean error in 12 bins of the belt phase, max - min (a phase that runs unevenly through the breath)
#   R2 belt     how much of the beat-to-beat belt value (sampled at R, 0.05-2 Hz) the ECG phase explains (cos / sin regression)
# Writes PREFIX_phase_eval.tsv (per beat: time, belt phase, ECG phases) and PREFIX_phase_eval.txt (the table).
# usage: perl heartbeat_resp_vs_belt.pl /path/RECORD [--hours 1,3] [--out PREFIX] [--belt ABD]
use strict; use warnings; use FindBin; use lib "$FindBin::Bin/../lib"; use PDL; use Getopt::Long;
use PDL::EEG::IO::WFDB qw(read_wfdb); use PDL::EEG::Heartbeat::BeatFeatures qw(beat_features);
use PDL::EEG::Filter qw(zp_butter analytic_signal); use PDL::EEG::Heartbeat qw(detect_r); use PDL::EEG::Heartbeat::Respiration qw(resp_phase phase_at phase_template phase_template_multi fourier_fit wrap_phase unwrap_phase circ_sd);
my %o = (belt => 'ABD', hours => ''); GetOptions(\%o, 'hours=s', 'out=s', 'belt=s') or die "bad option\n";
my $stem = shift @ARGV or die "usage: $0 /path/RECORD [options]\n"; $stem =~ s/\.(mat|hea)$//; (my $nm = $stem) =~ s{.*/}{}; my $P = $o{out} // $nm;
my $PI = 3.14159265358979; my $D = 180 / $PI;
open my $TXT, '>', "${P}_phase_eval.txt" or die; sub say_ { my $s = sprintf(shift, @_); print $s; print $TXT $s }
my $R = read_wfdb($stem); my $fs = $R->{fs}; my ($a0, $a1) = (0, $R->{nt});
if ($o{hours}) { my ($h0, $h1) = split /,/, $o{hours}; $a0 = int($h0 * 3600 * $fs); $a1 = int($h1 * 3600 * $fs); $a1 = $R->{nt} if $a1 > $R->{nt} }
my $nt = $a1 - $a0; my $chan = sub { $R->{chan}->($_[0], from => $a0, to => $a1) };

# ---- belt phase (scoring only) ----
my $belt = zp_butter($chan->($o{belt}), fs => $fs, low => 0.05, high => 2, order => 2);
my ($bre, $bim) = analytic_signal(zp_butter($chan->($o{belt}), fs => $fs, low => 0.1, high => 0.7, order => 2));
my $pu = unwrap_phase(atan2($bim, $bre)->slice('0:-1:10')->copy); my $tg = sequence(long, $pu->nelem) * 10; undef $bre; undef $bim;
my $kk = (($pu - $PI) / (2 * $PI))->floor; my $cr = which($kk->slice('1:-1') > $kk->slice('0:-2')); my $w6 = int(0.6 * $fs); my @onset;
for my $c ($cr->list) { my $s = $tg->at($c); my ($lo, $hi) = ($s - $w6 < 0 ? 0 : $s - $w6, $s + $w6 >= $nt ? $nt - 1 : $s + $w6); push @onset, $lo + $belt->slice("$lo:$hi")->minimum_ind->sclr }
my $bb = pdl(long, \@onset)->uniq; my $nbr = $bb->nelem - 1; my (@pk, @amp);
for my $i (0 .. $nbr - 1) { my ($s, $e) = ($bb->at($i), $bb->at($i + 1)); my $seg = $belt->slice("$s:$e"); push @pk, $s + $seg->maximum_ind->sclr; push @amp, ($seg->max - ($belt->at($s) + $belt->at($e)) / 2)->sclr }
my $amp = pdl(\@amp); my $len = ($bb->slice('1:-1') - $bb->slice('0:-2'))->double / $fs; my $tif = (pdl(long, \@pk) - $bb->slice('0:-2'))->double / $fs / $len;
my $medamp = pdl(map { my ($l, $h) = ($_ - 10 < 0 ? 0 : $_ - 10, $_ + 10 >= $nbr ? $nbr - 1 : $_ + 10); $amp->slice("$l:$h")->median } 0 .. $nbr - 1);
my $okbr = ($len >= 1.5) & ($len <= 8) & ($amp > 0.3 * $medamp);

# ---- ECG, exactly as for BN1-BN2 ----
my $ecg = zp_butter($chan->('ECG') * ($R->{units}[ { map { ($R->{labels}[$_] => $_) } 0 .. $#{ $R->{labels} } }->{ECG} ] =~ /^mV$/i ? 1000 : 1), fs => $fs, low => 0.5, high => 45, order => 4);
my ($r, $sg) = detect_r($ecg, $fs); my $B = beat_features($ecg * $sg, $fs, r => $r); $r = $B->{r}; my $tb = $B->{tb}; my $nb = $tb->nelem; my %F = %{ $B->{F} };
my $medw = sub { my ($v, $hw) = @_; my $n = $v->nelem; pdl(map { my ($l, $h) = ($_ - $hw < 0 ? 0 : $_ - $hw, $_ + $hw >= $n ? $n - 1 : $_ + $hw); $v->slice("$l:$h")->median } 0 .. $n - 1) };
my $mrr = $medw->($F{RRp}, 5); my $rmed = $F{R}->median; my $rmad = 1.4826 * abs($F{R} - $rmed)->median;
my $okb = (abs($F{RRp} - $mrr) < 0.25 * $mrr) & (abs($F{RRn} - $mrr) < 0.25 * $mrr) & (abs($F{R} - $rmed) < 6 * $rmad); my $use = zeroes(long, $nb); $use->index($B->{use}) .= 1; $okb &= $use;
my $ue = which($okb);                                                                   # beats the ECG methods use (ECG criteria only)
my $tk = $tb->index($ue); my $yR = $F{R}->index($ue);
my $hp = resp_phase($tk, $yR, band => [0.15, 0.5]); my $hph = phase_at($hp, $tk);
my %E; my @M; my $add = sub { my ($k, $p) = @_; my $c1 = fourier_fit($p, $yR, 1); $E{$k} = wrap_phase($p - atan2($c1->at(2, 0), $c1->at(1, 0))); push @M, $k };
$add->('H 0.15-0.5', $hph);
for my $cfg ([5, 3], [5, 6], [4, 3], [6, 3]) { my ($w, $it) = @$cfg; $add->("L w$w i$it", phase_template($tk, $yR, $E{'H 0.15-0.5'}, win => $w, K => 2, iter => $it)->{phi}) }
$add->('M R,Qwid,Tlat', phase_template_multi($tk, cat(map { $F{$_}->index($ue) } qw(R Qwid Tlat)), $E{'H 0.15-0.5'}, win => 5, K => 2, iter => 3, gate => 0.1)->{phi});

# ---- scoring on the beats that also lie in a kept breath ----
my ($cyc) = interpol($r->index($ue)->double, $bb->double, sequence($bb->nelem)); my $ci = $cyc->floor->long; my $cic = $ci->clip(0, $nbr - 1);
my $sc = which(($ci >= 0) & ($ci < $nbr) & $okbr->index($cic)); my $frac = ($cyc - $cyc->floor)->index($sc); my $pt = wrap_phase(2 * $PI * $frac); my $ti = $tif->index($cic->index($sc));
my $insp = $frac < $ti; my $timed = $ti->median->sclr; my $ts = $tk->index($sc);
my $bv = $belt->index($r->index($ue)->index($sc)); my $bvc = $bv - $bv->avg;
say_("%s %.2f-%.2f h, belt %s: %d breaths (%d kept, %.1f /min, inspiration share %.2f); ECG: %d beats, %d used, %d scored\n", $nm, $a0 / $fs / 3600, $a1 / $fs / 3600, $o{belt},
    $nbr, $okbr->sum, 60 / $len->where($okbr)->avg, $timed, $nb, $ue->nelem, $sc->nelem);
say_("%-15s %7s %6s | %6s %6s | %8s | %6s | %7s\n", 'ECG phase', 'offset', 'SD', '>45', '>90', 'off. SD', 'warp', 'R2 belt');
open my $T, '>', "${P}_phase_eval.tsv" or die; print $T join("\t", 't_s', 'belt', 'insp', @M), "\n";
for my $m (@M) { my $pe = $E{$m}->index($sc); my $d = wrap_phase($pe - $pt); my $off = atan2(sin($d)->avg, cos($d)->avg); my $e = wrap_phase($d - $off);
    my @wo; for (my $s = $ts->at(0); $s + 600 <= $ts->at(-1) + 1; $s += 600) { my $w = which(($ts >= $s) & ($ts < $s + 600)); next if $w->nelem < 100; my $dw = $d->index($w); push @wo, atan2(sin($dw)->avg, cos($dw)->avg)->sclr }
    my $osd = @wo > 1 ? circ_sd(pdl(@wo)) * $D : 0;
    my $bin = (($pt + $PI) / (2 * $PI) * 12)->floor->clip(0, 11); my @mb = map { my $w = which($bin == $_); $w->nelem ? atan2(sin($e->index($w))->avg, cos($e->index($w))->avg)->sclr * $D : 0 } 0 .. 11;
    my $X = cat(ones($pe->nelem), cos($pe), sin($pe))->transpose; my $bc = inv($X->transpose x $X) x ($X->transpose x $bvc->dummy(0)); my $f = ($X x $bc)->flat;
    say_("%-15s %+6.0f° %5.1f° | %5.1f%% %5.1f%% | %6.1f° | %5.0f° | %7.3f\n", $m, $off * $D, circ_sd($e) * $D, 100 * (abs($e) > $PI / 4)->avg, 100 * (abs($e) > $PI / 2)->avg,
        $osd, pdl(@mb)->max - pdl(@mb)->min, (sum($f ** 2) / sum($bvc ** 2))->sclr) }
say_("offset: ECG phase - belt phase; SD, >45/>90: error after the offset; off. SD: across 10-min windows;\nwarp: range of the mean error over 12 belt-phase bins; R2 belt: belt value explained by the ECG phase\n");
for my $i (0 .. $sc->nelem - 1) { print $T join("\t", sprintf('%.3f', $ts->at($i) + $a0 / $fs), sprintf('%.4f', $pt->at($i)), $insp->at($i) ? 1 : 0, map { sprintf '%.4f', $E{$_}->at($sc->at($i)) } @M), "\n" }
close $T; close $TXT;
