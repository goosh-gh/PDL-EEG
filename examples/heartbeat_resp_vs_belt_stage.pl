#!/usr/bin/env perl
# examples/heartbeat_resp_vs_belt_stage.pl -- as heartbeat_resp_vs_belt.pl (the ECG-only respiratory phase, made exactly as for BN1-BN2, scored against a breathing
# belt), but scored separately for each sleep stage (wake W, N1, N2, N3, REM R) of a Challenge 2018 training record, inside stretches
# where the heart's position is unlikely to have changed.
#
# Stages and events: <record>.arousal (WFDB annotations; stage changes are aux texts W N1 N2 N3 R, events "(name" ... "name)").  Beats
#   inside any event (arousals, apneas, hypopneas ...) or within 10 s after its end are not scored.
# Changes of the heart's position (turning in bed, sitting up): from the ECG alone, so that the same check can be run on BN1-BN2.
#   The average QRS (-60..+70 ms around R) of each 30 s window is compared with that of the previous window: d = RMS(difference) /
#   RMS(previous).  A change is flagged where d > max(--dmin 0.15, median + 6 x robust SD of d).  Flagged windows end a stretch.
# Stretches: runs of one stage, not crossing a flagged change, at least --minlen (120) s long.  The ECG phase is compared with the belt
#   phase in each stretch after removing that stretch's own offset (so a change of position, which moves the offset, does not count as
#   error); errors are pooled over the stretches of a stage.
# Methods (0 = largest R size): H  0.15-0.5 Hz band of the R size;  L  phase_template (R size, window 5 s, K 2, 3 rounds);
#   M  phase_template_multi (R size, R' stretch, T latency; window 5 s, K 2, 3 rounds, gate 0.1).  All three are computed on the
#   whole (or --hours) record from the ECG alone; the belt and the stages are used only for scoring.
# Scores per stage and method: stretches, beats, SD (circular SD of the error), >90 (share of beats with |error| > 90 deg), warp (range of
#   the mean error over 12 belt-phase bins), off. spread (circular SD of the stretch offsets: how much the ECG-belt relation moves
#   between stretches of the same stage).
# Output: PREFIX_stage_eval.txt (tables, the flagged position changes, the stretches), PREFIX_stage_eval.tsv (per scored beat).
# usage: perl heartbeat_resp_vs_belt_stage.pl /path/tr03-0005 [--hours 0,9] [--out PREFIX] [--belt ABD] [--dmin 0.15] [--minlen 120]
use strict; use warnings; use FindBin; use lib "$FindBin::Bin/../lib"; use PDL; use Getopt::Long;
use PDL::EEG::IO::WFDB qw(read_wfdb read_wfdb_annotations wfdb_stages_events); use PDL::EEG::Heartbeat::BeatFeatures qw(beat_features);
use PDL::EEG::Filter qw(zp_butter analytic_signal); use PDL::EEG::Heartbeat qw(detect_r beat_epochs); use PDL::EEG::Heartbeat::Respiration qw(resp_phase phase_at phase_template phase_template_multi fourier_fit wrap_phase unwrap_phase circ_sd);
my %o = (belt => 'ABD', hours => '', dmin => 0.15, minlen => 120); GetOptions(\%o, 'hours=s', 'out=s', 'belt=s', 'dmin=f', 'minlen=f', 'arousal=s') or die "bad option\n";
my $stem = shift @ARGV or die "usage: $0 /path/RECORD [options]\n"; $stem =~ s/\.(mat|hea|arousal)$//; (my $nm = $stem) =~ s{.*/}{}; my $P = $o{out} // $nm;
my $afile = $o{arousal} // "$stem.arousal"; die "no $afile (training records only)\n" unless -r $afile;
my $PI = 3.14159265358979; my $D = 180 / $PI;
open my $TXT, '>', "${P}_stage_eval.txt" or die; sub say_ { my $s = sprintf(shift, @_); print $s; print $TXT $s }
my $R = read_wfdb($stem); my $fs = $R->{fs}; my ($a0, $a1) = (0, $R->{nt});
if ($o{hours}) { my ($h0, $h1) = split /,/, $o{hours}; $a0 = int($h0 * 3600 * $fs); $a1 = int($h1 * 3600 * $fs); $a1 = $R->{nt} if $a1 > $R->{nt} }
my $nt = $a1 - $a0; my $chan = sub { $R->{chan}->($_[0], from => $a0, to => $a1) };
my ($stage_all, $sname, $events, $other) = wfdb_stages_events(read_wfdb_annotations($afile), $R->{nt}); my $stage = $stage_all->slice("$a0:" . ($a1 - 1));
my $bad = zeroes(byte, $nt); my %evn;                                                       # inside an event or 10 s after it
for my $e (@$events) { my ($n, $s, $t) = @$e; $evn{$n}++; $s -= $a0; $t = $t - $a0 + 10 * $fs; next if $t < 0 || $s >= $nt; $s = 0 if $s < 0; $t = $nt - 1 if $t >= $nt; $bad->slice("$s:$t") .= 1 }

# ---- belt phase (scoring only) ----
my $belt = zp_butter($chan->($o{belt}), fs => $fs, low => 0.05, high => 2, order => 2);
my ($bre, $bim) = analytic_signal(zp_butter($chan->($o{belt}), fs => $fs, low => 0.1, high => 0.7, order => 2));
my $pu = unwrap_phase(atan2($bim, $bre)->slice('0:-1:10')->copy); my $tg = sequence(long, $pu->nelem) * 10; undef $bre; undef $bim;
my $kk = (($pu - $PI) / (2 * $PI))->floor; my $cr = which($kk->slice('1:-1') > $kk->slice('0:-2')); my $w6 = int(0.6 * $fs); my @onset;
for my $c ($cr->list) { my $s = $tg->at($c); my ($lo, $hi) = ($s - $w6 < 0 ? 0 : $s - $w6, $s + $w6 >= $nt ? $nt - 1 : $s + $w6); push @onset, $lo + $belt->slice("$lo:$hi")->minimum_ind->sclr }
my $bb = pdl(long, \@onset)->uniq; my $nbr = $bb->nelem - 1; my @amp;
for my $i (0 .. $nbr - 1) { my ($s, $e) = ($bb->at($i), $bb->at($i + 1)); my $seg = $belt->slice("$s:$e"); push @amp, ($seg->max - ($belt->at($s) + $belt->at($e)) / 2)->sclr }
my $amp = pdl(\@amp); my $len = ($bb->slice('1:-1') - $bb->slice('0:-2'))->double / $fs;
my $medamp = pdl(map { my ($l, $h) = ($_ - 10 < 0 ? 0 : $_ - 10, $_ + 10 >= $nbr ? $nbr - 1 : $_ + 10); $amp->slice("$l:$h")->median } 0 .. $nbr - 1);
my $okbr = ($len >= 1.5) & ($len <= 8) & ($amp > 0.3 * $medamp);

# ---- ECG, exactly as for BN1-BN2 ----
my %ix; @ix{ @{ $R->{labels} } } = 0 .. $#{ $R->{labels} };
my $ecg = zp_butter($chan->('ECG') * ($R->{units}[ $ix{ECG} ] =~ /^mV$/i ? 1000 : 1), fs => $fs, low => 0.5, high => 45, order => 4);
my ($r, $sg) = detect_r($ecg, $fs); my $x = $ecg * $sg; my $B = beat_features($x, $fs, r => $r); $r = $B->{r}; my $tb = $B->{tb}; my $nb = $tb->nelem; my %F = %{ $B->{F} };
my $medw = sub { my ($v, $hw) = @_; my $n = $v->nelem; pdl(map { my ($l, $h) = ($_ - $hw < 0 ? 0 : $_ - $hw, $_ + $hw >= $n ? $n - 1 : $_ + $hw); $v->slice("$l:$h")->median } 0 .. $n - 1) };
my $mrr = $medw->($F{RRp}, 5); my $rmed = $F{R}->median; my $rmad = 1.4826 * abs($F{R} - $rmed)->median;
my $okb = (abs($F{RRp} - $mrr) < 0.25 * $mrr) & (abs($F{RRn} - $mrr) < 0.25 * $mrr) & (abs($F{R} - $rmed) < 6 * $rmad); my $use = zeroes(long, $nb); $use->index($B->{use}) .= 1; $okb &= $use;
my $ue = which($okb); my $tk = $tb->index($ue); my $yR = $F{R}->index($ue);

# ---- changes of the heart's position, from the average QRS of 30 s windows (ECG only) ----
my ($q0, $q1) = (int(-0.06 * $fs), int(0.07 * $fs)); my $wl = 30; my $nw = int($nt / $fs / $wl); my $ru = $r->index($ue);
my $Q = beat_epochs($x->dummy(0), $ru, $q0, $q1)->slice('(0)'); my $wi = ($ru / ($wl * $fs))->floor->long;
my (@T, @has); for my $w (0 .. $nw - 1) { my $s = which($wi == $w); if ($s->nelem >= 10) { push @T, $Q->dice_axis(1, $s)->mv(1, 0)->average; push @has, 1 } else { push @T, undef; push @has, 0 } }
my @d = (0) x $nw; for my $w (1 .. $nw - 1) { my $p = $w - 1; $p-- while $p >= 0 && !$has[$p]; next unless $has[$w] && $p >= 0; $d[$w] = (sqrt((($T[$w] - $T[$p]) ** 2)->avg) / sqrt(($T[$p] ** 2)->avg))->sclr }
my $dp = pdl(@d)->where(pdl(@has)); my $thr = $dp->median->sclr + 6 * 1.4826 * abs($dp - $dp->median)->median->sclr; $thr = $o{dmin} if $thr < $o{dmin};
my @chg = grep { $d[$_] > $thr } 1 .. $nw - 1;

# ---- ECG phases (ECG only, whole span) ----
my $hp = resp_phase($tk, $yR, band => [0.15, 0.5]); my $hph = phase_at($hp, $tk);
my %E; my @M; my $add = sub { my ($k, $p) = @_; my $c1 = fourier_fit($p, $yR, 1); $E{$k} = wrap_phase($p - atan2($c1->at(2, 0), $c1->at(1, 0))); push @M, $k };
$add->('H 0.15-0.5', $hph); $add->('L w5 i3', phase_template($tk, $yR, $E{'H 0.15-0.5'}, win => 5, K => 2, iter => 3)->{phi});
$add->('M R,Qwid,Tlat', phase_template_multi($tk, cat(map { $F{$_}->index($ue) } qw(R Qwid Tlat)), $E{'H 0.15-0.5'}, win => 5, K => 2, iter => 3, gate => 0.1)->{phi});

# ---- stretches and scoring ----
my ($cyc) = interpol($ru->double, $bb->double, sequence($bb->nelem)); my $ci = $cyc->floor->long; my $cic = $ci->clip(0, $nbr - 1);
my $inbr = ($ci >= 0) & ($ci < $nbr) & $okbr->index($cic); my $pt = wrap_phase(2 * $PI * ($cyc - $cyc->floor));
my $bst = $stage->index($ru); my $bbad = $bad->index($ru); my $bw = $wi;
my %isch = map { ($_ => 1) } @chg; my @run; my ($cur, $start) = (undef, undef);              # stretches: same stage, no change window inside
my $flush = sub { my ($s, $e) = @_; push @run, [ $s, $e ] if defined $s && ($tk->at($e) - $tk->at($s)) >= $o{minlen} };
for my $i (0 .. $ue->nelem - 1) { my $st = $bst->at($i); my $brk = !defined($cur) || $st != $cur || ($i > 0 && $bw->at($i) != $bw->at($i - 1) && $isch{ $bw->at($i) });
    if ($brk) { $flush->($start, $i - 1) if defined $start; ($cur, $start) = ($st, $i) } }
$flush->($start, $ue->nelem - 1);
say_("%s %.2f-%.2f h, belt %s; stages in the annotation: %s; events: %s\n", $nm, $a0 / $fs / 3600, $a1 / $fs / 3600, $o{belt},
    join(' ', map { my $c = $_; sprintf '%s %.0f min', $sname->{$c}, ($stage == $c)->sum / $fs / 60 } grep { ($stage == $_)->sum } 0 .. 5), join(', ', map { "$_ $evn{$_}" } sort keys %evn));
say_("annotation texts that are neither a stage nor an event: %s\n", join(', ', map { "'$_' $other->{$_}" } sort keys %$other) || '-');
say_("ECG: %d beats, %d used; position changes flagged from the 30 s average QRS (d > %.3f): %d at %s\n", $nb, $ue->nelem, $thr, scalar @chg,
    join(' ', map { sprintf '%.2fh', ($a0 / $fs + $_ * $wl) / 3600 } @chg) || '-');
my %bys; for my $s (@run) { my ($i0, $i1) = @$s; push @{ $bys{ $bst->at($i0) } }, $s }
say_("stretches (>= %g s, one stage, no flagged change): %s\n\n", $o{minlen}, join(', ', map { sprintf '%s %d (%.0f min)', $sname->{$_}, scalar @{ $bys{$_} }, sum(pdl(map { $tk->at($_->[1]) - $tk->at($_->[0]) } @{ $bys{$_} })) / 60 } sort keys %bys));
say_("%-6s %-15s %5s %6s | %6s %6s | %6s | %11s\n", 'stage', 'ECG phase', 'str.', 'beats', 'SD', '>90', 'warp', 'off. spread');
open my $T, '>', "${P}_stage_eval.tsv" or die; print $T join("\t", 't_s', 'stage', 'stretch', 'belt', @M), "\n"; my $sid = 0; my %done;
for my $st (sort keys %bys) { for my $m (@M) { my (@err, @ptb, @offs, $n);
        for my $s (@{ $bys{$st} }) { my ($i0, $i1) = @$s; my $sel = sequence(long, $i1 - $i0 + 1) + $i0; $sel = $sel->where($inbr->index($sel) & ($bbad->index($sel) == 0)); next if $sel->nelem < 30;
            my $d = wrap_phase($E{$m}->index($sel) - $pt->index($sel)); my $off = atan2(sin($d)->avg, cos($d)->avg)->sclr; push @offs, $off; push @err, wrap_phase($d - $off); push @ptb, $pt->index($sel);
            unless ($done{"$st:$i0"}++) { $sid++; for my $j ($sel->list) { print $T join("\t", sprintf('%.3f', $tk->at($j) + $a0 / $fs), $sname->{$st}, $sid, sprintf('%.4f', $pt->at($j)), map { sprintf '%.4f', $E{$_}->at($j) } @M), "\n" } } }
        next unless @err; my $e = null; $e = $e->append($_) for @err; my $pp = null; $pp = $pp->append($_) for @ptb;
        my $bin = (($pp + $PI) / (2 * $PI) * 12)->floor->clip(0, 11); my @mb = map { my $w = which($bin == $_); $w->nelem ? atan2(sin($e->index($w))->avg, cos($e->index($w))->avg)->sclr * $D : 0 } 0 .. 11;
        say_("%-6s %-15s %5d %6d | %5.1f° %5.1f%% | %5.0f° | %9.1f°\n", $sname->{$st}, $m, scalar @offs, $e->nelem, circ_sd($e) * $D, 100 * (abs($e) > $PI / 2)->avg, pdl(@mb)->max - pdl(@mb)->min,
            @offs > 1 ? circ_sd(pdl(@offs)) * $D : 0) } }
say_("SD, >90, warp: error of the ECG phase against the belt phase after each stretch's own offset; off. spread: circular SD of the\nstretch offsets within the stage. Beats in events (and 10 s after) and in rejected breaths are not scored.\n");
close $T; close $TXT;
