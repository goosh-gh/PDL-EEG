#!/usr/bin/env perl
# examples/heartbeat_resp_quality_vs_belt.pl -- can the ECG alone tell how good its own respiratory phase is?  For a Challenge 2018 PSG record, in windows of --win s
# (120, step --step 60), the true error of the ECG phase (phase_template on the R size, as for BN1-BN2; against the belt phase,
# offset removed per window) is set beside quality indices that use the ECG only:
#   q_fold   share of the R-size variance in the window explained by its fold on the ECG phase (cos / sin, 1x and 2x)
#   q_snr    amplitude of that fold (max - min of the fitted curve) / SD of the residual
#   q_band   share of the power of the R size (4 Hz series, window) in 0.12-0.6 Hz, out of 0.02-1.5 Hz
#   q_agree  circular SD between the H phase (0.15-0.5 Hz band) and the L phase in the window, offset removed (small = methods agree)
#   q_rr     share of the RR-interval variance explained by the ECG phase (cos / sin) -- RSA visible
# Windows with fewer than 40 scored beats are skipped.  If <record>.arousal exists, the stage of each window (majority) and the share of
# beats inside events are written too (events are NOT removed here: the indices should also flag disturbed windows).
# Output: PREFIX_quality.tsv (one row per window) and PREFIX_quality.txt (Spearman correlation of each index with the true error, and the
# true error by quartile of each index).  Pool the TSVs of several records with heartbeat_resp_quality_pool.pl.
# usage: perl heartbeat_resp_quality_vs_belt.pl /path/RECORD [--hours A,B] [--out PREFIX] [--belt ABD] [--win 120] [--step 60]
use strict; use warnings; use FindBin; use lib "$FindBin::Bin/../lib"; use PDL; use Getopt::Long;
use PDL::EEG::IO::WFDB qw(read_wfdb read_wfdb_annotations wfdb_stages_events); use PDL::EEG::Heartbeat::BeatFeatures qw(beat_features);
use PDL::EEG::Filter qw(zp_butter analytic_signal); use PDL::EEG::Heartbeat qw(detect_r); use PDL::EEG::Heartbeat::Respiration qw(resp_phase phase_at phase_template phase_template_multi fourier_fit wrap_phase unwrap_phase circ_sd);
my %o = (belt => 'ABD', hours => '', win => 120, step => 60); GetOptions(\%o, 'hours=s', 'out=s', 'belt=s', 'win=f', 'step=f') or die "bad option\n";
my $stem = shift @ARGV or die "usage: $0 /path/RECORD [options]\n"; $stem =~ s/\.(mat|hea|arousal)$//; (my $nm = $stem) =~ s{.*/}{}; my $P = $o{out} // $nm;
my $PI = 3.14159265358979; my $D = 180 / $PI;
my $R = read_wfdb($stem); my $fs = $R->{fs}; my ($a0, $a1) = (0, $R->{nt});
if ($o{hours}) { my ($h0, $h1) = split /,/, $o{hours}; $a0 = int($h0 * 3600 * $fs); $a1 = int($h1 * 3600 * $fs); $a1 = $R->{nt} if $a1 > $R->{nt} }
my $nt = $a1 - $a0; my $chan = sub { $R->{chan}->($_[0], from => $a0, to => $a1) };
my ($stage, $sname, $inev);                                                                   # optional annotations
if (-r "$stem.arousal") { my ($st, $nmap, $ev) = wfdb_stages_events(read_wfdb_annotations("$stem.arousal"), $R->{nt}); $stage = $st->slice("$a0:" . ($a1 - 1)); $sname = $nmap;
    $inev = zeroes(byte, $nt); for my $e (@$ev) { my ($n, $s, $t) = @$e; $s -= $a0; $t -= $a0; next if $t < 0 || $s >= $nt; $s = 0 if $s < 0; $t = $nt - 1 if $t >= $nt; $inev->slice("$s:$t") .= 1 } }
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
# ---- ECG as for BN1-BN2 ----
my %ix; @ix{ @{ $R->{labels} } } = 0 .. $#{ $R->{labels} };
my $ecg = zp_butter($chan->('ECG') * ($R->{units}[ $ix{ECG} ] =~ /^mV$/i ? 1000 : 1), fs => $fs, low => 0.5, high => 45, order => 4);
my ($r, $sg) = detect_r($ecg, $fs); my $B = beat_features($ecg * $sg, $fs, r => $r); $r = $B->{r}; my $tb = $B->{tb}; my $nb = $tb->nelem; my %F = %{ $B->{F} };
my $medw = sub { my ($v, $hw) = @_; my $n = $v->nelem; pdl(map { my ($l, $h) = ($_ - $hw < 0 ? 0 : $_ - $hw, $_ + $hw >= $n ? $n - 1 : $_ + $hw); $v->slice("$l:$h")->median } 0 .. $n - 1) };
my $mrr = $medw->($F{RRp}, 5); my $rmed = $F{R}->median; my $rmad = 1.4826 * abs($F{R} - $rmed)->median;
my $okb = (abs($F{RRp} - $mrr) < 0.25 * $mrr) & (abs($F{RRn} - $mrr) < 0.25 * $mrr) & (abs($F{R} - $rmed) < 6 * $rmad); my $use = zeroes(long, $nb); $use->index($B->{use}) .= 1; $okb &= $use;
my $ue = which($okb); my $tk = $tb->index($ue); my $yR = $F{R}->index($ue); my $rrn = $F{RRn}->index($ue); my $ru = $r->index($ue);
my $hp = resp_phase($tk, $yR, band => [0.15, 0.5]); my $hph = phase_at($hp, $tk);
my $pH = $hph; my $pL = phase_template($tk, $yR, $pH, win => 5, K => 2, iter => 3)->{phi};
my ($cyc) = interpol($ru->double, $bb->double, sequence($bb->nelem)); my $ci = $cyc->floor->long; my $cic = $ci->clip(0, $nbr - 1);
my $sc = ($ci >= 0) & ($ci < $nbr) & $okbr->index($cic); my $pt = wrap_phase(2 * $PI * ($cyc - $cyc->floor));
my $Rg = sub { my ($y, $X) = @_; my $beta = inv($X->transpose x $X) x ($X->transpose x $y->dummy(0)); my $f = ($X x $beta)->flat; my $yc = $y - $y->avg;
    return ((sum(($f - $y->avg) ** 2) / sum($yc ** 2))->sclr, $f) };
my $csd = sub { my $d = shift; my $off = atan2(sin($d)->avg, cos($d)->avg); circ_sd(wrap_phase($d - $off)) * $D };
# ---- windows ----
open my $T, '>', "${P}_quality.tsv" or die; print $T join("\t", qw(rec t0_s stage ev_share beats true_sd q_fold q_snr q_band q_agree q_rr)), "\n"; my @rows;
for (my $s = $tk->at(0); $s + $o{win} <= $tk->at(-1); $s += $o{step}) { my $w = which(($tk >= $s) & ($tk < $s + $o{win})); next if $w->nelem < 40;
    my $ws = $w->where($sc->index($w)); next if $ws->nelem < 40;
    my $true = $csd->(wrap_phase($pL->index($ws) - $pt->index($ws)));
    my $p = $pL->index($w); my $X = cat(ones($p->nelem), cos($p), sin($p), cos(2 * $p), sin(2 * $p))->transpose; my ($qf, $f) = $Rg->($yR->index($w), $X);
    my $res = $yR->index($w) - $f; my $qs = (($f->max - $f->min) / sqrt(($res ** 2)->avg))->sclr;
    my $t4 = $s + sequence(int($o{win} * 4)) / 4; my ($y4) = interpol($t4, $tk->index($w), $yR->index($w)); $y4 = ($y4 - $y4->avg) * (0.5 - 0.5 * cos(2 * $PI * sequence($y4->nelem) / ($y4->nelem - 1)));
    my $N = 1; $N *= 2 while $N < 4 * $y4->nelem; my ($xr, $xi) = (zeroes($N), zeroes($N)); $xr->slice('0:' . ($y4->nelem - 1)) .= $y4; PDL::FFT::fft($xr, $xi);
    my $fq = sequence($N / 2) * 4 / $N; my $pw = ($xr ** 2 + $xi ** 2)->slice('0:' . ($N / 2 - 1)); my $qb = ($pw->where(($fq >= 0.12) & ($fq <= 0.6))->sum / $pw->where(($fq >= 0.02) & ($fq <= 1.5))->sum)->sclr;
    my $qa = $csd->(wrap_phase($pH->index($w) - $pL->index($w)));
    my ($qr) = $Rg->($rrn->index($w), cat(ones($p->nelem), cos($p), sin($p))->transpose);
    my ($stl, $evs) = ('-', '');
    if (defined $stage) { my $st = $stage->index($ru->index($w)); my $best = (sort { ($st == $b)->sum <=> ($st == $a)->sum } 0 .. 5)[0]; $stl = $sname->{$best}; $evs = sprintf '%.2f', $inev->index($ru->index($w))->double->avg }   # byte: avg would be truncated
    my @v = ($nm, sprintf('%.0f', $s + $a0 / $fs), $stl, $evs, $ws->nelem, sprintf('%.1f', $true), sprintf('%.3f', $qf), sprintf('%.2f', $qs), sprintf('%.3f', $qb), sprintf('%.1f', $qa), sprintf('%.3f', $qr));
    print $T join("\t", @v), "\n"; push @rows, [ @v ] }
close $T;
open my $TXT, '>', "${P}_quality.txt" or die; my $say = sub { my $x = sprintf(shift, @_); print $x; print $TXT $x };
my $rank = sub { my @v = @_; my @i = sort { $v[$a] <=> $v[$b] } 0 .. $#v; my @r; $r[ $i[$_] ] = $_ for 0 .. $#i; pdl(@r) };
my $sp = sub { my ($x, $y) = @_; my ($rx, $ry) = ($rank->(@$x), $rank->(@$y)); my ($a1, $b1) = ($rx - $rx->avg, $ry - $ry->avg); (sum($a1 * $b1) / sqrt(sum($a1 ** 2) * sum($b1 ** 2)))->sclr };
my @true = map { $_->[5] } @rows; $say->("%s: %d windows of %g s; true error (L phase vs belt, per window): median %.1f deg\n", $nm, scalar @rows, $o{win}, pdl(@true)->median);
my %col = (q_fold => 6, q_snr => 7, q_band => 8, q_agree => 9, q_rr => 10);
for my $q (qw(q_fold q_snr q_band q_agree q_rr)) { my @x = map { $_->[ $col{$q} ] } @rows; my $rho = $sp->(\@x, \@true);
    my @o = sort { $x[$a] <=> $x[$b] } 0 .. $#x; my $n = @o; my @qm = map { my @s = @o[ int($_ * $n / 4) .. int(($_ + 1) * $n / 4) - 1 ]; @s ? pdl(@true[@s])->median->sclr : 0 } 0 .. 3;
    my @cut = map { $x[ $o[ int($_ * $n / 4) ] ] } 1 .. 3;
    $say->("  %-8s Spearman with true error %+.2f | true error by quartile of the index (low -> high): %s  (cuts %s)\n", $q, $rho, join(' ', map { sprintf '%.0f', $_ } @qm), join(' ', @cut)) }
close $TXT;
