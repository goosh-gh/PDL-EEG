#!/usr/bin/env perl
# examples/heartbeat_resp_phase.pl -- per-beat respiratory phase from BN1-BN2 (EDF) with the ECG-only quality index q_agree beside it.
#
#   perl examples/heartbeat_resp_phase.pl a.edf ... [--bn1 BN1] [--bn2 BN2] [--out DIR] [--win 120] [--min 40] [--cut 20,25,30]
#
# BN1-BN2 is band-passed 0.5-45 Hz (zero phase, order 4), R peaks are detected (detect_r), beats less than 1 s from either end of
# the record are dropped, and the per-beat features of PDL::EEG::Heartbeat::BeatFeatures are measured. For every usable beat:
#   phase_H  phase of the R size (R amplitude minus its 20 s running mean, %) in the 0.15-0.5 Hz band (resp_phase, phase_at)
#   phase_L  phase_template on the R size, started from phase_H (window 5 s, K = 2, 3 rounds) -- the phase to use
#   q_agree  phase_agreement of phase_H and phase_L: circular SD of their difference over the beats within +-win/2 s
#   class    from q_agree: good < 20 deg, ok 20-25, border 25-30, poor >= 30, na (fewer than --min beats in the window);
#            the cut-offs (--cut) are those where, on 10 Challenge 2018 training records (heartbeat_resp_quality_vs_belt.pl,
#            heartbeat_resp_quality_pool.pl), the true error of phase_L against the belt had a median of 17-27, 35, 45 and 64+ deg
# The phase origin (0 deg) is arbitrary: it is not the start of inspiration.
#
# Output per recording: OUTDIR/<name>_resp_phase.tsv (one row per usable beat: beat index among the detected beats, t_s, RR ending at
# the beat [ms], R size [%], phase_H, phase_L [deg], q_agree [deg], class) and a summary line on STDOUT.
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use PDL;
use PDL::EEG::IO::EDF                qw(read_edf);
use PDL::EEG::Filter                 qw(zp_butter);
use PDL::EEG::Heartbeat              qw(detect_r);
use PDL::EEG::Heartbeat::BeatFeatures qw(beat_features);
use PDL::EEG::Heartbeat::Respiration qw(resp_phase phase_at phase_template phase_agreement wrap_phase);

my %o = (bn1 => 'BN1', bn2 => 'BN2', out => '.', win => 120, min => 40, cut => '20,25,30');
GetOptions(\%o, 'bn1=s', 'bn2=s', 'out=s', 'win=f', 'min=i', 'cut=s') or die "bad options\n";
die "usage: $0 a.edf ... [options]\n" unless @ARGV;
my @cut = split /,/, $o{cut};
die "--cut needs 3 values\n" unless @cut == 3;
my @cname = qw(good ok border poor);
my $D = 180 / (4 * atan2(1, 1));

printf "%-14s %6s %6s | q_agree median | %6s %6s %6s %6s %6s  (share of beats)\n", 'recording', 'beats', 'min', @cname, 'na';
for my $f (@ARGV) {
    (my $nm = $f) =~ s{.*/}{};
    $nm =~ s/\.edf$//i;
    my $e  = read_edf($f);
    my $fs = $e->{fs};
    my %ix;
    @ix{ @{ $e->{labels} } } = 0 .. $#{ $e->{labels} };
    die "$f: no channel $o{bn1} / $o{bn2}\n" unless defined $ix{ $o{bn1} } && defined $ix{ $o{bn2} };
    my $x = zp_butter($e->{data}->slice("($ix{$o{bn1}})")->double - $e->{data}->slice("($ix{$o{bn2}})")->double,
                      fs => $fs, low => 0.5, high => 45, order => 4);
    my ($r) = detect_r($x, $fs);
    $r = $r->where(($r >= $fs) & ($r + $fs < $x->nelem));
    my $B  = beat_features($x, $fs, r => $r);
    my $u  = $B->{use};
    my $tk = $B->{tb}->index($u);
    my $y  = $B->{F}{R}->index($u);
    my $rr = $B->{F}{RRp}->index($u);
    my $pH = phase_at(resp_phase($tk, $y, band => [0.15, 0.5]), $tk);
    my $pL = phase_template($tk, $y, $pH, win => 5, K => 2, iter => 3)->{phi};
    my $q  = phase_agreement($tk, $pH, $pL, win => $o{win}, min => $o{min}) * $D;
    my @q  = $q->list;
    my @cl = map { $_ != $_ ? 'na' : $_ < $cut[0] ? 'good' : $_ < $cut[1] ? 'ok' : $_ < $cut[2] ? 'border' : 'poor' } @q;

    my $tsv = "$o{out}/${nm}_resp_phase.tsv";
    open my $T, '>', $tsv or die "$tsv: $!";
    print $T join("\t", qw(beat t_s rr_ms R_pct phase_H phase_L q_agree class)), "\n";
    my @ui = $u->list;
    my @t  = $tk->list;
    my @rl = $rr->list;
    my @yl = $y->list;
    my @hl = ($pH * $D)->list;
    my @ll = (wrap_phase($pL) * $D)->list;
    for my $k (0 .. $#t) {
        print $T join("\t", $ui[$k], sprintf('%.3f', $t[$k]), sprintf('%.0f', $rl[$k]), sprintf('%.3f', $yl[$k]),
            sprintf('%.1f', $hl[$k]), sprintf('%.1f', $ll[$k]), $q[$k] != $q[$k] ? 'nan' : sprintf('%.1f', $q[$k]), $cl[$k]), "\n";
    }
    close $T;
    my %cnt;
    $cnt{$_}++ for @cl;
    my $qf = $q->where($q->isfinite);
    printf "%-14s %6d %6.1f | %14s | %s\n", $nm, scalar @t, ($t[-1] - $t[0]) / 60, $qf->nelem ? sprintf('%.1f', $qf->median) : 'nan',
        join(' ', map { sprintf '%5.0f%%', 100 * ($cnt{$_} // 0) / @t } @cname, 'na');
}
print "q_agree cut-offs @cut deg, window $o{win} s centred on each beat, at least $o{min} beats; TSV per recording in $o{out}/\n";
