#!/usr/bin/env perl
# examples/heartbeat_hep_resp_grid.pl
#
#   perl examples/heartbeat_hep_resp_grid.pl REC.edf|REC.vhdr|REC.EEG [options]
#
#   REC.EEG (Nihon Kohden) is read with every waveform block (read_nk, all_blocks):
#   the blocks are joined end to end, so beats whose epoch comes within 2 s of a
#   recording break, and non-heartbeat time points as near, are left out.  The
#   breaks are the physical blocks, or for a file held in one block (EEG-1200A)
#   the .LOG "REC START" markers (PDL::EEG::IO::NihonKohden::block_ranges).  --nk-block N reads
#   block N only.  A BrainVision recording is given by its .vhdr.
#
#     --bn1 L --bn2 L   non-cephalic electrodes (cardiac regressors, R peaks on BN1-BN2)
#     --ecg L           instead of BN1, BN2: one ECG channel (a bipolar chest or limb
#                       lead), used for the R peaks, the respiratory variable and as
#                       the single cardiac regressor
#     --ref KIND        reference of the head channels: bne (BN = R*BN1 + (1-R)*BN2),
#                       sys (as recorded), a1, a2, ears ((A1+A2)/2), chans (the
#                       mean of the channels of --ref-chans, e.g. TP9,TP10), or avg
#                       (the mean of all head channels used, BN1 and BN2 not included)
#                       (bne; with --ecg: sys)
#     --ref-chans LIST  channels averaged for --ref chans
#     --head LIST       head channels, comma separated (the 10-20 channels found,
#                       old or new names, plus A1, A2)
#     --aux LIST        other channels (e.g. PPG) averaged over the same beats, as
#                       recorded, band-passed like the head; written to the TSVs
#                       as columns NAME* scaled to a peak-to-peak of 8 (arbitrary
#                       units, to be drawn with the uV channels), and the latency
#                       of their steepest rise after R is printed
#     --pulse CH        a pulse channel (e.g. PPG): for every beat, the time of the
#                       pulse upstroke gives the pulse transit time PTT (first guess:
#                       steepest rise of CH, low-passed at 15 Hz, 100..450 ms after
#                       R, the sign of CH taken from its average; then a template of
#                       the upstroke, -80..+80 ms, matched to each beat within +-40 ms,
#                       parabolic peak, three rounds; beats matching with r <= 0.8 are
#                       left out); then, for every head channel, the
#                       residual (cardiac field removed, without respiration) is
#                       compared aligned to R and aligned to the pulse (each beat
#                       shifted by its PTT minus the mean PTT), and the lag between
#                       the averages of the beats with the longest and the shortest
#                       third of PTT is measured.
#                       Writes PREFIX_avg_pulsealign.tsv (the pulse-aligned residual,
#                       on the R time axis shifted by the mean PTT)
#     --ratio R         R of the balanced non-cephalic reference              (0.673)
#     --leads KIND      cardiac regressors: diff (BN1-BN2 only, free of the
#                       head reference), both (BN1 and BN2 as recorded, against the
#                       system reference), or reref (BN1 and BN2 each re-referenced
#                       like the head channels; with --ref a1|a2|ears|chans two
#                       cardiac leads that bring no system reference into the
#                       residual: the residual is h - g1 b1 - g2 b2 - (1-g1-g2) r
#                       with h, b1, b2, r the potentials against any common point
#                       and r the chosen reference)
#                       (diff with --ref bne, both with --ref sys)
#     --cardiac LIST    more cardiac leads besides BN1 and BN2 (e.g. ClavL,ClavR,Th12),
#                       each re-referenced like the head channels; needs --leads reref.
#                       Leads in other directions (left-right) let the fit follow a
#                       T wave whose field points elsewhere than the QRS, and, with
#                       the respiration gains, a field that turns with breathing.
#                       They are not head channels: not shown, not used for the rejection.
#                       An entry X-Y is a bipolar lead X minus Y, free of any reference,
#                       and X, Y may be head channels: e.g. A1-A2 as a left-right lead
#                       for a recording without one (it also takes out any brain
#                       activity that differs between the two ears)
#     --skip S          ignore the first S seconds                              (10)
#     --band LO,HI      zero-phase Butterworth band of head and BN channels [Hz] (0.5,40)
#     --fit A,B         window [ms] in which the cardiac field is fitted         (-250,600)
#     --exclude LIST    head channels left out                                  (none)
#     --bn-channels     also show BN1 and BN2 themselves, in the chosen reference
#                       and with their own cardiac fit, as two more channels
#                       (not used for the beat rejection)
#     --reject UV       a beat is dropped when a head channel's peak-to-peak
#                       over -300..+700 ms exceeds this                          (200)
#     --nsurr N         non-heartbeat averages for the band                      (100)
#     --chunk N         beats per chunk when accumulating (memory)              (500)
#     --ylim U          y range of the grid [uV], +-U                            (6)
#     --show A,B        time range of the grid and joint plots [ms]             (-250,600)
#     --times LIST      latencies [ms] of the maps in the joint plots  (0,150,250,330,450,550)
#     --positive-up     positive up (default: negative up)
#     --out PREFIX      output prefix (REC without .edf + _hepresp)
#
# What is left of the R-locked average of the head channels once the cardiac
# field is taken out, with and without the respiratory variation of that field.
#
# Reference: by default the balanced non-cephalic reference BN = R*BN1 +
# (1-R)*BN2 (R = 0.673 minimises the QRS of the R-locked average of BN in
# task.edf, HBrBP 002), the closest to a monopolar recording.  BN1 and BN2 are
# recorded against the system reference on the head (the mean of C3 and C4).
#
# Cardiac field: every head channel is modelled, beat by beat, as a linear and
# instantaneous mixture of the cardiac leads of the same beat,
#   no respiration:    CFA_k(t) = sum_j g_j L_j,k(t)
#   with respiration:  CFA_k(t) = sum_j (g_j + h_j rho_k) L_j,k(t)
# The leads L are BN1-BN2 (--leads diff), which carries no head potential, so
# the fit cannot take any part of the head reference into the cardiac field and
# the reference stays the one chosen with --ref; or BN1 and BN2 as recorded (--leads both), which can
# fit any cardiac waveform the two electrodes see, but whose fit also decides,
# channel by channel, how much of the system reference enters the residual: with
# both leads the residual is the same for every reference that is a mixture of
# the system reference, BN1 and BN2, so --ref bne --leads both gives exactly the
# residual of --ref sys --leads both.
# rho_k is the respiratory variable of beat k from the ECG alone (R size ->
# resp_phase 0.15-0.5 Hz -> phase_template, window 5 s, K = 2, 3 rounds; the
# breath template at each beat's phase; mean 0, SD 1), as in heartbeat_bnref_resp.pl.
# The coefficients are fitted per channel on averages only, inside --fit:
#   mean_k E_k        = sum_j g_j mean(L_j) + h_j mean(rho L_j)
#   mean_k rho_k E_k  = sum_j g_j mean(rho L_j) + h_j mean(rho^2 L_j)
# (without respiration: the first line with the g_j only).  Single beats are not
# used for the fit: BN1 and BN2 are recorded against the system reference on
# the head, so a single beat of BN1 or BN2 carries the EEG at the reference,
# which is also in every head channel and would be fitted as if it were cardiac.
# In averages that EEG cancels.  The second line uses how the head channels and
# the two leads change with rho; it is what lets h1, h2 be estimated.
# The cardiac field in the average is the right-hand side of the first line,
# and the residual is the average minus it.  With the default window the fit
# takes whatever in the average has the time course of the two leads (an upper
# bound of the cardiac field); with --fit -45,55 only the QRS sets the gains.
#
# Non-heartbeat band: averages at --nsurr sets of time points drawn uniformly over the record,
# minus their cardiac fit with the no-respiration coefficients; pointwise 95th
# percentile of the absolute value.
#
# Text: share of the fitted averages explained by each model inside --fit
# (median and range over channels); per channel the residual over 200..300,
# 300..450 and 450..600 ms for both models with z against the non-heartbeat
# averages; the per-beat residual (with respiration) over the same windows
# regressed on rho (slope per SD, t).  PREFIX_avg_raw.tsv, PREFIX_avg_res.tsv,
# PREFIX_avg_resrho.tsv: the average and the residuals without and with
# respiration; PREFIX_avg_rho1..3.tsv: the residual (with respiration)
# averaged over the beats in the lowest, middle and highest third of rho; in the format of
# examples/avg_grid_scrub.pl (5x5 10-20 grid, scalp map by MAP2D, scrubbing in a
# giza-server window or --png).  PREFIX_joint_raw.png,
# PREFIX_joint_res.png, PREFIX_joint_resrho.png: PDL::EEG::MAP2D::plot_joint
# (waveforms of the scalp channels and maps at --times, one colour scale per
# figure) of the average, and of the residual without and with respiration.

use strict;
use warnings;
use PDL;
use PDL::MatrixOps ();
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use PDL::EEG::IO::EDF                 qw(read_edf);
use PDL::EEG::IO::BrainVision         qw(read_brainvision);
use PDL::EEG::Filter                  qw(zp_butter);
use PDL::EEG::Heartbeat               qw(detect_r beat_epochs );
use PDL::EEG::Heartbeat::BeatFeatures qw(beat_features);
use PDL::EEG::Heartbeat::Respiration  qw(resp_phase phase_at phase_template fourier_eval);
use PDL::EEG::MAP2D                   qw(plot_joint);
use PDL::EEG::IO::ASA                 ();

my ($bn1, $bn2, $skip, $band, $fitw, $exclude, $reject, $nsurr, $chunk, $ylim, $posup, $out)
    = ('BN1', 'BN2', 10, '0.5,40', '-250,600', '', 200, 100, 500, 6, 0);
my $show  = '-250,600';
my ($ref, $ratio, $leads, $withbn) = (undef, 0.673, undef, 0);
my ($ecgch, $refchans, $headlist, $auxlist, $pulsech, $cardlist) = (undef, '', '', '', undef, '');
my $nkblock;
my $mtimes = '0,150,250,330,450,550';
GetOptions('bn1=s' => \$bn1, 'bn2=s' => \$bn2, 'skip=f' => \$skip, 'band=s' => \$band,
           'fit=s' => \$fitw, 'exclude=s' => \$exclude, 'reject=f' => \$reject,
           'nsurr=i' => \$nsurr, 'chunk=i' => \$chunk, 'ylim=f' => \$ylim, 'show=s' => \$show, 'times=s' => \$mtimes,
           'positive-up' => \$posup, 'out=s' => \$out,
           'ref=s' => \$ref, 'ratio=f' => \$ratio, 'leads=s' => \$leads, 'bn-channels' => \$withbn,
           'ecg=s' => \$ecgch, 'ref-chans=s' => \$refchans, 'head=s' => \$headlist, 'aux=s' => \$auxlist, 'pulse=s' => \$pulsech, 'nk-block=i' => \$nkblock, 'cardiac=s' => \$cardlist) or die;
$ref //= defined $ecgch ? 'sys' : 'bne';
$ref =~ /^(bne|sys|a1|a2|ears|chans|avg)$/ or die "--ref must be bne, sys, a1, a2, ears, chans or avg\n";
if (defined $ecgch) {
    $ref ne 'bne' or die "--ref bne needs BN1 and BN2, not --ecg\n";
    !$withbn       or die "--bn-channels needs BN1 and BN2, not --ecg\n";
    !defined $leads || $leads eq 'diff' or die "with --ecg the only cardiac regressor is that channel\n";
    $leads = 'diff';
}
$leads //= $ref eq 'sys' ? 'both' : 'diff';
$leads =~ /^(diff|both|reref)$/ or die "--leads must be diff, both or reref\n";
$leads ne 'reref' || ($ref ne 'sys' && $ref ne 'bne') or die "--leads reref needs --ref a1, a2, ears, chans or avg\n";
$ref ne 'chans' || $refchans ne '' or die "--ref chans needs --ref-chans LIST\n";
my @xcard = grep { length } map { s/^\s+|\s+$//gr } split /,/, $cardlist;
my @xbip  = map { [ split /-/, $_, 2 ] } grep { /^[^-]+-[^-]+$/ } @xcard;   # bipolar X-Y
@xcard    = grep { !/-/ } @xcard;                                           # single channels
!(@xcard || @xbip) || ($leads eq 'reref' && !defined $ecgch) or die "--cardiac needs --leads reref (and BN1, BN2)\n";
my $file = shift @ARGV or die "usage: $0 REC.edf|REC.vhdr|REC.EEG [options]\n";
!@ARGV or die "unexpected argument(s) after $file: @ARGV (a space after a comma in a list option?)\n";
($out = $file) =~ s/\.(edf|vhdr|eeg)$//i, $out .= '_hepresp' unless defined $out;
my ($f_lo, $f_hi) = split /,/, $band;
my ($fw0, $fw1)   = split /,/, $fitw;
my %excl = map { $_ => 1 } split /,/, $exclude;

# BrainVision: read only the channels this run can use (a long recording of many
# channels does not fit in memory as double otherwise)
my @maybe = ($headlist ne '' ? split(/,/, $headlist)
                             : qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 T7 C3 Cz C4 T4 T8 T5 P7 P3 Pz P4 T6 P8 O1 O2 A1 A2),
             (defined $ecgch ? $ecgch : ($bn1, $bn2, @xcard, map { @$_ } @xbip)), grep({ length } split(/,/, $auxlist)),
             (defined $pulsech ? $pulsech : ()), map { s/^\s+|\s+$//gr } split(/,/, $refchans));
my %seen;
@maybe = grep { !$seen{$_}++ } @maybe;
my @joins;                                                   # samples where blocks were joined
my $e;
if ($file =~ /\.vhdr$/i) {
    $e = read_brainvision($file, channels => \@maybe, skip_missing => 1);
} elsif ($file =~ /\.eeg$/i) {
    require PDL::EEG::IO::NihonKohden;
    defined((PDL::EEG::IO::NihonKohden::nk_layout($file))[1])
        or die "$file: not a Nihon Kohden recording (a BrainVision recording is given by its .vhdr)\n";
    $e = defined $nkblock ? PDL::EEG::IO::NihonKohden::read_nk($file, block => $nkblock)
                          : PDL::EEG::IO::NihonKohden::read_nk($file, all_blocks => 1);
    # recording breaks: the physical blocks (EEG-1100 wfmblock files), or, for a file
    # held in one block (EEG-1200A extblock), the .LOG "REC START" markers
    @joins = grep { $_ > 0 } map { $_->{start} } @{ PDL::EEG::IO::NihonKohden::block_ranges($e) }
        unless defined $nkblock;
    printf "%s: Nihon Kohden, %s, %d block(s), %.1f min%s\n", $file,
        defined $nkblock ? "block $nkblock" : 'all blocks joined', $e->{n_blocks} // 1,
        $e->{data}->dim(1) / $e->{fs} / 60, @joins ? sprintf(', recording breaks at %s s', join(', ', map { sprintf '%.0f', $_ / $e->{fs} } @joins)) : '';
} else {
    $e = read_edf($file);
    # recording breaks of a joined recording: the "Block N" annotations that write_edf
    # (nk_to_edf.pl) puts at every block start after the first
    @joins = map { int($_->{onset} * $e->{fs} + 0.5) }
             grep { ($_->{label} // '') =~ /^Block \d+$/ && $_->{onset} > 0 } @{ $e->{events} // [] };
    printf "%s: EDF, %.1f min%s\n", $file, $e->{data}->dim(1) / $e->{fs} / 60,
        @joins ? sprintf(', recording breaks (Block annotations) at %s s', join(', ', map { sprintf '%.0f', $_ / $e->{fs} } @joins))
               : ', no Block annotations (no recording breaks excluded)';
}
my $fs = $e->{fs};
my %ix; @ix{ @{ $e->{labels} } } = 0 .. $#{ $e->{labels} };
my @cardiac = defined $ecgch ? ($ecgch) : ($bn1, $bn2, @xcard);   # channels after the head
my %headset = map { $_ => 1 } grep { exists $ix{$_} && !$excl{$_} }
              ($headlist ne '' ? split(/,/, $headlist)
                               : qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 T7 C3 Cz C4 T4 T8 T5 P7 P3 Pz P4 T6 P8 O1 O2 A1 A2));
my %incard = map { $_ => 1 } @cardiac;
for my $nm (map { @$_ } @xbip) {                             # bipolar parts that are neither: load them too
    push @cardiac, $nm unless $headset{$nm} || $incard{$nm}++;
}
exists $ix{$_} or die "channel '$_' not in $file\n" for @cardiac;
my @head = grep { exists $ix{$_} && !$excl{$_} }
           ($headlist ne '' ? split(/,/, $headlist)
                            : qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 T7 C3 Cz C4 T4 T8 T5 P7 P3 Pz P4 T6 P8 O1 O2 A1 A2));
my %ishead = map { $_ => 1 } @head;
$ishead{$_} and die "--cardiac: $_ is a head channel (leave it out with --exclude $_)\n" for @xcard;
my $nh = @head;
my $raw = $e->{data};

# ---- R peaks and the respiratory variable (ECG only) --------------------------------
my $ecg = zp_butter(defined $ecgch ? $raw->slice("($ix{$ecgch})")->double
                                  : $raw->slice("($ix{$bn1})")->double - $raw->slice("($ix{$bn2})")->double,
                    fs => $fs, low => 0.5, high => 45, order => 4);
my ($rall, $sign) = detect_r($ecg, $fs);
$ecg  = $ecg * $sign;
$rall = $rall->where(($rall >= $fs) & ($rall + $fs < $ecg->nelem));
my $BF  = beat_features($ecg, $fs, r => $rall);
my $u   = $BF->{use};
my $tk  = $BF->{tb}->index($u);
my $ysz = $BF->{F}{R}->index($u);
my $pL  = phase_template($tk, $ysz, phase_at(resp_phase($tk, $ysz, band => [0.15, 0.5]), $tk),
                         win => 5, K => 2, iter => 3);
my $rho = fourier_eval($pL->{c}, $pL->{phi})->flat;
my $phi = $pL->{phi}->flat->copy;                           # respiratory phase of each beat [rad]
my $rb  = $BF->{r}->index($u);

# ---- head (as recorded) and the two leads, same band --------------------------------
my $Yr = zp_butter($raw->dice_axis(0, pdl(long, [ @ix{@head, @cardiac} ]))->double,
                   fs => $fs, low => $f_lo, high => $f_hi, order => 4);   # (nh+1 or nh+2, nt), as recorded
my ($B1, $B2) = defined $ecgch ? () : map { $Yr->slice("($_),:") } $nh, $nh + 1;
my $nrej = $nh;                                              # channels used for the beat rejection
my $Hd = $Yr->slice("0:" . ($nh - 1) . ",:");
my %ypos; @ypos{ @head, @cardiac } = 0 .. $nh + $#cardiac;    # rows of $Yr, as recorded
my %hpos; @hpos{@head} = 0 .. $nh - 1;
if ($withbn) {                                               # BN1, BN2 as two more channels
    $Hd = $Yr->slice("0:" . ($nh + 1) . ",:");
    push @head, $bn1, $bn2;
    $nh += 2;
}
my @rc = grep { length } map { s/^\s+|\s+$//gr } split /,/, $refchans;   # 'TP9, TP10' as well
my $avgref = $ref eq 'avg';
($ref, @rc) = ('chans', @head[0 .. $nrej - 1]) if $avgref;    # every head channel, not BN1/BN2
for my $need ($ref eq 'a1' ? 'A1' : $ref eq 'a2' ? 'A2' : $ref eq 'ears' ? ('A1', 'A2') : $ref eq 'chans' ? @rc : ()) {
    exists $hpos{$need} or die "--ref $ref needs channel $need among the head channels\n";
}
my $refsig = $ref eq 'bne'  ? $ratio * $B1 + (1 - $ratio) * $B2
           : $ref eq 'a1'   ? $Hd->slice("($hpos{A1}),:")
           : $ref eq 'a2'   ? $Hd->slice("($hpos{A2}),:")
           : $ref eq 'ears' ? 0.5 * ($Hd->slice("($hpos{A1}),:") + $Hd->slice("($hpos{A2}),:"))
           : $ref eq 'chans' ? $Hd->dice_axis(0, pdl(long, [ @hpos{@rc} ]))->average
           : undef;
$Hd = $Hd - $refsig->dummy(0) if defined $refsig;
my $Ld = defined $ecgch     ? $Yr->slice("($nrej),:")->dummy(0)  # the ECG channel (after the head)
       : $leads eq 'diff'  ? ($B1 - $B2)->dummy(0)
       : $leads eq 'reref' ? cat((map { $Yr->slice("($_),:") - $refsig } $nrej .. $nrej + 1 + @xcard),
                                 (map { $Yr->slice("($ypos{$_->[0]}),:") - $Yr->slice("($ypos{$_->[1]}),:") } @xbip)
                                )->transpose                                              # (2+nx+nb, nt)
       :                     $Yr->slice("$nrej:" . ($nrej + 1) . ",:");
my $nl = $Ld->dim(0);
my $Y  = $Hd->glue(0, $Ld);                                  # (nh+nl, nt): head, then leads
my $nt = $Y->dim(1);
my $s  = sub { my $v = shift() * 1e-3 * $fs; $v >= 0 ? int($v + 0.5) : -int(-$v + 0.5) };
my ($n0, $n1) = ($s->(-300), $s->(700));
my $nep = $n1 - $n0;
my @bw  = ($s->(-75) - $n0, $s->(-50) - $n0 - 1);
my @fw  = ($s->($fw0) - $n0, $s->($fw1) - $n0 - 1);
my $tms = (sequence($nep) + $n0) * (1000 / $fs);

my $ok = ($rb >= $s->(1000 * $skip) - $n0) & ($rb + $n1 <= $nt) & ($rb + $n0 >= 0) & _off_joins($rb);
my @count = ($rall->nelem, $u->nelem, $ok->sum);
$rb  = $rb->where($ok);
$rho = $rho->where($ok);
$phi = $phi->where($ok);
my ($keep, $nbad) = _accept($rb);
printf "%s: beats: %d R peaks, %d with beat features, %d away from the ends and breaks, %d within %g uV\n",
    $file, @count, $keep->sum, $reject;
if ($keep->sum < 0.5 * $rb->nelem) {                          # which channels throw the beats out
    my @w = grep { $nbad->at($_) > 0 } 0 .. $nrej - 1;
    @w = sort { $nbad->at($b) <=> $nbad->at($a) } @w;
    printf "  beats over %g uV peak-to-peak, by channel: %s\n", $reject,
        join(', ', map { sprintf '%s %d', $head[$_], $nbad->at($_) } @w) || '-';
}
$keep->sum >= 30 or die sprintf "%s: only %d beat(s) left; %s\n", $file, $keep->sum,
    $rb->nelem == 0 ? 'no R peak usable (check the R detection on BN1-BN2 and the break positions)'
                    : 'leave out the channels above with --exclude, or raise --reject';
$rb  = $rb->where($keep);
$rho = $rho->where($keep);
$phi = $phi->where($keep);
$rho = $rho - $rho->avg;
$rho = $rho / sqrt(($rho**2)->avg);
my $nb = $rb->nelem;
my $refdesc = $ref eq 'bne'  ? sprintf('BN = %.3f %s + %.3f %s', $ratio, $bn1, 1 - $ratio, $bn2)
            : $ref eq 'a1'   ? 'A1'
            : $ref eq 'a2'   ? 'A2'
            : $ref eq 'ears' ? '(A1+A2)/2'
            : $avgref         ? sprintf('average of the %d head channels', scalar @rc)
            : $ref eq 'chans' ? 'mean of ' . join(',', @rc)
            :                  'system reference (as recorded)';
my $leaddesc = defined $ecgch ? $ecgch : $leads eq 'diff' ? "$bn1-$bn2"
             : $leads eq 'reref' ? join(', ', $bn1, $bn2, @xcard) . " against the head reference" . join('', map { ", $_->[0]-$_->[1] (bipolar)" } @xbip) : "$bn1, $bn2 as recorded";
printf "%s: fs %g Hz, %d head channels, reference %s, cardiac leads %s, %d beats (reject %g uV), fit window %s ms\n",
    $file, $fs, $nh, $refdesc, $leaddesc, $nb, $reject, $fitw;

# ---- averages: head, rho x head, and the moments of the leads ---------------------------
# Mbar rows: mean L (nl), mean rho L (nl), mean rho^2 L (nl)
my ($Ebar, $Erho, $Mbar) = (zeroes(double, $nh, $nep), zeroes(double, $nh, $nep), zeroes(double, 3 * $nl, $nep));
_chunks($rb, sub {
    my ($E, $xs) = @_;                                       # (nh+2, nep, nb), (nb)
    my $x3 = $xs->dummy(0)->dummy(0);                        # (1, 1, nb)
    my $H  = $E->slice("0:" . ($nh - 1) . ",:,:");
    my $L  = $E->slice("$nh:" . ($nh + $nl - 1) . ",:,:");     # (nl, nep, nb)
    $Ebar += $H->mv(2, 0)->sumover;
    $Erho += ($H * $x3)->mv(2, 0)->sumover;
    $Mbar += $L->glue(0, $L * $x3, $L * $x3**2)->mv(2, 0)->sumover;
});
$_ /= $nb for $Ebar, $Erho, $Mbar;

# ---- fit of the gains on the averages inside the fit window ----------------------------
my ($f0, $f1) = @fw;
my %G;                                                       # model => (nh, p): row j = regressor j
my %r2;
{
    # no respiration: Ebar = sum_j g_j L_j
    my $X  = $Mbar->slice("0:" . ($nl - 1) . ",$f0:$f1");    # (nl, nw): nw rows x nl cols
    my $Yt = $Ebar->slice(":,$f0:$f1");                      # (nh, nw)
    $G{nores} = PDL::MatrixOps::inv($X->transpose x $X) x ($X->transpose x $Yt);
    $r2{nores} = 1 - ((($Yt - $X x $G{nores})**2)->mv(1, 0)->sumover / ($Yt**2)->mv(1, 0)->sumover);
    # with respiration: both moment equations stacked over the window
    my $X1 = $Mbar->slice("0:" . (2 * $nl - 1) . ",$f0:$f1");        # Ebar on L, rho L
    my $X2 = $Mbar->slice("$nl:" . (3 * $nl - 1) . ",$f0:$f1");      # Erho on rho L, rho^2 L
    my $Xs = $X1->glue(1, $X2);                              # (2 nl, 2 nw)
    my $Ys = $Yt->glue(1, $Erho->slice(":,$f0:$f1"));         # (nh, 2 nw)
    $G{resp} = PDL::MatrixOps::inv($Xs->transpose x $Xs) x ($Xs->transpose x $Ys);
    my $fitE = $X1 x $G{resp};
    $r2{resp} = 1 - ((($Yt - $fitE)**2)->mv(1, 0)->sumover / ($Yt**2)->mv(1, 0)->sumover);
    my $Er = $Erho->slice(":,$f0:$f1");
    my $r2rho_n = 1 - ((($Er - $Mbar->slice("$nl:" . (2 * $nl - 1) . ",$f0:$f1") x $G{nores})**2)->mv(1, 0)->sumover / ($Er**2)->mv(1, 0)->sumover);
    my $r2rho_r = 1 - ((($Er - $X2 x $G{resp})**2)->mv(1, 0)->sumover / ($Er**2)->mv(1, 0)->sumover);
    printf "share of the average explained inside the fit window: without rho median %.3f (%.3f..%.3f), with rho %.3f (%.3f..%.3f)\n",
        map { my $v = $r2{$_}->where($r2{$_}->isfinite); ($v->median, $v->minmax) } qw(nores resp);   # the reference channel itself is 0/0
    printf "share of the rho-weighted average (how the head changes with breathing) explained: without rho %.3f, with rho %.3f (medians)\n",
        map { $_->where($_->isfinite)->median } $r2rho_n, $r2rho_r;
}

# cardiac field in the average and the residuals
my %res;
for my $nm (qw(nores resp)) {
    my $p = $G{$nm}->dim(1);
    my $C = $Mbar->slice("0:" . ($p - 1) . ",:") x $G{$nm};  # (nep x p)(p x nh) -> (nh, nep)
    $res{$nm} = $Ebar - $C;
}

# ---- non-heartbeat averages, fitted with the no-respiration coefficients -----------------
my @sres;
for my $k (1 .. $nsurr) {
    # time points drawn uniformly over the usable record, as many as the beats:
    # shifting every R by 0.25..0.60 s (surrogate_r) keeps the points at a
    # similar place in a ~1 s cardiac cycle, so those averages share the
    # smeared heartbeat-locked signal and their spread underestimates the noise
    CORE::srand($k);
    my ($lo, $hi) = ($s->(1000 * $skip) - $n0, $nt - $n1);
    my @pts;
    while (@pts < $nb) {
        my $p = $lo + int(CORE::rand() * ($hi - $lo));
        push @pts, $p if _off_joins(pdl(long, [$p]))->sclr;
    }
    my $sr = pdl(long, \@pts)->qsort;
    my $Es = zeroes(double, $nh + $nl, $nep);
    _chunks($sr, sub { $Es += $_[0]->mv(2, 0)->sumover }, 1);
    $Es /= $sr->nelem;
    my $C = $Es->slice("$nh:" . ($nh + $nl - 1) . ",:") x $G{nores};   # (nep x nl)(nl x nh) -> (nh, nep)
    push @sres, $Es->slice("0:" . ($nh - 1) . ",:") - $C;
}
my $S    = cat(@sres);                                       # (nh, nep, nsurr)
my $band95 = $S->abs->mv(2, 0)->pctover(0.95);               # (nh, nep)
# channels the fit takes out exactly (BN1, BN2 with --bn-channels --leads reref are the
# regressors themselves): residuals at rounding level, so z and t mean nothing there
my $exact = $S->clump(1, 2)->mv(1, 0)->abs->maximum < 1e-6;   # (nh)

# ---- text --------------------------------------------------------------------------------
my @win = ([200, 300], [300, 450], [450, 600]);
print "\nresidual mean [uV] (z against the non-heartbeat averages): without rho | with rho\n";
printf "  %-4s %s\n", 'ch', join('', map { sprintf '%22s', "$_->[0]..$_->[1] ms" } @win);
for my $c (0 .. $nh - 1) {
    my @cells;
    for my $w (@win) {
        my ($i0, $i1) = ($s->($w->[0]) - $n0, $s->($w->[1]) - $n0 - 1);
        my $sm = $S->slice("($c),$i0:$i1,:")->average;            # (nsurr) window means
        my $sd = sqrt((($sm - $sm->avg)**2)->sum / ($sm->nelem - 1));
        my ($a1, $a2) = map { $res{$_}->slice("($c),$i0:$i1")->avg } qw(nores resp);
        push @cells, $exact->at($c) ? sprintf('%6.2f(%4s) %6.2f(%4s)', $a1, '-', $a2, '-')
                                    : sprintf('%6.2f(%4.1f) %6.2f(%4.1f)', $a1, $a1 / $sd, $a2, $a2 / $sd);
    }
    printf "  %-4s %s\n", $head[$c], join(' ', @cells);
}

# ---- respiration and what is left -------------------------------------------------------
# Per beat, the residual after the cardiac field with respiration,
#   R_k(t) = E_k(t) - sum_j (g_j + h_j rho_k) L_j,k(t),
# its means over the latency windows regressed on rho_k (slope per SD of rho,
# standard error from the beats), and the R-locked averages of R in thirds of rho.
my $gq = $G{resp}->slice(":,0:" . ($nl - 1));                # (nh, nl)
my $hq = $G{resp}->slice(":,$nl:" . (2 * $nl - 1));
my @wi = map { [ $s->($_->[0]) - $n0, $s->($_->[1]) - $n0 - 1 ] } @win;
my (@wm, @rk);                                               # per-beat window means, rho in the same order
my $tert = zeroes(double, $nh, $nep, 3);
my $q = $rho->qsorti;
my $grp = zeroes(long, $nb);
$grp->index($q->slice(int($nb / 3) . ":" . (int(2 * $nb / 3) - 1))) .= 1;
$grp->index($q->slice(int(2 * $nb / 3) . ":-1")) .= 2;
my $ntert = pdl(map { ($grp == $_)->sum } 0 .. 2);
my $bi = 0;
_chunks($rb, sub {
    my ($E, $xs) = @_;
    my $L  = $E->slice("$nh:" . ($nh + $nl - 1) . ",:,:");     # (nl, nep, nb)
    my $Rk = $E->slice("0:" . ($nh - 1) . ",:,:") - ($L x $gq) - ($L x $hq) * $xs->dummy(0)->dummy(0);
    push @wm, cat(map { $Rk->slice(":,$_->[0]:$_->[1],:")->mv(1, 0)->average } @wi);   # (nh, nb, 3)
    my $g = $grp->slice("$bi:" . ($bi + $xs->nelem - 1));
    for my $t (0 .. 2) {
        my $sel = which($g == $t);
        $tert->slice(":,:,($t)") += $Rk->dice_axis(2, $sel)->mv(2, 0)->sumover if $sel->nelem;
    }
    $bi += $xs->nelem;
});
$tert /= $ntert->dummy(0)->dummy(0);
my $WM = $wm[0]->glue(1, @wm[1 .. $#wm]);                    # (nh, nb, nwin)
my $xc = $rho - $rho->avg;
my $sxx = ($xc**2)->sum;
my $bsl = ($WM * $xc->dummy(0))->mv(1, 0)->sumover / $sxx;  # (nh, nwin) slope per SD
my $res_ = $WM - $WM->mv(1, 0)->average->dummy(1) - $bsl->dummy(1) * $xc->dummy(0);
my $se  = sqrt(($res_**2)->mv(1, 0)->sumover / ($nb - 2) / $sxx);
my $tv  = $bsl / $se->setvaltobad(0);                        # a reference channel has se 0: t is bad
$tv = $tv->setbadif($exact->dummy(1, scalar @win));          # and so has a channel fitted exactly
printf "\nwhat is left after the cardiac field (with rho), against rho: slope [uV per SD of rho] (t)\n";
printf "  rho thirds: %d / %d / %d beats\n", $ntert->list;
printf "  %-4s %s\n", 'ch', join('', map { sprintf '%16s', "$_->[0]..$_->[1] ms" } @win);
for my $c (0 .. $nh - 1) {
    printf "  %-4s %s\n", $head[$c], join('', map {
        $tv->slice("($c),($_)")->isbad->sclr ? sprintf('%9.2f(%5s)', $bsl->at($c, $_), '-')
                                              : sprintf('%9.2f(%5.1f)', $bsl->at($c, $_), $tv->at($c, $_))
    } 0 .. $#win);
}
{
    my $tt = $tv->flat->where($tv->flat->isgood);
    printf "  |t| > 2 in %d of %d channel-windows (about %.1f expected by chance; '-' channels left out)\n",
        ($tt->abs > 2)->sum, $tt->nelem, 0.0455 * $tt->nelem;
}

# ---- auxiliary channels (e.g. PPG): R-locked average over the same beats ------------------
my @aux = grep { length } split /,/, $auxlist;
my ($Aux, @auxlab);
if (@aux) {
    exists $ix{$_} or die "--aux: channel '$_' not in $file\n" for @aux;
    my $X = zp_butter($raw->dice_axis(0, pdl(long, [ @ix{@aux} ]))->double,
                      fs => $fs, low => $f_lo, high => $f_hi, order => 4);
    $Aux = zeroes(double, scalar @aux, $nep);
    for (my $i = 0; $i < $nb; $i += $chunk) {
        my $j = ($i + $chunk < $nb ? $i + $chunk : $nb) - 1;
        my $E = beat_epochs($X, $rb->slice("$i:$j"), $n0, $n1)->double;
        $E -= $E->slice(":,$bw[0]:$bw[1]")->mv(1, 0)->average->dummy(1);
        $Aux += $E->mv(2, 0)->sumover;
    }
    $Aux /= $nb;
    my $i0 = $s->(0) - $n0;
    print "\nauxiliary channels (R-locked average, as recorded):\n";
    for my $k (0 .. $#aux) {
        my $a  = $Aux->slice("($k),:");
        my $d  = $a->slice("1:-1") - $a->slice("0:-2");          # slope per sample
        my $ds = $d->slice("$i0:-1");
        my $im = $ds->maximum_ind->sclr + $i0;                  # steepest rise after R
        my $in = $ds->minimum_ind->sclr + $i0;                  # steepest fall after R
        printf "  %-6s peak-to-peak %.3g; steepest rise at %.0f ms, steepest fall at %.0f ms, maximum at %.0f ms, minimum at %.0f ms\n",
            $aux[$k], ($a->max - $a->min)->sclr, $tms->at($im), $tms->at($in),
            $tms->at($a->slice("$i0:-1")->maximum_ind->sclr + $i0), $tms->at($a->slice("$i0:-1")->minimum_ind->sclr + $i0);
    }
    # scaled copies for drawing with the uV channels: centred, peak-to-peak 8
    my $pp = ($Aux->mv(1, 0)->maximum - $Aux->mv(1, 0)->minimum)->clip(1e-30, undef);
    $Aux = ($Aux - $Aux->mv(1, 0)->average->dummy(1)) * (8 / $pp)->dummy(1);   # mean over time, per channel
    @auxlab = map { "$_*" } @aux;
}

# ---- pulse-locked or R-locked? ----------------------------------------------------------------
# A wave that follows the arrival of the pulse at the head moves, beat by beat,
# with the pulse transit time; a response locked to the heartbeat at a fixed
# latency after R does not.  Three readings, every head channel:
#  1. the residual averaged aligned to R and aligned to the pulse upstroke
#     (beats shifted by PTT_k - mean PTT): power of the average in 150..550 ms;
#     jitter smears a wave, so a pulse-locked wave keeps more power when
#     aligned to the pulse (ratio > 1), an R-locked one when aligned to R
#     (ratio < 1); noise gives ratios scattered around 1.  Controls: the
#     cardiac lead (R-locked) and the pulse channel itself (pulse-locked).
#  2. beats in thirds of PTT, residual aligned to R: the lag (+-100 ms) that best
#     matches the average of the longest-PTT third to that of the shortest, in
#     150..550 ms; a pulse-locked wave is delayed by about the PTT difference
#     (lag / PTT difference near 1), an R-locked one is not (near 0).
#  3. all beats at once: the beat-to-beat deviations from the R-aligned average,
#     regressed on -d_k m'(t) (d_k = PTT_k - mean PTT, m' the slope of that
#     average), alone and jointly with -e_k m'(t) (e_k = the preceding RR minus
#     its mean), since PTT and the timing of RR-dependent R-locked waves both
#     change with the heart rate; beta relative to that of the pulse channel:
#     near 1 if the wave moves with the pulse, near 0 if it stays at its latency
#     after R; 95% intervals from resampling the beats.
my $Apulse;
if (defined $pulsech) {
    exists $ix{$pulsech} or die "--pulse: channel '$pulsech' not in $file\n";
    my $P  = zp_butter($raw->slice("($ix{$pulsech})")->double, fs => $fs, low => $f_lo, high => 15, order => 4);
    my $Ep = beat_epochs($P->dummy(0), $rb, $n0, $n1)->double->clump(0, 1);   # (nep, nb)
    my ($w0, $w1) = ($s->(100) - $n0, $s->(450) - $n0);
    my $sm = sub { my $x = shift; my $k = _rnd2(0.010 * $fs); ($x->slice("$k:-1") - $x->slice("0:-" . ($k + 1))) / $k };  # slope over 10 ms
    my $dA = $sm->($Ep->mv(1, 0)->average)->slice("$w0:$w1");
    my $sg = $dA->max->sclr >= -$dA->min->sclr ? 1 : -1;      # upstroke = steepest rise of sg*CH
    my $dE = $sm->($Ep * $sg)->slice("$w0:$w1");               # (nw, nb)
    my $up0 = ($dE->maximum_ind + $w0)->double;                # first guess: steepest rise, sample in the epoch
    # refine every beat by matching a template of the upstroke (-80..+80 ms around it):
    # correlation over lags of +-40 ms around the guess, parabolic peak (sub-sample),
    # template rebuilt from the refined positions, three rounds
    my ($hw, $lg) = (_rnd2(0.080 * $fs), _rnd2(0.040 * $fs));
    my $js  = sequence(long, 2 * $hw + 1) - $hw;               # (ns) offsets in the template
    my $Ps  = $P * $sg;
    my $seg = sub {                                            # (ns, nb) segments centred at R + n0 + pos
        my $c = shift;                                         # (nb) long, sample in the record
        return $Ps->index(($js->dummy(1) + $c->dummy(0))->clip(0, $nt - 1));
    };
    my $dm  = sub { my $x = shift; my $y = $x - $x->average->dummy(0); $y / sqrt(($y**2)->sumover)->clip(1e-30, undef)->dummy(0) };
    my $base = $rb + $n0;                                      # record sample of epoch sample 0
    my ($up, $rmax) = ($up0, undef);
    for my $round (1 .. 3) {
        my $T  = $dm->($seg->(($base + $up)->rint->long)->mv(1, 0)->average->dummy(1))->flat;   # (ns) template
        my $c0 = ($base + $up)->rint->long;
        my $R  = cat(map { ($dm->($seg->($c0 + $_)) * $T->dummy(1))->sumover } -$lg .. $lg);   # (nb, nlag)
        my $im = $R->mv(1, 0)->maximum_ind->clip(1, 2 * $lg - 1);                            # (nb)
        my ($rm, $r0, $rp) = map { $R->index2d(sequence(long, $nb), $im + $_) } -1, 0, 1;
        my $den = $rm - 2 * $r0 + $rp;
        my $dl  = ($den < 0) * 0.5 * ($rm - $rp) / ($den + ($den >= 0));                        # parabolic, 0 if not a peak
        $up   = ($c0 - $base)->double + ($im - $lg)->double + $dl->clip(-1, 1);
        $rmax = $r0;
    }
    my $okb = $rmax > 0.8;                                     # beats whose upstroke matches the template
    my $ptt = ($up + $n0) * (1000 / $fs);                      # ms after R, sub-sample
    my $ptg = $ptt->where($okb);
    my $p0g = (($up0 + $n0) * (1000 / $fs))->where($okb);
    printf "\npulse channel %s%s: upstroke matched (r > 0.8) in %d of %d beats; PTT (R to the upstroke) mean %.0f ms, SD %.1f ms, 5..95%% %.0f..%.0f ms\n",
        $pulsech, $sg < 0 ? ' (inverted)' : '', $okb->sum, $nb, $ptg->avg, sqrt((($ptg - $ptg->avg)**2)->avg), $ptg->pct(0.05), $ptg->pct(0.95);
    printf "  PTT by the steepest rise alone: SD %.1f ms; correlation of the two estimates %.2f\n",
        sqrt((($p0g - $p0g->avg)**2)->avg), _corr2($p0g, $ptg);
    my $ib = which($okb);
    my $rbp = $rb->index($ib);
    my $shift = ($ptt->index($ib) - $ptg->avg) * ($fs / 1000);
    $shift = $shift->rint->long;                               # nearest sample
    my $gq0 = $G{nores};                                       # (nh, nl)
    # residual epochs of the head (cardiac field removed) and the pulse channel at the
    # samples $r + $sh, baselined at -75..-50 ms of the R-aligned epoch
    my $resid = sub {
        my ($rp, $sh) = @_;
        my $Rs = zeroes(double, $nh + 2, $nep);                # head channels, the lead, the pulse
        for (my $i = 0; $i < $rp->nelem; $i += $chunk) {
            my $j  = ($i + $chunk < $rp->nelem ? $i + $chunk : $rp->nelem) - 1;
            my $rr = $rp->slice("$i:$j");
            my $E0 = beat_epochs($Y, $rr, $n0, $n1)->double;
            my $E  = beat_epochs($Y, $rr + $sh->slice("$i:$j"), $n0, $n1)->double;
            my $bl = $E0->slice(":,$bw[0]:$bw[1]")->mv(1, 0)->average->dummy(1);
            $E -= $bl;
            my $L  = $E->slice("$nh:" . ($nh + $nl - 1) . ",:,:");
            my $H  = $E->slice("0:" . ($nh - 1) . ",:,:") - ($L x $gq0);
            my $Pe = beat_epochs($P->dummy(0), $rr + $sh->slice("$i:$j"), $n0, $n1)->double * $sg;
            $Rs += $H->glue(0, $L->slice("(0)")->dummy(0), $Pe)->mv(2, 0)->sumover;
        }
        return $Rs / $rp->nelem;
    };
    my $okshift = ($rbp + $shift + $n0 >= 0) & ($rbp + $shift + $n1 <= $nt);
    $rbp   = $rbp->where($okshift);
    $shift = $shift->where($okshift);
    my $AR = $resid->($rbp, zeroes(long, $rbp->nelem));        # aligned to R
    my $AP = $resid->($rbp, $shift);                           # aligned to the pulse
    my @lab = (@head, "$leaddesc (lead)", "$pulsech (pulse)");
    my ($s0, $s1) = ($s->(150) - $n0, $s->(550) - $n0 - 1);
    # power in 150..550 ms of the average (its window mean removed), per channel
    my $pw  = sub { my $W = shift->slice(":,$s0:$s1"); (($W - $W->mv(1, 0)->average->dummy(1))**2)->mv(1, 0)->average };
    my ($pR, $pP) = map { $pw->($_) } $AR, $AP;
    # noise: the same power in the non-heartbeat averages (median over them), head channels only;
    # a channel's numbers mean something only where the heartbeat average stands well above it
    my $SW = $S->slice(":,$s0:$s1,:");
    my $pN = (($SW - $SW->mv(1, 0)->average->dummy(1))**2)->mv(1, 0)->average->mv(1, 0)->medover;   # (nh)
    print "  1. power of the average in 150..550 ms [uV^2]: aligned to R | aligned to the pulse | ratio pulse/R | R-aligned power / non-heartbeat power\n";
    for my $c (0 .. $#lab) {
        next if $pR->at($c) == 0;
        printf "     %-18s %10.3f %10.3f   %5.2f   %s\n", $lab[$c], $pR->at($c), $pP->at($c), $pP->at($c) / $pR->at($c),
            $c < $nh ? sprintf('%6.1f%s', $pR->at($c) / $pN->at($c), $pR->at($c) / $pN->at($c) >= 4 ? '' : '  (near noise: ignore)') : '     -';
    }
    # 2. thirds of PTT: lag of the latest third against the earliest, by cross-correlation
    my $pq = $ptt->index($ib)->where($okshift);
    my $qo = $pq->qsorti;
    my $nq = $pq->nelem;
    my (@pm, @Ag);
    for my $g (0 .. 2) {
        my $sel = $qo->slice(int($g * $nq / 3) . ":" . (int(($g + 1) * $nq / 3) - 1));
        push @Ag, $resid->($rbp->index($sel), zeroes(long, $sel->nelem));
        push @pm, $pq->index($sel)->avg;
    }
    my $mx = _rnd2(0.100 * $fs);                               # lags up to +-100 ms
    printf "  2. thirds of PTT (mean %.0f / %.0f / %.0f ms; PTT difference third 3 - third 1: %.0f ms):\n", @pm, $pm[2] - $pm[0];
    print "     lag [ms] of third 3 against third 1 that best matches them in 150..550 ms (correlation), R-aligned\n";
    for my $c (0 .. $#lab) {
        my ($a1, $a3) = map { $_->slice("($c),:") } $Ag[0], $Ag[2];
        my $x = $a1->slice("$s0:$s1"); $x = $x - $x->avg;
        my ($best, $bl) = (-2, 0);
        for my $l (-$mx .. $mx) {
            my $y = $a3->slice(($s0 + $l) . ":" . ($s1 + $l)); $y = $y - $y->avg;
            my $den = sqrt(($x**2)->sum * ($y**2)->sum)->sclr;
            next unless $den > 0;
            my $rr = ($x * $y)->sum->sclr / $den;
            ($best, $bl) = ($rr, $l) if $rr > $best;
        }
        next if $best < -1;
        printf "     %-18s lag %+5.0f ms (r %.2f); lag / PTT difference %5.2f%s\n", $lab[$c], $bl * 1000 / $fs, $best,
            ($bl * 1000 / $fs) / ($pm[2] - $pm[0]),
            ($c < $nh && $pR->at($c) / $pN->at($c) < 4) ? '  (near noise: ignore)' : '';
    }
    # 3. every beat at once: if a wave m(t) is locked to the pulse, the beat whose pulse
    # comes d_k later carries m(t - d_k) ~ m(t) - d_k m'(t).  With a_k = -sum_t (y_k - m) m'
    # over 150..550 ms and M = sum_t m'^2, a_k ~ M beta d_k.  PTT changes with the heart
    # rate, and so does the timing of R-locked waves that depend on the RR interval (the
    # T wave, for one), so beta alone mixes the two; the joint fit
    #     a_k ~ M (beta d_k + gamma e_k),  e_k = RR before beat k minus its mean (samples)
    # gives beta with the RR-related timing taken out.  95% intervals: beats resampled.
    {
        my %pos;
        my @ra = $rall->list;
        @pos{@ra} = 0 .. $#ra;
        my $rrp = pdl(double, [ map { my $q = $pos{$_}; defined $q && $q > 0 ? $_ - $ra[$q - 1] : 'NaN' } $rbp->list ]);
        my $dk  = ($pq - $pq->avg) * ($fs / 1000);              # d_k in samples, beats as in $rbp
        my $okr = $rrp->isfinite;
        my $ek  = ($rrp - $rrp->where($okr)->avg)->setnantobad->setbadtoval(0);   # RR deviation, samples
        my $k   = _rnd2(0.010 * $fs);
        my $Mr  = $AR;                                         # (nh+2, nep) R-aligned mean
        my $mp  = ($Mr->slice(":,$k:-1") - $Mr->slice(":,0:-" . ($k + 1))) / $k;   # slope, centred k/2 later
        my $hk  = int($k / 2);
        my $mpw = $mp->slice(":," . ($s0 - $hk) . ":" . ($s1 - $hk));            # (nh+2, nw)
        my $Mm  = ($mpw**2)->mv(1, 0)->sumover;               # (nh+2)
        my @ak;
        for (my $i = 0; $i < $rbp->nelem; $i += $chunk) {
            my $j  = ($i + $chunk < $rbp->nelem ? $i + $chunk : $rbp->nelem) - 1;
            my $rr = $rbp->slice("$i:$j");
            my $E  = beat_epochs($Y, $rr, $n0, $n1)->double;
            $E    -= $E->slice(":,$bw[0]:$bw[1]")->mv(1, 0)->average->dummy(1);
            my $L  = $E->slice("$nh:" . ($nh + $nl - 1) . ",:,:");
            my $H  = $E->slice("0:" . ($nh - 1) . ",:,:") - ($L x $gq0);
            my $Pe = beat_epochs($P->dummy(0), $rr, $n0, $n1)->double * $sg;
            my $Yk = $H->glue(0, $L->slice("(0)")->dummy(0), $Pe)->slice(":,$s0:$s1,:") - $Mr->slice(":,$s0:$s1")->dummy(2);
            push @ak, -($Yk * $mpw->dummy(2))->mv(1, 0)->sumover;   # (nh+2, nb)
        }
        my $A   = $ak[0]->glue(1, @ak[1 .. $#ak]);             # (nh+2, nb)
        my $sel = which($okr);
        ($A, $dk, $ek) = ($A->dice_axis(1, $sel), $dk->index($sel), $ek->index($sel));
        my $nbb = $A->dim(1);
        printf "  3. all beats (%d with a preceding RR): correlation of PTT and RR %.2f; RR SD %.1f ms\n",
            $nbb, _corr2($dk, $ek), sqrt(($ek**2)->avg) * 1000 / $fs;
        # (beta alone, beta and gamma jointly) from sums over a set of beats
        my $fit = sub {
            my ($Ab, $d, $e) = @_;
            my ($Sdd, $See, $Sde) = map { $_->sum } $d**2, $e**2, $d * $e;
            my $Sad = ($Ab * $d->dummy(0))->mv(1, 0)->sumover;
            my $Sae = ($Ab * $e->dummy(0))->mv(1, 0)->sumover;
            my $det = $Sdd * $See - $Sde**2;
            return ($Sad / $Sdd / $Mm, ($Sad * $See - $Sae * $Sde) / $det / $Mm, ($Sae * $Sdd - $Sad * $Sde) / $det / $Mm);
        };
        my ($b1, $bj, $gj) = $fit->($A, $dk, $ek);
        CORE::srand(11);
        my (@B1, @BJ, @GJ);
        for (1 .. 1000) {
            my $ix = pdl(long, [ map { int(CORE::rand() * $nbb) } 1 .. $nbb ]);
            my @f = $fit->($A->dice_axis(1, $ix), $dk->index($ix), $ek->index($ix));
            push @B1, $f[0]; push @BJ, $f[1]; push @GJ, $f[2];
        }
        my $ip = $#lab;
        my $ci = sub { my ($X, $rel) = @_; my $Z = cat(@$X); $Z = $Z / $Z->slice("($ip),:")->dummy(0) if $rel;
                       map { $Z->mv(1, 0)->pctover($_) } 0.025, 0.975 };
        my ($r1lo, $r1hi) = $ci->(\@B1, 1);
        my ($rjlo, $rjhi) = $ci->(\@BJ, 1);
        my ($gjlo, $gjhi) = $ci->(\@GJ, 0);
        print "     beta relative to that of the pulse channel (near 1: moves with the pulse; near 0: stays at its latency after R),\n";
        print "     alone and with the RR interval in the fit; gamma: shift per sample of RR (95% intervals, beats resampled)\n";
        printf "     %-18s %-24s %-24s %s\n", '', 'beta / pulse, alone', 'beta / pulse, with RR', 'gamma (RR)';
        for my $c (0 .. $#lab) {
            next if $pR->at($c) == 0;
            printf "     %-18s %5.2f (%5.2f .. %5.2f)   %5.2f (%5.2f .. %5.2f)   %5.2f (%5.2f .. %5.2f)%s\n", $lab[$c],
                $b1->at($c) / $b1->at($ip), $r1lo->at($c), $r1hi->at($c),
                $bj->at($c) / $bj->at($ip), $rjlo->at($c), $rjhi->at($c),
                $gj->at($c), $gjlo->at($c), $gjhi->at($c),
                ($c < $nh && $pR->at($c) / $pN->at($c) < 4) ? '  (near noise: ignore)' : '';
        }
        printf "     (beta of the pulse channel itself: alone %.2f, with RR %.2f)\n", $b1->at($ip), $bj->at($ip);
    }
    $Apulse = $AP->slice("0:" . ($nh - 1) . ",:");
}

# ---- averages for examples/avg_grid_scrub.pl (5x5 10-20 grid with a MAP2D map) --------------
for my $tv (['raw', $Ebar], ['res', $res{nores}], ['resrho', $res{resp}],
            (map { [ "rho$_", $tert->slice(":,:,(" . ($_ - 1) . ")") ] } 1 .. 3),
            (defined $Apulse ? ([ 'pulsealign', $Apulse ]) : ())) {
    my ($tag, $M) = @$tv;
    my $fn = "${out}_avg_$tag.tsv";
    open my $fh, '>', $fn or die "$fn: $!";
    print $fh join("\t", 'time_ms', @head, @auxlab), "\n";
    $M = $M->glue(0, $Aux) if defined $Aux;                  # the same scaled aux columns in every file
    my @rows = $M->dog;                                     # one (nh + naux) piddle per sample
    for my $i (0 .. $nep - 1) {
        print $fh join("\t", sprintf('%.1f', $tms->at($i)), map { sprintf '%.4f', $_ } $rows[$i]->list), "\n";
    }
    close $fh;
    print "wrote $fn\n";
}

if (eval { require PDL::Graphics::Cairo; 1 }) {
    # joint plots (waveforms + maps) of the scalp channels that have a position
    # in the 10-20 montage (not A1, A2, BN1, BN2, nor channels such as TP9), MAP2D
    my %elc  = (T3 => 'T7', T4 => 'T8', T5 => 'P7', T6 => 'P8');
    my $mon  = PDL::EEG::IO::ASA::read_elc("$FindBin::Bin/../t/data/standard_1020_subset.elc");
    my @sc   = grep { my $n = $elc{ $head[$_] } // $head[$_];
                      $head[$_] !~ /^A[12]$/ && $head[$_] ne $bn1 && $head[$_] ne $bn2 && exists $mon->{pos}{$n} } 0 .. $nh - 1;
    my $sci  = pdl(long, \@sc);
    my @tlab = map { $elc{ $head[$_] } // $head[$_] } @sc;
    my @mt   = split /,/, $mtimes;
    my ($x0, $x1) = split /,/, $show;
    for my $jp (['raw', $Ebar, 'average as recorded'],
                ['res', $res{nores}, 'minus cardiac field'],
                ['resrho', $res{resp}, 'minus cardiac field with rho']) {
        my ($tag, $M, $what) = @$jp;
        plot_joint(avg => $M->dice_axis(0, $sci), t => $tms, labels => \@tlab,
                   montage => "$FindBin::Bin/../t/data/standard_1020_subset.elc",
                   times => \@mt, xlim => [ $x0, $x1 ], ylim => [ -$ylim, $ylim ],
                   negative_up => ($posup ? 0 : 1), nave => $nb,
                   title => sprintf('%s: %s (%s; leads %s, gains fitted %s ms)', $file =~ s{.*/}{}r, $what, $refdesc, $leaddesc, $fitw),
                   xlabel => 'time from R peak [ms]', outfile => "${out}_joint_$tag.png");
        print "wrote ${out}_joint_$tag.png\n";
    }
}

# ---------------------------------------------------------------------------
# epochs (nh+nl, nep, nb) of the head and the leads, baselined, in chunks;
# the callback gets (epochs, rho of those beats)
sub _chunks {
    my ($rp, $cb, $norho) = @_;
    for (my $i = 0; $i < $rp->nelem; $i += $chunk) {
        my $j = ($i + $chunk < $rp->nelem ? $i + $chunk : $rp->nelem) - 1;
        my $E = beat_epochs($Y, $rp->slice("$i:$j"), $n0, $n1)->double;
        $E -= $E->slice(":,$bw[0]:$bw[1]")->mv(1, 0)->average->dummy(1);
        $cb->($E, $norho ? undef : $rho->slice("$i:$j"));
    }
}
# (keep flag per beat, number of beats over the limit per rejection channel)
sub _accept {
    my $rp   = shift;
    my $keep = zeroes(byte, $rp->nelem);
    my $nbad = zeroes(long, $nrej);
    for (my $i = 0; $i < $rp->nelem; $i += $chunk) {
        my $j  = ($i + $chunk < $rp->nelem ? $i + $chunk : $rp->nelem) - 1;
        my $Em = beat_epochs($Y, $rp->slice("$i:$j"), $n0, $n1)->slice("0:" . ($nrej - 1) . ",:,:")->mv(1, 0);
        my $over = ($Em->maximum - $Em->minimum) > $reject;          # (nrej, nb chunk)
        $keep->slice("$i:$j") .= !$over->orover;
        $nbad += $over->mv(1, 0)->sumover;
    }
    return ($keep, $nbad);
}

sub _rnd2 { my $v = shift; $v >= 0 ? int($v + 0.5) : -int(-$v + 0.5) }
sub _corr2 { my ($x, $y) = @_; $x = $x - $x->avg; $y = $y - $y->avg; (($x * $y)->sum / sqrt(($x**2)->sum * ($y**2)->sum))->sclr }

# 1 where the epoch (n0..n1 around each sample) stays 2 s away from every block join
sub _off_joins {
    my $p  = shift;
    my $ok = ones(byte, $p->nelem);
    my $m  = 2 * $fs;
    for my $j (@joins) {
        $ok &= !(($p + $n1 + $m > $j) & ($p + $n0 - $m < $j));
    }
    return $ok;
}
