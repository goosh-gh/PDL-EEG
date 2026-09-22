#!/usr/bin/env perl
# examples/raw_plot.pl
#
# A raw.plot()-style span annotator for EEG recordings (EDF or Nihon Kohden),
# built on the
# same P:G:C giza-server interactive viewer as read_nihonkohden.pl (which this
# references for its trace layout and gain/scroll sliders). It adds MNE-like
# mouse span selection so you can mark blink / horizontal-saccade / bad segments
# by eye and save them; the spans then feed GED (S from marked segments, R from
# clean) or ICA segment selection via PDL::EEG::Spans.
#
# The viewer already forwards mouse events (P:G:C Driver::GS on_pick/on_cursor,
# = giza-server GSP_MSG_PICK/CURSOR), and Axes->image_frac_to_data() converts the
# click's image fraction to the plot's data-coordinate seconds exactly, so no
# giza-server/protocol change is needed -- this is pure Perl on top of the
# existing viewer.
#
# MOUSE (no keyboard events are delivered to the window, so everything is mouse):
#   left click   : place a span endpoint. First left = start, second = end ->
#                  the span is committed with --kind and auto-saved.
#   right click  : delete the span under the cursor (or cancel a pending start).
#   (the two sliders keep their read_nihonkohden.pl meaning: time scroll / gain)
#
# The span kind is set per run by --kind (default hsaccade) and is free-form:
# any name works (blink, hsaccade, bad, or any GED class you later declare with
# artifact_remove.pl --ged NAME, e.g. slowA, vert_EOG, sweat_like). To label more
# than one kind, re-run with a different --kind: the new spans are appended to the same
# sidecar (loaded on start), so you mark all hsaccade, then all clean, etc.
#
# Spans are written to a sidecar TSV (--out, default <file>.spans.tsv) on every
# change, and re-loaded on start (--load, default = --out if it exists).
#
# Usage:
#   perl -Ilib -I/path/to/PDL-EEG-ICA/lib \
#        examples/raw_plot.pl FILE.edf|subject.eeg \
#        [--block N] [--sec S] [--nch N] [--uv U] [--chans LIST|all] [--aux MODE] \
#   --chans all : show every channel.  --aux MODE: same | auto | <uV/div>
#                 (default 3000000 = 3 V/div, so 2-3 V DC triggers fit one
#                  division; EEG/EOG/ECG are never rescaled)
#        [--kind KIND] [--out FILE] [--load FILE]
#
# Then, to build a horizontal GED filter from what you marked:
#   S_horiz = covariance of samples in the 'hsaccade' mask (vEOG quiet),
#   R       = covariance of 'clean' samples,   via PDL::EEG::Spans::spans_to_mask.

use strict;
use warnings;
use PDL;
use PDL::EEG::IO::NihonKohden qw(read_nk);
use PDL::EEG::IO::EDF         qw(read_edf clean_edf_label);
use PDL::Graphics::Cairo qw(subplots);
use PDL::Graphics::Cairo::Driver::GS;
use PDL::EEG::Spans qw(read_spans write_spans);
use Getopt::Long;

binmode(STDOUT, ':encoding(UTF-8)');
binmode(STDERR, ':encoding(UTF-8)');

my ($blk, $chans, $out, $load);
my $nsec = 10; my $nch = 8; my $uv_init = 100; my $aux_mode = '3000000';   # aux (non-EEG/EOG/ECG) default: 3 V per division
my $cur_kind = 'hsaccade';
GetOptions(
    'block=i'  => \$blk,   'sec=f'   => \$nsec,   'nch=i'  => \$nch,
    'uv=f'     => \$uv_init,'chans=s' => \$chans,  'aux=s'  => \$aux_mode,
    'kind=s'   => \$cur_kind,'out=s'  => \$out,    'load=s' => \$load,
) or die "bad options\n";
my $file = shift @ARGV or die "Usage: $0 FILE.edf|subject.eeg [options]  (--chans all shows every channel)\n";
$out  //= "$file.spans.tsv";
$load //= (-e $out ? $out : undef);

# ---- EOG label variants + aux classifier -----------------------------------
# EOG naming varies: lhEOG/hEOG = horizontal; lvEOG/vEOG/rvEOG = vertical.
# X1 is commonly wired at the rvEOG position, so treat X1 as vertical EOG too.
sub is_eog { my $l = shift // ''; $l =~ /eog/i || $l =~ /^X1$/i }
sub eog_kind {
    my $l = shift // '';
    return 'h' if $l =~ /h\s*eog/i;                  # lhEOG, hEOG
    return 'v' if $l =~ /v\s*eog/i || $l =~ /^X1$/i; # lvEOG, vEOG, rvEOG, X1
    return '?';
}
# 'aux' = anything that should get the fixed aux gain. EOG and ECG are
# physiological and are NEVER aux (they stay at the EEG scale).
sub is_aux {
    my $l = shift // '';
    return 0 if is_eog($l) || $l =~ /^(?:e?k?cg|ECG|EKG)$/i;
    return $l =~ /^(?:DC\d|STIM|PAD|COM|BN\d|Pulse|CO2|SpO2|Mark|Events|\$)/i
        || $l =~ /EtCO2|Trigger|Photo|Resp/i;
}

# ---- load: EDF (.edf) via read_edf, else Nihon Kohden raw via read_nk --------
my ($rec, $abs0, @labels, $device);
$abs0 = 0;
if ($file =~ /\.edf$/i) {
    warn "--block ignored for EDF (EDF has no Nihon Kohden blocks)\n" if defined $blk;
    $rec    = read_edf($file);
    @labels = map { clean_edf_label($_) } @{ $rec->{labels} };   # "EEG Fp1-Ref" -> "Fp1"
    $device = 'EDF' . (defined $rec->{edf_type} ? " ($rec->{edf_type})" : '');
}
else {
    $rec    = defined $blk ? read_nk($file, block => $blk)
                           : read_nk($file, all_blocks => 1);
    @labels = @{ $rec->{labels} };
    $device = $rec->{device} // 'NK';
}
my $data       = $rec->{data};                       # (n_ch, n_samp) uV
my $fs         = $rec->{fs};
my $n_ch_valid = $rec->{n_ch_valid} // scalar @labels;
my $n_total    = $data->dim(1);
printf "Loaded : %s  %d ch, %d samp = %.1f s @ %g Hz\n",
    $device, $data->dim(0), $n_total, $n_total/$fs, $fs;

# ---- channel selection ------------------------------------------------------
my @sel;
if (length $chans) {
    if (lc $chans eq 'all') { @sel = 0 .. $#labels; }
    else {
    my %idx; $idx{$labels[$_]} = $_ for 0..$#labels;
    for my $name (split /\s*,\s*/, $chans) {
        defined $idx{$name} or die "no such channel: $name\n";
        push @sel, $idx{$name};
    }
    }
} else {
    my $np = $nch < $n_ch_valid ? $nch : $n_ch_valid;
    @sel = (0 .. $np - 1);
}
my $n_plot     = scalar @sel;
my @sel_labels = @labels[@sel];
printf "Channels: %s\n", join(' ', @sel_labels);

# ---- span state (shared between render and mouse callbacks) ------------------
my @spans;
@spans = @{ read_spans($load) } if $load && -e $load;
printf "Loaded %d existing spans from %s\n", scalar @spans, $load if $load;
my $pending;                       # data-coord seconds of a placed start, or undef
my $cursor_t;                      # last cursor time (for live preview line)
my $shared_ax;                     # most recently drawn Axes (for image_frac_to_data)
my $ylo = 0; my $yhi = 1;          # current y-range (set in render)

# --kind is free-form: any name may be marked (blink, hsaccade, bad, and any
# GED artifact class you declare later with artifact_remove.pl --ged NAME, e.g.
# slowA, vert_EOG, sweat_like). 'bad' has a fixed meaning (excluded at
# averaging, never GED-filtered). Known kinds get fixed colours; anything else
# is assigned a stable colour from a palette so it is distinguishable on screen.
$cur_kind =~ /^\S+$/ or die "--kind must be a single non-empty token\n";
my %KCOL  = (hsaccade=>'#00BCD4', blink=>'#FF9800', bad=>'#F44336', clean=>'#4CAF50');
my @PALETTE = ('#9C27B0','#3F51B5','#009688','#795548','#E91E63','#607D8B','#8BC34A');
sub kind_color {
    my ($k) = @_; return $KCOL{$k} if exists $KCOL{$k};
    my $h = 0; $h = ($h*31 + ord $_) % scalar(@PALETTE) for split //, $k;   # stable
    return $PALETTE[$h];
}
$KCOL{$cur_kind} //= kind_color($cur_kind);

sub save_spans {
    write_spans($out, \@spans);
    printf STDERR "spans: %d saved to %s\n", scalar @spans, $out;
}
sub span_at {                      # index of span whose [t0,t1] contains $t, else undef
    my $t = shift;
    for my $i (0..$#spans) { return $i if $t >= $spans[$i]{t0} && $t <= $spans[$i]{t1} }
    return undef;
}

# ---- viewer -----------------------------------------------------------------
my $gs = PDL::Graphics::Cairo::Driver::GS->new(
    width  => 1100,
    height => 80 + $n_plot * 70,
);
my $sv_init = 1.0 - log(($uv_init > 0 ? $uv_init : 100) / 10) / log(50);
$sv_init = 0 if $sv_init < 0; $sv_init = 1 if $sv_init > 1;

$gs->show_interactive(
    init           => { 0 => 0.0, 1 => $sv_init },
    cursor_overlay => 1,           # makes clicks/moves trigger a redraw
    on_cursor => sub {
        my ($fx, $fy) = @_;
        return unless $shared_ax;
        my ($xd) = $shared_ax->image_frac_to_data($fx, $fy);
        $cursor_t = $xd;           # may be undef outside the plot frame
    },
    on_pick => sub {
        my ($fx, $fy, $btn) = @_;
        return unless $shared_ax;
        my ($xd) = $shared_ax->image_frac_to_data($fx, $fy);
        return unless defined $xd;                                    # ignore clicks off-plot
        if ($btn == 3) {                                             # right: delete / cancel
            if (defined $pending) { $pending = undef; return; }
            my $i = span_at($xd);
            if (defined $i) { splice @spans, $i, 1; save_spans(); }
            return;
        }
        # left: place endpoint
        if (!defined $pending) { $pending = $xd; return; }
        my ($a, $b) = sort { $a <=> $b } ($pending, $xd);
        push @spans, { kind => $cur_kind, t0 => $a, t1 => $b };
        $pending = undef;
        save_spans();
    },
    render => sub {
        my ($state, $w, $h) = @_;

        # slider 0: time offset; slider 1: gain (µV/div) -- as read_nihonkohden.pl
        my $n_show = int($nsec * $fs); $n_show = 1 if $n_show < 1;
        $n_show    = $n_total if $n_show > $n_total;
        my $t_off  = ($state->{0} // 0.0) * (($n_total - $n_show) / $fs);
        my $s0     = int($t_off * $fs); $s0 = $n_total - $n_show if $s0 + $n_show > $n_total;
        $s0        = 0 if $s0 < 0;
        my $s1     = $s0 + $n_show - 1; $s1 = $n_total - 1 if $s1 >= $n_total;

        my $sv         = $state->{1} // 0.5;
        my $uv_per_div = 10 * (50 ** (1.0 - $sv));
        my $spacing    = $uv_per_div * 2.0;

        my ($fig, $ax) = subplots(1, 1, width => $w, height => $h);
        $shared_ax = $ax;                          # for image_frac_to_data in callbacks

        my $plot_w = $w - 120; $plot_w = 200 if $plot_w < 200;
        my $n_raw  = $s1 - $s0 + 1;
        my $step   = int($n_raw / $plot_w) || 1; $step = 5 if $step > 5;
        my $t;
        for my $i (0 .. $n_plot - 1) {
            my $ci  = $sel[$i];
            my $raw = $data->slice("($ci),$s0:$s1");
            my $sig = ($step > 1) ? $raw->slice("0:-1:$step") : $raw;
            $t //= (sequence($sig->dim(0)) * $step + $s0 + $abs0) / $fs;
            my $aux    = is_aux($sel_labels[$i]);
            my $offset = ($n_plot - 1 - $i) * $spacing;
            my $y;
            if (!$aux || $aux_mode eq 'same') {
                $y = $sig + $offset;                       # physiological scale (raw uV)
            } elsif ($aux_mode eq 'auto') {
                my $mean = $sig->avg;
                my $amp  = ($sig - $mean)->abs->max; $amp = 1 if $amp <= 0;
                $y = ($sig - $mean) / $amp * ($spacing * 0.4) + $offset;
            } else {                                       # fixed N uV per division (default 3 V)
                my $N = ($aux_mode + 0) || 1;
                $y = ($sig - $sig->avg) * ($uv_per_div / $N) + $offset;
            }
            $ax->line($t, $y, color => ($aux ? '#CC0000' : '#1565C0'), lw => 0.7);
            $ax->axhline($offset, color => '#DDDDDD', lw => 0.3);
        }

        my $x_lo   = ($s0 + $abs0) / $fs;
        my $x_hi   = $t->max;
        my $gutter = ($x_hi - $x_lo) * 0.05 + 1e-9;
        $ylo = -$spacing * 1.0;
        $yhi = ($n_plot - 0.2) * $spacing;

        # committed spans as shaded full-height bands
        for my $s (@spans) {
            next if $s->{t1} < $x_lo || $s->{t0} > $x_hi + $gutter;
            $ax->axvspan($s->{t0}, $s->{t1},
                         color => (kind_color($s->{kind})), alpha => 0.25);
        }
        # pending start (solid) and live cursor (thin) markers
        if (defined $pending && $pending >= $x_lo && $pending <= $x_hi + $gutter) {
            $ax->line(PDL->new([$pending, $pending]), PDL->new([$ylo, $yhi]),
                      color => (kind_color($cur_kind)), lw => 1.6);
        }
        if (defined $cursor_t && $cursor_t >= $x_lo && $cursor_t <= $x_hi + $gutter) {
            $ax->line(PDL->new([$cursor_t, $cursor_t]), PDL->new([$ylo, $yhi]),
                      color => '#607D8B', lw => 0.6);
        }

        # EEG calibration bar in the right gutter
        my $bar_t = $x_hi + $gutter * 0.5;
        $ax->line(PDL->new([$bar_t, $bar_t]), PDL->new([0, $uv_per_div]),
                  color => '#CC0000', lw => 2.0);

        my @tick_vals = map { ($n_plot - 1 - $_) * $spacing } 0 .. $n_plot - 1;
        $ax->yticks(PDL->new(\@tick_vals), [@sel_labels]);
        $ax->xlim($x_lo, $x_hi + $gutter);
        $ax->ylim($ylo, $yhi);
        $ax->xlabel(sprintf("Time (s)  [%.1f-%.1fs step=%d]", $x_lo, $x_hi, $step));

        my %count; $count{$_->{kind}}++ for @spans;
        my %seen = (%count, $cur_kind => ($count{$cur_kind}//0));
        my @show = sort keys %seen;
        $fig->suptitle(sprintf(
            "%s  %gHz  aux=%s  | KIND=%s (--kind)  L=mark R=delete  | %s%s",
            $device, $fs, ($aux_mode eq 'same'||$aux_mode eq 'auto' ? $aux_mode : sprintf('%gV/div',($aux_mode+0)/1e6)), uc $cur_kind,
            (join(' ', map { "$_:".($count{$_}//0) } @show)),
            (defined $pending ? "  [start set-click end]" : '')),
            fontsize => 9);
        $fig->tight_layout;
        return $fig;
    },
);
