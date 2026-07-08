#!/usr/bin/env perl
# examples/read_nihonkohden.pl
# Usage:
#   perl -Ilib examples/read_nihonkohden.pl patient.eeg [--block N] [--plot]
#              [--sec S] [--nch N] [--uv U] [--chans LIST] [--aux MODE]
#
# Reads either Nihon Kohden layout (wfmblock / extblock) — dispatch is entirely
# inside read_nk(), so this viewer needs no format knowledge.
#
# Display options:
#   --sec S     : seconds per screen (default 10)
#   --nch N     : number of channels to show (default 8; ignored if --chans)
#   --uv U      : µV per channel division for EEG traces (default 100)
#   --chans L   : comma-separated channel NAMES to show, in order
#                 e.g. --chans Fp1,Cz,Pz,DC01,DC02,STIM
#                 (any label from the file, incl. DC*/STIM/$A1/… ; overrides --nch)
#   --aux MODE  : how to scale "aux" channels (DC*/STIM/PAD/COM/$*/BN*/Pulse/CO2…)
#                   same  : draw as-is, same µV scale as EEG (default).
#                           Triggers may overlap neighbours — often desirable
#                           for lining a TTL up against the EEG.
#                   auto  : auto-scale each aux channel to fit its own slot
#                           (keeps TTL square waves visible, no overlap).
#                   <N>   : give aux channels a fixed N µV/div (e.g. --aux 2000)
#                           — a middle ground: reduced but still overlapping.
#
# EEG channels always use the slider-controlled µV/div. Only aux channels are
# affected by --aux.

use strict;
use warnings;
use PDL;
use PDL::EEG::IO::NihonKohden qw(read_nk);

my @args   = @ARGV;
my $file   = shift @args or die "Usage: $0 subject.eeg [--block N] [--allblocks] [--plot] [--sec S] [--nch N] [--uv U] [--chans LIST] [--aux MODE]\n";
my $plot   = grep { $_ eq '--plot' } @args;
my $allblk = grep { $_ eq '--allblocks' } @args;

sub _arg {
    my ($name, $default, @a) = @_;
    for my $i (0..$#a) {
        return $1       if $a[$i] =~ /^--\Q$name\E=(.+)$/;
        return $a[$i+1] if $a[$i] eq "--$name" && defined $a[$i+1] && $a[$i+1] !~ /^--/;
    }
    return $default;
}
my $blk      = _arg('block', 0,      @args);
my $nsec     = _arg('sec',   10,     @args);
my $nch      = _arg('nch',   8,      @args);
my $uv_init  = _arg('uv',    100,    @args);  # µV per channel division (EEG)
my $chans    = _arg('chans', '',     @args);  # comma list of names
my $aux_mode = _arg('aux',   'same', @args);  # same | auto | <µV/div>

my $rec        = read_nk($file, block => $blk, ($allblk ? (all_blocks => 1) : ()));
my $data       = $rec->{data};
my $fs         = $rec->{fs};
my @labels     = @{ $rec->{labels} };
my @events     = @{ $rec->{events} };
my $n_ch_valid = $rec->{n_ch_valid};
my $n_total    = $data->dim(1);
my $duration   = $n_total / $fs;
my @gap_bounds = @{ $rec->{gap_bounds} // [] };   # [[start_samp,end_samp],...]

printf "Device : %s\n",   $rec->{device};
printf "Start  : %s\n",   $rec->{t_start};
printf "Loaded : %d ch (%d valid), %d samples @ %g Hz  [%s]\n",
    $data->dim(0), $n_ch_valid, $n_total, $fs,
    ($allblk ? sprintf("all %d blocks concat, %d gaps", $rec->{n_blocks}, scalar @gap_bounds)
             : sprintf("block %d/%d", $blk, $rec->{n_blocks} - 1));
printf "Events : %s\n",
    join(', ', map { "$_->{label}(\@$_->{t}s)" } @events) || 'none';

# --- aux channel classification (by label) ---
sub is_aux {
    my $l = shift // '';
    return 1 if $l =~ /^DC/i;
    return 1 if $l eq 'STIM' || $l eq 'PAD' || $l eq 'COM';
    return 1 if $l =~ /^\$/;                 # reference channels $A1, $Cz, ...
    return 1 if $l =~ /^BN/i;
    return 1 if $l =~ /^(Pulse|CO2|SpO2|EtCO2)/i;
    return 0;
}

# --- resolve which channel indices to display ---
my @sel;   # 0-based indices into @labels / $data dim0
if (length $chans) {
    my %by_name;
    for my $i (0 .. $#labels) { $by_name{lc $labels[$i]} //= $i }
    for my $nm (split /\s*,\s*/, $chans) {
        my $i = $by_name{lc $nm};
        if (defined $i) { push @sel, $i }
        else { warn "channel '$nm' not found (have: @labels)\n" }
    }
    die "no valid channels selected\n" unless @sel;
} else {
    my $n_plot = ($n_ch_valid < $nch) ? $n_ch_valid : $nch;
    @sel = (0 .. $n_plot - 1);
}
my $n_plot = scalar @sel;
my @sel_labels = @labels[@sel];

exit unless $plot;

# --- Load Cairo once, outside render loop ---
require PDL::Graphics::Cairo;
PDL::Graphics::Cairo->import(qw(subplots));
require PDL::Graphics::Cairo::Driver::GS;

my $gs = PDL::Graphics::Cairo::Driver::GS->new(
    width  => 1100,
    height => 80 + $n_plot * 70,
);

$gs->show_interactive(
    init   => { 0 => 0.0, 1 => 0.5 },
    render => sub {
        my ($state, $w, $h) = @_;

        # --- Slider 0: time offset ---
        my $max_offset = $duration - $nsec;
        $max_offset    = 0 if $max_offset < 0;
        my $t_off      = ($state->{0} // 0.0) * $max_offset;
        my $n_show     = int($nsec * $fs);
        $n_show        = 1 if $n_show < 1;
        my $s0         = int($t_off * $fs);
        $s0            = $n_total - $n_show if $s0 + $n_show > $n_total;
        $s0            = 0 if $s0 < 0;
        my $s1         = $s0 + $n_show - 1;
        $s1            = $n_total - 1 if $s1 >= $n_total;
        $s0            = $s1 if $s0 > $s1;

        # --- Slider 1: gain (µV/div) for EEG, logarithmic 10..500 ---
        my $sv         = $state->{1} // 0.5;
        my $uv_per_div = 10 * (50 ** (1.0 - $sv));   # 10..500 µV/div
        my $spacing    = $uv_per_div * 2.0;          # channel slot height (µV)

        my ($fig, $ax) = subplots(1, 1, width => $w, height => $h);

        # --- decimation to display width ---
        my $plot_w   = $w - 120; $plot_w = 200 if $plot_w < 200;
        my $n_raw    = $s1 - $s0 + 1;
        my $step_max = 5;
        my $step     = int($n_raw / $plot_w) || 1;
        $step        = $step_max if $step > $step_max;
        my $t;

        # --- draw each selected channel ---
        for my $i (0 .. $n_plot - 1) {
            my $ci  = $sel[$i];
            my $raw = $data->slice("($ci),$s0:$s1");
            my $sig = ($step > 1) ? $raw->slice("0:-1:$step") : $raw;
            $t //= sequence($sig->dim(0)) * ($step / $fs);

            my $offset = ($n_plot - 1 - $i) * $spacing;
            my $lab    = $sel_labels[$i];
            my $aux    = is_aux($lab);

            my $y;
            if ($aux && $aux_mode ne 'same') {
                my $mean = $sig->avg;
                if ($aux_mode eq 'auto') {
                    my $amp = ($sig - $mean)->abs->max; $amp = 1 if $amp <= 0;
                    $y = ($sig - $mean) / $amp * ($spacing * 0.4) + $offset;
                } else {                                  # fixed µV/div for aux
                    my $val = $aux_mode + 0; $val = 1 if $val <= 0;
                    $y = ($sig - $mean) * ($uv_per_div / $val) + $offset;
                }
            } else {
                $y = $sig + $offset;                      # as-is (µV), may overlap
            }

            my $col = $aux ? '#CC0000' : '#1565C0';
            $ax->line($t, $y, color => $col, lw => 0.7);
            $ax->axhline($offset, color => '#DDDDDD', lw => 0.3);
        }

        # --- recording-break (gap) markers: grey verticals in this window ---
        my $y_lo = -$spacing * 0.8;
        my $y_hi = ($n_plot - 0.2) * $spacing;
        for my $g (@gap_bounds) {
            my ($gs, $ge) = @$g;
            next if $ge < $s0 || $gs > $s1;             # gap not in this window
            my $gc = ($gs + $ge) / 2;                   # gap centre (sample)
            my $gx = ($gc - $s0) / $fs;                 # local x (s)
            $ax->line(PDL->new([$gx, $gx]), PDL->new([$y_lo, $y_hi]),
                      color => '#FF9800', lw => 1.5);   # orange = recording break
        }

        # --- EEG scale bar: uv_per_div vertical, in a right gutter (not over data) ---
        my $xmax   = $t->max;
        my $gutter = $xmax * 0.05 + 1e-9;      # empty margin right of the last sample
        my $bar_t  = $xmax + $gutter * 0.5;    # calibration bar sits here, clear of traces
        $ax->line(
            PDL->new([$bar_t, $bar_t]),
            PDL->new([0, $uv_per_div]),
            color => '#CC0000', lw => 2.0
        );

        my @tick_vals = map { ($n_plot - 1 - $_) * $spacing } 0 .. $n_plot - 1;
        $ax->yticks(PDL->new(\@tick_vals), [@sel_labels]);
        $ax->xlim(0, $xmax + $gutter);         # include the gutter so the bar is inside the frame but off the data
        $ax->ylim(-$spacing * 0.8, ($n_plot - 0.2) * $spacing);
        $ax->xlabel(sprintf("Time (s)  [%.1fs + %.1fs  step=%d]",
            $t_off, $nsec, $step));

        $fig->suptitle(
            sprintf("%s  blk=%d  %s  %gHz  | EEG %.0fµV/div  aux=%s  %gs/screen",
                $rec->{device}, $blk, $rec->{t_start}, $fs,
                $uv_per_div, $aux_mode, $nsec),
            fontsize => 9
        );
        $fig->tight_layout;
        return $fig;
    },
);
