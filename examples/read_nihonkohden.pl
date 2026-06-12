#!/usr/bin/env perl
# examples/read_nihonkohden.pl
# Usage:
#   perl -Ilib examples/read_nihonkohden.pl patient.eeg [--block N] [--plot]
#              [--sec S] [--nch N] [--uv U]
#
# --sec S : seconds per screen (default 10)
# --nch N : channels to show (default 8)
# --uv U  : µV per channel height (default 100; clinical standard ~50-100)

use strict;
use warnings;
use PDL;
use PDL::EEG::IO::NihonKohden qw(read_nk);

my @args   = @ARGV;
my $file   = shift @args or die "Usage: $0 patient.eeg [--block N] [--plot] [--sec S] [--nch N] [--uv U]\n";
my $plot   = grep { $_ eq '--plot' } @args;

sub _arg {
    my ($name, $default, @a) = @_;
    for my $i (0..$#a) {
        return $1       if $a[$i] =~ /^--$name=(\S+)$/;
        return $a[$i+1] if $a[$i] eq "--$name" && defined $a[$i+1] && $a[$i+1] =~ /^\d/;
    }
    return $default;
}
my $blk      = _arg('block', 0,   @args);
my $nsec     = _arg('sec',   10,  @args);
my $nch      = _arg('nch',   8,   @args);
my $uv_init  = _arg('uv',    100, @args);  # µV per channel division

my $rec        = read_nk($file, block => $blk);
my $data       = $rec->{data};
my $fs         = $rec->{fs};
my @labels     = @{ $rec->{labels} };
my @events     = @{ $rec->{events} };
my $n_ch_valid = $rec->{n_ch_valid};
my $n_total    = $data->dim(1);
my $duration   = $n_total / $fs;

printf "Device : %s\n",   $rec->{device};
printf "Start  : %s\n",   $rec->{t_start};
printf "Loaded : %d ch (%d valid), %d samples @ %g Hz  [block %d/%d]\n",
    $data->dim(0), $n_ch_valid, $n_total, $fs,
    $blk, $rec->{n_blocks} - 1;
printf "Events : %s\n",
    join(', ', map { "$_->{label}(\@$_->{t}s)" } @events) || 'none';

exit unless $plot;

# --- Load Cairo once, outside render loop ---
require PDL::Graphics::Cairo;
PDL::Graphics::Cairo->import(qw(subplots));
require PDL::Graphics::Cairo::Driver::GS;

my $n_plot   = ($n_ch_valid < $nch) ? $n_ch_valid : $nch;
# spacing = 2 × uv_per_div  (channel center-to-center distance)
my $spacing0 = $uv_init * 2.0;

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

        # --- Slider 1: gain (µV/div) ---
        # slider=1.0 (up)   → small µV/div → sensitive (large display)
        # slider=0.0 (down) → large µV/div → compressed
        # Range: 500 µV/div (down) to 10 µV/div (up)
        my $sv       = $state->{1} // 0.5;
        # Logarithmic scale: 10µV at top, 500µV at bottom
        my $uv_per_div = 10 * (50 ** (1.0 - $sv));   # 10..500 µV/div
        my $spacing    = $uv_per_div * 2.0;

        my ($fig, $ax) = subplots(1, 1, width => $w, height => $h);

        # --- Downsample to display pixel width (1 point/pixel) ---
        # Note: SVG/PDF save also uses this render callback, so step is kept
        # conservative (max 5x) to preserve waveform quality in saved files.
        # For pure display speed, step could be set to n_raw/plot_w.
        my $plot_w   = $w - 120;
        $plot_w      = 200 if $plot_w < 200;
        my $n_raw    = $s1 - $s0 + 1;
        my $step_max = 5;   # max decimation: preserves ~200Hz resolution at 1000Hz
        my $step     = int($n_raw / $plot_w) || 1;
        $step        = $step_max if $step > $step_max;

        # Slice all channels at once: [n_plot, n_raw]
        my $block    = $data->slice("0:${\($n_plot-1)},$s0:$s1");

        # Decimate: [n_plot, n_decimated]
        my $decimated = ($step > 1)
            ? $block->slice(":,0:-1:$step")
            : $block;
        my $n_dec    = $decimated->dim(1);
        my $t        = sequence($n_dec) * ($step / $fs);

        # Draw all channels with offset
        for my $i (0 .. $n_plot - 1) {
            my $ch     = $decimated->slice("($i),:");
            my $offset = ($n_plot - 1 - $i) * $spacing;
            $ax->line($t, $ch + $offset, color => '#1565C0', lw => 0.7);
            $ax->axhline($offset, color => '#DDDDDD', lw => 0.3);
        }

        # Scale bar: 50µV vertical line at right edge
        my $bar_t    = $nsec * 0.97;
        my $bar_base = 0;
        $ax->line(
            PDL->new([$bar_t, $bar_t]),
            PDL->new([$bar_base, $bar_base + $uv_per_div]),
            color => '#CC0000', lw => 2.0
        );

        my @tick_vals  = map { ($n_plot - 1 - $_) * $spacing } 0 .. $n_plot - 1;
        $ax->yticks(PDL->new(\@tick_vals), [@labels[0 .. $n_plot - 1]]);
        $ax->xlim(0, $t->max);
        $ax->ylim(-$spacing * 0.8, ($n_plot - 0.2) * $spacing);
        $ax->xlabel(sprintf("Time (s)  [%.1fs + %.1fs  step=%d]",
            $t_off, $nsec, $step));

        $fig->suptitle(
            sprintf("%s  blk=%d  %s  %gHz  | %.0fµV/div  %gs/screen",
                $rec->{device}, $blk, $rec->{t_start}, $fs,
                $uv_per_div, $nsec),
            fontsize => 9
        );
        $fig->tight_layout;
        return $fig;
    },
);
