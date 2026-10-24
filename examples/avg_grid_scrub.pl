#!/usr/bin/env perl
# examples/avg_grid_scrub.pl
#
#   perl examples/avg_grid_scrub.pl AVG.tsv [MORE.tsv ...] [options]
#
#     --xlim A,B       time range [ms]                       (-250,800)
#     --ylim A,B       amplitude range [uV], same in every panel   (auto)
#     --positive-up    positive up (default: negative up)
#     --mark LIST      latencies [ms] drawn as vertical lines in every panel
#     --clim V         colour scale of the scalp map, +-V uV
#                      (default: the largest |amplitude| of the map at each
#                      latency, rounded up to a multiple of --cstep)
#     --cstep S        spacing of the numbers on the colour bar [uV], e.g. 0.5   (1)
#     --montage ELC    ASA .elc with the electrode positions
#                      (t/data/standard_1020_subset.elc next to this script's directory)
#     --no-map         waveforms only
#     --png FILE       write the grid to FILE and exit (no window)
#     --size W,H       window / image size in pixels          (1100,850)
#     --cell LIST      extra channels placed in free cells, NAME:ROW:COL[,...]
#                      (rows and columns from 0, e.g. BN1:4:0,BN2:4:4); a channel
#                      without a position in the montage is left out of the map
#     --names LIST     names of the files in the key, comma separated (the file names)
#     --colors LIST    line colours, one per file, comma separated
#                      (#1565C0, #9E9E9E, #C62828, #2E7D32, #6A1B9A, #EF6C00)
#
# AVG.tsv is a tab-separated average with a header row: time_ms, then one column
# per channel (e.g. the PREFIX_avg_<ref>.tsv files of heartbeat_bne_survey.pl).
# With more than one file, the waveforms of every file are drawn in the same
# panels (the first file on top of the others is not guaranteed: it is drawn
# first, so the read-out of an enlarged panel follows it), with a key in the top
# right cell. The files must have the same time axis; a channel missing from a
# later file is left out of that file's lines. The scalp map is that of the first
# file.
#
# Shows the channels as a 5x5 grid laid out like the 10-20 system
# (PDL::EEG::MAP2D::grid_layout_1020)
#
#          Fp1 [map] Fp2
#     F7   F3   Fz  F4   F8
#     T3   C3   Cz  C4   T4
#     T5   P3   Pz  P4   T6
#          O1       O2
#
# All panels have the same plot area. The axes are the same in every panel; the
# tick numbers are on one panel only (the leftmost of the lowest row, O1). The
# scalp map has its colour bar on its right and the latency under the head.
#
# The grid is shown in a giza-server window
# (PDL::Graphics::Cairo::Driver::GS::show_grid_scrub):
#
#   click or drag in a      chooses that latency: a vertical line in every panel
#     panel of the grid     of the grid, and the scalp map (top centre) redrawn
#                           for it
#   double-click a panel    opens that channel enlarged in a new window
#   move in that window     reads out the nearest sample (latency, amplitude)
#   click or drag there     chooses the latency in the same way and marks the sample
#   double-click the map    opens the map enlarged
#
# The scalp map is PDL::EEG::MAP2D::plot_topomap of the channels shown, at the
# sample nearest to the chosen latency. Without --clim its colour scale follows
# the map, so it changes with the latency; the colour bar gives the values, with
# numbers at every --cstep uV. Before a latency is chosen it shows the
# first --mark value, or stays empty. With --png the map is drawn for the first
# --mark value.
#
# Needs a running giza-server that reports double-click and drag (macOS Cocoa
# build), and a PDL::Graphics::Cairo that has show_grid_scrub with main_scrub
# and _scrub_skip, add_gridspec(cell_margins), label_only and colorbar(cmap).
# With --png nothing is displayed and giza-server is not needed.

use strict;
use warnings;
use PDL;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use File::Basename qw(basename);
use POSIX ();
use PDL::Graphics::Cairo qw(figure);
use PDL::EEG::MAP2D qw(plot_topomap grid_layout_1020);

# ---- settings ---------------------------------------------------------------
# Grid position (row, column) of each channel, and the free cell of the map:
# the 10-20 layout of PDL::EEG::MAP2D. To place other channels or move some,
# add or change entries here, e.g.  $CELL{Oz} = [4, 2];
my $GRID = grid_layout_1020();
my %CELL = %{ $GRID->{cell} };
my ($NROWS, $NCOLS, @MAPCELL) = ($GRID->{nrows}, $GRID->{ncols}, @{ $GRID->{map} });
# -----------------------------------------------------------------------------

my ($xlim, $ylim, $posup, $mark, $png, $size) = ('-250,800', undef, 0, '', undef, '1100,850');
my ($clim, $montage, $nomap, $cstep) = (undef, undef, 0, 1);
my ($names, $colors, $cells) = ('', '', '');
GetOptions('xlim=s' => \$xlim, 'ylim=s' => \$ylim, 'positive-up' => \$posup,
           'mark=s' => \$mark, 'png=s' => \$png, 'size=s' => \$size,
           'clim=f' => \$clim, 'montage=s' => \$montage, 'no-map' => \$nomap,
           'cstep=f' => \$cstep, 'names=s' => \$names, 'colors=s' => \$colors, 'cell=s' => \$cells) or die;
for my $c (split /,/, $cells) {
    my ($nm, $r, $cc) = split /:/, $c;
    defined $cc && $r =~ /^\d+$/ && $cc =~ /^\d+$/ && $r < $NROWS && $cc < $NCOLS
        or die "--cell: '$c' is not NAME:ROW:COL inside the ${NROWS}x$NCOLS grid\n";
    $CELL{$nm} = [$r, $cc];
}
$cstep > 0 or die "--cstep must be positive\n";
$montage //= "$FindBin::Bin/../t/data/standard_1020_subset.elc";
$nomap or -r $montage or die "montage $montage not readable: give --montage ELC or --no-map\n";
my @files = @ARGV or die "usage: $0 AVG.tsv [MORE.tsv ...] [options]\n";
my $file  = $files[0];
my @fname = $names ne '' ? split(/,/, $names) : map { basename($_) } @files;
my @fcol  = $colors ne '' ? split(/,/, $colors) : ('#1565C0', '#9E9E9E', '#C62828', '#2E7D32', '#6A1B9A', '#EF6C00');
@fname >= @files or die "--names: fewer names than files\n";
@fcol  >= @files or die "--colors: fewer colours than files\n";
my ($W, $H)   = split /,/, $size;
my ($x0, $x1) = split /,/, $xlim;
my @mark      = split /,/, $mark;

# ($M (ncol, ntime), @labels) of one average
sub read_avg {
    my $f = shift;
    open my $fh, '<', $f or die "$f: $!";
    chomp(my $head = <$fh>);
    my (undef, @lab) = split /\t/, $head;
    my @rows;
    while (<$fh>) { chomp; next unless length; push @rows, [ split /\t/ ] }
    close $fh;
    return (pdl(\@rows), @lab);
}
my ($M, @labels) = read_avg($file);              # (ncol, ntime)
my $t = $M->slice("(0)");

# channels that have a grid cell, in the order their panels are created; the
# index in this list is the panel index the driver reports on a double-click
my @show = grep { exists $CELL{ $labels[$_] } } 0 .. $#labels;
@show or die "no channel of $file has a grid position\n";
my @y    = map { $M->slice("(" . ($_ + 1) . ")") } @show;
my @name = @labels[@show];
# the other files: $yo[$f][$k] is the waveform of channel $name[$k] (undef if missing)
my @yo;
for my $f (1 .. $#files) {
    my ($Mf, @lf) = read_avg($files[$f]);
    $Mf->dim(1) == $t->nelem && ($Mf->slice("(0)") - $t)->abs->max < 1e-6
        or die "$files[$f]: time axis differs from $file\n";
    my %col; @col{@lf} = 1 .. @lf;
    $yo[$f] = [ map { exists $col{$_} ? $Mf->slice("($col{$_})") : undef } @name ];
}
my ($y0, $y1) = $ylim ? (split /,/, $ylim) : do {
    my $in = which(($t >= $x0) & ($t <= $x1));
    my $m  = 1.05 * pdl(map { $_->index($in)->abs->max } @y, grep { defined } map { @{ $yo[$_] } } 1 .. $#files)->max->sclr;
    (-$m, $m);
};
# all lines of channel $k into $ax (the first file first)
sub lines_of {
    my ($ax, $k, $lw) = @_;
    $ax->line($t, $y[$k], color => $fcol[0], lw => $lw);
    for my $f (1 .. $#files) {
        $ax->line($t, $yo[$f][$k], color => $fcol[$f], lw => $lw) if defined $yo[$f][$k];
    }
}
# key of the files (two or more), drawn into $ax
sub draw_key {
    my $ax = shift;
    $ax->{_scrub_skip} = 1;
    $ax->axis('off');
    $ax->xlim(0, 10);
    $ax->set_ylim(0, @files + 1);
    for my $f (0 .. $#files) {
        my $yk = @files - $f;
        $ax->line(pdl(0.3, 2.0), pdl($yk, $yk), color => $fcol[$f], lw => 2);
        $ax->text(2.4, $yk, $fname[$f], ha => 'left', halign => 'left', va => 'center',
                  valign => 'center', fontsize => 10, color => 'black');
    }
}
my $KEYCELL = [0, $NCOLS - 1];
my $key_index = @files > 1 ? @show + ($nomap ? 0 : 1) : -1;   # panel index of the key

my $mon  = $nomap ? undef : PDL::EEG::IO::ASA::read_elc($montage);
# the channels of the map: those with a position in the montage (aliases T3..T6
# as plot_topomap resolves them)
my %ALIAS_ELC = (T3 => 'T7', T4 => 'T8', T5 => 'P7', T6 => 'P8');
my @mapk = $nomap ? () : grep { my $n = $name[$_];
    exists $mon->{pos}{$n} || (exists $ALIAS_ELC{$n} && exists $mon->{pos}{ $ALIAS_ELC{$n} }) } 0 .. $#name;
my $Y    = cat(@y[@mapk ? @mapk : 0 .. $#y]);         # (ntime, nmap)
# figure title: the data file and, when a map is drawn, the montage file
my $title = $file . ($nomap ? '' : '    montage: ' . basename($montage));

# scalp map of the shown channels at the sample nearest to $lat, into $ax
# %cbar: options of the colour bar ($ax->colorbar of PDL::Graphics::Cairo)
sub draw_map {
    my ($fig, $ax, $lat, $res, %cbar) = @_;
    my $i = ($t - $lat)->abs->minimum_ind->sclr;
    $ax->{_scrub_skip} = 1;          # not a time axis: no latency line, not a click target
    my $v = $Y->slice("($i),:");
    # automatic colour scale: the bar ends on one of its numbers, and is never
    # shorter than one step
    my $lim = $clim // do { my $n = POSIX::ceil($v->abs->max->sclr / $cstep - 1e-9);
                            ($n < 1 ? 1 : $n) * $cstep };
    plot_topomap(values => $v, labels => [ @name[@mapk] ], montage => $mon,
                 clim => $lim, contours => 0, names => 0, sensor_size => 3, res => $res,
                 title => '', _ax => $ax, _fig => $fig,
                 colorbar => { location => 'right', tick_step => $cstep, %cbar });
    # latency, just under the head outline (the image has a blank band below it)
    $fig->add_text(text => sprintf('%.0f ms', $t->at($i)), ax => $ax,
                   xy => [0.5, 0], coords => 'axes', ha => 'center', va => 'bottom',
                   offset => [0, 3], fontsize => 10);
}

# The one panel that carries the tick numbers (the axes are the same in every
# panel): the leftmost panel of the lowest row shown, O1 in the full layout.
my ($numbered) = sort { $CELL{ $name[$b] }[0] <=> $CELL{ $name[$a] }[0]
                     || $CELL{ $name[$a] }[1] <=> $CELL{ $name[$b] }[1] } 0 .. $#name;

sub style {
    my ($ax, $big) = @_;
    $ax->axhline(0, color => '#999999', lw => 0.5);
    $ax->axvline(0, color => '#999999', lw => 0.5);
    $ax->axvline($_, color => '#C62828', lw => 1.0) for @mark;
    $ax->xlim($x0, $x1);
    $ax->set_ylim($y0, $y1);
    $ax->invert_yaxis unless $posup;
    $ax->clip_pad_top(0);
}

my $grid_render = sub {
    my ($state, $w, $h) = @_;
    my $fig = figure(width => $w, height => $h);
    # cell_margins => 0: every cell is a plot rectangle, so all panels have the
    # same size; the gaps between them are wspace / hspace
    my $gs  = $fig->add_gridspec($NROWS, $NCOLS, left => 0.02, right => 0.98, top => 0.955,
                                 bottom => 0.045, wspace => 0.12, hspace => 0.17,
                                 cell_margins => 0);
    for my $k (0 .. $#show) {
        my ($row, $col) = @{ $CELL{ $name[$k] } };
        my $ax = $fig->add_subplot($gs->at($row, $col));
        lines_of($ax, $k, 0.9);
        style($ax);
        # channel name inside the panel, at the top-left corner
        $ax->text($x0 + 0.03 * ($x1 - $x0), ($posup ? $y1 : $y0) + ($posup ? -1 : 1) * 0.04 * ($y1 - $y0),
                  $name[$k], ha => 'left', halign => 'left', va => 'top', valign => 'top',
                  fontsize => 11, color => 'black');
    }
    # the map panel comes after the channel panels, so the panel index the driver
    # reports for a channel is still its index in @show
    unless ($nomap) {
        my $lat = $state->{_scrub_lat} // $mark[0];
        # reserve => 0: the bar stands beside the head without moving it off
        # the centre of the column
        draw_map($fig, $fig->add_subplot($gs->at(@MAPCELL)), $lat, 110,
                 width => 7, pad => 2, length => 0.8, reserve => 0) if defined $lat;
    }
    # the key comes last, after the map panel
    draw_key($fig->add_subplot($gs->at(@$KEYCELL))) if @files > 1;
    $fig->label_only($fig->axes_list->[$numbered]);
    $fig->suptitle($title, fontsize => 11);
    return $fig;
};

my $cell_render = sub {
    my ($k, $state, $w, $h) = @_;
    my $fig = figure(width => $w, height => $h);
    if ($k == $key_index) {                    # the key was double-clicked
        draw_key($fig->axes);
        return $fig;
    }
    if ($k > $#show) {                         # the map panel was double-clicked
        # The driver puts its read-out and marker on the lines of the first
        # panel; an empty first panel keeps them off the head outline.
        my $gs = $fig->add_gridspec(1, 2, width_ratios => [1, 60], left => 0.0, right => 0.95,
                                    top => 0.88, bottom => 0.05);
        my $dummy = $fig->add_subplot($gs->at(0, 0));
        $dummy->axis('off');
        my $lat = $state->{_scrub_lat} // $mark[0] // 0;
        draw_map($fig, $fig->add_subplot($gs->at(0, 1)), $lat, 220, length => 0.6, label => 'uV');
        return $fig;
    }
    my $ax  = $fig->axes;
    lines_of($ax, $k, 1.3);
    style($ax);
    $ax->set_title($name[$k]);
    $ax->set_xlabel('time [ms]');
    $ax->set_ylabel('uV' . ($posup ? '' : ' (negative up)'));
    return $fig;
};

if (defined $png) {
    $grid_render->({}, $W, $H)->save($png);
    print "wrote $png\n";
    exit 0;
}

require PDL::Graphics::Cairo::Driver::GS;
my $drv = PDL::Graphics::Cairo::Driver::GS->new(width => $W, height => $H, title => "grid: $file");
$drv->show_grid_scrub(
    grid_render => $grid_render,
    cell_render => $cell_render,
    main_scrub  => 1,
    # everything but the scalp map (a _scrub_skip panel) is the same for every
    # latency: draw it once and keep it (redrawn only when the window size changes)
    static_key  => sub { '' },
    zoom_size   => [760, 560],
    vline       => { color => [0.78, 0.16, 0.16], lw => 1.2 },
    marker      => { color => [0.78, 0.16, 0.16], ms => 8, marker => 'o' },
    readout     => sub { my ($x, $yv) = @_; sprintf("%.0f ms\n%+.2f uV", $x, $yv) },
    cell_title  => sub { $_[0] == $key_index ? 'key' : $_[0] > $#show ? 'scalp map' : $name[ $_[0] ] },
);
