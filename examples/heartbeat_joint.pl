#!/usr/bin/env perl
# examples/heartbeat_joint.pl
#
#   perl examples/heartbeat_joint.pl AVG.tsv [options]
#
#     --times LIST     latencies of the scalp maps [ms], e.g. --times 230,320
#                      (default: @TIMES below)
#     --xlim A,B       time range of the waveform panel   (-250,800)
#     --ylim A,B       amplitude range [uV]               (auto)
#     --clim V         colour scale of the maps, +-V uV   (largest |amplitude| at --times)
#     --positive-up    positive up (default: negative up)
#     --nave N         number of beats, printed in the panel
#     --title TEXT
#     --montage ELC    ASA .elc              (t/data/standard_1020_subset.elc)
#     --out FILE       PNG                   (AVG without .tsv + _joint.png)
#
# AVG.tsv is a tab-separated average with a header row: time_ms, then one column
# per channel (the PREFIX_avg_<ref>.tsv files of heartbeat_bne_survey.pl). Draws
# the figure with PDL::EEG::MAP2D::plot_joint:
#
#   lower panel   every channel as one trace, coloured by the position of its
#                 electrode (right -> red, anterior -> blue, superior -> green)
#   small head    inside the lower panel, top left: the key to those colours.
#                 One dot per electrode at its scalp position (nose up), drawn
#                 in the colour of that electrode's trace
#   maps          scalp distribution at each --times latency, one colour scale
#                 for all (bar at the right), each joined by a line to its
#                 latency on the time axis
#   N = ...       above the top right corner of the lower panel, from --nave

use strict;
use warnings;
use PDL;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use PDL::EEG::MAP2D qw(plot_joint);

# ---- settings ---------------------------------------------------------------
# Latencies [ms] of the scalp maps, used when --times is not given. Any number
# of values, in any spacing; each map is taken at the sample nearest to its
# value and labelled with the value written here.
my @TIMES = (200, 250, 300, 350, 400, 450, 500, 550, 600, 650);
# -----------------------------------------------------------------------------

my ($times, $xlim, $ylim, $clim, $posup, $nave, $title, $montage, $out) = (undef, '-250,800');
GetOptions('times=s' => \$times, 'xlim=s' => \$xlim, 'ylim=s' => \$ylim, 'clim=f' => \$clim,
           'positive-up' => \$posup, 'nave=i' => \$nave, 'title=s' => \$title,
           'montage=s' => \$montage, 'out=s' => \$out) or die;
my $file = shift @ARGV or die "usage: $0 AVG.tsv [options]\n";
($out = $file) =~ s/\.tsv$//i, $out .= '_joint.png' unless defined $out;
$montage //= "$FindBin::Bin/../t/data/standard_1020_subset.elc";
my @times = defined $times ? (split /,/, $times) : @TIMES;
@times or die "no map latency: give --times or fill \@TIMES\n";

open my $fh, '<', $file or die "$file: $!";
chomp(my $head = <$fh>);
my ($tcol, @labels) = split /\t/, $head;
my @rows;
while (<$fh>) { chomp; next unless length; push @rows, [ split /\t/ ] }
close $fh;
my $M   = pdl(\@rows);                         # ncol columns x ntime rows -> dims (ncol, ntime)
my $t   = $M->slice("(0)");
my $avg = $M->slice("1:-1");                   # (nch, ntime)

plot_joint(
    avg => $avg, t => $t, labels => \@labels, montage => $montage,
    times => \@times, xlim => [ split /,/, $xlim ],
    ($ylim ? (ylim => [ split /,/, $ylim ]) : ()),
    (defined $clim ? (clim => $clim) : ()),
    negative_up => ($posup ? 0 : 1),
    (defined $nave ? (nave => $nave) : ()),
    (defined $title ? (title => $title) : ()),
    xlabel => 'time from R peak [ms]', outfile => $out,
);
print "wrote $out\n";
