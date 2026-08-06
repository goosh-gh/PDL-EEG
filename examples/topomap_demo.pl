#!/usr/bin/env perl

use eeg; use mat; use File::Basename; use PDL; use PDL::MatrixOps qw(inv);
use lib '/Users/goosh/src/PDL_Graphics_Cairo/lib';
use lib '/Users/goosh/src/PDL_EEG/lib'; use PDL::EEG::IO::ASA ();
use lib '/Users/goosh/src/PDL_EEG_MAP2D/lib';
use PDL::EEG::MAP2D qw(plot_topomap);

my $infile  =$ARGV[0];
my $lat     =$ARGV[1];
my $vol     =$ARGV[2];
my $verbose =$ARGV[3];

my ($rd, $rh) = &eeg::NK_READ_ONEAVGFILE2("$infile", "trgCH", 0); # &eeg::NK_heatmap_gnuplot($infile, $lat, $vol, $verbose);
my @channel_names4topo; for($i=0;$i<$rh->{general}->{ChannelN};$i++) { $channel_names4topo[$i] = $rh->{ChannelName}->{$i}; }
print "Got channels names for topo:"; &mat::DISPLAY_1DY(\@channel_names4topo);
my $map_time = $lat + $rh->{average}->{pretrigger};
my $basename = basename($infile);
my $outfile = "$basename"."_"."$lat"."ms"."$vol"."uV".".png"; print "outfilename: $outfile\n";

# --- montage: standard_1020.elc (subset ships with PDL-EEG t/data) ----------
# my $elc = shift @ARGV || 'standard_1020_subset.elc';
# my $elc = '/Users/goosh/src/PDL_IO_NYHead/standard_1020_plus_4eog.7.elc' || '/Users/goosh/src/NYHead/standard_1020.elc'; 
my $elc = '/Users/goosh/src/PDL_EEG/t/data/standard_1020_plus_4eog.7.elc' || '/Users/goosh/src/NYHead/standard_1020.elc'; 
# --- channel order of the recording (old 10-20 names) -----------------------
# my @labels = qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 C3 Cz C4 T4 T5 P3 Pz P4 T6 O1 O2);


my ($fig, $ax) = plot_topomap(
    avg     => $rd,                    # PDL (nchan,ntime) でも AoA でも可
    time    => $map_time,           # このサンプルの断面を描く
    labels  => \@channel_names4topo,         # $avg の行順に一致
    montage => $elc,
    clim    => $vol,                        # ±µV（省略で max|v| 自動対称）
    contours=> 10,
    names   => 1,
    margin  => 0.15,             # electrodes sit inside the head circle by this
    overshoot => 0.2,           # color extends this far *outside* the head circle
    eog_interp => 1,           # let EOG channels feed the interpolation (default off)
    cbar_label_x => 1.8,       # push the colour-bar numbers right (default 3.4), smaller more to right
    title   => "$basename, $lat ms, $vol uV",
    outfile => $outfile,              # 省略時は ($fig,$ax) を返し、pngは保存されない
    # device => 'gs',            # giza-server window ('osx'/'aqua' are aliases for 'gs')
    # device => 'osx',           # native Cocoa window (giza /osx)  "退役しました", should be aliased to gs
    # device defaults to 'png' -> writes 'outfile'
);
if($fig  =~ /png/) { print "wrote $outfile \n"; }

sub make_sep_data {
# --- synthetic average: $avg is (nchan, ntime) -----------------------------
#     SEP-like left centro-parietal positivity peaking mid-epoch.
my $ntime = 100;
my %peak = (Fp1=>-0.5,Fp2=>-0.6,F7=>-0.3,F3=>1.0,Fz=>1.5,F4=>0.4,F8=>-0.2, T3=>0.2,C3=>4.5,Cz=>3.0,C4=>0.6,T4=>-0.4,T5=>0.3,P3=>2.0,Pz=>0.8,P4=>0.0, T6=>-0.6,O1=>-1.2,O2=>-1.5);
my $t = sequence($ntime);
my $env = exp(-(($t-40)/8)**2);                       # Gaussian time course
my $avg = zeroes(scalar(@labels), $ntime);
$avg->slice("($_),:") .= $peak{$labels[$_]} * $env for 0..$#labels;
}

