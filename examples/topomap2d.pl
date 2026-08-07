#!/usr/bin/env perl
use strict; use warnings; use eeg; use mat; use File::Basename; use PDL; use PDL::MatrixOps qw(inv);
use lib '/Users/goosh/src/PDL_Graphics_Cairo/lib';
use lib '/Users/goosh/src/PDL_EEG/lib'; use PDL::EEG::IO::ASA ();
use lib '/Users/goosh/src/PDL_EEG_MAP2D/lib/'; use PDL::EEG::MAP2D qw(plot_topomap  plot_topomap_panels);

my $infile      =$ARGV[0];
my $lat         =$ARGV[1];
my $vol         =$ARGV[2];
my $orientation =$ARGV[3];
my $device      =$ARGV[4];
my $verbose     =$ARGV[5];

my ($rd, $rh) = &eeg::NK_READ_ONEAVGFILE2("$infile", "trgCH", 0); # &eeg::NK_heatmap_gnuplot($infile, $lat, $vol, $verbose);
my @channel_names4topo; for(my $i=0;$i<$rh->{general}->{ChannelN};$i++) { $channel_names4topo[$i] = $rh->{ChannelName}->{$i}; }
print "Got channels names for topo:"; &mat::DISPLAY_1DY(\@channel_names4topo);
my $map_time = $lat + $rh->{average}->{pretrigger};
my $basename = basename($infile); my $outfile; my $how2cut; my ($fig, $ax);
my $elc = '/Users/goosh/src/PDL_EEG_MAP2D/t/data/standard_1020_eog_nose.elc' || '/Users/goosh/src/PDL_IO_NYHead/standard_1020_plus_4eog.7.elc' || '/Users/goosh/src/NYHead/standard_1020.elc'; 
my $silhouette = '/Users/goosh/src/PDL_EEG/t/data/sagittal_silhouette.poly' || 'Users/goosh/src/PDL_EEG_MAP2D/t/data/sagittal_silhouette.poly';

if($orientation =~ /[Aa]x/) { $how2cut = "axial"; } elsif($orientation =~ /[Ll]t/) { $how2cut = "sagittal-left"; } elsif($orientation =~ /[Rr]t/) { $how2cut = "sagittal-right"; } else {$how2cut = "panel3"; }

if($how2cut eq "panel3") {
 $outfile = "$basename"."_"."$lat"."ms"."$vol"."uV"."_"."$how2cut".".png"  ; print "outfilename: $outfile \n";
 ($fig, $ax) = &topo_panel(); if($fig  =~ /png/) { print "wrote $outfile \n"; }
}
elsif( $how2cut eq "sagittal-left" ) {
 $outfile = "$basename"."_"."$lat"."ms"."$vol"."uV"."_"."$how2cut".".png"  ; print "outfilename: $outfile \n";
 ($fig, $ax) = &topo_1(); if($fig  =~ /png/) { print "wrote $outfile \n"; }
}
elsif( $how2cut eq "sagittal-right" ) {
 $outfile = "$basename"."_"."$lat"."ms"."$vol"."uV"."_"."$how2cut".".png"  ; print "outfilename: $outfile \n";
 ($fig, $ax) = &topo_1(); if($fig  =~ /png/) { print "wrote $outfile \n"; }
}
else {
 $outfile = "$basename"."_"."$lat"."ms"."$vol"."uV"."_"."$how2cut".".png"  ; print "outfilename: $outfile \n";
 ($fig, $ax) = &topo_1(); if($fig  =~ /png/) { print "wrote $outfile \n"; }
}
 

sub topo_1 { ### どれかの断面を1枚だけ
 ($fig, $ax) = 
 plot_topomap(
  avg     => $rd,                    # PDL (nchan,ntime) でも AoA でも可
  time    => $map_time,              # このサンプルのtopomapを描く
  labels  => \@channel_names4topo,   # $avg の行順に一致
  montage => $elc,
  clim    => $vol,                   # ±µV（省略で max|v| 自動対称）
  orientation=>$how2cut,
  silhouette=>$silhouette,
  cbar_label_x => 1.8,       # push the colour-bar numbers right (default 3.4), smaller more to right
  margin  => 0.20,           # electrodes sit inside the head circle by this   start:0.15
  overshoot => 0.15,         # color extends this far *outside* the head circle    start:0.1
  eog_interp => 1,           # let EOG channels feed the interpolation (1:on, 0: off)
  names=>1,
  title      => "$basename, $lat ms, $vol uV",
  # outfile=>"$outfile",
  device => 'gs',            # giza-server window ('osx'/'aqua' are aliases for 'gs')
 )
};


sub topo_panel { ### 3断面
 ($fig, $ax) = 
 plot_topomap_panels(
  avg        => $rd,
  time       => $map_time,
  labels     => \@channel_names4topo,
  montage    => $elc,
  clim       => $vol,                        # ±µV（省略で max|v| 自動対称）
  silhouette => $silhouette,
  cbar_label_x => 1.8,       # push the colour-bar numbers right (default 3.4), smaller more to right
  margin  => 0.20,           # electrodes sit inside the head circle by this
  overshoot => 0.15,          # color extends this far *outside* the head circle
  eog_interp => 1,           # let EOG channels feed the interpolation (1:on, 0: off)
  names      => 1,
  title      => "$basename, $lat ms, $vol uV",
  # outfile    => $outfile,
  device => 'gs',            # giza-server window ('osx'/'aqua' are aliases for 'gs')
 );
}


sub old_plot_topomap{
# my ($fig, $ax) = plot_topomap(
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
# );
# if($fig  =~ /png/) { print "wrote $outfile \n"; }
}

sub make_sep_data {
# --- synthetic average: $avg is (nchan, ntime) -----------------------------
#     SEP-like left centro-parietal positivity peaking mid-epoch.
my $ntime = 100; my @labels;
my %peak = (Fp1=>-0.5,Fp2=>-0.6,F7=>-0.3,F3=>1.0,Fz=>1.5,F4=>0.4,F8=>-0.2, T3=>0.2,C3=>4.5,Cz=>3.0,C4=>0.6,T4=>-0.4,T5=>0.3,P3=>2.0,Pz=>0.8,P4=>0.0, T6=>-0.6,O1=>-1.2,O2=>-1.5);
my $t = sequence($ntime);
my $env = exp(-(($t-40)/8)**2);                       # Gaussian time course
my $avg = zeroes(scalar(@labels), $ntime);
$avg->slice("($_),:") .= $peak{$labels[$_]} * $env for 0..$#labels;
}

