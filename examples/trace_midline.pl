#!/usr/bin/env perl
#
# trace_midline.pl — the midline electrode positions (10-20 or 10-10) along a line drawn on an
# image: a white line on a black background (or --invert) traced from the nasion end to the inion
# end, e.g. the skin of the midline traced by hand over a sagittal view.  The line is thinned,
# ordered from end to end, smoothed, and the positions are put at fractions of its length
# (Nasion 0, Fpz 10 %, Fz 30 %, Cz 50 %, Pz 70 %, Oz 90 %, Inion 100 %; with --system 10-10 also
# AFz FCz CPz POz).
#
#   perl examples/trace_midline.pl midline.png
#   perl examples/trace_midline.pl midline.png --out head1 --background OpenAnatomyAtlas.png --length-mm 360
#   perl examples/trace_midline.pl --selftest           # synthetic, no image needed
#
#   IMAGE            PNG (read with PDL::IO::PNG, else PDL::IO::Pic)
#   --out PREFIX     PREFIX_midline.tsv and PREFIX_midline.png                       (the image name)
#   --system S       10-20 | 10-10                                                         (10-20)
#   --start END      which end of the line is the nasion: left | right | top | bottom       (left)
#   --smooth W       moving average over W points before measuring (odd; 1 = the pixel chain) (5)
#   --threshold T    on the mean of the colour channels, 0 .. 1          ((min + max) / 2)
#   --invert         a dark line on a light background
#   --length-mm L    the measured nasion-inion length along the skin: positions also in mm
#   --background IMG draw the result over this image (same size; e.g. the picture the line was
#                    traced on) instead of over the line image
#   --no-labels      the PNG without the electrode names
#   --no-marks       the PNG without the position marks
#   --no-line        the PNG without the smoothed line (yellow)
#
# Printed: the size of the image, the skeleton, the length of the pixel chain and of the
# smoothed line, and for every position its fraction, image coordinates (x right, y down, origin
# at the top left, sub-pixel), length from the nasion end (pixels, and mm with --length-mm), and
# how far the position moves when measured along the pixel chain instead (the chain counts
# slanted stretches up to 8 % long).  The fractions are of the length of the traced line; the
# line is the skin as drawn, not a measured head.

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use PDL;
use PDL::EEG::ElectrodeLocation::TraceLine qw(binarize centerline);
use PDL::EEG::ElectrodeLocation::Polyline qw(arc_length point_at_fraction midline_labels);

my %o = (out => undef, system => '10-20', start => 'left', smooth => 5, threshold => undef, invert => 0,
         length_mm => undef, background => undef, labels => 1, marks => 1, line => 1, selftest => 0);
GetOptions('out=s' => \$o{out}, 'system=s' => \$o{system}, 'start=s' => \$o{start}, 'smooth=i' => \$o{smooth},
           'threshold=f' => \$o{threshold}, 'invert' => \$o{invert}, 'length-mm=f' => \$o{length_mm},
           'background=s' => \$o{background}, 'labels!' => \$o{labels}, 'marks!' => \$o{marks}, 'line!' => \$o{line},
           'selftest' => \$o{selftest})
    or die "bad options (see the header of $0)\n";
if ($o{selftest}) { exit(selftest() ? 0 : 1) }
my $file = shift @ARGV;
defined $file or die "usage: $0 IMAGE [options]  (see the header)\n";
($o{out} //= $file) =~ s/\.png$//i;

my $img = read_png($file);                                       # (c,nx,ny), 0 .. 1, row 0 at the top
my $res = trace({ %o, file => $file }, $img);
report(\%o, $res);
write_tsv(\%o, $res);
plot(\%o, $res, defined $o{background} ? read_png($o{background}) : $img);
exit 0;

# =====================================================================================
sub read_png {
    my ($f) = @_;
    -f $f or die "not found: $f\n";
    if (eval { require PDL::IO::PNG; 1 }) { return PDL::IO::PNG::rpng($f) }   # (c,nx,ny) float 0..1, row 0 at the top
    require PDL::IO::Pic;
    my $p = PDL::IO::Pic::rpic($f);                              # (c,nx,ny), row 0 at the bottom
    return ($p->slice(':,:,-1:0')->float / ($p->max > 1 ? 255 : 1))->copy;
}

sub trace {
    my ($opt, $img) = @_;
    my $mask = binarize($img, (defined $opt->{threshold} ? (threshold => $opt->{threshold}) : ()), invert => $opt->{invert});
    my $cl = centerline($mask, start => $opt->{start}, smooth => $opt->{smooth});
    my ($lab, $f) = midline_labels($opt->{system});
    my $cum  = arc_length($cl->{smooth});
    my $xy   = point_at_fraction($cl->{smooth}, $f, $cum);       # (k,2)
    my $xy_c = point_at_fraction($cl->{raw}, $f);                # the same from the pixel chain
    my $shift = sqrt((($xy - $xy_c)**2)->mv(1, 0)->sumover);     # (k): over x, y (dim 1)
    return { %$cl, labels => $lab, frac => $f, xy => $xy, shift => $shift, nx => $mask->dim(0), ny => $mask->dim(1),
             mask_px => $mask->sum, file => $opt->{file} // 'image' };
}

sub report {
    my ($opt, $r) = @_;
    printf "%s: %d x %d, %d line pixels; skeleton %d pixels (%d iterations), %d end points, %d connected part%s (the largest used)\n",
        $r->{file} // 'image', $r->{nx}, $r->{ny}, $r->{mask_px}, $r->{skeleton}->sum, $r->{iterations}, $r->{endpoints}->dim(0),
        $r->{n_components}, $r->{n_components} == 1 ? '' : 's';
    printf "length: pixel chain %.1f px, smoothed (%d points) %.1f px (chain %+.1f %%)\n", $r->{length_raw}, $opt->{smooth},
        $r->{length}, 100 * ($r->{length_raw} / $r->{length} - 1);
    my $mm = defined $opt->{length_mm} ? $opt->{length_mm} / $r->{length} : undef;
    printf "scale: %.4f mm per pixel (--length-mm %g)\n", $mm, $opt->{length_mm} if defined $mm;
    printf "%-7s %6s %8s %8s %9s%s %13s\n", 'label', '%', 'x', 'y', 'arc px', defined $mm ? '   arc mm' : '', 'chain shift';
    for my $i (0 .. $#{ $r->{labels} }) {
        my $s = $r->{frac}->at($i) * $r->{length};
        printf "%-7s %6.1f %8.1f %8.1f %9.1f%s %10.1f px\n", $r->{labels}[$i], 100 * $r->{frac}->at($i),
            $r->{xy}->at($i, 0), $r->{xy}->at($i, 1), $s, defined $mm ? sprintf('%9.1f', $s * $mm) : '', $r->{shift}->at($i);
    }
}

sub write_tsv {
    my ($opt, $r) = @_;
    my $f = "$opt->{out}_midline.tsv";
    open my $fo, '>', $f or die "$f: $!";
    print $fo "# trace_midline.pl $r->{file}: image coordinates (x right, y down, origin at the top left); smoothed over $opt->{smooth} points\n";
    my $mm = defined $opt->{length_mm} ? $opt->{length_mm} / $r->{length} : undef;
    print $fo join("\t", qw(label percent x y arc_px), defined $mm ? 'arc_mm' : ()), "\n";
    for my $i (0 .. $#{ $r->{labels} }) {
        my $s = $r->{frac}->at($i) * $r->{length};
        print $fo join("\t", $r->{labels}[$i], sprintf('%.1f', 100 * $r->{frac}->at($i)), sprintf('%.2f', $r->{xy}->at($i, 0)),
                       sprintf('%.2f', $r->{xy}->at($i, 1)), sprintf('%.2f', $s), defined $mm ? sprintf('%.1f', $s * $mm) : ()), "\n";
    }
    close $fo;
    print "wrote $f\n";
}

sub plot {
    my ($opt, $r, $bg) = @_;
    if (!eval { require PDL::Graphics::Cairo; 1 }) { print "PDL::Graphics::Cairo not found: no PNG\n"; return }
    my ($nx, $ny) = ($r->{nx}, $r->{ny});
    if ($bg->dim(1) != $nx || $bg->dim(2) != $ny) {
        printf "--background is %d x %d, the line image %d x %d: drawn over the line image\n", $bg->dim(1), $bg->dim(2), $nx, $ny;
        $bg = undef;
    }
    my $rgb = defined $bg ? $bg->slice('0:2,:,:') : do { my $s = $r->{skeleton}->float; cat($s, $s, $s)->reorder(2, 0, 1) };
    my $fig = PDL::Graphics::Cairo::figure(width => 2 * $nx + 120, height => 2 * $ny + 100);
    my $ax  = $fig->axes;
    $ax->imshow($rgb->reorder(2, 1, 0));                         # (rows, cols, 3), row 0 drawn at the top
    my $Y = sub { $ny - $_[0] };                                 # image y (down) -> plot y (up)
    my $sm = $r->{smooth};
    $ax->line($sm->slice(':,(0)'), $Y->($sm->slice(':,(1)')), color => 'yellow', lw => 1) if $opt->{line} // 1;
    for my $i (0 .. $#{ $r->{labels} }) {
        my ($x, $y) = ($r->{xy}->at($i, 0), $r->{xy}->at($i, 1));
        my $lab = $r->{labels}[$i];
        my $col = $lab =~ /^(Nasion|Inion)$/ ? 'red' : 'cyan';
        $ax->scatter(pdl($x), pdl($Y->($y)), color => $col, s => 6) if $opt->{marks} // 1;
        $ax->text($x, $Y->($y), " $lab", color => $col, fontsize => 12) if $opt->{labels} // 1;
    }
    $ax->title(sprintf '%s, %s; length %.1f px%s', $r->{file}, $opt->{system}, $r->{length},
               defined $opt->{length_mm} ? sprintf(' = %g mm', $opt->{length_mm}) : '');
    my $f = "$opt->{out}_midline.png";
    $fig->save($f);
    print "wrote $f\n";
}

# =====================================================================================
sub selftest {
    my $ok = 1;
    my $chk = sub { my ($c, $m) = @_; printf "  [%s] %s\n", $c ? 'PASS' : 'FAIL', $m; $ok &&= $c };
    my $PI = 3.14159265358979;
    # a semicircle (centre (250, 260), R 180, 3 px wide), the open side down, with a side branch
    my ($nx, $ny, $R) = (500, 300, 180);
    my $X = xvals(double, $nx, $ny) - 250; my $Yc = yvals(double, $nx, $ny) - 260;
    my $rr = sqrt($X**2 + $Yc**2);
    my $line = (abs($rr - $R) <= 1.5) * ($Yc <= 0);
    $line |= (abs($X - 0) <= 1) * ($Yc >= -$R + 1) * ($Yc <= -$R + 15);  # a 14 px spur hanging from the top inside
    my $img = cat(map { $line->float } 1 .. 3)->reorder(2, 0, 1)->copy;   # (3,nx,ny)
    my %opt = (%o, system => '10-20', start => 'left', smooth => 5);
    my $r = trace(\%opt, $img);
    my $true = $PI * $R;
    $chk->(abs($r->{length} / $true - 1) < 0.01, sprintf 'semicircle length %.1f px (true %.1f; pixel chain %.1f)', $r->{length}, $true, $r->{length_raw});
    $chk->($r->{length_raw} > 1.02 * $true, 'the pixel chain over-counts the length');
    $chk->($r->{endpoints}->dim(0) >= 3, sprintf 'the spur makes %d end points; the longest path still runs end to end', $r->{endpoints}->dim(0));
    my $ang = atan2(-($r->{xy}->slice(':,(1)') - 260), -($r->{xy}->slice(':,(0)') - 250)) * 180 / $PI;   # 0 at the left end, 180 at the right
    my $want = pdl(0, 18, 54, 90, 126, 162, 180);
    $chk->(($ang - $want)->abs->max < 1, sprintf 'positions at the right angles (max error %.2f deg)', ($ang - $want)->abs->max);
    $chk->(abs($r->{xy}->at(0, 0) - 70) < 2 && abs($r->{xy}->at(6, 0) - 430) < 2, 'Nasion at the left end, Inion at the right (start left)');
    my $rr2 = trace({ %opt, start => 'right' }, $img);
    $chk->(abs($rr2->{xy}->at(0, 0) - 430) < 2, 'start right: Nasion at the right end');
    my $r10 = trace({ %opt, system => '10-10' }, $img);
    $chk->(@{ $r10->{labels} } == 11 && $r10->{labels}[4] eq 'FCz', '10-10: eleven positions');
    # a dark line on white
    my $inv = trace({ %opt, invert => 1 }, 1 - $img);
    $chk->(abs($inv->{length} - $r->{length}) < 1e-9, '--invert reads a dark line on white the same');
    # output files
    require File::Temp;
    my $dir = File::Temp::tempdir(CLEANUP => 1);
    my $out = '';
    { open my $fh, '>', \$out or die; my $old = select $fh; my %oo = (%opt, out => "$dir/t", length_mm => 360);
      report(\%oo, $r); write_tsv(\%oo, $r); plot(\%oo, $r, $img); select $old; close $fh }
    print $out;
    my ($cz) = grep { /^Cz\t/ } do { open my $ft, '<', "$dir/t_midline.tsv" or die; <$ft> };
    my @czf = split /\t/, ($cz // ''); chomp @czf;
    $chk->(scalar(@czf == 6 && abs($czf[2] - 250) < 1 && abs($czf[3] - 80) < 1 && $czf[5] == 180),
           sprintf 'TSV: Cz at the top (%s, %s), %s mm of 360', @czf[2, 3, 5]);
    $chk->(scalar(-s "$dir/t_midline.png" || $out =~ /no PNG/), 'PNG written (or PDL::Graphics::Cairo absent)');
    if (-s "$dir/t_midline.png") {                                # the same picture without names and marks is smaller
        { open my $fh, '>', \my $o2 or die; my $old = select $fh; plot({ %opt, out => "$dir/n", labels => 0, marks => 0 }, $r, $img); select $old; close $fh }
        $chk->(-s "$dir/n_midline.png" && -s "$dir/n_midline.png" < -s "$dir/t_midline.png", '--no-labels --no-marks: a PNG without them');
    }
    print $ok ? "selftest passed\n" : "selftest FAILED\n";
    return $ok;
}
