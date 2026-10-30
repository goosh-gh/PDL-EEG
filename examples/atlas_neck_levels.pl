#!/usr/bin/env perl
#
# atlas_neck_levels.pl — the length along the skin of the back of the neck from the inion to
# each spinous process (C1 .. T3), measured on two side views of an anatomy atlas taken with the
# same camera: one with the skin (SKIN, e.g. Open Anatomy Atlas / BodyParts3D with the skin
# shown), one with the bones only (BONE, cropped differently).  The two are aligned by phase
# correlation (a shift only: same zoom, no rotation), the back outline of the body is taken from
# SKIN, the inion and the tips of the spinous processes from BONE, and every bony point is put
# on the outline at its nearest point.  The lengths are scaled so that inion -> C7 is --c7-mm
# (a length measured on a subject; or --mm-per-px).  The outline (in mm, inion at 0) and the
# levels are written to a TSV that PDL-IO-NYHead nyhead_neck_electrodes.pl --extend-with reads.
#
#   perl examples/atlas_neck_levels.pl OpenAnatomyAtlas_wide.png OpenAtlasBone.png
#   perl examples/atlas_neck_levels.pl SKIN.png BONE.png --c7-mm 165 --mark 'Cv12=85,Cv45=125,BN1=165,Th12=200'
#   perl examples/atlas_neck_levels.pl --selftest
#
#   SKIN BONE        PNG (PDL::IO::PNG, else PDL::IO::Pic); the face to the left
#   --out PREFIX     PREFIX_neck_profile.tsv and PREFIX_neck_levels.png       (atlas)
#   --inion-row Y    row of SKIN (y down) at which the inion is taken: the most posterior bone
#                    pixel of that row, put on the outline                     (266)
#   --tips LIST      NAME=X,Y;...  guesses of the spinous-process tips in SKIN coordinates; each
#                    is moved to the most posterior bone pixel within --box of it
#                    (C1..T3 of the Open Anatomy Atlas views of 2026-10-28, bp3d43; in the side
#                    view T2 and T3 lie over the ribs, so the tips found there are less certain,
#                    and the deep C1 tip moves between two parts of the hollow of the neck with
#                    --smooth: 59 .. 81 mm for H = 2 .. 8, the other levels within 3 mm)
#   --box DX,DY      half size of the search box around a guess, pixels       (15,8)
#   --smooth H       the outline x(y) is averaged over 2H+1 rows               (4)
#   --c7-mm MM       inion -> C7 along the skin, mm (sets the scale)           (165)
#   --mm-per-px S    the scale given directly (instead of --c7-mm)
#   --mark LIST      NAME=MM,...  lengths from the inion to show against the levels
#   --no-png         no picture
#
# Printed: the shift between the two views (and how distinct the correlation peak is), the
# scale, the inion, and for every level the tip, its nearest point on the outline, the length
# from the inion (px, mm) and the depth of the tip under the skin (mm); with --mark the level
# interval of every mark.  The outline is the silhouette of the body (the outermost skin, a
# little off the midline where the neck muscles bulge), one person's model in one posture.
#
# The TSV: name, s_mm (length from the inion), y_mm (forward +), z_mm (up +), kind (profile |
# level), and for the levels tip_y_mm, tip_z_mm, depth_mm; y and z are relative to the inion
# point on the outline, in the plane of the picture.

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use PDL;
use PDL::EEG::ElectrodeLocation::Polyline qw(arc_length point_at_length);
use PDL::EEG::ElectrodeLocation::AtlasView qw(read_png atlas_levels level_of);

my $TIPS = 'C1=415,331;C2=422,388;C3=419,414;C4=415,437;C5=422,458;C6=435,474;C7=449,502;'
         . 'T1=456,536;T2=463,571;T3=482,610';
my %o = (out => 'atlas', inion_row => 266, tips => $TIPS, box => '15,8', smooth => 4, c7_mm => 165,
         mm_per_px => undef, mark => '', png => 1, selftest => 0);
GetOptions('out=s' => \$o{out}, 'inion-row=i' => \$o{inion_row}, 'tips=s' => \$o{tips}, 'box=s' => \$o{box},
           'smooth=i' => \$o{smooth}, 'c7-mm=f' => \$o{c7_mm}, 'mm-per-px=f' => \$o{mm_per_px},
           'mark=s' => \$o{mark}, 'png!' => \$o{png}, 'selftest' => \$o{selftest})
    or die "bad options (see the header of $0)\n";
if ($o{selftest}) { exit(selftest() ? 0 : 1) }
@ARGV == 2 or die "usage: $0 SKIN.png BONE.png [options]  (see the header)\n";
my ($fskin, $fbone) = @ARGV;

my $skin = read_png($fskin);                                      # (3,nx,ny), 0 .. 1, row 0 at the top
my $bone = read_png($fbone);
my $r = atlas_levels(\%o, $skin, $bone);
report(\%o, $r, $fskin, $fbone);
write_tsv(\%o, $r, $fskin, $fbone);
plot(\%o, $r, $skin) if $o{png};
exit 0;

# =====================================================================================

sub report {
    my ($opt, $r, $fs, $fb) = @_;
    printf "%s and %s: BONE is shifted by (%d, %d) px in SKIN (correlation peak %.1f x the next highest)\n",
        $fs, $fb, $r->{dx}, $r->{dy}, $r->{peak};
    printf "back outline: %d rows (y %d .. %d), averaged over %d rows\n", $r->{P}->dim(0), $r->{P}->at(0, 1),
        $r->{P}->at(-1, 1), 2 * $opt->{smooth} + 1;
    printf "inion: bone (%.0f, %.0f) -> outline (%.1f, %.1f), %.1f mm under the skin\n", $r->{inion}->list,
        $r->{inion_skin}->list, $r->{inion_depth} * $r->{mm};
    printf "scale: %.4f mm per px (%s)\n", $r->{mm}, $r->{by_scale};
    printf "  %-4s %13s %15s %9s %8s %9s\n", 'name', 'guess (x,y)', 'tip (x,y)', 'from inion', '', 'depth';
    for my $l (@{ $r->{levels} }) {
        printf "  %-4s %6.0f,%-6.0f %7.1f,%-7.1f %7.1f px %6.1f mm %6.1f mm\n", $l->{name}, @{ $l->{guess} },
            $l->{tip}->list, $l->{s}, $l->{s} * $r->{mm}, $l->{depth} * $r->{mm};
    }
    return if $opt->{mark} eq '';
    print "marks (mm from the inion along the outline):\n";
    for my $kv (split /,/, $opt->{mark}) {
        my ($n, $L) = $kv =~ /^([^=]+)=([\d.]+)$/ or die "--mark NAME=MM, not '$kv'\n";
        printf "  %-6s %6.1f mm  %s\n", $n, $L, level_of($r, $L);
    }
}

sub write_tsv {
    my ($opt, $r, $fs, $fb) = @_;
    my $f = "$opt->{out}_neck_profile.tsv";
    my ($P, $mm) = ($r->{P}, $r->{mm});
    my ($xI, $yI) = $r->{inion_skin}->list;
    my $s  = ($r->{cum} - $r->{sI}) * $mm;                        # (n) mm from the inion
    my $Y  = -($P->slice(':,(0)') - $xI) * $mm;                   # forward + (the face is to the left)
    my $Z  = -($P->slice(':,(1)') - $yI) * $mm;                   # up +
    open my $fo, '>', $f or die "$f: $!";
    print $fo "# atlas_neck_levels.pl $fs $fb: shift ($r->{dx}, $r->{dy}) px, $mm mm/px ($r->{by_scale}), inion row $opt->{inion_row}\n";
    print $fo "# s_mm: length along the back outline from the inion; y_mm forward +, z_mm up +, from the inion point on the outline\n";
    print $fo join("\t", qw(name s_mm y_mm z_mm kind tip_y_mm tip_z_mm depth_mm)), "\n";
    for my $l (@{ $r->{levels} }) {
        my $q = point_at_length($P, pdl($l->{s} + $r->{sI}), $r->{cum})->slice('(0),');
        printf $fo "%s\t%.2f\t%.2f\t%.2f\tlevel\t%.2f\t%.2f\t%.2f\n", $l->{name}, $l->{s} * $mm,
            -($q->at(0) - $xI) * $mm, -($q->at(1) - $yI) * $mm,
            -($l->{tip}->at(0) - $xI) * $mm, -($l->{tip}->at(1) - $yI) * $mm, $l->{depth} * $mm;
    }
    my $keep = which($s >= -10);                                  # the outline from 10 mm above the inion down
    printf $fo "p%04d\t%.2f\t%.2f\t%.2f\tprofile\n", $_, $s->at($_), $Y->at($_), $Z->at($_) for $keep->list;
    close $fo;
    print "wrote $f\n";
}

sub plot {
    my ($opt, $r, $skin) = @_;
    if (!eval { require PDL::Graphics::Cairo; 1 }) { print "PDL::Graphics::Cairo not found: no PNG\n"; return }
    my ($nx, $ny) = ($skin->dim(1), $skin->dim(2));
    my $fig = PDL::Graphics::Cairo::figure(width => 2 * $nx + 120, height => 2 * $ny + 100);
    my $ax  = $fig->axes;
    $ax->imshow($skin->reorder(2, 1, 0));                         # (rows, cols, 3), row 0 drawn at the top
    my $Y = sub { $ny - $_[0] };                                  # image y (down) -> plot y (up)
    my ($P, $mm) = ($r->{P}, $r->{mm});
    my $s = ($r->{cum} - $r->{sI}) * $mm;
    my $w = which(($s >= 0) & ($s <= 260));
    $ax->line($P->slice(':,(0)')->index($w), $Y->($P->slice(':,(1)')->index($w)), color => 'lime', lw => 1.5);
    for my $l (@{ $r->{levels} }) {
        my ($x, $y) = $l->{tip}->list;
        $ax->scatter(pdl($x), pdl($Y->($y)), color => 'red', s => 5);
        $ax->text($x, $Y->($y), " $l->{name}", color => 'red', fontsize => 11);
    }
    my @mk = ([ 'Inion', 0 ], map { [ split /=/ ] } grep { $_ ne '' } split /,/, $opt->{mark});
    for my $m (@mk) {
        my $q = point_at_length($P, pdl($m->[1] / $mm + $r->{sI}), $r->{cum})->slice('(0),');
        $ax->scatter(pdl($q->at(0)), pdl($Y->($q->at(1))), color => 'gold', s => 7);
        $ax->text($q->at(0), $Y->($q->at(1)), " $m->[0] $m->[1]", color => 'gold', fontsize => 12);
    }
    $ax->title(sprintf 'back outline from the inion (lime), spinous processes (red); %.3f mm/px (%s)', $mm, $r->{by_scale});
    my $f = "$opt->{out}_neck_levels.png";
    $fig->save($f);
    print "wrote $f\n";
}

# =====================================================================================
# synthetic views: a body whose back outline is x(y) = 300 + 30 sin(2 pi y / 400), and five
# "spinous processes" (bars ending at known tips) drawn in a bone layer.  SKIN = body + bones,
# BONE = the bones only, cropped with a known offset.
sub selftest {
    my $ok = 1;
    my $chk = sub { my ($c, $m) = @_; printf "  [%s] %s\n", $c ? 'PASS' : 'FAIL', $m; $ok &&= $c };
    my ($nx, $ny) = (420, 700);
    my $PI = 3.14159265358979;
    my $X = xvals(double, $nx, $ny); my $Yv = yvals(double, $nx, $ny);
    my $edge = 300 + 30 * sin(2 * $PI * $Yv / 400);
    my $body = ($X <= $edge) & ($X >= 60);                       # background at the top-left pixel, as in the real views
    my @tips = map { [ 240 + 6 * $_, 200 + 70 * $_ ] } 0 .. 4;      # tips (x, y)
    my $bonel = zeroes(byte, $nx, $ny);
    for my $t (@tips) {                                            # a bar 60 px long, 7 px high, ending at the tip
        $bonel |= (($X <= $t->[0]) & ($X >= $t->[0] - 60) & (abs($Yv - $t->[1]) <= 3));
    }
    $bonel |= (($X - 230)**2 + ($Yv - 120)**2 <= 50**2) & ($X <= 280);   # a "skull" with its back at x 280, row 120
    my $bgc = pdl(30, 37, 46) / 255;
    my $mk = sub {
        my ($m, $col) = @_;
        my $img = zeroes(3, $nx, $ny) + $bgc;
        $img = $img * (1 - $m->dummy(0)) + $m->dummy(0) * pdl(@$col)->dummy(1)->dummy(2);
        return $img;
    };
    my $skin = $mk->($body, [0.80, 0.65, 0.55]);
    $skin = $skin * (1 - $bonel->dummy(0)) + $bonel->dummy(0) * pdl(0.90, 0.90, 0.85)->dummy(1)->dummy(2);
    my $bonefull = $mk->($bonel, [0.90, 0.90, 0.85]);
    # textures so that the correlation has something besides the bones
    my ($ox, $oy) = (37, 21);
    my $bone = $bonefull->slice(sprintf('%s,%d:%d,%d:%d', ':', $ox, $nx - 1, $oy, $ny - 1))->copy;   # BONE(x,y) = full(x+ox, y+oy)
    my $tipstr = join ';', map { sprintf 'V%d=%d,%d', $_ + 1, $tips[$_][0] + 2, $tips[$_][1] - 3 } 0 .. $#tips;
    my %opt = (%o, inion_row => 120, tips => $tipstr, box => '15,8', smooth => 4, mm_per_px => 0.5);
    my $r = atlas_levels(\%opt, $skin, $bone);
    $chk->($r->{dx} == $ox && $r->{dy} == $oy, sprintf 'shift (%d, %d), expected (%d, %d); peak %.1f x the next', $r->{dx}, $r->{dy}, $ox, $oy, $r->{peak});
    my $terr = pdl(map { my $l = $r->{levels}[$_]; sqrt((($l->{tip} - pdl(@{ $tips[$_] }))**2)->sum) } 0 .. $#tips);
    $chk->($terr->max < 0.6, sprintf 'tips found from guesses 2-3 px off (max error %.2f px)', $terr->max);
    $chk->(abs($r->{inion}->at(0) - 280) <= 1, sprintf 'inion bone point at x %.0f (expected 280)', $r->{inion}->at(0));
    # exact lengths along x(y) from the inion row 120 to the row of each tip's nearest outline point
    my $yy = sequence(70001) * 0.01;                              # 0 .. 700
    my $dxy = 30 * 2 * $PI / 400 * cos(2 * $PI * $yy / 400);
    my $cl = pdl(0)->append((sqrt(1 + $dxy**2)->slice('1:-1') * 0.01)->cumusumover);
    my $xe = 300 + 30 * sin(2 * $PI * $yy / 400);
    my $len_at = sub {                                            # nearest point of the exact outline to (x, y)
        my ($x, $y) = @_;
        my $k = (($xe - $x)**2 + ($yy - $y)**2)->minimum_ind->sclr;
        return $cl->at($k);
    };
    # the outline itself: its length between rows 150 and 550 (the pixel edge is floor(x), a
    # staircase; the row average takes most of the steps out)
    my $P   = $r->{P};
    my $cum = arc_length($P);
    my $sel = which(($P->slice(':,(1)') >= 150) & ($P->slice(':,(1)') <= 550));
    my $lo  = $cum->index($sel)->max - $cum->index($sel)->min;
    my $le  = $cl->at(55000) - $cl->at(15000);
    $chk->(abs($lo / $le - 1) < 0.005, sprintf 'outline length rows 150..550: %.1f px (exact %.1f, %+.2f %%)', $lo, $le, 100 * ($lo / $le - 1));
    # the tips are 40-60 px inside the outline; their nearest outline point moves along a curved
    # outline with small changes of it, so these lengths are looser than the outline itself
    my $s0 = $len_at->(280, 120);
    my $lerr = pdl(map { $r->{levels}[$_]{s} - ($len_at->(@{ $tips[$_] }) - $s0) } 0 .. $#tips);
    $chk->($lerr->abs->max < 4, sprintf 'tip lengths from the inion within 4 px of those on the exact outline (max %.2f px over %.0f .. %.0f px)',
        $lerr->abs->max, $r->{levels}[0]{s}, $r->{levels}[-1]{s});
    my %o2 = (%opt, mm_per_px => undef, c7_mm => 100, tips => $tipstr =~ s/V3=/C7=/r);
    my $r2 = atlas_levels(\%o2, $skin, $bone);
    $chk->(abs($r2->{levels}[2]{s} * $r2->{mm} - 100) < 1e-9, '--c7-mm: the level named C7 is at that length');
    my $mid = $r->{mm} * 0.5 * ($r->{levels}[0]{s} + $r->{levels}[1]{s});
    my $lv  = level_of($r, $mid);
    $chk->($lv eq 'V1-V2', "a length between V1 and V2 is reported as V1-V2 ($lv)");
    $lv = level_of($r, $r->{mm} * $r->{levels}[3]{s} + 1);
    $chk->($lv eq 'V4', "a length 1 mm from V4 is reported as V4 ($lv)");
    print $ok ? "selftest passed\n" : "selftest FAILED\n";
    return $ok;
}
