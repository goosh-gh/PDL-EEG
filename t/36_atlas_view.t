use strict;
use warnings;
use Test::More;
use PDL;
use File::Temp qw(tempdir);
use PDL::EEG::ElectrodeLocation::AtlasView qw(read_png foreground shift_between back_outline body_outline refine_tip atlas_levels level_of);

# Synthetic side views with known answers (the background colour of the Open Anatomy Atlas views,
# a body whose back edge is a known line, bars ending at known tips).
my $bg = pdl(30, 37, 46) / 255;
my $paint = sub {                       # (3,nx,ny) image: background, then masks painted in colours
    my ($nx, $ny, @mc) = @_;
    my $img = zeroes(3, $nx, $ny) + $bg;
    while (my ($m, $c) = splice(@mc, 0, 2)) {
        $img = $img * (1 - $m->dummy(0)) + $m->dummy(0) * pdl(@$c)->dummy(1)->dummy(2);
    }
    return $img;
};
my ($nx, $ny) = (160, 220);
my $X = xvals(double, $nx, $ny); my $Y = yvals(double, $nx, $ny);

# ---- foreground and body_outline: a filled rectangle -----------------------------------------
my $rect = ($X >= 40) & ($X <= 59) & ($Y >= 100) & ($Y <= 109);   # 20 x 10
my $fg = foreground($paint->($nx, $ny, $rect, [0.8, 0.6, 0.5]));
is($fg->sum, 200, 'foreground: the 200 pixels of a 20 x 10 rectangle');
my $bo = body_outline($fg);
is($bo->dim(0), 2 * (20 + 10) - 4, 'body_outline: 56 boundary pixels of the rectangle');
my ($bx, $by) = ($bo->slice(':,(0)'), $bo->slice(':,(1)'));
ok(all((($bx == 40) | ($bx == 59) | ($by == 100) | ($by == 109))), 'body_outline: every pixel on the rectangle edge');

# ---- back_outline: back edge x = 100 + y/4 (rows 20 .. 199) --------------------------------------
my $edge = 100 + $Y / 4;
my $body = ($X <= $edge) & ($X >= 20) & ($Y >= 20) & ($Y <= 199);
my $P = back_outline(foreground($paint->($nx, $ny, $body, [0.8, 0.6, 0.5])), 2);
is($P->dim(1), 2, 'back_outline: (n,2)');
is($P->at(0, 1), 22, 'back_outline: the first row with a full window of 5 rows (h = 2)');
is($P->at(-1, 1), 197, 'back_outline: the last such row');
my $want = (100 + $P->slice(':,(1)') / 4)->floor;          # the last body column of a row is floor(edge)
# averaged over 5 rows the floor steps average out: within 0.5 px of the mean of floor over the window
my $wavg = pdl(map { my $y = $_; (100 + (sequence(5) + $y - 2) / 4)->floor->avg } $P->slice(':,(1)')->list);
ok(($P->slice(':,(0)') - $wavg)->abs->max < 1e-9, 'back_outline: x = the mean of the last body column over 2h+1 rows');

# ---- refine_tip: a bar ending at (90, 60) -------------------------------------------------------
my $bar = ($X >= 50) & ($X <= 90) & (abs($Y - 60) <= 3);
my $fb = foreground($paint->($nx, $ny, $bar, [0.9, 0.9, 0.85]));
my ($tx, $ty) = refine_tip($fb, 85, 64, 10, 6);
ok($tx == 90 && abs($ty - 60) < 1e-9, "refine_tip: from a guess (85, 64) to the tip (90, 60): ($tx, $ty)");
ok(!defined((refine_tip($fb, 140, 150, 5, 5))[0]), 'refine_tip: undef when there is no bone in the box');

# ---- shift_between: a textured image and a crop of it ------------------------------------------------
my $tex = (($X * 7919 + $Y * 104729 + $X * $Y * 31) % 1009) / 1009;   # not periodic (a periodic texture gives side peaks)
my ($ox, $oy) = (13, 7);
my $crop = $tex->slice("$ox:-1,$oy:-1")->copy;                # crop(x, y) = tex(x + ox, y + oy)
my ($dx, $dy, $pk) = shift_between($tex, $crop);
ok($dx == $ox && $dy == $oy, "shift_between: ($dx, $dy), expected ($ox, $oy)");
ok($pk > 3, sprintf 'shift_between: a distinct peak (%.1f x the next)', $pk);

# ---- atlas_levels and level_of: skin view + bone view (cropped by a known shift) ------------
my $skinb = ($X <= $edge) & ($X >= 20) & ($Y >= 10) & ($Y <= 215);
my @tips  = ([96, 70], [110, 120], [120, 170]);                # x < edge everywhere (edge 117.5 .. 142.5)
my $bones = zeroes(byte, $nx, $ny);
$bones |= (($X <= $_->[0]) & ($X >= $_->[0] - 30) & (abs($Y - $_->[1]) <= 2)) for @tips;
$bones |= (($X - 60)**2 + ($Y - 40)**2 <= 30**2) & ($X <= 85);    # a "skull", back at x 85, row 40
my $skin = $paint->($nx, $ny, $skinb, [0.8, 0.65, 0.55], $bones, [0.9, 0.9, 0.85]);
my $bonev = $paint->($nx, $ny, $bones, [0.9, 0.9, 0.85])->slice(':,5:-1,3:-1')->copy;   # shift (5, 3)
my %opt = (inion_row => 40, tips => join(';', map { sprintf 'L%d=%d,%d', $_ + 1, $tips[$_][0] - 2, $tips[$_][1] + 2 } 0 .. 2),
           box => '8,5', smooth => 4, c7_mm => 100, mm_per_px => 0.5);
my $r = atlas_levels(\%opt, $skin, $bonev);
is("$r->{dx},$r->{dy}", '5,3', 'atlas_levels: the shift of the bone view');
is(scalar @{ $r->{levels} }, 3, 'atlas_levels: three levels');
ok(abs($r->{inion}->at(0) - 85) < 1e-9, 'atlas_levels: inion bone point at the back of the skull row');
# the length along the back line x = 100 + y/4 (slope 1/4: ds = sqrt(1 + 1/16) dy) between the
# nearest outline points of the inion and of a tip; nearest point of a line from (x0, y0):
# y* = (y0 + (x0 - 100)/4) / (1 + 1/16)
my $ys = sub { my ($x0, $y0) = @_; ($y0 + ($x0 - 100) / 4) / (1 + 1 / 16) };
my $yI = $ys->(85, 40);
my @exp = map { ($ys->(@$_) - $yI) * sqrt(1 + 1 / 16) } @tips;
my $err = pdl(map { $r->{levels}[$_]{s} - $exp[$_] } 0 .. 2)->abs->max;
# the outline is a staircase of whole pixels (averaged over 2h+1 rows); the nearest point of a
# staircase to a point 20-25 px inside moves by 1-2 px against that of the true line (h = 2: 2.4,
# h = 4: 0.9, h = 6: 1.6 px here), hence 2 px
ok($err < 2.0, sprintf 'atlas_levels: lengths from the inion along the back within 2 px of the line (max %.2f px)', $err);
is(level_of($r, 0.5 * 0.5 * ($r->{levels}[0]{s} + $r->{levels}[1]{s})), 'L1-L2', 'level_of: between L1 and L2');
is(level_of($r, 0.5 * $r->{levels}[2]{s} + 1), 'L3', 'level_of: within 2 mm of L3');
is(level_of($r, 0.5 * $r->{levels}[2]{s} + 10), 'below L3', 'level_of: below the last level');
my $r2 = atlas_levels({ %opt, mm_per_px => undef, tips => $opt{tips} =~ s/L2=/C7=/r }, $skin, $bonev);
ok(abs($r2->{levels}[1]{s} * $r2->{mm} - 100) < 1e-9, '--c7-mm: the level named C7 at 100 mm');

# ---- read_png: a written PNG read back (PDL::IO::PNG) -----------------------------------------
SKIP: {
    skip 'PDL::IO::PNG not installed', 2 unless eval { require PDL::IO::PNG; 1 };
    my $dir = tempdir(CLEANUP => 1);
    my $small = $skin->slice(':,0:39,0:29')->copy;            # (3,40,30)
    PDL::IO::PNG::wpng($small->float, "$dir/s.png");
    my $back = read_png("$dir/s.png");
    is(join(',', $back->dims), '3,40,30', 'read_png: (3,nx,ny)');
    ok(($back - $small)->abs->max < 1 / 255 + 1e-6, 'read_png: the pixels back to 8 bits');
}

done_testing;
