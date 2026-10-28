use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::ElectrodeLocation::Polyline qw(arc_length point_at_length point_at_fraction project_point midline_labels);
use PDL::EEG::ElectrodeLocation::TraceLine qw(binarize thin_zhang_suen skeleton_endpoints centerline);

my $PI = 3.14159265358979;

# ---- Polyline: a 3-D polyline with a repeated point (a zero-length segment) -------------------
my $P = pdl([[0, 0, 0], [3, 4, 0], [3, 4, 0], [3, 4, 12]])->transpose->copy;   # (4,3): lengths 5, 0, 12
is(join(',', $P->dims), '4,3', 'polyline (n,d)');
my $cum = arc_length($P);
ok(($cum - pdl(0, 5, 5, 17))->abs->max < 1e-12, 'arc_length 0 5 5 17');
my $X = point_at_length($P, pdl(0, 2.5, 5, 11, 17, 99), $cum);
my $want = pdl([[0, 0, 0], [1.5, 2, 0], [3, 4, 0], [3, 4, 6], [3, 4, 12], [3, 4, 12]])->transpose;
ok(($X - $want)->abs->max < 1e-12, 'point_at_length: interpolated, past the zero-length segment, clipped at the end');
my $Xf = point_at_fraction($P, pdl(0.5));
ok(($Xf->slice('(0),') - pdl(3, 4, 3.5))->abs->max < 1e-12, 'point_at_fraction 0.5 (8.5 of 17)');
my ($s, $off) = project_point($P, pdl(4, 4, 6), $cum);
ok(abs($s - 11) < 1e-12 && abs($off - 1) < 1e-12, 'project_point: 11 along, 1 off');
my ($l20, $f20) = midline_labels('10-20');
is("@$l20", 'Nasion Fpz Fz Cz Pz Oz Inion', '10-20 labels');
ok(($f20 - pdl(0, .1, .3, .5, .7, .9, 1))->abs->max < 1e-12, '10-20 fractions');
my ($l10, $f10) = midline_labels('10-10');
is(scalar @$l10, 11, '10-10: eleven');
eval { midline_labels('10-5') }; like($@, qr/10-20.*10-10/, 'unknown system croaks');

# ---- TraceLine: thinning of a thick bar keeps one connected line, end to end ------------
my $bar = zeroes(byte, 60, 20);
$bar->slice('5:54,7:12') .= 1;                                   # 50 x 6
my ($sk, $it) = thin_zhang_suen($bar);
ok($sk->sum > 30 && $sk->sum < 60, sprintf 'bar thinned to %d pixels in %d iterations', $sk->sum, $it);
ok(all($sk <= $bar), 'skeleton inside the shape');
my $ep = skeleton_endpoints($sk);
is($ep->dim(0), 2, 'two end points');
# the same thinning as an explicit per-pixel Zhang-Suen (the published rules written out)
my $ref = _zs_reference($bar);
ok(all($ref == $sk), 'vectorised thinning = per-pixel reference');
my $ring = do { my $x = xvals(80, 80) - 40; my $y = yvals(80, 80) - 40; my $r = sqrt($x**2 + $y**2); byte((abs($r - 30) <= 2) * ($y < 20)) };
ok(all(_zs_reference($ring) == (thin_zhang_suen($ring))[0]), 'vectorised thinning = reference on an arc');

# ---- centerline on a semicircle: length, order, the spur ignored ------------------------
my ($nx, $ny, $R) = (420, 260, 150);
my $x = xvals(double, $nx, $ny) - 210; my $y = yvals(double, $nx, $ny) - 220;
my $r = sqrt($x**2 + $y**2);
my $line = byte((abs($r - $R) <= 1.5) * ($y <= 0));
$line |= byte((abs($x) <= 1) * ($y >= -$R + 1) * ($y <= -$R + 12));   # a spur
my $cl = centerline($line, start => 'left');
ok(abs($cl->{length} / ($PI * $R) - 1) < 0.01, sprintf 'smoothed length %.1f (true %.1f)', $cl->{length}, $PI * $R);
ok($cl->{length_raw} / ($PI * $R) > 1.03, sprintf 'pixel chain %.1f is longer', $cl->{length_raw});
ok($cl->{raw}->at(0, 0) < $cl->{raw}->at(-1, 0), 'starts at the left end');
my $mid = point_at_fraction($cl->{smooth}, pdl(0.5));
ok(abs($mid->at(0, 0) - 210) < 1.5 && abs($mid->at(0, 1) - (220 - $R)) < 1.5, 'the half-way point is the top of the arc');
my $clb = centerline($line, start => 'right');
ok($clb->{raw}->at(0, 0) > $clb->{raw}->at(-1, 0), 'start right');
# ---- curves whose curvature is not even: an ellipse, and a head-like profile (a flat front,
#      a round back), drawn 3 px wide from a dense parametric curve; the truth is that curve
for my $shape (['ellipse', sub { my $t = shift; (200 * cos($t), -110 * sin($t)) }],
               ['head-like', sub { my $t = shift; (190 * cos($t) * (1 - 0.18 * sin($t)), -150 * sin($t) * (1 + 0.22 * cos($t))) }]) {
    my ($name, $f) = @$shape;
    my $t  = sequence(double, 8001) / 8000 * $PI;                 # 0 = front (+x, nasion side drawn on the right), pi = back
    my ($cx, $cy) = $f->($t);
    my $D = cat($cx + 250, $cy + 240)->copy;                      # (8001,2) image coordinates
    my $img = zeroes(byte, 500, 300);
    for my $dx (-1 .. 1) { for my $dy (-1 .. 1) {                 # 3 x 3 brush along the dense curve
        my $ix = ($D->slice(':,(0)') + $dx)->rint->long; my $iy = ($D->slice(':,(1)') + $dy)->rint->long;
        $img->indexND(cat($ix, $iy)->transpose) .= 1;
    } }
    my $truth = arc_length($D);
    my ($lab, $fr) = midline_labels('10-20');
    my $want  = point_at_fraction($D, $fr, $truth);               # (7,2)
    my $c = centerline($img, start => 'right');                   # the front (t = 0) is on the right
    my $ec = sqrt(((point_at_fraction($c->{raw}, $fr) - $want)**2)->mv(1, 0)->sumover);      # (7): over x, y (dim 1)
    my $es = sqrt(((point_at_fraction($c->{smooth}, $fr) - $want)**2)->mv(1, 0)->sumover);
    my $L  = $truth->at(-1);
    my ($ends, $inner) = ($es->index(pdl(long, 0, 6))->max->sclr, $es->slice('1:5')->max->sclr);
    my $inner_c = $ec->slice('1:5')->max->sclr;
    note(sprintf "%s: length %.1f px; pixel chain %+.2f %%, smoothed %+.2f %%; Fpz..Oz error: chain %.2f px, smoothed %.2f px; ends %.2f px",
         $name, $L, 100 * ($c->{length_raw} / $L - 1), 100 * ($c->{length} / $L - 1), $inner_c, $inner, $ends);
    ok(abs($c->{length} / $L - 1) < 0.01, sprintf '%s: smoothed length within 1 %% (%+.2f %%)', $name, 100 * ($c->{length} / $L - 1));
    ok($ends < 3, sprintf '%s: Nasion and Inion within 3 px of the ends of the drawn curve (%.2f px; thinning moves the ends)', $name, $ends);
    ok($inner < 0.004 * $L, sprintf '%s: Fpz..Oz within 0.4 %% of the length (%.2f px of %.1f)', $name, $inner, 0.004 * $L);
    if ($name eq 'head-like') {
        ok($inner < $inner_c, sprintf 'head-like: smoothing brings Fpz..Oz closer (%.2f px, pixel chain %.2f px)', $inner, $inner_c);
    }
}

# binarize: RGB, threshold, invert
my $rgb = cat(map { $line->float } 1 .. 3)->reorder(2, 0, 1);
ok(all(binarize($rgb) == $line), 'binarize RGB (c,nx,ny)');
ok(all(binarize(1 - $rgb, invert => 1) == $line), 'binarize --invert');
eval { centerline(zeroes(byte, 10, 10)) }; like($@, qr/no line/, 'empty mask croaks');

done_testing;

# Zhang & Suen (1984), written pixel by pixel from the paper: P2..P9 = N NE E SE S SW W NW
sub _zs_reference {
    my ($m) = @_;
    my ($w, $h) = $m->dims;
    my @p = map { [ (0) x ($w + 2) ] } 1 .. $h + 2;
    for my $j (0 .. $h - 1) { for my $i (0 .. $w - 1) { $p[ $j + 1 ][ $i + 1 ] = $m->at($i, $j) ? 1 : 0 } }
    while (1) {
        my $chg = 0;
        for my $step (0, 1) {
            my @del;
            for my $j (1 .. $h) {
                for my $i (1 .. $w) {
                    next unless $p[$j][$i];
                    my @n = ($p[ $j - 1 ][$i], $p[ $j - 1 ][ $i + 1 ], $p[$j][ $i + 1 ], $p[ $j + 1 ][ $i + 1 ],
                             $p[ $j + 1 ][$i], $p[ $j + 1 ][ $i - 1 ], $p[$j][ $i - 1 ], $p[ $j - 1 ][ $i - 1 ]);
                    my $B = 0; $B += $_ for @n;
                    next if $B < 2 || $B > 6;
                    my $A = grep { !$n[$_] && $n[ ($_ + 1) % 8 ] } 0 .. 7;
                    next unless $A == 1;
                    if ($step == 0) { next if $n[0] * $n[2] * $n[4] || $n[2] * $n[4] * $n[6] }
                    else            { next if $n[0] * $n[2] * $n[6] || $n[0] * $n[4] * $n[6] }
                    push @del, [ $j, $i ];
                }
            }
            $p[ $_->[0] ][ $_->[1] ] = 0 for @del;
            $chg += @del;
        }
        last unless $chg;
    }
    my $o = zeroes(byte, $w, $h);
    for my $j (0 .. $h - 1) { for my $i (0 .. $w - 1) { $o->set($i, $j, $p[ $j + 1 ][ $i + 1 ]) } }
    return $o;
}
