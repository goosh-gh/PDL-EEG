package PDL::EEG::ElectrodeLocation::AtlasView;
# ============================================================================
# PDL::EEG::ElectrodeLocation::AtlasView — side views of an anatomy atlas (screenshots of a 3-D
# model with the face to the left, e.g. Open Anatomy Atlas / BodyParts3D): read a PNG, the body
# against the background, the shift between two views taken with the same camera (phase
# correlation), the back outline of the body, the outline of the whole body, the tips of the
# spinous processes and the inion on a bone view, and their lengths along the back outline.
#
# Image arrays are (3,nx,ny) double 0 .. 1 with row 0 at the top (x right, y down); masks are
# (nx,ny) byte; outlines and point sets are (n,2) as (x, y) in pixels.
#
# 由来: PDL::IO::NYHead06(アトラスの体表の線から頸部電極の椎骨の高さを読む、NYHead との重ね合わせ)
# ============================================================================
use strict;
use warnings;
use Carp;
use PDL;
use PDL::EEG::FFT::Dispatcher qw(fft_inplace ifft_inplace);
use PDL::EEG::ElectrodeLocation::Polyline qw(arc_length point_at_length project_point);
use Exporter 'import';

our @EXPORT_OK = qw(read_png foreground shift_between back_outline body_outline refine_tip atlas_levels level_of);
our $VERSION   = '0.01';

sub read_png {
    my ($f) = @_;
    -f $f or die "not found: $f\n";
    my $img;
    if (eval { require PDL::IO::PNG; 1 }) { $img = PDL::IO::PNG::rpng($f) }   # (c,nx,ny) float 0..1, row 0 at the top
    else {
        require PDL::IO::Pic;
        my $p = PDL::IO::Pic::rpic($f);                           # (c,nx,ny), row 0 at the bottom
        $img = ($p->slice(':,:,-1:0')->float / ($p->max > 1 ? 255 : 1))->copy;
    }
    return $img->slice('0:2,:,:')->double->copy;
}

# foreground: differs from the background colour (the top-left pixel) by more than 0.05
# summed over r, g, b.  Returns (nx,ny) byte.
sub foreground {
    my ($img) = @_;
    my $bg = $img->slice(':,(0),(0)');
    return (abs($img - $bg)->sumover > 0.05)->byte;
}

# fft2($re, $im): 2-D transform of (nx,ny) in place of copies (the 1-D transform runs along
# dim 0 and broadcasts over dim 1; the two passes are separated by a transpose)
sub fft2 {
    my ($re, $im, $inv) = @_;
    ($re, $im) = ($re->copy, $im->copy);
    $inv ? ifft_inplace($re, $im) : fft_inplace($re, $im);
    my ($rt, $it) = ($re->xchg(0, 1)->copy, $im->xchg(0, 1)->copy);
    $inv ? ifft_inplace($rt, $it) : fft_inplace($rt, $it);
    return ($rt->xchg(0, 1)->copy, $it->xchg(0, 1)->copy);
}

# shift_between($a (nx,ny), $b (mx,my)): integer (dx, dy) with a(x, y) ~ b(x - dx, y - dy), by
# phase correlation of the two grey images zero-padded to (nx+mx, ny+my).  Also returns the
# peak height over the next highest value outside 5 px of the peak.
sub shift_between {
    my ($ga, $gb) = @_;
    my ($W, $H) = ($ga->dim(0) + $gb->dim(0), $ga->dim(1) + $gb->dim(1));
    my ($pa, $pb) = (zeroes($W, $H), zeroes($W, $H));
    $pa->slice('0:' . ($ga->dim(0) - 1) . ',0:' . ($ga->dim(1) - 1)) .= $ga - $ga->avg;
    $pb->slice('0:' . ($gb->dim(0) - 1) . ',0:' . ($gb->dim(1) - 1)) .= $gb - $gb->avg;
    my ($ar, $ai) = fft2($pa, zeroes($W, $H));
    my ($br, $bi) = fft2($pb, zeroes($W, $H));
    my $cr = $ar * $br + $ai * $bi;                               # A conj(B)
    my $ci = $ai * $br - $ar * $bi;
    my $m  = sqrt($cr**2 + $ci**2)->clip(1e-12, undef);
    my ($rr) = fft2($cr / $m, $ci / $m, 1);
    my $k  = $rr->flat->maximum_ind->sclr;
    my ($dx, $dy) = ($k % $W, int($k / $W));
    my $peak = $rr->flat->at($k);
    my $mask = ((abs(xvals($W, $H) - $dx) <= 5) & (abs(yvals($W, $H) - $dy) <= 5));
    my $next = ($rr * (1 - $mask) + $mask * $rr->min)->max;
    $dx -= $W if $dx > $W / 2;
    $dy -= $H if $dy > $H / 2;
    return ($dx, $dy, $peak / $next);
}

# back_outline($fg (nx,ny), $h): the most posterior (largest x) foreground pixel of every row,
# averaged over 2h+1 rows.  Returns the outline (n,2) as (x, y) for the rows that have a full
# window of foreground rows, top first.
sub back_outline {
    my ($fg, $h) = @_;
    my ($nx, $ny) = $fg->dims;
    my $xr = ($fg * (xvals(long, $nx, $ny) + 1))->maximum - 1;    # (ny): last foreground column, -1 if none
    my $ok = ($xr >= 0);
    my $cs = pdl(0)->append($xr->double->cumusumover);            # windowed sums by differences of the cumulative sum
    my $cn = pdl(0)->append($ok->long->cumusumover);
    my $y  = sequence(long, $ny - 2 * $h) + $h;
    my $n  = $cn->index($y + $h + 1) - $cn->index($y - $h);
    my $sel = which($n == 2 * $h + 1);
    $y = $y->index($sel);
    my $xs = ($cs->index($y + $h + 1) - $cs->index($y - $h)) / (2 * $h + 1);
    return cat($xs, $y->double);                                  # (n,2): cat stacks along a new last dim
}

# refine_tip($fg (nx,ny), $x, $y, $bx, $by): the most posterior foreground pixel in the box
# x-bx .. x+bx, y-by .. y+by (the mean row of the pixels in the most posterior column)
sub refine_tip {
    my ($fg, $x, $y, $bx, $by) = @_;
    my ($nx, $ny) = $fg->dims;
    my ($x0, $x1, $y0, $y1) = (max2(0, $x - $bx), min2($nx - 1, $x + $bx), max2(0, $y - $by), min2($ny - 1, $y + $by));
    my $box = $fg->slice("$x0:$x1,$y0:$y1");
    $box->max > 0 or return;
    my $xv = ($box * (xvals(long, $box->dims) + 1))->max - 1;     # most posterior column in the box
    my $rows = which($box->slice("($xv),"));
    return ($x0 + $xv, $y0 + $rows->avg);
}
sub max2 { $_[0] > $_[1] ? $_[0] : $_[1] }
sub min2 { $_[0] < $_[1] ? $_[0] : $_[1] }

sub atlas_levels {
    my ($opt, $skin, $bone) = @_;
    my ($dx, $dy, $pk) = shift_between($skin->average, $bone->average);   # grey: mean over r, g, b (dim 0)
    my $fs = foreground($skin);
    my $fb = foreground($bone);
    my $P  = back_outline($fs, $opt->{smooth});                   # (n,2) in SKIN pixels
    my $cum = arc_length($P);
    # the inion: the most posterior bone pixel of --inion-row, on the outline
    my $yb = $opt->{inion_row} - $dy;
    $yb >= 0 && $yb < $fb->dim(1) or die "--inion-row $opt->{inion_row} is outside the bone view (shift $dx, $dy)\n";
    my $row = $fb->slice(":,($yb)");
    $row->max > 0 or die "--inion-row $opt->{inion_row}: no bone in that row\n";
    my $inion = pdl(which($row)->max + $dx, $opt->{inion_row});
    my ($sI, $dI) = project_point($P, $inion, $cum);
    my ($bx, $by) = split /,/, $opt->{box};
    my @lev;
    for my $kv (split /;/, $opt->{tips}) {
        my ($n, $x, $y) = $kv =~ /^\s*([^=]+)=([\d.]+),([\d.]+)\s*$/ or die "--tips NAME=X,Y;..., not '$kv'\n";
        my ($tx, $ty) = refine_tip($fb, int($x - $dx + 0.5), int($y - $dy + 0.5), $bx, $by);
        defined $tx or do { print "  $n: no bone within the box around ($x, $y); left out\n"; next };
        my $t = pdl($tx + $dx, $ty + $dy);
        my ($s, $d) = project_point($P, $t, $cum);
        push @lev, { name => $n, guess => [$x, $y], tip => $t, s => $s - $sI, depth => $d };
    }
    my %at = map { $lev[$_]{name} => $_ } 0 .. $#lev;
    my $mm = $opt->{mm_per_px};
    if (!defined $mm) {
        exists $at{C7} or die "--c7-mm needs a level named C7 in --tips (or give --mm-per-px)\n";
        $mm = $opt->{c7_mm} / $lev[ $at{C7} ]{s};
    }
    my $pI = (point_at_length($P, pdl($sI), $cum))->slice('(0),');
    return { dx => $dx, dy => $dy, peak => $pk, P => $P, cum => $cum, sI => $sI, inion => $inion, inion_skin => $pI,
             inion_depth => $dI, levels => \@lev, mm => $mm, by_scale => defined $opt->{mm_per_px} ? 'given' : "C7 = $opt->{c7_mm} mm" };
}

# the level interval of a length (mm from the inion)
sub level_of {
    my ($r, $L) = @_;
    my @lv = sort { $a->{s} <=> $b->{s} } @{ $r->{levels} };
    my @s  = map { $_->{s} * $r->{mm} } @lv;
    return "above $lv[0]{name}" if $L < $s[0] - 2;
    for my $i (0 .. $#lv) {
        return $lv[$i]{name} if abs($L - $s[$i]) <= 2;
        return "$lv[$i]{name}-$lv[$i + 1]{name}" if $i < $#lv && $L < $s[$i + 1] - 2;
    }
    return "below $lv[-1]{name}";
}

# body_outline($fg (nx,ny)): the boundary pixels of the body (foreground pixels with a
# background 4-neighbour, or on the image edge are left out) as (n,2) (x, y), not ordered
sub body_outline {
    my ($fg) = @_;
    my ($nx, $ny) = $fg->dims;
    my $c  = $fg->slice('1:-2,1:-2');
    my $in = $fg->slice('0:-3,1:-2') & $fg->slice('2:-1,1:-2') & $fg->slice('1:-2,0:-3') & $fg->slice('1:-2,2:-1');
    my $b  = whichND($c & !$in) + 1;                              # (2,n): x, y of the inner frame
    return $b->xchg(0, 1)->double->copy;                          # (n,2)
}

1;
