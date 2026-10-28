package PDL::EEG::ElectrodeLocation::TraceLine;
# ============================================================================
# PDL::EEG::ElectrodeLocation::TraceLine — 画像に描いた線(黒地に白い線など)から、端から端までの
# 順序付きの中心線(点の列)を取り出す。
#
#   1. 2 値化(しきい値)
#   2. Zhang-Suen の細線化。近傍 8 個を slice でずらした配列の論理演算で、画像全体を一度に
#      処理する(画素ごとの Perl のループを使わない)
#   3. 骨格の画素を 8 近傍でつないだグラフで、最も画素の多い連結成分の中の最長の測地線
#      (二度掃き: 任意の点から最も遠い点 A、A から最も遠い点 B)。小さな枝(ひげ)は自然に外れる
#   4. 向きをそろえる(start => 'left' なら x の小さい端から)
#   5. 移動平均でならす(端点は動かさない)。8 近傍の画素の鎖の長さ(縦横 1、斜め √2)は、
#      斜めに近い曲線で階段を数えるので長めに出る。ならした線で長さを測る
#
# 座標は画像の座標(x 右、y 下、行 0 が上)。マスクは (nx, ny)。
# 由来: PDL::IO::NYHead05(ChatGPT の trace_midline2.pl の手順を PDL で書き直したもの)。
# ============================================================================
use strict;
use warnings;
use Carp;
use PDL;
use Exporter 'import';
use PDL::EEG::ElectrodeLocation::Polyline qw(arc_length);

our @EXPORT_OK = qw(binarize thin_zhang_suen skeleton_endpoints centerline);
our $VERSION   = '0.01';

# binarize($img, threshold => T) -> (nx,ny) byte
#   $img: (nx,ny) の濃淡、または (c,nx,ny)(c = 3 / 4。RGB の平均、アルファは使わない)。
#   T の既定は (最小 + 最大) / 2。invert => 1 で暗い線(白地に黒い線)
sub binarize {
    my ($img, %opt) = @_;
    my $g = $img->ndims == 3 ? $img->slice('0:' . (($img->dim(0) >= 3 ? 3 : 1) - 1) . ',:,:')->double->average : $img->double;
    my $t = $opt{threshold} // ($g->min + $g->max) / 2;
    return byte($opt{invert} ? $g < $t : $g > $t);
}

# thin_zhang_suen($mask (nx,ny)) -> (skeleton (nx,ny) byte, iterations)
#   Only the line pixels are visited: the padded image is flattened, the line pixels are kept as
#   flat indices, and their 8 neighbours are read with fixed index offsets, so one iteration costs
#   the number of line pixels, not the size of the image.
sub thin_zhang_suen {
    my ($mask, %opt) = @_;
    my ($nx, $ny) = $mask->dims;
    my $w2 = $nx + 2;
    my $p = zeroes(byte, $w2, $ny + 2);                          # one pixel of padding
    $p->slice('1:-2,1:-2') .= $mask != 0;
    my $flat = $p->flat;                                         # a view: writes go into $p
    my $fg = which($flat);                                       # the line pixels
    # P2..P9 = N NE E SE S SW W NW (y down: N is the row above)
    my @off = (-$w2, -$w2 + 1, 1, $w2 + 1, $w2, $w2 - 1, -1, -$w2 - 1);
    my $it = 0;
    while (1) {
        my $removed = 0;
        for my $step (0, 1) {
            last unless $fg->nelem;
            my @n = map { $flat->index($fg + $_)->long } @off;
            my $B = zeroes(long, $fg->nelem); $B += $_ for @n;
            my $A = zeroes(long, $fg->nelem); $A += (!$n[$_]) * $n[ ($_ + 1) % 8 ] for 0 .. 7;
            my ($P2, $P4, $P6, $P8) = @n[0, 2, 4, 6];
            my $del = ($B >= 2) * ($B <= 6) * ($A == 1);
            $del *= $step == 0 ? (($P2 * $P4 * $P6) == 0) * (($P4 * $P6 * $P8) == 0)
                               : (($P2 * $P4 * $P8) == 0) * (($P2 * $P6 * $P8) == 0);
            my $d = which($del);
            next unless $d->nelem;
            $flat->index($fg->index($d)) .= 0;
            $fg = $fg->where(!$del);
            $removed += $d->nelem;
        }
        $it++;
        last unless $removed;
        croak "thin_zhang_suen: did not converge" if $it > ($opt{max_iter} // 10000);
    }
    return ($p->slice('1:-2,1:-2')->copy, $it);
}

# skeleton_endpoints($skel) -> (k,2) long: 8 近傍がちょうど 1 個の画素 (x, y)
sub skeleton_endpoints {
    my ($skel) = @_;
    my ($nx, $ny) = $skel->dims;
    my $p = zeroes(long, $nx + 2, $ny + 2);
    $p->slice('1:-2,1:-2') .= $skel != 0;
    my $nb = zeroes(long, $nx, $ny);
    $nb += $p->slice($_) for ('1:-2,0:-3', '2:-1,0:-3', '2:-1,1:-2', '2:-1,2:-1', '1:-2,2:-1', '0:-3,2:-1', '0:-3,1:-2', '0:-3,0:-3');
    my $w = whichND(($skel != 0) * ($nb == 1));                  # (2,k): x, y
    return $w->nelem ? $w->transpose : zeroes(long, 0, 2);
}

# _graph: skeleton pixels -> index, coordinates, 8-neighbour lists
sub _graph {
    my ($skel) = @_;
    my ($nx, $ny) = $skel->dims;
    my $w = whichND($skel != 0);                                 # (2,n)
    my $n = $w->dim(1);
    croak "centerline: no line found (empty mask after thresholding)" unless $n;
    my @x = $w->slice('(0),')->list; my @y = $w->slice('(1),')->list;
    my %id; $id{ $y[$_] * $nx + $x[$_] } = $_ for 0 .. $n - 1;
    my @adj;
    for my $i (0 .. $n - 1) {
        for my $d ([-1, -1], [0, -1], [1, -1], [-1, 0], [1, 0], [-1, 1], [0, 1], [1, 1]) {
            my ($xx, $yy) = ($x[$i] + $d->[0], $y[$i] + $d->[1]);
            next if $xx < 0 || $yy < 0 || $xx >= $nx || $yy >= $ny;
            my $j = $id{ $yy * $nx + $xx };
            push @{ $adj[$i] }, [ $j, ($d->[0] && $d->[1]) ? sqrt(2) : 1 ] if defined $j;
        }
    }
    return (\@x, \@y, \@adj, $n);
}

# _dijkstra($adj, $n, $src) -> (\@dist, \@prev): binary heap
sub _dijkstra {
    my ($adj, $n, $src) = @_;
    my @dist = (9**9**9) x $n; my @prev = (-1) x $n;
    $dist[$src] = 0;
    my @h = ([0, $src]);
    while (@h) {
        my $top = $h[0];
        my $last = pop @h;
        if (@h) {                                                # sift down
            $h[0] = $last;
            my $i = 0;
            while (1) {
                my ($l, $r, $m) = (2 * $i + 1, 2 * $i + 2, $i);
                $m = $l if $l < @h && $h[$l][0] < $h[$m][0];
                $m = $r if $r < @h && $h[$r][0] < $h[$m][0];
                last if $m == $i;
                @h[ $i, $m ] = @h[ $m, $i ]; $i = $m;
            }
        }
        my ($du, $u) = @$top;
        next if $du > $dist[$u];
        for my $e (@{ $adj->[$u] || [] }) {
            my ($v, $w) = @$e;
            my $alt = $du + $w;
            next unless $alt < $dist[$v];
            $dist[$v] = $alt; $prev[$v] = $u;
            push @h, [ $alt, $v ];
            my $i = $#h;                                         # sift up
            while ($i > 0) { my $pa = int(($i - 1) / 2); last if $h[$pa][0] <= $h[$i][0]; @h[ $i, $pa ] = @h[ $pa, $i ]; $i = $pa }
        }
    }
    return (\@dist, \@prev);
}

# centerline($mask, %opt) -> {
#     raw => (n,2) the pixel chain (x, y), smooth => (n,2) after the moving average,
#     length_raw, length (of the smoothed line), skeleton, endpoints (k,2), n_components, iterations }
#   start  => 'left' | 'right' | 'top' | 'bottom'  which end comes first                   ('left')
#   smooth => W  moving average over W points (odd; 1 = none)                              (5)
sub centerline {
    my ($mask, %opt) = @_;
    my ($skel, $it) = thin_zhang_suen($mask);
    my ($x, $y, $adj, $n) = _graph($skel);
    # the largest connected component
    my @comp = (-1) x $n; my ($nc, $best, $bsize) = (0, 0, 0);
    for my $s (0 .. $n - 1) {
        next if $comp[$s] >= 0;
        my @q = ($s); $comp[$s] = $nc; my $size = 0;
        while (@q) { my $u = shift @q; $size++; for my $e (@{ $adj->[$u] || [] }) { next if $comp[ $e->[0] ] >= 0; $comp[ $e->[0] ] = $nc; push @q, $e->[0] } }
        ($best, $bsize) = ($nc, $size) if $size > $bsize;
        $nc++;
    }
    my ($s0) = grep { $comp[$_] == $best } 0 .. $n - 1;
    my ($d0) = _dijkstra($adj, $n, $s0);
    my @in = grep { $comp[$_] == $best } 0 .. $n - 1;
    my $end1 = $in[ pdl(@{$d0}[@in])->maximum_ind->sclr ];      # farthest from any point: one end
    my ($d1, $pr) = _dijkstra($adj, $n, $end1);
    my $end2 = $in[ pdl(@{$d1}[@in])->maximum_ind->sclr ];      # farthest from that end: the other end
    my @path = ($end2);
    push @path, $pr->[ $path[-1] ] while $path[-1] != $end1;
    @path = reverse @path;                                       # end1 .. end2
    my $start = $opt{start} // 'left';
    my ($f, $l) = ($path[0], $path[-1]);
    my $flip = $start eq 'left'   ? $x->[$f] > $x->[$l]
             : $start eq 'right'  ? $x->[$f] < $x->[$l]
             : $start eq 'top'    ? $y->[$f] > $y->[$l]
             : $start eq 'bottom' ? $y->[$f] < $y->[$l]
             : croak "centerline: start must be left, right, top or bottom";
    @path = reverse @path if $flip;
    my $raw = pdl([ map { [ $x->[$_], $y->[$_] ] } @path ])->transpose->double->copy;   # (n,2)
    my $sm = _smooth($raw, $opt{smooth} // 5);
    return { raw => $raw, smooth => $sm, length_raw => arc_length($raw)->at(-1), length => arc_length($sm)->at(-1),
             skeleton => $skel, endpoints => skeleton_endpoints($skel), n_components => $nc, iterations => $it };
}

# _smooth($P (n,2), W): centred moving average, the window shrinking near the ends (end points fixed)
sub _smooth {
    my ($P, $w) = @_;
    my $n = $P->dim(0);
    return $P->copy if $w <= 1 || $n < 3;
    my $h  = int($w / 2);
    my $i  = sequence(long, $n);
    my $hh = cat($i, $n - 1 - $i, pdl(long, $h) * ones(long, $n))->mv(1, 0)->minimum;   # (n): the half window at each point
    my $cs = zeroes(double, $n + 1, 2);
    $cs->slice('1:-1,:') .= $P->cumusumover;                     # (n+1,2): cumulative sums along the points (dim 0)
    my $lo = $i - $hh; my $hi = $i + $hh + 1;
    return ($cs->dice_axis(0, $hi) - $cs->dice_axis(0, $lo)) / ($hi - $lo)->dummy(1);
}

1;

__END__

=head1 NAME

PDL::EEG::ElectrodeLocation::TraceLine - the ordered centre line of a line drawn on an image

=head1 SYNOPSIS

  use PDL::IO::PNG;
  use PDL::EEG::ElectrodeLocation::TraceLine qw(binarize centerline);
  use PDL::EEG::ElectrodeLocation::Polyline qw(midline_labels point_at_fraction);
  my $mask = binarize(PDL::IO::PNG::rpng('midline.png'));       # (nx,ny), row 0 at the top
  my $cl   = centerline($mask, start => 'left', smooth => 5);
  my ($lab, $f) = midline_labels('10-20');
  my $xy   = point_at_fraction($cl->{smooth}, $f);              # (7,2) image coordinates

=head1 DESCRIPTION

Zhang-Suen thinning on whole arrays, the longest geodesic path of the largest skeleton
component (side branches drop out), oriented from the chosen end, smoothed by a moving average
(the end points stay). The pixel chain over-counts the length of slanted curves (unit and sqrt(2)
steps); length the smoothed line.

=cut
