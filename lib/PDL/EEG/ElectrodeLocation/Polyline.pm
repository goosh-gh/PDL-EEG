package PDL::EEG::ElectrodeLocation::Polyline;
# ============================================================================
# PDL::EEG::ElectrodeLocation::Polyline — 点の列(折れ線)の弧長と、長さ・割合で点を置く、点を
# 折れ線に投影する。2 次元(画像の上の線)でも 3 次元(頭の面の上の線)でも同じ。
# 10-20 / 10-10 法の正中のラベルと割合も持つ。
#
# 次元の約束: 折れ線 $P は (n, d)(dim 0 = 点、dim 1 = 座標)。
#
# 由来: PDL::IO::NYHead05(画像の上の正中の線から 10-20 の位置を決める。
# PDL-IO-NYHead nyhead_neck_electrodes.pl の --arc、PDL-EEG nyhead_alpha_sim.pl の
# 正中の長さと同じ計算)。
# ============================================================================
use strict;
use warnings;
use Carp;
use PDL;
use Exporter 'import';

our @EXPORT_OK = qw(arc_length point_at_length point_at_fraction project_point midline_labels);
our $VERSION   = '0.01';

# arc_length($P (n,d)) -> (n): 始点からの弧長(各点まで)
sub arc_length {
    my ($P) = @_;
    $P->ndims == 2 && $P->dim(0) >= 2 or croak "arc_length: the polyline must be (n,d) with n >= 2, got (" . join(',', $P->dims) . ")";
    my $seg = sqrt((($P->slice('1:-1,:') - $P->slice('0:-2,:'))**2)->mv(1, 0)->sumover);   # (n-1): over the coordinates (dim 1)
    return pdl(0)->append($seg->cumusumover);
}

# point_at_length($P, $s (k), $cum) -> (k,d): 弧長 $s の点(線形補間、0 .. 全長に切る)
sub point_at_length {
    my ($P, $s, $cum) = @_;
    $cum //= arc_length($P);
    $s = pdl($s)->double->flat;
    my $L = $cum->at(-1);
    $L > 0 or croak "point_at_length: the polyline has zero length";
    $s = $s->clip(0, $L);
    my $n = $cum->nelem;
    # 区間の添字 k: cum(k) <= s < cum(k+1)。長さ 0 の区間は飛ばす(vsearch は右端の同値を選ぶ)
    my $k = vsearch_insert_rightmost($s, $cum) - 1;               # (k)
    $k = $k->clip(0, $n - 2);
    my $c0 = $cum->index($k); my $c1 = $cum->index($k + 1);
    my $t  = (($s - $c0) / ($c1 - $c0)->clip(1e-300, undef))->clip(0, 1);
    my $A  = $P->dice_axis(0, $k); my $B = $P->dice_axis(0, $k + 1);   # (k,d)
    return $A + $t->dummy(1) * ($B - $A);
}

# point_at_fraction($P, $f (k), $cum) -> (k,d): 全長の割合 $f(0 .. 1)の点
sub point_at_fraction {
    my ($P, $f, $cum) = @_;
    $cum //= arc_length($P);
    return point_at_length($P, pdl($f)->double * $cum->at(-1), $cum);
}

# project_point($P, $q (d), $cum) -> (弧長, 折れ線からの距離): 折れ線の上で $q に最も近い点
sub project_point {
    my ($P, $q, $cum) = @_;
    $cum //= arc_length($P);
    my $A  = $P->slice('0:-2,:'); my $B = $P->slice('1:-1,:');
    my $AB = $B - $A;                                             # (n-1,d)
    my $l2 = ($AB**2)->mv(1, 0)->sumover;                         # (n-1): over the coordinates (dim 1)
    my $t  = ((($q->dummy(0) - $A) * $AB)->mv(1, 0)->sumover / $l2->clip(1e-300, undef))->clip(0, 1);
    my $X  = $A + $t->dummy(1) * $AB;
    my $dd = sqrt((($X - $q->dummy(0))**2)->mv(1, 0)->sumover);
    my $j  = $dd->minimum_ind->sclr;
    return ($cum->at($j) + $t->at($j) * sqrt($l2->at($j)), $dd->at($j));
}

# midline_labels('10-20' | '10-10') -> ([labels], pdl(fractions)): 鼻根から後頭結節までの正中
sub midline_labels {
    my ($sys) = @_;
    $sys //= '10-20';
    my @l = $sys eq '10-20' ? ([Nasion => 0], [Fpz => 0.1], [Fz => 0.3], [Cz => 0.5], [Pz => 0.7], [Oz => 0.9], [Inion => 1])
          : $sys eq '10-10' ? ([Nasion => 0], [Fpz => 0.1], [AFz => 0.2], [Fz => 0.3], [FCz => 0.4], [Cz => 0.5],
                               [CPz => 0.6], [Pz => 0.7], [POz => 0.8], [Oz => 0.9], [Inion => 1])
          : croak "midline_labels: '10-20' or '10-10', not '$sys'";
    return ([ map { $_->[0] } @l ], pdl(map { $_->[1] } @l));
}

1;

__END__

=head1 NAME

PDL::EEG::ElectrodeLocation::Polyline - arc length along a polyline, points at given lengths or fractions, projection onto it

=head1 SYNOPSIS

  use PDL::EEG::ElectrodeLocation::Polyline qw(arc_length point_at_fraction project_point midline_labels);
  my $cum = arc_length($P);                       # $P (n,d), 2-D or 3-D
  my ($lab, $f) = midline_labels('10-20');        # Nasion Fpz Fz Cz Pz Oz Inion at 0 .1 .3 .5 .7 .9 1
  my $X = point_at_fraction($P, $f, $cum);        # (7,d), linear interpolation along the polyline
  my ($s, $off) = project_point($P, $q, $cum);    # length along it and distance from it

=head1 DESCRIPTION

Dimension-agnostic polyline helpers for electrode placement: the same functions serve a line
traced on an image (PDL::EEG::ElectrodeLocation::TraceLine) and a profile on a head surface.

=cut
