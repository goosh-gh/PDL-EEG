#!/usr/bin/env perl
# ===========================================================================
# silhouette_from_nyhead.pl
#
# NYHead skin サーフェスから「正中矢状面の横顔シルエット」を一度だけ抽出し、
# 共有 10-20 電極19点で NYHead(MNI)→対象モンタージュ(.elc)フレームへ相似
# Procrustes 合わせ込みして、montage 座標系の (Y,Z) 閉ポリラインを書き出す。
# MAP2D の sagittal view はこのポリラインを読むだけ(描画時 NYHead 非依存)。
#
#   perl silhouette_from_nyhead.pl sa_nyhead.mat montage.elc [out.poly]
#
# 出力 out.poly: "# sagittal silhouette (montage frame, mm)" + "Y Z" 行。
# 左 view は screen_x=-Y, 右 view は screen_x=+Y で MAP2D 側が向きを決める。
# ===========================================================================
use strict; use warnings;
use PDL;
use PDL::MatrixOps qw(svd det identity);
use lib '/Users/goosh/src/PDL_IO_NYHead/lib'; use PDL::IO::NYHead;
use lib '/Users/goosh/src/PDL_EEG//lib'; use PDL::EEG::IO::ASA ();

my $mat = "/Users/goosh/src/NYHead/sa_nyhead.mat";
my $elc = "/Users/goosh/src/NYHead/standard_1020_eog_nose.elc";
# my ($mat,$elc,$out)=@ARGV;
die "usage: $0 sa_nyhead.mat montage.elc [out.poly]\n" unless $mat && $elc;
my $out //= 'sagittal_silhouette.poly';

my %ALI=(T3=>'T7',T4=>'T8',T5=>'P7',T6=>'P8',T7=>'T3',T8=>'T4',P7=>'T5',P8=>'T6');
sub canon { my ($h,$n)=@_; return $n if $h->{$n};
            my $u=uc $n; for ($u,$ALI{$u}//()) { return $_ if defined && $h->{$_} } undef }
sub to_N3 { my $p=shift; ($p->dim(1)==3)?$p:($p->dim(0)==3)?$p->xchg(0,1)->sever:$p }
sub xf   { my($R,$P)=@_; $R x $P }
sub mean3{ my $P=shift; pdl(map { $P->slice(":,($_)")->avg } 0..2) }
sub fit_similarity {
    my ($src,$dst)=@_; $src=to_N3($src); $dst=to_N3($dst); my $N=$src->dim(0);
    my $mu_s=mean3($src); my $mu_d=mean3($dst);
    my $Xs=$src-$mu_s->dummy(0); my $Xd=$dst-$mu_d->dummy(0);
    my $Cov=($Xd->dummy(2)*$Xs->dummy(1))->sumover/$N;
    my ($U,$D,$V)=svd($Cov);
    my $S=identity(3); my $g=($U x $V->xchg(0,1))->det; $S->set(2,2,$g<0?-1:1);
    my $R=$V x $S x $U->xchg(0,1);
    my $s=($D*$S->diagonal(0,1))->sum/(($Xs**2)->sum/$N);
    my $t=$mu_d - $s*xf($R,$mu_s->dummy(0))->flat;
    my $rms=sqrt(((($s*xf($R,$src)+$t->dummy(0))-$dst)**2)->sumover->avg);
    return ($s,$R,$t,$rms);
}
# 2D convex hull (Andrew monotone chain) -> CCW outline of (Y,Z)
sub convex_hull {
    my ($y,$z)=@_;
    my @p = map { [$y->at($_),$z->at($_)] } 0..$y->nelem-1;
    @p = sort { $a->[0]<=>$b->[0] or $a->[1]<=>$b->[1] } @p;
    my $x = sub { my($o,$a,$b)=@_; ($a->[0]-$o->[0])*($b->[1]-$o->[1])-($a->[1]-$o->[1])*($b->[0]-$o->[0]) };
    my (@lo,@up);
    for my $q (@p){ pop @lo while @lo>=2 && $x->($lo[-2],$lo[-1],$q)<=0; push @lo,$q }
    for my $q (reverse @p){ pop @up while @up>=2 && $x->($up[-2],$up[-1],$q)<=0; push @up,$q }
    pop @lo; pop @up; return (@lo,@up);
}

my $ny  = PDL::IO::NYHead->new($mat);
my $mon = PDL::EEG::IO::ASA::read_elc($elc);
my $nl  = $ny->electrode_labels_19;
my $np  = to_N3($ny->electrode_pos_19);
my (@si,@sn);
for my $i (0..$#$nl){ my $k=canon($mon->{pos},$nl->[$i]); next unless $k; push @si,$i; push @sn,$k }
die "only ".scalar(@si)." shared electrodes - need >=8\n" if @si<8;
my $src=$np->dice_axis(0,pdl(long,@si));
my $dst=to_N3(pdl(map { my $p=$mon->{pos}{$_}; [$p->{x},$p->{y},$p->{z}] } @sn));
my ($s,$R,$t,$rms)=fit_similarity($src,$dst);

# skin vertices -> montage frame, project to sagittal (Y,Z), hull.
# The outline (nose bump included) comes purely from the NYHead skin surface -
# it does NOT depend on which electrodes (nose/LM/RM ...) the montage happens to
# carry.  A small outward buffer keeps on-scalp electrodes comfortably inside;
# MAP2D grows the outline further at draw time if a montage still pokes out.
my $vc  = to_N3($ny->surface('head')->{vc});
my $vm  = $s*xf($R,$vc) + $t->dummy(0);            # (Nv,3) montage mm
my @hull = convex_hull($vm->slice(":,(1)"), $vm->slice(":,(2)"));
{ # 3% outward buffer about the hull centroid
    my ($cy,$cz) = (0,0); $cy+=$_->[0] for @hull; $cz+=$_->[1] for @hull;
    $cy/=@hull; $cz/=@hull;
    $_ = [ $cy+($_->[0]-$cy)*1.03, $cz+($_->[1]-$cz)*1.03 ] for @hull;
}

open my $fh,'>',$out or die "$out: $!";
print $fh "# sagittal silhouette (montage frame, mm)  co-reg RMS=".sprintf('%.2f',$rms)." mm\n";
print $fh "# Y Z\n";
printf $fh "%.3f %.3f\n", @$_ for @hull;
close $fh;
printf STDERR "shared electrodes : %d\n", scalar @sn;
printf STDERR "co-reg RMS        : %.2f mm %s\n", $rms, $rms<5?"(good)":"(!! check units/frame)";
printf STDERR "silhouette        : %d vertices -> %s\n", scalar @hull, $out;
