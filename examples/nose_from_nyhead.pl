#!/usr/bin/env perl
# ===========================================================================
# nose_from_nyhead.pl
#
# New York Head の skin サーフェスから鼻尖(pronasale)頂点を取り、共有する
# 10-20 電極19点で NYHead フレーム(MNI, AC 原点)を対象モンタージュ(.elc)の
# フレームへ相似 Procrustes で合わせ込み、モンタージュ座標系での "nose"
# 位置を .elc 追記行として出力する。
#
#   perl nose_from_nyhead.pl sa_nyhead.mat your_montage.elc
#
# 出力の最後の1行(x y z)を .elc の Positions 末尾へ、"nose" を Labels 末尾へ
# 追記すれば、その montage の枠内で nose が使える。co-reg RMS も表示するので
# 妥当性(数 mm 以内か)を確認すること。座標系が違う NYHead mm を生で使うと
# ずれるため、この合わせ込みは必須。
# ===========================================================================
use strict; use warnings;
use PDL;
use PDL::MatrixOps qw(svd det identity);
use PDL::IO::NYHead;
use PDL::EEG::IO::ASA ();

my ($mat,$elc)=@ARGV;
die "usage: $0 sa_nyhead.mat montage.elc\n" unless $mat && $elc;

# name normalisation so T3/T7 etc. and case don't split the shared set
my %ALI=(T3=>'T7',T4=>'T8',T5=>'P7',T6=>'P8',T7=>'T3',T8=>'T4',P7=>'T5',P8=>'T6');
sub canon { my ($h,$n)=@_; return $n if $h->{$n};
            my $u=uc $n; for ($u,$ALI{$u}//()) { return $_ if defined && $h->{$_} } undef }

# ---- transform primitives (validated: xf(R,P)=R x P for P=(N,3)) ----------
# normalise any (N,3)/(3,N) point array to (N,3); leave (3,3) alone (assume N,3)
sub to_N3 { my $p=shift; ($p->dim(1)==3) ? $p
                        : ($p->dim(0)==3) ? $p->xchg(0,1)->sever : $p }
sub xf   { my($R,$P)=@_; $R x $P }
sub mean3{ my $P=shift; pdl(map { $P->slice(":,($_)")->avg } 0..2) }
sub fit_similarity {                       # (src,dst each (N,3)) -> (s,R,t,rms)
    my ($src,$dst)=@_; my $N=$src->dim(0);
    my $mu_s=mean3($src); my $mu_d=mean3($dst);
    my $Xs=$src-$mu_s->dummy(0); my $Xd=$dst-$mu_d->dummy(0);
    my $Cov=($Xd->dummy(2)*$Xs->dummy(1))->sumover/$N;
    my ($U,$D,$V)=svd($Cov);
    my $S=identity(3); my $g=($U x $V->xchg(0,1))->det; $S->set(2,2,$g<0?-1:1);
    my $R=$V x $S x $U->xchg(0,1);          # V S U^T
    my $s=($D*$S->diagonal(0,1))->sum/(($Xs**2)->sum/$N);
    my $t=$mu_d - $s*xf($R,$mu_s->dummy(0))->flat;
    my $rms=sqrt(((($s*xf($R,$src)+$t->dummy(0))-$dst)**2)->sumover->avg);
    return ($s,$R,$t,$rms);
}
sub nose_tip {                             # head vc (Nv,3) MNI mm -> (3)
    my ($vc,$xtol)=@_; $xtol//=15;
    my $mid=which($vc->slice(":,(0)")->abs < $xtol);
    $mid=sequence($vc->dim(0)) if $mid->nelem<1;
    my $ys=$vc->slice(":,(1)")->index($mid);
    return $vc->slice("(".$mid->at($ys->maximum_ind)."),")->flat;
}

# ---- load ------------------------------------------------------------------
my $ny  = PDL::IO::NYHead->new($mat);
my $mon = PDL::EEG::IO::ASA::read_elc($elc);
my $nl  = $ny->electrode_labels_19;        # 19 names
my $np  = to_N3($ny->electrode_pos_19);    # (19,3) MNI mm  (orientation-normalised)

# shared subset (NYHead <-> montage) by canonical name
my (@si,@sn);
for my $i (0..$#$nl) { my $k=canon($mon->{pos},$nl->[$i]); next unless $k; push @si,$i; push @sn,$k }
die "only ".scalar(@si)." shared electrodes - need >=8 for a stable fit\n" if @si<8;
my $src = $np->dice_axis(0, pdl(long,@si));                 # NY (M,3)
my $dst = to_N3(pdl(map { my $p=$mon->{pos}{$_}; [$p->{x},$p->{y},$p->{z}] } @sn)); # (M,3)

my ($s,$R,$t,$rms)=fit_similarity($src,$dst);
my $nose_ny  = nose_tip(to_N3($ny->surface('head')->{vc}));
my $nose_mon = ($s*xf($R,$nose_ny->dummy(0))+$t->dummy(0))->flat;

printf STDERR "shared electrodes : %d (%s)\n", scalar(@sn), join(",",@sn);
printf STDERR "co-reg scale=%.4f  RMS=%.2f mm  %s\n",
    $s,$rms, $rms<5 ? "(good)" : "(!! large - check montage frame/units)";
printf STDERR "nose tip (NY MNI) : %.1f %.1f %.1f\n", list $nose_ny;
printf STDERR "--- append to %s : Positions line below, and 'nose' to Labels ---\n",$elc;
printf "%9.4f %9.4f %9.4f\n", list $nose_mon;
