#!/usr/bin/env perl
#
# dump_nyhead19.pl — sa_nyhead.mat から標準10-20の19ch＋fiducialの座標を抜き、
#                    overlay_nyhead.pl 用の nyhead19.txt を書き出す。
#
#   電極座標は /sa/locs_3D_orig(MNI mm)。137MB の V_fem_normal は読まない。
#
#   前提(P:EEG14 実測で確定):
#     - v7.3 / HDF5、PDL::IO::HDF5 で読む
#     - PDL は locs_3D_orig を (231,3)=(電極,xyz) で落とす。ただし dataset ごとに
#       h5ls 表示と PDL の落ち方が食い違うので、実測 dims から xyz 軸を判定する
#     - @IDX19 / %NYFID は clab_electrodes に対する 0-based 添字(凍結済み)
#
#   組み込みサニティ: 出力 fiducial が既知の NY 値と一致するか照合して表示。
#     LPA(-83.5,-17.5,-38.5) RPA(83.5,-17.5,-38.5) Nz(0,83.5,-41)
#     一致すれば「フィールド名・転置・添字」がすべて正しいことの証明。
#
# 使い方: perl dump_nyhead19.pl [sa_nyhead.mat] [nyhead19.txt]
#
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;

my $MAT = shift // 'sa_nyhead.mat';
my $OUT = shift // 'nyhead19.txt';

# --- 凍結: 10-20 標準19ch と NY Head clab の 0-based 添字 -------------------
my @CLAB  = qw(Fp1 Fp2 F7 F3 Fz F4 F8 T7 C3 Cz C4 T8 P7 P3 Pz P4 P8 O1 O2);
my @IDX19 = (156,157,119,115,159,116,120,222, 19,48,20,223,
             182,178,220,179,183,168,169);
my %NYFID = (LPA => 162, RPA => 221, Nz => 167);

# 照合用の既知 NY fiducial (P:EEG14)
my %EXPECT = (LPA => [-83.5,-17.5,-38.5], RPA => [83.5,-17.5,-38.5], Nz => [0,83.5,-41]);

# --- read locs_3D_orig -----------------------------------------------------
my $h = PDL::IO::HDF5->new($MAT) or die "open $MAT: $!";
my $locs = $h->dataset('/sa/locs_3D_orig')->get;
printf "locs_3D_orig dims = %s   type = %s\n", join('x', $locs->dims), $locs->type;

# xyz 軸(サイズ3)を dim1、電極軸を dim0 に揃える(推測せず実測 dims で分岐)
if    ($locs->dim(1) == 3) { }                       # 既に (Nelec, 3)
elsif ($locs->dim(0) == 3) { $locs = $locs->xchg(0,1)->sever }  # (3,Nelec) -> (Nelec,3)
else  { die "locs_3D_orig has no size-3 axis: dims=", join('x',$locs->dims), "\n" }
my $nelec = $locs->dim(0);
printf "  -> electrode axis = dim0 (size %d), xyz axis = dim1\n", $nelec;
my $need = 0; for (@IDX19, values %NYFID) { $need = $_ if $_ > $need }
die "index $need out of range (only $nelec electrodes)\n" if $need >= $nelec;

sub row { my ($i) = @_; $locs->slice("($i),")->flat }   # 電極 i の (x,y,z)

# --- fiducial サニティ照合 -------------------------------------------------
print "\n=== fiducial sanity vs P:EEG14 known values (mm) ===\n";
printf "  %-4s %26s %26s %8s\n", '', 'from .mat', 'expected', 'd(mm)';
my $ok = 1;
for my $f (qw(LPA RPA Nz)) {
    my $got = row($NYFID{$f});
    my $exp = pdl(double, @{$EXPECT{$f}});
    my $d   = sqrt((($got - $exp)**2)->sum)->sclr;
    $ok &&= ($d < 1.0);
    printf "  %-4s (%8.1f %8.1f %8.1f)   (%8.1f %8.1f %8.1f) %8.2f\n",
        $f, $got->list, $exp->list, $d;
}
print $ok ? "  OK: fiducials match -> field/transpose/indices all correct\n"
          : "  !! MISMATCH: check field name, transpose, or IDX/NYFID indices\n";

# --- write nyhead19.txt ----------------------------------------------------
open my $o, '>', $OUT or die "write $OUT: $!";
print $o "# name x y z   (NY Head, MNI mm, from /sa/locs_3D_orig)\n";
for my $i (0 .. $#CLAB) {
    printf $o "%-4s %9.3f %9.3f %9.3f\n", $CLAB[$i], row($IDX19[$i])->list;
}
for my $f (qw(LPA RPA Nz)) {
    printf $o "%-4s %9.3f %9.3f %9.3f\n", $f, row($NYFID{$f})->list;
}
close $o;
printf "\nwrote %s  (19ch + LPA/RPA/Nz)\n", $OUT;
print  "next: perl overlay_nyhead.pl --elc standard_1020.elc --ny $OUT\n";
