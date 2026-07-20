#!/usr/bin/env perl
#
# overlay_nyhead.pl — standard_1020.elc を New York Head の 19ch に重畳し、
#   電極対応・座標処理の正しさを残差(mm)で独立検証する。
#
#   .elc 側は PDL::EEG::IO::ASA::read_elc で読む(P:EEG15 で統合したパーサ)。
#   read_elc は LPA/RPA/Nz を fiducials として自動抽出するので、そのまま使う。
#
#   検証の考え方(P:EEG14 の筋を踏襲):
#     fiducial 3点 (LPA/RPA/Nz) だけで rigid + uniform-scale の座標枠を合わせ、
#     19ch 記録電極の残差を「独立」に測る。fiducial は合わせに使い、19ch は
#     検証に使う — 電極対応が正しければ 19ch が勝手に合う。
#
#   出力:
#     (a) 生残差 (アラインなし)     — 両者が既に同一 MNI 枠かの直接確認
#     (b) fiducial 枠アライン後の残差 — worst ch も名指し(対応ミス炙り出し)
#     electrodes_overlay.xyz         — GS3D 用点群 (NY実点 + .elc アライン点)
#
# 使い方:
#   perl overlay_nyhead.pl --elc standard_1020.elc --ny nyhead19.txt
#       nyhead19.txt: "name x y z" を1行1電極。19ch(Fp1..O2)+ LPA/RPA/Nz。
#       (read_nyhead_19.pl の出力をこの形で書き出せばよい)
#
#   perl overlay_nyhead.pl --elc standard_1020.elc --selftest
#       .elc から既知変換で疑似NYを合成し、(a)大 / (b)≈0 を確認(数式の自己検証)。
#
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PDL;
use PDL::EEG::IO::ASA qw(read_elc);

my %opt = (elc => 'standard_1020.elc', ny => undef, selftest => 0);
while (@ARGV) {
    my $a = shift @ARGV;
    if    ($a eq '--elc')      { $opt{elc}      = shift @ARGV }
    elsif ($a eq '--ny')       { $opt{ny}       = shift @ARGV }
    elsif ($a eq '--selftest') { $opt{selftest} = 1 }
    else { die "unknown option: $a\n" }
}

# 10-20 標準19ch(NY Head clab と同綴り: T7/T8/P7/P8)
my @CLAB = qw(Fp1 Fp2 F7 F3 Fz F4 F8 T7 C3 Cz C4 T8 P7 P3 Pz P4 P8 O1 O2);

# ---- helpers --------------------------------------------------------------
sub v3   { pdl(double, @_[0..2]) }
sub unit { my ($v) = @_; $v / sqrt(($v**2)->sum) }
sub vmag { my ($v) = @_; sqrt(($v**2)->sum) }
sub cross {
    my ($a, $b) = @_;
    pdl(double,
        $a->at(1)*$b->at(2) - $a->at(2)*$b->at(1),
        $a->at(2)*$b->at(0) - $a->at(0)*$b->at(2),
        $a->at(0)*$b->at(1) - $a->at(1)*$b->at(0));
}

# fiducial 3点 -> (origin, 正規直交基底 ex/ey/ez, 特性スケール)
#   o  = LPA-RPA 中点 / ex = 右向き / ey = 前方(exに直交) / ez = 上 = ex×ey
#   s  = 3 fiducial の重心からの RMS 距離(uniform scale の代表長)
sub fid_frame {
    my ($L, $R, $N) = @_;
    my $o  = 0.5 * ($L + $R);
    my $ex = unit($R - $L);
    my $vy = $N - $o;
    my $ey = unit($vy - (($vy*$ex)->sum) * $ex);
    my $ez = cross($ex, $ey);
    my $c  = ($L + $R + $N) / 3;
    my $s  = sqrt((vmag($L-$c)**2 + vmag($R-$c)**2 + vmag($N-$c)**2) / 3);
    return { o => $o, ex => $ex, ey => $ey, ez => $ez, s => $s };
}

# 点 p を srcフレームの座標にほどき、dstフレームへ同一 rigid+scale で置き直す
sub map_frame {
    my ($p, $src, $dst) = @_;
    my $d = $p - $src->{o};
    my $c0 = ($d * $src->{ex})->sum / $src->{s};   # 無次元座標
    my $c1 = ($d * $src->{ey})->sum / $src->{s};
    my $c2 = ($d * $src->{ez})->sum / $src->{s};
    return $dst->{o}
         + $dst->{s} * ($c0*$dst->{ex} + $c1*$dst->{ey} + $c2*$dst->{ez});
}

# ---- load .elc via read_elc ----------------------------------------------
my $mon = read_elc($opt{elc});
for my $nm (@CLAB) { die "elc missing $nm\n" unless $mon->{pos}{$nm} }
for my $f (qw(LPA RPA Nz)) { die "elc missing fiducial $f\n" unless $mon->{fiducials}{$f} }
my %elc = map { $_ => v3(@{$mon->{pos}{$_}}{qw(x y z)}) } @CLAB;
my %elc_fid = map { $_ => v3(@{$mon->{fiducials}{$_}}) } qw(LPA RPA Nz);
printf "elc  : %s (%d electrodes; 19ch + LPA/RPA/Nz found)\n", $opt{elc}, $mon->{n};

# ---- obtain NY Head 19ch + fiducials -------------------------------------
my (%ny, %ny_fid, $mode);
if ($opt{selftest}) {
    $mode = 'SELFTEST (synthetic NY = known transform of .elc)';
    # 既知の rigid+scale: z軸まわり7deg回転 + 平行移動 + scale 0.94、微小ノイズ無し
    my $th = 7 * 4*atan2(1,1) / 180;
    my ($ct, $st) = (cos($th), sin($th));
    my $t = v3(5, -3, 8);
    my $sc = 0.94;
    my $xform = sub {
        my ($p) = @_;
        my $rp = v3($ct*$p->at(0) - $st*$p->at(1),
                    $st*$p->at(0) + $ct*$p->at(1),
                    $p->at(2));
        return $sc*$rp + $t;
    };
    %ny     = map { $_ => $xform->($elc{$_}) } @CLAB;
    %ny_fid = map { $_ => $xform->($elc_fid{$_}) } qw(LPA RPA Nz);
}
elsif (defined $opt{ny}) {
    $mode = "NY from $opt{ny}";
    open my $fh, '<', $opt{ny} or die "open $opt{ny}: $!";
    my %raw;
    while (<$fh>) {
        next if /^\s*#/ || /^\s*$/;
        my ($nm, $x, $y, $z) = split;
        $raw{$nm} = v3($x, $y, $z);
    }
    close $fh;
    for my $nm (@CLAB)        { die "NY file missing $nm\n"        unless exists $raw{$nm} }
    for my $f (qw(LPA RPA Nz)){ die "NY file missing fiducial $f\n" unless exists $raw{$f} }
    %ny     = map { $_ => $raw{$_} } @CLAB;
    %ny_fid = map { $_ => $raw{$_} } qw(LPA RPA Nz);
}
else {
    die "need --ny nyhead19.txt (name x y z: 19ch + LPA/RPA/Nz)  or  --selftest\n";
}
print "mode : $mode\n\n";

# ---- fiducial coords side-by-side ----------------------------------------
printf "=== fiducial coords (MNI mm) ===\n  %-4s %26s %26s\n", '', 'NY Head', '.elc';
for my $f (qw(LPA RPA Nz)) {
    printf "  %-4s (%8.1f %8.1f %8.1f)   (%8.1f %8.1f %8.1f)\n",
        $f, $ny_fid{$f}->list, $elc_fid{$f}->list;
}
print "\n";

# ---- (a) raw residual (no alignment) -------------------------------------
print "=== (a) raw residual  ||NY - .elc||  (no alignment) ===\n";
my @raw;
for my $nm (@CLAB) {
    my $r = vmag($ny{$nm} - $elc{$nm})->sclr;
    push @raw, $r;
}
report(\@raw);

# ---- (b) residual after fiducial-frame alignment -------------------------
my $ny_frame  = fid_frame(@ny_fid{qw(LPA RPA Nz)});
my $elc_frame = fid_frame(@elc_fid{qw(LPA RPA Nz)});
printf "=== (b) residual after fiducial-frame alignment (scale NY/elc = %.3f) ===\n",
    $ny_frame->{s} / $elc_frame->{s};
my (@ali, @dump);
push @dump, "# name x y z set   (NY = NY Head 19ch, ELC = .elc RAW native MNI — both already coregistered)\n";
for my $nm (@CLAB) {
    my $ny_p = $ny{$nm};
    my $el_al = map_frame($elc{$nm}, $elc_frame, $ny_frame);  # aligned: for (b) residual only
    push @ali, vmag($ny_p - $el_al)->sclr;
    push @dump, sprintf("%-4s %8.2f %8.2f %8.2f NY\n",  $nm, $ny_p->list);
    push @dump, sprintf("%-4s %8.2f %8.2f %8.2f ELC\n", $nm, $elc{$nm}->list);  # RAW, not aligned
}
report(\@ali, 1);

# ---- GS3D-ready point cloud ----------------------------------------------
my $OUT = 'electrodes_overlay.xyz';
open my $o, '>', $OUT or die "write $OUT: $!";
print $o @dump;
close $o;
printf "wrote %s  (%d points: NY 19ch + .elc 19ch RAW; color by NY/ELC set)\n",
    $OUT, scalar(@CLAB) * 2;

# ---- residual reporter ----------------------------------------------------
sub report {
    my ($r, $name_worst) = @_;
    printf "  %-4s %8s\n", 'ch', 'mm';
    printf "  %-4s %8.1f\n", $CLAB[$_], $r->[$_] for 0 .. $#$r;
    my $p = pdl($r);
    printf "  -> mean=%.1f  median=%.1f  max=%.1f  (mm)\n",
        $p->avg, $p->qsort->at(int(@$r/2)), $p->max;
    if ($name_worst) {
        my $wi = $p->maximum_ind->sclr;
        printf "  worst: %s = %.1f mm\n", $CLAB[$wi], $p->at($wi);
    }
    print "\n";
}
