#!/usr/bin/env perl
# ---------------------------------------------------------------------------
# powers_to_usda.pl  (P:EEG22)
#
#   小さいラッパー: <indir> の powers.<lat>.dat(sep_n20_sweep.pl が書いたもの,
#   74382 行=cortex75K 順)を読み、cortex(10K など)の頂点色にして、
#   潜時ごとに **静的 USDA を 1 枚** <outdir>/<res>.<lat>.usda に書き出す。
#   源計算は持たない(dat を読むだけ)。単発で潜時を残す/配る用途。
#   (スクラブ再生したいなら 1 本アニメの sweep_usda_anim.pl の方)
#
#   例:
#     # まず dat を作る(範囲を絞ると軽い)
#     perl ... sep_n20_sweep.pl ... --lat-min 18.4 --lat-max 40 --lat-step 0.1 --outdir LORETA_output
#     # それを USDA 連番に
#     perl -I<PDL_IO_NYHead/lib> powers_to_usda.pl --mat sa_nyhead.mat \
#          --indir LORETA_output --surf cortex10K --lat-min 18.4 --lat-max 40 \
#          --cmap inferno --axes 90 --outdir usda_seq
# ---------------------------------------------------------------------------
use strict; use warnings;
use PDL; use PDL::IO::Misc; use PDL::IO::NYHead; use File::Path qw(make_path); use Getopt::Long; ## using rcol

my ($mat, $indir);
my $surf = "cortex10K";
my $cmap = "inferno";
my $overlay = 1; my $threshold = 0.25; my $base_grey = 0.55;
my $norm = "global";               # global(全ファイル最大) | frame(各ファイル最大) | fixed(--vmax)
my $vmax;                          # --norm fixed のときの絶対最大
my $axes = 90;
my ($lat_min, $lat_max);
my $outdir = "usda_seq";
GetOptions(
    "mat=s"=>\$mat, "indir=s"=>\$indir, "surf=s"=>\$surf, "cmap=s"=>\$cmap,
    "overlay!"=>\$overlay, "threshold=f"=>\$threshold, "base-grey=f"=>\$base_grey,
    "norm=s"=>\$norm, "vmax=f"=>\$vmax, "axes=f"=>\$axes,
    "lat-min=f"=>\$lat_min, "lat-max=f"=>\$lat_max, "outdir=s"=>\$outdir,
) or die "bad options\n";
$mat && $indir or die "need --mat and --indir\n";
(my $res = $surf) =~ s{^/?sa/}{}; $res =~ s{/$}{};

sub _cmap_rgb {
    my ($t, $name) = @_; $t = 0 if $t < 0; $t = 1 if $t > 1;
    my %STOPS = (
        inferno => [[0.0,0.001,0.000,0.014],[0.25,0.258,0.039,0.406],
                    [0.5,0.578,0.148,0.404],[0.75,0.865,0.317,0.226],
                    [0.9,0.988,0.645,0.140],[1.0,0.988,0.998,0.645]],
        viridis => [[0.0,0.267,0.005,0.329],[0.25,0.229,0.322,0.545],
                    [0.5,0.128,0.567,0.551],[0.75,0.369,0.789,0.383],[1.0,0.993,0.906,0.144]],
        hot     => [[0.0,0,0,0],[0.375,1,0,0],[0.75,1,1,0],[1.0,1,1,1]],
    );
    my $s = $STOPS{$name} || $STOPS{inferno};
    for my $i (0 .. $#$s-1) {
        my ($t0,@c0) = @{$s->[$i]}; my ($t1,@c1) = @{$s->[$i+1]};
        if ($t <= $t1 || $i == $#$s-1) {
            my $f = ($t1 > $t0) ? ($t-$t0)/($t1-$t0) : 0; $f=0 if $f<0; $f=1 if $f>1;
            return [ $c0[0]+($c1[0]-$c0[0])*$f, $c0[1]+($c1[1]-$c0[1])*$f, $c0[2]+($c1[2]-$c0[2])*$f ];
        }
    }
    return [1,1,1];
}
sub _overlay_rgb {
    my ($t, $cm, $grey, $thr) = @_;
    my $rgb = _cmap_rgb($t, $cm);
    return $rgb unless defined $grey;
    my $a = ($thr < 1) ? ($t-$thr)/(1-$thr) : 0; $a=0 if $a<0; $a=1 if $a>1;
    $a = $a*$a*(3-2*$a);
    return [ map { $grey*(1-$a) + $rgb->[$_]*$a } 0..2 ];
}

# ---- cortex メッシュ + in_from(頂点→75K 0-based) ----
my $ny  = PDL::IO::NYHead->new($mat);
my $cx  = $ny->cortex($res);
my $vc  = $cx->{vc};  my $tri = $cx->{tri};  my $inf = $cx->{in_from};
my $Nv  = $vc->dim(0);  my $Nf = $tri->dim(0);
printf "surface %s: %d verts, %d faces\n", $res, $Nv, $Nf;

# ---- 対象ファイルを列挙(powers.<lat>.dat)、潜時で絞る/並べる ----
my @files;
for my $f (glob "$indir/powers.*.dat") {
    next unless $f =~ /powers\.([-0-9.]+)\.dat$/;
    my $lat = $1 + 0;
    next if defined $lat_min && $lat < $lat_min - 1e-9;
    next if defined $lat_max && $lat > $lat_max + 1e-9;
    push @files, [$lat, $f];
}
@files = sort { $a->[0] <=> $b->[0] } @files;
@files or die "no powers.<lat>.dat in $indir (range) — run sep_n20_sweep.pl first\n";
printf "found %d powers files\n", scalar @files;

# ---- pass1: 読み込み+10K へ写像(+global 最大) ----
my (@pw, @lats);
my $gmax = 0;
for my $e (@files) {
    my ($p) = rcols($e->[1]);                 # (74382)
    my $v10 = $p->index($inf)->sever;         # (Nv)
    push @pw, $v10; push @lats, $e->[0];
    my $m = $v10->max; $gmax = $m if $m > $gmax;
}
$gmax = 1e-300 if $gmax <= 0;

# 色 LUT
my $NLUT = 1024; my $grey = $overlay ? $base_grey : undef;
my @LUT = map { _overlay_rgb($_/($NLUT-1), $cmap, $grey, $threshold) } 0..$NLUT-1;

# 静的メッシュ部品(全ファイル共通)を一度だけ作る
my $counts = join(",", (3) x $Nf);
my $idxstr = join(",", (long($tri)-1)->xchg(0,1)->flat->list);
my @vx = $vc->slice(":,(0)")->list; my @vy = $vc->slice(":,(1)")->list; my @vz = $vc->slice(":,(2)")->list;
my $ptstr  = join(",", map { sprintf("(%.3f,%.3f,%.3f)", $vx[$_],$vy[$_],$vz[$_]) } 0..$Nv-1);
my $axesblk = "";
if ($axes > 0) {
    my ($L,$W) = ($axes, $axes/50);
    $axesblk = <<"USDA";
    def BasisCurves "Axes"
    {
        uniform token type = "linear"
        int[] curveVertexCounts = [2, 2, 2]
        point3f[] points = [(0,0,0),($L,0,0), (0,0,0),(0,$L,0), (0,0,0),(0,0,$L)]
        color3f[] primvars:displayColor = [(1,0,0),(0,1,0),(0,0,1)] (interpolation = "uniform")
        float[] widths = [$W] (interpolation = "constant")
    }
USDA
}

# ---- pass2: 潜時ごとに静的 USDA を書く ----
make_path($outdir) unless -d $outdir;
my $n = 0;
for my $k (0 .. $#pw) {
    my $lat = $lats[$k];
    my $den = $norm eq 'frame' ? ($pw[$k]->max || 1e-300)
            : $norm eq 'fixed' ? ($vmax // $gmax)
            :                     $gmax;
    my @li  = ( ($pw[$k]/$den)->clip(0,1) * ($NLUT-1) + 0.5 )->long->list;
    my $col = join(",", map { my $c=$LUT[$_]; sprintf("(%.3f,%.3f,%.3f)",@$c) } @li);

    my $out = sprintf("%s/%s.%05.1f.usda", $outdir, $res, $lat);
    open my $fh, '>', $out or die "open $out: $!";
    print $fh "#usda 1.0\n(\n    upAxis = \"Z\"\n    metersPerUnit = 0.001\n)\n\n";
    print $fh "def Xform \"NYHead\"\n{\n    def Mesh \"$res\"\n    {\n";
    print $fh "        int[] faceVertexCounts = [$counts]\n";
    print $fh "        int[] faceVertexIndices = [$idxstr]\n";
    print $fh "        point3f[] points = [$ptstr]\n";
    print $fh "        color3f[] primvars:displayColor = [$col] (interpolation = \"vertex\")\n";
    print $fh "    }\n$axesblk}\n";
    close $fh;
    $n++;
}
printf "wrote %d USDA files to %s/  (%s.<lat>.usda, norm=%s%s)\n",
    $n, $outdir, $res, $norm, ($norm eq 'fixed' ? sprintf(" vmax=%.3g", $vmax//$gmax) : "");
