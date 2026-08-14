#!/usr/bin/env perl
# ---------------------------------------------------------------------------
# sweep_usda_anim.pl  (P:EEG22)
#
#   潜時窓 [lat-min, lat-max] を lat-step おきに走査し、各潜時の源パワーを
#   cortex(10K など)の頂点色にして、**1 本のアニメ USD** に
#   color3f[] primvars:displayColor.timeSamples として書き出す。
#   メッシュ・軸は静的なので、usdview/Keynote で止めてグリグリ回せる。
#
#   源計算まで内蔵(中間 powers.dat 不要)。逆作用素は b 非依存なので 1 回だけ
#   構築し、窓を 1 回の行列積で一括 solve。頂点色は cmap＋灰オーバレイ(ビューアと
#   同一)。既定は窓全体で正規化(global)=N20 が明るく立ち上がる時間経過が見える。
#
#   例:
#     perl -I<PDL_EEG/lib> -I<PDL_IO_NYHead/lib> sweep_usda_anim.pl \
#       --mat sa_nyhead.mat --data sep_ave.txt --labels sep_labels.txt \
#       --sfreq 10000 --tmin -0.05 --method eloreta \
#       --lat-min 15 --lat-max 30 --lat-step 0.1 \
#       --surf cortex10K --cmap inferno --axes 90 --out n20_anim.usda
#
#   時間コード = 潜時 ×10 (例 201 = 20.1 ms)。startTimeCode..endTimeCode で
#   スクラブ。fps=timeCodesPerSecond。
# ---------------------------------------------------------------------------
use strict;
use warnings;
use PDL;
use PDL::IO::NYHead;
use Getopt::Long;

my @INV_FN = qw(avg_reference inverse_operator source_power);
my $HAVE_INV = eval { require PDL::EEG::Inverse::MinimumNorm;
                      PDL::EEG::Inverse::MinimumNorm->import(@INV_FN); 1 };

my ($data_file, $labels_file, $times_file, $mat);
my $sfreq = 0; my $tmin = 0; my $avg_ms = 0;
my $method = "eloreta"; my $reg = 0.05;
my $lat_min = 15.0; my $lat_max = 30.0; my $lat_step = 0.1;
my $surf = "cortex10K";
my $cmap = "inferno";
my $overlay = 1; my $threshold = 0.25; my $base_grey = 0.55;
my $norm = "global";               # global | frame
my $axes = 90;                     # 0=軸なし
my $cbar = 1;                      # カラーバーを独立プリムとして出力(--no-colorbar で無し)
my $fps  = 100;                    # timeCodesPerSecond
my $out  = "n20_anim.usda";
GetOptions(
    "data=s"=>\$data_file, "labels=s"=>\$labels_file, "times=s"=>\$times_file,
    "sfreq=f"=>\$sfreq, "tmin=f"=>\$tmin, "mat=s"=>\$mat, "avg-ms=f"=>\$avg_ms,
    "method=s"=>\$method, "reg-frac=f"=>\$reg,
    "lat-min=f"=>\$lat_min, "lat-max=f"=>\$lat_max, "lat-step=f"=>\$lat_step,
    "surf=s"=>\$surf, "cmap=s"=>\$cmap,
    "overlay!"=>\$overlay, "threshold=f"=>\$threshold, "base-grey=f"=>\$base_grey,
    "norm=s"=>\$norm, "axes=f"=>\$axes, "colorbar!"=>\$cbar, "fps=f"=>\$fps, "out=s"=>\$out,
) or die "bad options\n";
$data_file && $labels_file && $mat or die "need --data --labels --mat\n";
$HAVE_INV or die "PDL::EEG::Inverse::MinimumNorm not found (\@INC に PDL-EEG)\n";
$lat_step > 0 or die "--lat-step must be > 0\n";
(my $res = $surf) =~ s{^/?sa/}{};  $res =~ s{/$}{};   # '/sa/cortex10K' -> 'cortex10K'

# --- cmap + 灰オーバレイ(demo_gs3d_nyhead.pl と同一) ---
sub _cmap_rgb {
    my ($t, $name) = @_;  $t = 0 if $t < 0; $t = 1 if $t > 1;
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

# ---------------------------------------------------------- read waveform ---
my @rows;
open my $df, '<', $data_file or die "open $data_file: $!";
while (<$df>) { next unless /\S/; push @rows, [ split ' ' ]; }
close $df;
my $data = pdl(\@rows);                 # (n_time, n_ch)
my ($ntime, $nch) = ($data->dim(0), $data->dim(1));
open my $lf, '<', $labels_file or die "open $labels_file: $!";
my @labels = map { s/\s+$//r } grep { /\S/ } <$lf>; close $lf;
@labels == $nch or die "labels != data rows\n";
my $times;
if ($times_file) {
    open my $tf, '<', $times_file or die "open $times_file: $!";
    $times = pdl( map { chomp; $_ } grep { /\S/ } <$tf> ); close $tf;
} else { $sfreq > 0 or die "need --sfreq (or --times)\n"; $times = $tmin + sequence($ntime)/$sfreq; }

# -------------------------------------------- map electrodes -> NYHead ---
my $ny  = PDL::IO::NYHead->new($mat);
my $ecl = $ny->electrode_labels;
my %ny_idx; $ny_idx{ lc $ecl->[$_] } = $_ for 0 .. $#$ecl;
my %alias = (t3=>'t7', t4=>'t8', t5=>'p7', t6=>'p8', a1=>'', a2=>'');
my (@file_i, @ny_i);
for my $c (0 .. $nch-1) {
    my $k = lc $labels[$c]; $k =~ s/\s+//g; $k = $alias{$k} if exists $alias{$k};
    if ($k ne '' && exists $ny_idx{$k}) { push @file_i, $c; push @ny_i, $ny_idx{$k}; }
}
my $ne = @ny_i;  $ne >= 4 or die "only $ne electrodes matched NYHead\n";
printf "matched %d/%d electrodes\n", $ne, $nch;

my $sel  = pdl(long, \@ny_i);
my $fidx = pdl(long, \@file_i);
my $K    = $ny->leadfield->dice_axis(0, $sel)->transpose->sever;   # (Ns, Ne)
my $Ns   = $K->dim(0);
my $Kc   = avg_reference($K);
my $op   = inverse_operator($Kc, method=>$method, ref=>'car', reg_frac=>$reg);

# ---- cortex(10K など): vc / tri / in_from(10K頂点→75K 0-based) ----
my $cx = $ny->cortex($res);
my $vc = $cx->{vc};                       # (Nv,3) MNI mm
my $tri= $cx->{tri};                      # (Nf,3) 1-based(低解像内)
my $inf= $cx->{in_from};                  # (Nv) 0-based 75K index
my $Nv = $vc->dim(0);  my $Nf = $tri->dim(0);
printf "surface %s: %d verts, %d faces\n", $res, $Nv, $Nf;

# ------------------------------------- 窓を一括 solve → 10K 頂点パワー ---
my $wsamp = ($avg_ms>0 && $sfreq>0) ? int($avg_ms/1000*$sfreq+0.5) : 0;
my $lo10 = int($lat_min*10+0.5); my $hi10 = int($lat_max*10+0.5);
my $st10 = int($lat_step*10+0.5); $st10 = 1 if $st10 < 1;
my @codes; for (my $l=$lo10; $l<=$hi10; $l+=$st10) { push @codes, $l; }  # 時間コード=lat*10
my $Nlat = @codes;

my $B = zeroes($ne, $Nlat);
for my $k (0..$Nlat-1) {
    my $lat = $codes[$k]/10;
    my $ti  = int(($lat/1000-$tmin)*$sfreq+0.5); $ti=0 if $ti<0; $ti=$ntime-1 if $ti>$ntime-1;
    my $topo;
    if ($wsamp>0) { my $a=$ti-$wsamp;$a=0 if$a<0;my $bb=$ti+$wsamp;$bb=$ntime-1 if$bb>$ntime-1;
                    $topo=$data->slice("$a:$bb,:")->mv(0,0)->average; }
    else { $topo=$data->slice("($ti),:"); }
    $B->slice(":,($k)") .= $topo->index($fidx);
}
my $Bc = $B - $B->average->dummy(0);

my $POW;
my $batched = eval { my $p = source_power($op, $Bc);
    ($p->ndims==2 && $p->dim(0)==$Ns && $p->dim(1)==$Nlat) ? do{$POW=$p;1} : 0 };
if (!$batched) {
    $POW = zeroes($Ns, $Nlat);
    $POW->slice(":,($_)") .= source_power($op, $Bc->slice(":,($_)")->sever)->flat for 0..$Nlat-1;
}
my $P10 = $POW->dice_axis(0, $inf)->sever;         # (Nv, Nlat) 表示メッシュ上のパワー
printf "solved %d latencies %s; coloring %d verts...\n",
    $Nlat, ($batched?"batched":"per-latency"), $Nv;

# 正規化基準
my $gmax = ($norm eq 'frame') ? undef : $P10->max;
$gmax = 1e-300 if defined $gmax && $gmax <= 0;

# 色 LUT(cmap＋オーバレイを 0..1 で 1024 段プリ計算)
my $NLUT = 1024; my $grey = $overlay ? $base_grey : undef;
my @LUT = map { _overlay_rgb($_/($NLUT-1), $cmap, $grey, $threshold) } 0..$NLUT-1;

# ---------------------------------------------------------- write USDA ---
open my $fh, '>', $out or die "open $out: $!";
printf $fh "#usda 1.0\n";
printf $fh "# SEP source animation: threshold=%.3g cmap=%s overlay=%s base_grey=%.3g norm=%s lat=%.1f..%.1f ms\n",
    $threshold, $cmap, ($overlay?"on":"off"), $base_grey, $norm, $codes[0]/10, $codes[-1]/10;
printf $fh "(\n    upAxis = \"Z\"\n    metersPerUnit = 0.001\n";
printf $fh "    startTimeCode = %d\n    endTimeCode = %d\n    timeCodesPerSecond = %g\n", $codes[0],$codes[-1],$fps;
printf $fh "    customLayerData = {\n";
printf $fh "        double sep_threshold = %g\n", $threshold;
printf $fh "        string sep_cmap = \"%s\"\n", $cmap;
printf $fh "        string sep_overlay = \"%s\"\n", ($overlay?"on":"off");
printf $fh "        double sep_base_grey = %g\n", $base_grey;
printf $fh "        string sep_norm = \"%s\"\n", $norm;
printf $fh "        string sep_latency_ms = \"%.1f..%.1f\"\n", $codes[0]/10, $codes[-1]/10;
printf $fh "    }\n)\n\n";
print  $fh "def Xform \"NYHead\"\n{\n    def Mesh \"$res\"\n    {\n";
print  $fh "        int[] faceVertexCounts = [", join(",", (3) x $Nf), "]\n";
my @idx = (long($tri)-1)->xchg(0,1)->flat->list;                 # 0-based
print  $fh "        int[] faceVertexIndices = [", join(",", @idx), "]\n";
my @vx = $vc->slice(":,(0)")->list; my @vy = $vc->slice(":,(1)")->list; my @vz = $vc->slice(":,(2)")->list;
print  $fh "        point3f[] points = [",
           join(",", map { sprintf("(%.3f,%.3f,%.3f)", $vx[$_],$vy[$_],$vz[$_]) } 0..$Nv-1), "]\n";
print  $fh "        color3f[] primvars:displayColor (\n            interpolation = \"vertex\"\n        )\n";
print  $fh "        color3f[] primvars:displayColor.timeSamples = {\n";
for my $k (0..$Nlat-1) {
    my $den = defined $gmax ? $gmax : ($P10->slice(":,($k)")->max || 1e-300);
    my @li  = ( ($P10->slice(":,($k)") / $den)->clip(0,1) * ($NLUT-1) + 0.5 )->long->list;
    my $str = join(",", map { my $c=$LUT[$_]; sprintf("(%.3f,%.3f,%.3f)",@$c) } @li);
    print $fh "            $codes[$k]: [", $str, "],\n";
}
print  $fh "        }\n    }\n";

if ($axes > 0) {
    my ($L,$W) = ($axes, $axes/50);
    print $fh <<"USDA";
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
if ($cbar) {
    # 独立プリム: usdview の outliner で Axes と別に表示/非表示。頭の右(+X)に立てる
    # (ワールド固定=シーンと一緒に回る)。勾配は cortex と同じ cmap＋灰オーバレイ。
    my $NBAR = 64;
    my ($xl,$xr,$yy) = (95, 108, 0);   # バー左右x, y平面
    my ($z0,$z1)     = (-10, 70);       # 下(t=0)→上(t=1)
    my (@bpts,@bcol);
    for my $i (0..$NBAR) {
        my $tt = $i/$NBAR;  my $z = $z0 + ($z1-$z0)*$tt;
        push @bpts, sprintf("(%.2f,%.2f,%.2f)",$xl,$yy,$z), sprintf("(%.2f,%.2f,%.2f)",$xr,$yy,$z);
        my $c = _overlay_rgb($tt, $cmap, $grey, $threshold);
        my $cs = sprintf("(%.3f,%.3f,%.3f)", @$c);
        push @bcol, $cs, $cs;
    }
    my @bidx; for my $i (0..$NBAR-1) { push @bidx, 2*$i, 2*$i+1, 2*$i+3, 2*$i+2; }  # quad
    # 目盛りノッチ: 0/0.5/1.0=白、threshold=シアン(長め=灰↔カラー境界を明示)
    my @allt = ([0,0.9,0.9,0.95,6],[0.5,0.9,0.9,0.95,6],[1,0.9,0.9,0.95,6]);
    push @allt, [$threshold, 0.1,0.95,0.95, 11] if defined $grey;
    my (@tp,@tc,@tcol);
    for my $a (@allt) {
        my ($frac,$r,$g,$b,$len) = @$a;  my $z = $z0 + ($z1-$z0)*$frac;
        push @tp, sprintf("(%.2f,%.2f,%.2f)",$xr,$yy,$z), sprintf("(%.2f,%.2f,%.2f)",$xr+$len,$yy,$z);
        push @tc, 2;  push @tcol, sprintf("(%.2f,%.2f,%.2f)",$r,$g,$b);
    }
    print $fh "    def Xform \"ColorBar\"\n    {\n";
    print $fh "        def Mesh \"Bar\"\n        {\n";
    print $fh "            uniform bool doubleSided = 1\n";
    print $fh "            int[] faceVertexCounts = [", join(",",(4) x $NBAR), "]\n";
    print $fh "            int[] faceVertexIndices = [", join(",",@bidx), "]\n";
    print $fh "            point3f[] points = [", join(",",@bpts), "]\n";
    print $fh "            color3f[] primvars:displayColor = [", join(",",@bcol), "] (interpolation = \"vertex\")\n";
    print $fh "        }\n";
    print $fh "        def BasisCurves \"Ticks\"\n        {\n";
    print $fh "            uniform token type = \"linear\"\n";
    print $fh "            int[] curveVertexCounts = [", join(",",@tc), "]\n";
    print $fh "            point3f[] points = [", join(",",@tp), "]\n";
    print $fh "            color3f[] primvars:displayColor = [", join(",",@tcol), "] (interpolation = \"uniform\")\n";
    print $fh "            float[] widths = [1.5] (interpolation = \"constant\")\n";
    print $fh "        }\n";
    print $fh "    }\n";
}
print  $fh "}\n";
close $fh;

printf "wrote %s  (%d frames, codes %d..%d = %.1f..%.1f ms, norm=%s)\n",
    $out, $Nlat, $codes[0], $codes[-1], $codes[0]/10, $codes[-1]/10, $norm;
printf "open in usdview / Keynote; scrub time %d..%d, pause & rotate.\n", $codes[0], $codes[-1];
