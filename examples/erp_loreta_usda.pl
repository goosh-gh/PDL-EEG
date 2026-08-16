#!/usr/bin/env perl
#
# erp_loreta_usda.pl
#   多チャンネル・多潜時 LORETA -> NYHead -> 潜時アニメ USDA(1本)
#
#   通常 ERP(1000Hz, 10-20 等)の加算波形を eeg.pm(NK_READ_ONEAVGFILE2)で
#   直接読み込み、[start_latency, end_latency] の各潜時で皮質源パワーを解いて、
#   色だけが潜時で変化する単一 .usda(cortex mesh + displayColor.timeSamples)を書く。
#   ※MNE 経由の text ダンプは経由しない(eeg.pm は既に Perl 構造なので round-trip 不要)。
#
# 使い方:
#   perl erp_loreta_usda.pl INFILE START_MS END_MS [options]
#     INFILE    加算ファイル(eeg.pm が読む)
#     START_MS  開始潜時(ms, 刺激=0)
#     END_MS    終了潜時(ms)
#   出力: "${INFILE}${START_MS}_${END_MS}_loreta.usda"
#
# options(既定で 3 引数だけで走る):
#   --mat PATH        NYHead リードフィールド sa_nyhead.mat
#   --method NAME     sloreta(既定)|eloreta|mne     ERP は分布源(sloreta/mne)が無難
#   --res NAME        cortex10K(既定)|cortex5K|cortex2K  ※75K は面数が u16 を超えるので不可
#   --step MS         潜時ステップ(既定 0 = 1 サンプル毎)。窓が広い時は間引きに使う
#   --cmap NAME       inferno(既定)|viridis|hot
#   --threshold F     [0,1] これ未満は灰色(既定 0)
#   --norm global|frame  global(既定)=窓最大で規格化 / frame=各潜時で規格化
#   --no-axes         原点の XYZ 軸トライアドを出さない
#   --no-colorbar     カラーバーを出さない
#
# ---- 実機で合わせる箇所(このコンテナでは検証不可)-------------------------
#   (1) eeg.pm の @INC 追加と load 方法(下の "use eeg;")
#   (2) $rd の中身: (ch x time) の AoA を仮定(ref なら PDL としても扱う)。
#       実データ 1 本で $data->info を見て dim0=電極/dim1=時間 を確認。
#   (3) 参照合わせ(avg_reference / centering)と inverse_operator の引数名は
#       sep_n20_inverse.pl と一致させる(下記コメント参照)。
#   (4) USDA writer は sweep_usda_anim.pl と二重管理になる。理想は両者から呼ぶ
#       共通 sub に括り出すこと。ここでは仕様どおり自前で書いている。
# --------------------------------------------------------------------------

use strict; use warnings; use PDL; use File::Basename; use Getopt::Long;

# (1) eeg.pm を読む。実機では topomap2d.pl と同じ流儀で:
#   use lib '/Users/goosh/src/...';   # eeg.pm / mat.pm のあるパス
#   use eeg;
# ここでは private module 依存なので require を実行時に。
BEGIN { eval { require eeg; eeg->import if eeg->can('import'); }; }

use lib '/Users/goosh/src/PDL_IO_NYHead/lib/'; use PDL::IO::NYHead;
use lib '/Users/goosh/src/PDL_EEG/lib/'; use PDL::EEG::Inverse::MinimumNorm
    qw(inverse_operator source_power avg_reference);

# ---- defaults -------------------------------------------------------------
my $MAT       = '/Users/goosh/src/NYHead/sa_nyhead.mat';
my $method    = 'sloreta';
my $res       = 'cortex10K';
my $step_ms   = 0;          # 0 = 1 サンプル毎
my $cmap      = 'inferno';
my $threshold = 0.0;
my $norm      = 'global';
my $fps       = 0;          # framesPerSecond。0=自動(=sfreq/stride) → Frame を 1 サンプル刻みに。
                            # 明示値を渡すと動画的再生になるが 1ms 刻みに止まれなくなる
my $axes      = 1;
my $colorbar  = 1;
my $base_grey = 0.45;       # 閾値下の灰色レベル

GetOptions(
    'mat=s'       => \$MAT,
    'method=s'    => \$method,
    'res=s'       => \$res,
    'step=f'      => \$step_ms,
    'cmap=s'      => \$cmap,
    'threshold=f' => \$threshold,
    'norm=s'      => \$norm,
    'fps=f'       => \$fps,
    'grey=f'      => \$base_grey,
    'axes!'       => \$axes,
    'colorbar!'   => \$colorbar,
) or die "bad options\n";

my ($infile, $start_lat, $end_lat) = @ARGV;
die "usage: $0 INFILE START_MS END_MS [--mat .. --method .. --res .. --step ..]\n"
    unless defined $end_lat;

# 出力名は指定どおり(拡張子は剥がさない)
my $outfile = "$infile" ."."."$start_lat" . "_" . "$end_lat" . "_loreta.usda";

# ---- 1. eeg.pm で加算読込 -------------------------------------------------
my ($rd, $rh) = eeg::NK_READ_ONEAVGFILE2("$infile", "trgCH", 0);

my $nch   = $rh->{general}{ChannelN};
my $sfreq = $rh->{average}{sampling_freq};
my $pre   = $rh->{average}{pretrigger};        # baseline サンプル数(t=0 の index)
my @labels = map { $rh->{ChannelName}{$_} } 0 .. $nch - 1;

# データを (Ne, Nt) piddle に正規化(dim0=電極, dim1=時間)
my $data = (ref($rd) eq 'PDL') ? $rd->copy : pdl($rd);
if ($data->dim(0) != $nch && $data->dim(1) == $nch) {
    $data = $data->xchg(0, 1)->sever;          # (Nt,Nch) -> (Nch,Nt)
}
die "channel dim ($nch) not found in data dims [".join('x',$data->dims)."]\n"
    unless $data->dim(0) == $nch;
my $nt = $data->dim(1);

printf STDERR "loaded: %d ch x %d samp @ %g Hz, pretrigger=%d samp (t=0)\n",
    $nch, $nt, $sfreq, $pre;

# ---- 2. 潜時窓 -> サンプル index ------------------------------------------
my $spms = $sfreq / 1000.0;                     # samples per ms
my $i0 = int($pre + $start_lat * $spms + 0.5);
my $i1 = int($pre + $end_lat   * $spms + 0.5);
($i0, $i1) = ($i1, $i0) if $i0 > $i1;
$i0 = 0        if $i0 < 0;
$i1 = $nt - 1  if $i1 > $nt - 1;
my $stride = $step_ms > 0 ? int($step_ms * $spms + 0.5) : 1;
$stride = 1 if $stride < 1;
# framesPerSecond: 自動なら sfreq/stride。timeCodesPerSecond(1000)/fps = 1フレーム=1サンプル。
my $fps_out = ($fps > 0) ? $fps : $sfreq / $stride;

my @frame_idx;
for (my $i = $i0; $i <= $i1; $i += $stride) { push @frame_idx, $i; }
printf STDERR "latency window: %g..%g ms -> samp %d..%d, stride %d, %d frames\n",
    $start_lat, $end_lat, $i0, $i1, $stride, scalar(@frame_idx);
printf STDERR "framesPerSecond = %g (timeCodesPerSecond 1000 / fps = %g ms per step)\n",
    $fps_out, 1000.0 / $fps_out;

# ---- 3. 電極 -> NYHead リードフィールド対応(ラベル一致)-------------------
# LORETA では .elc の座標は不要。leadfield 行(NYHead の電極集合)にデータの
# チャンネルをラベルで対応づけるだけ。旧 10-20 <-> 10-10 別名を吸収。
my $ny = PDL::IO::NYHead->new($MAT);

my $ny_labels = $ny->electrode_labels;          # arrayref(231 labels)
die "NYHead electrode_labels empty (h5dump 不在?)\n"
    unless ref($ny_labels) eq 'ARRAY' && @$ny_labels;

my %ALIAS = (T3 => 'T7', T4 => 'T8', T5 => 'P7', T6 => 'P8',
             LM => 'M1', RM => 'M2');
my %ny_idx;
for my $j (0 .. $#$ny_labels) {
    my $u = uc $ny_labels->[$j];
    $u =~ s/\s+//g;
    $ny_idx{$u} = $j unless exists $ny_idx{$u};
}

my (@data_rows, @ny_rows, @used, @skipped);
for my $i (0 .. $#labels) {
    my $u = uc $labels[$i];
    $u =~ s/\s+//g;
    $u = $ALIAS{$u} if exists $ALIAS{$u} && !exists $ny_idx{$u};
    if (exists $ny_idx{$u}) {
        push @data_rows, $i;
        push @ny_rows,   $ny_idx{$u};
        push @used,      $labels[$i];
    } else {
        push @skipped, $labels[$i];
    }
}
die "no electrodes matched NYHead\n" unless @ny_rows;
printf STDERR "matched %d electrodes: %s\n", scalar(@used), join(',', @used);
warn "skipped (no NYHead match): @skipped\n" if @skipped;

my $drows = pdl(long, @data_rows);
my $nrows = pdl(long, @ny_rows);

# ---- 4. 逆作用素を 1 回だけ構築 -------------------------------------------
# leadfield は NYHead で (Ne=231, Ns=74382)=(電極,源)。inverse_operator は数学 K=(Ns,Ne)
# を要求する(Ns>Ne を assert)。電極を subset して転置し (Ns, Ne_sub) で渡す。
my $K = $ny->leadfield->dice_axis(0, $nrows)->transpose->sever;   # (Ns, Ne_sub)
$K = avg_reference($K);   # subset 上で再 CAR: 各源の電極平均を引き K と data の参照を揃える
my $op = inverse_operator($K, method => $method);   # ref=>'car' 既定, sloreta は標準化込み
my $Ns = $K->dim(0);

# データ側も同じ電極・同じ順に切り出し(以降 (Ne_sub, Nt))
my $dsub = $data->dice_axis(0, $drows)->sever;

# ---- 5. 潜時ごとに源パワー -----------------------------------------------
my @frames;         # 各要素 = 源パワー (Ns)
my @lat_ms;
my $gmax = 0;
for my $i (@frame_idx) {
    my $b = $dsub->slice(":,($i)")->sever;    # (Ne_sub) その潜時のトポ
    $b = $b - $b->avg;                        # トポを CAR(電極平均を引く)
                                              #   ※centering() は n×n 行列生成子なので不可
    my $pw = source_power($op, $b);           # (Ns) 標準化源パワー(>=0)
    push @frames, $pw;
    push @lat_ms, ($i - $pre) / $spms;
    my $m = $pw->max->sclr;
    $gmax = $m if $m > $gmax;
}
$gmax = 1 if $gmax <= 0;

# ---- 6. 皮質メッシュ & 正規化 & 頂点色 -----------------------------------
my $cx  = $ny->cortex($res);                  # {vc=>(Nv,3), tri=>(Nf,3), in_from=>(Nv)}
my $vc  = $cx->{vc};
my $tri = $cx->{tri};
my $inf = $cx->{in_from};                      # local -> 75K(0-based)。75K なら identity
my $Nv  = $vc->dim(0);

my @frame_rgb;                                 # 各要素 = (Nv,3) in [0,1]
my @peak_lines;                                # 診断: 潜時ごとの peak
for my $f (0 .. $#frames) {
    my $pw = $frames[$f];
    # 75K -> res 頂点へ
    my $pv = (defined $inf) ? $pw->index($inf) : $pw;   # (Nv)
    my $norm_max = ($norm eq 'frame') ? ($pv->max->sclr || 1) : $gmax;
    my $t = ($pv / $norm_max)->clip(0, 1);              # (Nv) in [0,1]
    push @frame_rgb, overlay_rgb($t, $cmap, $threshold, $base_grey);

    # 診断: peak 頂点(75K)と HO ラベル
    my $vmax = $pw->maximum_ind->sclr;
    my $area = eval { $ny->can('area_of_vertex') ? $ny->area_of_vertex($vmax) : '' } // '';
    push @peak_lines, sprintf("%.1f\t%d\t%s\t%.4g",
        $lat_ms[$f], $vmax, $area, $pw->max->sclr);
}
print STDERR "lat_ms\tpeak_v75k\tHO_area\tpeak_pow\n";
print STDERR "$_\n" for @peak_lines;

# ---- 7. USDA 書き出し -----------------------------------------------------
write_anim_usda(
    outfile  => $outfile,
    vc       => $vc,
    tri      => $tri,
    frames   => \@frame_rgb,
    lat_ms   => \@lat_ms,
    axes     => $axes,
    colorbar => $colorbar,
    cmap     => $cmap,
    threshold=> $threshold,
    grey     => $base_grey,
    norm     => $norm,
    fps      => $fps_out,
    infile   => $infile,
    method   => $method,
    res      => $res,
);
printf STDERR "wrote %s (%d frames, %d verts, %s)\n",
    $outfile, scalar(@frame_rgb), $Nv, $res;

# ==========================================================================
#  helpers
# ==========================================================================

# 値[0,1] -> RGB(Nv,3): 閾値下は灰色, 以上は cmap
sub overlay_rgb {
    my ($t, $name, $thr, $grey) = @_;
    my $rgb = cmap_rgb($t, $name);                 # (Nv,3)
    my $a   = ($t >= $thr)->double;                # 0/1 alpha (Nv)
    my $g   = zeroes($rgb->dims) + $grey;
    return $g * (1 - $a->dummy(1)) + $rgb * $a->dummy(1);
}

# 線形補間 colormap。値[0,1](Nv) -> (Nv,3)。ビューア一致は各自の _cmap_rgb と要調整。
sub cmap_rgb {
    my ($t, $name) = @_;
    my %A = (
        inferno => [
            [0.001,0.000,0.014],[0.259,0.039,0.408],[0.576,0.149,0.404],
            [0.867,0.318,0.227],[0.988,0.647,0.039],[0.988,0.998,0.645],
        ],
        viridis => [
            [0.267,0.005,0.329],[0.283,0.141,0.458],[0.254,0.265,0.530],
            [0.164,0.471,0.558],[0.478,0.821,0.318],[0.993,0.906,0.144],
        ],
        hot => [
            [0,0,0],[0.4,0,0],[0.8,0,0],[1,0.4,0],[1,0.8,0],[1,1,1],
        ],
    );
    my $anch = $A{$name} || $A{inferno};
    my $n = scalar @$anch;
    my $x = ($t->clip(0,1)) * ($n - 1);            # (Nv) in [0,n-1]
    my $lo = $x->floor->clip(0, $n - 2)->long;     # (Nv)
    my $fr = ($x - $lo)->dummy(1);                 # (Nv,1)
    # anchor 行列 (n,3)
    my $M = pdl($anch);                            # (3,n)? -> ensure (n,3)
    $M = $M->xchg(0,1)->sever if $M->dim(0) == 3 && $M->dim(1) == $n;
    my $c0 = $M->dice_axis(0, $lo);                # (Nv,3)
    my $c1 = $M->dice_axis(0, $lo + 1);            # (Nv,3)
    return $c0 * (1 - $fr) + $c1 * $fr;            # (Nv,3)
}

# (Nv,3) piddle -> "(x, y, z), (x, y, z), ..." 文字列(小数6桁)
sub triples_str {
    my ($p, $fmt) = @_;                            # $p: (N,3)
    $fmt ||= "%.6f";
    my $N = $p->dim(0);
    my @out;
    for my $k (0 .. $N - 1) {
        push @out, sprintf("($fmt, $fmt, $fmt)",
            $p->at($k,0), $p->at($k,1), $p->at($k,2));
    }
    return join(", ", @out);
}

sub write_anim_usda {
    my %o = @_;
    my $vc  = $o{vc};                              # (Nv,3) mm
    my $tri = $o{tri};                             # (Nf,3) 1-based(NYHead 規約)
    my $frames = $o{frames};                       # arrayref of (Nv,3)
    my $lat = $o{lat_ms};
    my $Nv  = $vc->dim(0);
    my $Nf  = $tri->dim(0);

    # 面: 1-based -> 0-based, per-face 3 連
    my $tri0 = ($tri - 1)->long;
    my $fvi  = join(", ", $tri0->xchg(0,1)->flat->list);   # f0a,f0b,f0c,f1a,...
    my $fvc  = join(", ", (3) x $Nf);
    my $pts  = triples_str($vc);

    # 時間コード = 潜時(ms)そのまま。timeCodesPerSecond=1000 → 1 tc = 1ms、
    # タイムラインが ms で読める(1000Hz なら整数; %g で非整数 sfreq も保持)。
    my @tc = map { sprintf("%g", $_) } @$lat;
    my $tc0 = $tc[0]; my $tc1 = $tc[-1];

    # bbox
    my $mn = $vc->minimum; my $mx = $vc->maximum;   # (3)
    my $ext = sprintf("(%.6f, %.6f, %.6f), (%.6f, %.6f, %.6f)",
        $mn->at(0),$mn->at(1),$mn->at(2), $mx->at(0),$mx->at(1),$mx->at(2));

    # doc / provenance(USD 文字列は \ と " をエスケープ)
    my $esc  = sub { my $s = shift // ''; $s =~ s/\\/\\\\/g; $s =~ s/"/\\"/g; $s };
    my $base = basename($o{infile} // 'unknown');
    my $in_e = $esc->($base);
    my $doc  = $esc->(sprintf('New York Head ERP LORETA | input=%s | %g-%g ms | method=%s | res=%s',
                              $base, $lat->[0], $lat->[-1], $o{method}//'', $o{res}//''));

    open my $fh, '>', $o{outfile} or die "open $o{outfile}: $!";
    print $fh <<"HDR";
#usda 1.0
(
    defaultPrim = "NYHead"
    upAxis = "Z"
    metersPerUnit = 0.001
    startTimeCode = $tc0
    endTimeCode = $tc1
    timeCodesPerSecond = 1000
    framesPerSecond = $o{fps}
    doc = "$doc"
    customLayerData = {
        string sep_input = "$in_e"
        string sep_method = "$o{method}"
        string sep_res = "$o{res}"
        string sep_cmap = "$o{cmap}"
        double sep_threshold = $o{threshold}
        double sep_base_grey = $o{grey}
        string sep_norm = "$o{norm}"
    }
)

def Xform "NYHead"
{
    def Mesh "cortex"
    {
        int[] faceVertexCounts = [$fvc]
        int[] faceVertexIndices = [$fvi]
        point3f[] points = [$pts]
        uniform token subdivisionScheme = "none"
        float3[] extent = [$ext]
HDR

    # displayColor 既定値 = 先頭フレーム(メタデータは値の後ろ)
    my $c0 = triples_str($frames->[0]);
    print $fh "        color3f[] primvars:displayColor = [$c0] (\n";
    print $fh "            interpolation = \"vertex\"\n";
    print $fh "        )\n";
    print $fh "        color3f[] primvars:displayColor.timeSamples = {\n";
    for my $f (0 .. $#$frames) {
        my $cs = triples_str($frames->[$f]);
        print $fh "            $tc[$f]: [$cs],\n";
    }
    print $fh "        }\n";
    print $fh "    }\n";   # end cortex

    write_axes($fh, $vc)     if $o{axes};
    write_colorbar($fh, $vc, $o{cmap}, $o{threshold}) if $o{colorbar};

    print $fh "}\n";        # end NYHead
    close $fh;
}

# 原点の XYZ 軸(赤/緑/青)BasisCurves
sub write_axes {
    my ($fh, $vc) = @_;
    my $L = ($vc->maximum - $vc->minimum)->max->sclr * 0.6;  # 軸長
    my @ax = ([[$L,0,0],[1,0,0]], [[0,$L,0],[0,1,0]], [[0,0,$L],[0,0,1]]);
    print $fh "    def Scope \"Axes\"\n    {\n";
    for my $k (0..2) {
        my ($tip, $col) = @{$ax[$k]};
        my $pts = sprintf("(0,0,0), (%.3f, %.3f, %.3f)", @$tip);
        my $c   = sprintf("(%.3f, %.3f, %.3f)", @$col);
        print $fh <<"AX";
        def BasisCurves "axis$k"
        {
            uniform token type = "linear"
            int[] curveVertexCounts = [2]
            point3f[] points = [$pts]
            color3f[] primvars:displayColor = [$c] (
                interpolation = "constant"
            )
            float[] widths = [1.5] (
                interpolation = "constant"
            )
        }
AX
    }
    print $fh "    }\n";
}

# ワールド固定カラーバー(回すと一緒に回る。画面固定は usdview 拡張の宿題)
sub write_colorbar {
    my ($fh, $vc, $name, $thr) = @_;
    my $mx = $vc->maximum; my $mn = $vc->minimum;
    my $x  = $mx->at(0) + 15;                       # 頭の右
    my $z0 = $mn->at(2); my $z1 = $mx->at(2);
    my $w  = 8; my $ns = 32;                        # リボン幅・分割
    my (@pts, @idx, @cols); my $vi = 0;
    for my $s (0 .. $ns - 1) {
        my $t0 = $s / $ns; my $t1 = ($s + 1) / $ns;
        my $za = $z0 + ($z1 - $z0) * $t0;
        my $zb = $z0 + ($z1 - $z0) * $t1;
        push @pts, "($x, 0, $za)", "(".($x+$w).", 0, $za)",
                   "(".($x+$w).", 0, $zb)", "($x, 0, $zb)";
        push @idx, $vi, $vi+1, $vi+2, $vi+3; $vi += 4;
        my $c = cmap_rgb(pdl([($t0+$t1)/2]), $name)->slice("(0),:");
        my $cc = sprintf("(%.3f, %.3f, %.3f)", $c->at(0),$c->at(1),$c->at(2));
        push @cols, ($cc) x 4;
    }
    my $P = join(", ", @pts);
    my $I = join(", ", @idx);
    my $C = join(", ", @cols);
    my $FC = join(", ", (4) x $ns);
    print $fh <<"CB";
    def Xform "ColorBar"
    {
        def Mesh "Bar"
        {
            int[] faceVertexCounts = [$FC]
            int[] faceVertexIndices = [$I]
            point3f[] points = [$P]
            color3f[] primvars:displayColor = [$C] (
                interpolation = "vertex"
            )
            uniform bool doubleSided = true
            uniform token subdivisionScheme = "none"
        }
    }
CB
}
