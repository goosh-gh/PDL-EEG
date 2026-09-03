#!/usr/bin/env perl
# ---------------------------------------------------------------------------
# sep_latency_movie.pl
#
# 加算平均 SEP を「潜時 1 フレーム」で走る 3 パネル動画にする:
#
#     [ 左 ] F4 / C4 の波形（縦カーソルが現在潜時）
#     [ 中 ] その潜時の 2D 頭皮トポマップ（PDL::EEG::MAP2D）
#     [ 右 ] ECD Goodness of Fit(潜時) 曲線（縦カーソルが現在潜時）
#
# 各潜時で 1 枚 PNG を描き、末尾で ffmpeg が mp4/gif に束ねる。トポマップは
# MAP2D の _ax/_fig 埋め込みフック（plot_topomap_panels と同じ経路）で中央
# パネルへ直接描く＝PNG 往復も imread も不要。カラーバーは MAP2D と同じ
# RdBu(_colorize)で右端に共有 1 本。波形と GoF は P:G:C の line + axvline。
#
# スタックは前スレのまま: PDL::EEG::MAP2D + PDL::Graphics::Cairo。
#
# 入力 = P:EEG22/23 と同じ MNE エクスポート:
#   --ave      sep_ave.txt    : np.savetxt(ev.data) の行列。ev.data=(n_ch, n_time) なので
#                               1 行 = 1 チャネル(その全時間サンプル) = [channel, sampling-point]。
#   --names    sep_labels.txt : 1 行 1 チャネル名。
#   --sfreq    10000          : サンプリング周波数 [Hz]。
#   --tmin     -0.05          : 先頭サンプルの時刻 [s]。
#   --montage  xxx.elc        : ASA .elc（MAP2D の電極位置。例 standard_1020_*.elc）。
#
# GoF は「自前で計算済みのものを読む」= SEP HFO 01 / NYHead 単一ダイポール
# スキャンで出した曲線と一致させる:
#   --gof gof.tsv : 2 列 "latency_ms<TAB>gof"。任意の潜時格子で可。
#
# 例:
#   perl sep_latency_movie.pl \
#     --ave sep_ave.txt --names sep_labels.txt --sfreq 10000 --tmin -0.05 \
#     --montage standard_1020_eog_nose.elc --gof gof.tsv \
#     --wave-chans F4,C4 \
#     --ctx-min-ms -5 --ctx-max-ms 40 \
#     --anim-min-ms 10 --anim-max-ms 32 --step-ms 0.2 \
#     --out sep_latency.mp4 --fps 12
# ---------------------------------------------------------------------------
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Graphics::Cairo qw(figure);
use PDL::EEG::MAP2D qw(plot_topomap);
use Getopt::Long;
use File::Temp qw(tempdir);
use File::Path qw(make_path);

# --------------------------- args ---------------------------
my %o = (
    'wave-chans' => 'F4,C4',
    'ctx-min-ms' => -5,   'ctx-max-ms' => 40,
    'anim-min-ms'=> 10,   'anim-max-ms'=> 32,  'step-ms' => 0.2,
    'clim-pct'   => 99,   'clim'       => undef,
    'wave-clim'  => 'auto','wave-blank-ms' => 2.0, 'neg-up' => 1,
    'wave-color' => 'polarity','polarity-ms' => 20,
    'outdir'     => undef,'out'        => 'sep_latency.mp4',
    'fps'        => 12,   'unit'       => 'uV',
    'contours'   => 6,    'keep-frames'=> 0,
    'figw'       => 900,  'figh'       => 560,  'wspace' => 0.012,
    'cbar-width' => 11,   'cbar-pad'   => 6,    'cbar-attach' => 'map-bottom',
    'gof-min-ms' => undef,'gof-max-ms' => undef,
    'left'       => 0.045,'right'      => 0.955,'top' => 0.90, 'bottom' => 0.15,
);
GetOptions(\%o,
    'ave=s','names=s','sfreq=f','tmin=f','montage=s','gof=s',
    'wave-chans=s','ctx-min-ms=f','ctx-max-ms=f',
    'anim-min-ms=f','anim-max-ms=f','step-ms=f',
    'clim=f','clim-pct=f','unit=s','contours=i',
    'wave-clim=s','wave-blank-ms=f','neg-up!','wave-color=s','polarity-ms=f',
    'outdir=s','out=s','fps=f','keep-frames!','figw=i','figh=i','wspace=f',
    'cbar-width=f','cbar-pad=f','cbar-attach=s',
    'gof-min-ms=f','gof-max-ms=f',
    'left=f','right=f','top=f','bottom=f',
) or die "bad options\n";
for my $req (qw(ave names sfreq tmin montage gof)) {
    defined $o{$req} or die "--$req is required\n";
}

# --------------------------- load evoked ---------------------------
# sep_ave.txt: n_ch 行 × n_time 列（np.savetxt of (n_ch,n_time)）
my @rows;
open my $fh, '<', $o{ave} or die "open $o{ave}: $!";
while (<$fh>) { next unless /\S/; push @rows, [split] }
close $fh;
my $avg = pdl(\@rows)->xchg(0,1)->sever;   # (Ne, Nt), V   (pdl(AoA) は内側=dim0=時間)
my ($Ne, $Nt) = $avg->dims;

my @labels;
open my $lh, '<', $o{names} or die "open $o{names}: $!";
while (<$lh>) { s/\s+$//; next unless /\S/; push @labels, $_ }
close $lh;
@labels == $Ne or die "labels (".scalar(@labels).") != channels ($Ne)\n";

my $sf   = $o{sfreq};
my $tmin = $o{tmin};
my $tms  = ($tmin + sequence($Nt) / $sf) * 1000;   # 各サンプルの潜時[ms]

sub ms2idx {                                       # 潜時[ms] → サンプル index
    my $ms = shift;
    my $i  = sclr(rint(($ms/1000 - $tmin) * $sf));
    $i = 0        if $i < 0;
    $i = $Nt - 1  if $i > $Nt - 1;
    return $i;
}

# --------------------------- GoF ---------------------------
# 2 列 "ms  gof"。ヘッダ/コメント行(数字で始まらない行)は読み飛ばす
my (@gm, @gv);
open my $gh, '<', $o{gof} or die "open $o{gof}: $!";
while (<$gh>) {
    next if /^\s*#/;
    my @f = split;
    next unless @f >= 2 && $f[0] =~ /^[-+]?\d/ && $f[1] =~ /^[-+]?[\d.]/;
    push @gm, $f[0]; push @gv, $f[1];
}
close $gh;
@gm >= 2 or die "no numeric rows in $o{gof}\n";
my $gms = pdl(\@gm); my $gof = pdl(\@gv);          # ms, gof
my $gi = $gms->qsorti;
$gms = $gms->index($gi)->sever; $gof = $gof->index($gi)->sever;

# --------------------------- frames + 固定 clim ---------------------------
my $i0 = ms2idx($o{'anim-min-ms'});
my $i1 = ms2idx($o{'anim-max-ms'});
$i1 > $i0 or die "anim window empty\n";
my $stride = int(rint($o{'step-ms'}/1000 * $sf)); $stride = 1 if $stride < 1;
my @frames;
for (my $i = $i0; $i <= $i1; $i += $stride) { push @frames, $i }
@frames >= 2 or die "too few frames\n";

my $vmax;                                          # µV
if (defined $o{clim}) { $vmax = $o{clim} }
else {
    my $s = $avg(:, $i0:$i1)->abs->flat->qsort;    # V, sorted
    my $k = int(($o{'clim-pct'}/100) * ($s->nelem - 1));
    $vmax = sclr($s($k)) * 1e6;                    # µV
    $vmax = sclr($avg(:,$i0:$i1)->abs->max)*1e6 if $vmax <= 0;
}
printf "frames: %d  (%.1f..%.1f ms, step %.2f ms)   clim=+-%.2f %s\n",
    scalar(@frames), $tms->at($i0), $tms->at($i1), $o{'step-ms'}, $vmax, $o{unit};

# --------------------------- 波形チャネル ---------------------------
my @wch = split /,/, $o{'wave-chans'};
my %lab2i = map { $labels[$_] => $_ } 0..$#labels;
for my $c (@wch) { exists $lab2i{$c} or die "wave channel '$c' not found\n" }
my @wcolor = ('#d62728','#1f77b4','#2ca02c','#9467bd','#8c564b');   # --wave-color cycle 用

# 波形の色分け(既定 --wave-color polarity): N20 潜時(--polarity-ms 既定20)での符号で、
# 陰性=寒色、陽性=暖色。同符号内は離散パレット順(濃淡でなく別色)。白背景で判別しやすく黄色回避。
my @COOL = ('#1f77b4', '#2ca02c', '#17a2b8', '#555555');   # 青・緑・シアン・濃グレー
my @WARM = ('#c0392b', '#6a3d9a', '#d95f02', '#c71585');   # 濃赤・濃紫・濃橙・マゼンタ
my @wave_colors;
{
    my $pol_idx = ms2idx($o{'polarity-ms'});
    my ($nc, $nw) = (0, 0);
    for my $k (0..$#wch) {
        if (lc $o{'wave-color'} eq 'cycle') {
            push @wave_colors, $wcolor[$k % @wcolor];
        } else {
            my $val = sclr($avg(($lab2i{$wch[$k]}), ($pol_idx)));   # 符号だけ見る(単位無関係)
            push @wave_colors, ($val < 0) ? $COOL[$nc++ % @COOL] : $WARM[$nw++ % @WARM];
        }
    }
    printf "wave colors @ %.1f ms: %s\n", $o{'polarity-ms'},
        join(' ', map { "$wch[$_]=$wave_colors[$_]" } 0..$#wch);
}

# 波形 y スケール: トポの clim(±vmax)と切り離す。
#   auto  = 表示窓内で描画チャネルの実ピークに合わせる（刺激アーチ窓 |t|<blank は
#           除外＝アーチは振り切れるが F2 等の実信号は切れない）。既定。
#   match = トポと同じ ±vmax(前回の挙動)。
#   <数値>= 明示 ±値。
my $wave_vmax;
my $wc = lc $o{'wave-clim'};
if    ($wc eq 'match')            { $wave_vmax = $vmax }
elsif ($wc =~ /^[-+]?[\d.]+$/)    { $wave_vmax = abs($wc) + 0 }
else {                              # auto
    my $blank = $o{'wave-blank-ms'};
    my $inwin = ($tms >= $o{'ctx-min-ms'}) & ($tms <= $o{'ctx-max-ms'}) & (abs($tms) >= $blank);
    my $mm = 0;
    for my $c (@wch) {
        my $y   = ($avg(($lab2i{$c}), :) * 1e6);
        my $sel = $y->where($inwin);
        my $m   = $sel->nelem ? sclr($sel->abs->max) : 0;
        $mm = $m if $m > $mm;
    }
    $wave_vmax = $mm > 0 ? $mm * 1.10 : $vmax;
}
printf "wave y-scale: +-%.2f %s (%s)\n", $wave_vmax, $o{unit}, $wc;

# context 範囲（波形・GoF の x 軸）
my $cx0  = ms2idx($o{'ctx-min-ms'});
my $cx1  = ms2idx($o{'ctx-max-ms'});
my $tctx = $tms($cx0:$cx1)->sever;

# 波形 ylim はトポマップと同じ ±vmax(µV)。刺激アーチは振り切れてよい(=クリップ)。

# GoF x 範囲（既定=ctx。--gof-min-ms/--gof-max-ms で時間軸を狭めて拡大できる）
my $gof_lo = defined $o{'gof-min-ms'} ? $o{'gof-min-ms'} : $o{'ctx-min-ms'};
my $gof_hi = defined $o{'gof-max-ms'} ? $o{'gof-max-ms'} : $o{'ctx-max-ms'};

# GoF ylim（GoF 表示窓内）
my $gmask = ($gms >= $gof_lo) & ($gms <= $gof_hi);
my ($glo, $ghi);
if ($gmask->any) { $glo = sclr($gof->where($gmask)->min); $ghi = sclr($gof->where($gmask)->max) }
else             { ($glo,$ghi) = (sclr($gof->min), sclr($gof->max)) }
my $gpad = 0.05 * (($ghi - $glo) + 1e-9); $glo -= $gpad; $ghi += $gpad;
my $gmx  = $gms->where($gmask)->sever;
my $gvy  = $gof->where($gmask)->sever;

# GoF y 目盛りは 0.1 刻み
my @gyt; { my $k = int($glo*10); $k++ while $k/10 < $glo - 1e-9;
           for (; $k/10 <= $ghi + 1e-9; $k++) { push @gyt, $k/10 } }

# --------------------------- 出力ディレクトリ ---------------------------
my $framedir = $o{outdir} // tempdir('sepmov_XXXX', TMPDIR => 1, CLEANUP => 0);
make_path($framedir);

# カラーバー用 cmap: MAP2D のトポマップと同一の RdBu(青=低→白→赤=高)を
# P:G:C の ColorMap オブジェクトとして構築し、ネイティブ colorbar に渡す。
# (色をトポマップと厳密に一致させるため。anchors は MAP2D の @RDBU_R を写したもの)
my @RDBU = (
    [0.019,0.188,0.380],[0.129,0.400,0.674],[0.262,0.576,0.764],
    [0.572,0.772,0.870],[0.819,0.898,0.941],[0.968,0.968,0.968],
    [0.992,0.858,0.780],[0.956,0.647,0.509],[0.839,0.376,0.302],
    [0.698,0.094,0.168],[0.403,0.000,0.121],
);
my @cmap_pts = map { [ $_/$#RDBU, @{$RDBU[$_]} ] } 0..$#RDBU;
my $rdbu_cmap = bless { name=>'map2d_rdbu', pts=>\@cmap_pts },
                'PDL::Graphics::Cairo::ColorMap';

# --------------------------- フレームループ ---------------------------
my $fi = 0;
for my $s (@frames) {
    my $t   = $tms->at($s);
    my $gat = sclr(pdl($t)->interpol($gms, $gof));   # 現在潜時の GoF

    my $fig = figure(width => $o{figw}, height => $o{figh});
    my $gs  = $fig->add_gridspec(1, 3, wspace => 0);   # 位置は下で set_position 上書き
    my $axw = $fig->add_subplot($gs->at(0,0));   # 波形
    my $axt = $fig->add_subplot($gs->at(0,1));   # トポマップ
    my $axg = $fig->add_subplot($gs->at(0,2));   # GoF
    # 明示配置。gridspec の隙間・set_aspect の中央寄せ(円の左右の空帯)を回避するため、
    # マップの描画域を正方形にして円で埋める。カラーバー分は箱の高さに足す(正方形は削らない)。
    {
        my ($W,$H)  = ($o{figw}, $o{figh});
        my $py = $o{bottom}; my $ph = $o{top} - $o{bottom};
        my $cbar_res = ($o{'cbar-attach'} eq 'map-bottom') ? 54 : 0;  # 下カラーバーの確保px
        my $Dpx  = $ph*$H - $cbar_res;                 # 円の直径(=正方形描画域の一辺)px
        my $mw   = ($Dpx + 10) / $W;                   # マップ箱の幅(frac, 軸マージン余裕込み)
        my $gap  = $o{wspace};                          # パネル間の隙間(frac)

        my $avail = $o{right} - $o{left};
        my $col   = ($avail - 2*$gap) / 3;        # 3 等分の 1 列幅
        my $wf_w  = $col;
        my $gof_w = $col;
        my $wf_x  = $o{left};
        my $mp_x  = $wf_x + $col + $gap;          # 中央列の左端
        # 中央列(幅 $col)の中でマップ正方形($mw)を中央寄せ。列より広ければ列いっぱい。
        if ($mw < $col) { $mp_x += ($col - $mw) / 2; }
        else            { $mw = $col; }           # 列に収まらなければ列幅に合わせる
        my $gof_x = $o{left} + 2*($col + $gap);   # 右列の左端(3等分の3番目)

#        my $avail = $o{right} - $o{left};
#        my $rest  = $avail - $mw - 2*$gap;
#        my $wf_w  = $rest/2; my $gof_w = $rest/2;
#        my $wf_x  = $o{left};
#        my $mp_x  = $wf_x + $wf_w + $gap;
#        my $gof_x = $mp_x + $mw + $gap;

        $axw->set_position([$wf_x,  $py, $wf_w,  $ph]);
        $axt->set_position([$mp_x,  $py, $mw,    $ph]);
        $axg->set_position([$gof_x, $py, $gof_w, $ph]);
    }

    # 左: 波形 + カーソル
    for my $k (0..$#wch) {
        my $ci = $lab2i{$wch[$k]};
        my $y  = ($avg(($ci), $cx0:$cx1) * 1e6)->sever;
        $axw->line($tctx, $y, color => $wave_colors[$k], lw => 1.4);
        my $label_top = $o{'neg-up'} ? -$wave_vmax : $wave_vmax;   # neg-up なら上端=負
        $axw->text($o{'ctx-min-ms'} + 1, $label_top * (0.92 - 0.13*$k),
                   $wch[$k], color => $wave_colors[$k],
                   ha => 'left', halign => 'left', fontsize => 11);
    }
    $axw->axhline(0,  color => '#cccccc', lw => 0.6);
    $axw->axvline(0,  color => '#999999', lw => 0.8);
    $axw->axvline($t, color => 'black',   lw => 1.4);
    $axw->xlim($o{'ctx-min-ms'}, $o{'ctx-max-ms'});
    $axw->ylim(-$wave_vmax, $wave_vmax);        # 波形は自前スケール(--wave-clim)
    $axw->invert_yaxis() if $o{'neg-up'};        # 生理/ERP 慣習: 陰性を上(既定)。--no-neg-up で正を上
    $axw->set_title('Waveforms');
    $axw->xlabel('Time (ms)');
    $axw->ylabel('Amplitude ('.$o{unit}.')');

    # 中: トポマップ（MAP2D を軸へ埋め込み）
    my $vec = ($avg(:, ($s)) * 1e6)->sever;          # (Ne) µV
    plot_topomap(
        values   => $vec,
        labels   => \@labels,
        montage  => $o{montage},
        clim     => [ -$vmax, $vmax ],
        contours => $o{contours},
        _ax      => $axt,
        _fig     => $fig,
        colorbar => 0,
    );

    # 右: GoF + カーソル
    $axg->line($gmx, $gvy, color => '#333333', lw => 1.6);
    $axg->axvline(0,  color => '#999999', lw => 0.8);
    $axg->axvline($t, color => 'black',   lw => 1.4);
    $axg->xlim($gof_lo, $gof_hi);
    $axg->ylim($glo, $ghi);
    $axg->set_yticks(pdl(\@gyt)) if @gyt;            # 0.1 刻み
    $axg->set_title('ECD Goodness of Fit');
    $axg->xlabel('Time (ms)');
    $axg->ylabel('GoF');

    # カラーバー = P:G:C ネイティブ(数値目盛りはビルド非依存)。
    if ($o{'cbar-attach'} eq 'gof') {
        # フォールバック(パッチ不要): GoF 軸(表示ON)の左マージンに縦カラーバー。
        $axg->_colorbar({ cmap => $rdbu_cmap, vmin => -$vmax, vmax => $vmax });
        $axg->colorbar(location => 'left', label => $o{unit},
                       width => $o{'cbar-width'},
                       pad => ($o{'cbar-pad'} < 20 ? 26 : $o{'cbar-pad'}));
    } else {
        # 既定 map-bottom: 地図の下に横カラーバー。横幅を食わず、GoF ラベルと衝突しない。
        # トポは軸OFF なので _colorbar_force を立てて描かせる(要 pgc_colorbar_force.patch)。
        $axt->{_colorbar_force} = 1;
        $axt->_colorbar({ cmap => $rdbu_cmap, vmin => -$vmax, vmax => $vmax });
        $axt->colorbar(location => 'bottom', label => $o{unit},
                       width => $o{'cbar-width'}, pad => $o{'cbar-pad'});
    }

    $fig->suptitle(sprintf("t = %.1f ms      GoF = %.3f", $t, $gat));
    $fig->save(sprintf("%s/frame_%04d.png", $framedir, $fi));

    $fi++;
    printf "  frame %3d/%d  t=%6.1f ms  GoF=%.3f\n", $fi, scalar(@frames), $t, $gat
        if $fi == 1 || $fi % 20 == 0;
}

# --------------------------- ffmpeg で束ねる ---------------------------
my $out = $o{out};
my $cmd;
if ($out =~ /\.gif$/i) {
    my $pal = "$framedir/palette.png";
    system("ffmpeg -y -loglevel error -framerate $o{fps} -i $framedir/frame_%04d.png -vf palettegen $pal");
    $cmd = "ffmpeg -y -loglevel error -framerate $o{fps} -i $framedir/frame_%04d.png -i $pal -lavfi paletteuse $out";
} else {
    $cmd = "ffmpeg -y -loglevel error -framerate $o{fps} -i $framedir/frame_%04d.png -c:v libx264 -pix_fmt yuv420p $out";
}
my $rc = system($cmd);
if ($rc == 0) { print "[done] saved: $out\n" }
else          { warn "[warn] ffmpeg failed; frames kept in $framedir\n" }

unless ($o{'keep-frames'} || $rc != 0) {
    unlink glob("$framedir/frame_*.png");
    unlink "$framedir/palette.png";
    rmdir $framedir unless defined $o{outdir};
} else {
    print "[info] frames kept in $framedir\n";
}
