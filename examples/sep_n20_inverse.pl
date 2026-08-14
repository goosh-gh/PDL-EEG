#!/usr/bin/env perl
# ---------------------------------------------------------------------------
# sep_n20_inverse.pl  (P:EEG22)
#
#   実 SEP 加算波形の N20 潜時トポを NYHead リードフィールドに載せ、
#   MinimumNorm 逆問題で皮質源パワーを求め powers.dat に書き出す glue。
#
#   パイプライン:
#     [a] 加算波形(ch×time)を読み、N20 潜時でトポ b(電極ベクトル)を切り出す
#     [b] 記録電極ラベル → NYHead 231 リードフィールド行に label 対応付け → K_sub
#     [c] subset で CAR 再参照 → inverse_operator → source_power
#     [d] powers.dat(74382 行 = cortex75K 頂点順)+ ピーク頂点の HO 野名を診断出力
#
#   表示(グリグリ)は既存の GS3D ビューアにそのまま渡す:
#     demo_gs3d_nyhead.pl --mat sa_nyhead.mat --surf /sa/cortex10K \
#        --no-electrodes --source-power powers.dat --cmap inferno
#
#   ※ ピーク頂点が刺激対側の Postcentral Gyrus(S1)に載れば検証成功。
#     area 3b は Brodmann、HO は肉眼解剖なので直接は出ない(目標=postcentral)。
# ---------------------------------------------------------------------------
use strict;
use warnings;
use PDL;
use PDL::IO::NYHead;
use Getopt::Long;

# --- MinimumNorm は Mac 側 PDL-EEG に置いてある。無い環境でも読込診断まで走るよう遅延ロード ---
#     関数は @EXPORT_OK なので明示 import する(自動 export されない)。
my @INV_FN = qw(avg_reference inverse_operator source_power forward_project);
my $HAVE_INV = eval { require PDL::EEG::Inverse::MinimumNorm;
                      PDL::EEG::Inverse::MinimumNorm->import(@INV_FN); 1 };

# ---------------------------------------------------------------- options ---
my ($data_file, $labels_file, $times_file, $mat);
my $sfreq   = 0;        # Hz   (--times を使わないとき必須)
my $tmin    = 0;        # s    最初のサンプルの時刻(例 -0.05)
my $latency;            # ms   明示 N20 潜時。省略時は GFP ピークで自動
my $search  = "16,26";  # ms   自動ピックの探索窓
my $avg_ms  = 0;        # ms   トポを ±avg_ms 平均(0=単サンプル)
my $method  = "eloreta";
my $reg     = 0.05;     # Tikhonov reg_frac
my $out     = "powers.dat";
my $side;               # L|R  刺激手(対側判定のヒント。任意)
my $expect;             # "x,y,z" 期待 S1 の MNI(省略時は side から自動)
my $topn    = 12;       # 診断で出す上位頂点数
my $diag    = 1;        # 診断ブロック(--no-diag で抑止)
GetOptions(
    "data=s"     => \$data_file,
    "labels=s"   => \$labels_file,
    "times=s"    => \$times_file,
    "sfreq=f"    => \$sfreq,
    "tmin=f"     => \$tmin,
    "mat=s"      => \$mat,
    "latency=f"  => \$latency,
    "search=s"   => \$search,
    "avg-ms=f"   => \$avg_ms,
    "method=s"   => \$method,
    "reg-frac=f" => \$reg,
    "out=s"      => \$out,
    "side=s"     => \$side,
    "expect=s"   => \$expect,
    "topn=i"     => \$topn,
    "diag!"      => \$diag,
) or die "bad options\n";
$data_file   or die "need --data ch-by-time text matrix\n";
$labels_file or die "need --labels one-label-per-line\n";
$mat         or die "need --mat sa_nyhead.mat\n";

# ---------------------------------------------------------- read waveform ---
# --data: n_ch 行 × n_time 列(空白区切り)。np.savetxt(ev.data) がそのまま。
my @rows;
open my $df, '<', $data_file or die "open $data_file: $!";
while (<$df>) { next unless /\S/; push @rows, [ split ' ' ]; }
close $df;
my $data = pdl(\@rows);                 # dims = (n_time, n_ch)  ※行=ch, 列=time
my ($ntime, $nch) = ($data->dim(0), $data->dim(1));

open my $lf, '<', $labels_file or die "open $labels_file: $!";
my @labels = map { s/\s+$//r } grep { /\S/ } <$lf>;   # strip CR/末尾空白
close $lf;
@labels == $nch or die "labels($#labels+1) != data rows($nch)\n";

# 時間軸(秒)
my $times;
if ($times_file) {
    open my $tf, '<', $times_file or die "open $times_file: $!";
    $times = pdl( map { chomp; $_ } grep { /\S/ } <$tf> );
    close $tf;
    $times->nelem == $ntime or die "times != n_time\n";
} else {
    $sfreq > 0 or die "need --sfreq (or --times)\n";
    $times = $tmin + sequence($ntime) / $sfreq;
}
printf "loaded %s: %d ch x %d samp,  %.1f..%.1f ms\n",
    $data_file, $nch, $ntime,
    1000*$times->at(0), 1000*$times->at($ntime-1);

# ----------------------------------------------------- pick N20 latency ---
my $ti;   # 時間インデックス
if (defined $latency) {
    $ti = ( ($times*1000 - $latency)->abs )->minimum_ind;
} else {
    my ($lo, $hi) = split /,/, $search;
    my $ms   = $times * 1000;
    my $mask = ($ms >= $lo) & ($ms <= $hi);
    $mask->any or die "search window $search ms has no samples\n";
    # GFP = 電極方向の空間 SD(符号非依存でN20/P20双極子の最大をとる)
    my $mean = $data->mv(1,0)->average;               # (n_time)  各時刻の電極平均
    my $gfp  = (($data - $mean->dummy(1))**2)->mv(1,0)->average->sqrt;  # (n_time)
    my $gfpm = $gfp->where($mask);
    my $idxm = which($mask);
    $ti = $idxm->at( $gfpm->maximum_ind );
}
my $lat_ms = 1000 * $times->at($ti);
printf "N20 latency = %.2f ms  (sample %d)%s\n",
    $lat_ms, $ti, (defined $latency ? " [given]" : " [GFP auto]");

# トポ b(±avg_ms 平均 or 単サンプル)
my $topo;
if ($avg_ms > 0 && $sfreq > 0) {
    my $w = int($avg_ms/1000 * $sfreq + 0.5);
    my $a = $ti - $w; $a = 0        if $a < 0;
    my $b = $ti + $w; $b = $ntime-1 if $b > $ntime-1;
    $topo = $data->slice("$a:$b,:")->mv(0,0)->average;   # (n_ch) 時間平均
    printf "topo averaged over samples %d..%d (+/-%.1f ms)\n", $a, $b, $avg_ms;
} else {
    $topo = $data->slice("($ti),:")->sever;              # (n_ch)
}

# -------------------------------------------- map electrodes -> NYHead ---
my $ny  = PDL::IO::NYHead->new($mat);
my $ecl = $ny->electrode_labels;                 # arrayref, 231
my %ny_idx;                                      # lc(label) -> 0-based row
for my $i (0 .. $#$ecl) { $ny_idx{ lc $ecl->[$i] } = $i; }

# 旧10-20 → 新10-05 別名(NYHead は新命名)
my %alias = (t3=>'t7', t4=>'t8', t5=>'p7', t6=>'p8', a1=>'', a2=>'');
my (@file_i, @ny_i, @unmatched);
for my $c (0 .. $nch-1) {
    my $key = lc $labels[$c];
    $key =~ s/\s+//g;
    $key = $alias{$key} if exists $alias{$key};
    if ($key ne '' && exists $ny_idx{$key}) {
        push @file_i, $c;  push @ny_i, $ny_idx{$key};
    } else {
        push @unmatched, $labels[$c];
    }
}
my $ne = @ny_i;
$ne >= 4 or die "only $ne electrodes matched NYHead — check labels\n";
printf "matched %d/%d scalp electrodes to NYHead (dropped %d: %s)\n",
    $ne, $nch, scalar(@unmatched),
    (@unmatched > 12 ? join(',',@unmatched[0..11]).",..." : join(',',@unmatched));

my $sel = pdl(long, \@ny_i);
my $K   = $ny->leadfield->dice_axis(0, $sel)->transpose->sever;  # (Ns, Ne)  ← 転置必須
my $b   = $topo->index( pdl(long, \@file_i) )->sever;            # (Ne) 同順

my $Ns = $K->dim(0);
$Ns > $ne or die "sanity: expected Ns($Ns) >> Ne($ne)\n";

# ------------------------------------------------------ inverse solve ---
# ★ ここが Mac 側 PDL::EEG::Inverse::MinimumNorm の実シグネチャに依存する唯一の箇所。
#   前スレの nyhead_inverse.pl と同じ呼び方に合わせること(下は申し送りからの再構成)。
#   - K は (Ns,Ne) 格納。leadfield は (Ne,Ns) なので上で transpose 済。
#   - モンタージュ部分選択なので subset 上で再平均参照(avg_reference)してから解く。
$HAVE_INV or die
    "PDL::EEG::Inverse::MinimumNorm が見つかりません(Mac 側 PDL-EEG を \@INC に)。\n";

my $Kc = avg_reference($K);                 # 各源列を電極軸で中心化(subset CAR)
my $bc = $b - $b->avg;                       # トポも subset CAR
my $op  = inverse_operator($Kc, method => $method, ref => 'car', reg_frac => $reg);
my $pow = source_power($op, $bc);            # (Ns) 標準化源パワー(無次元)

# ----------------------------------------------------- diagnostics ------
my $pv   = $pow->flat;
my $vc   = $ny->cortex('cortex75K')->{vc};   # (Ns,3) MNI mm : dim0=頂点, dim1=xyz
my $gmax = $pv->max;
my $peak = $pv->maximum_ind;                 # 0-based 75K 頂点

sub _mni { my $v=shift; ($vc->at($v,0), $vc->at($v,1), $vc->at($v,2)) }

if ($diag) {
    # (1) N20 トポ b の双極子パターン(CAR済 bc を標準電極で。左手刺激なら対側=右で
    #     C4/CP4 が負(N20)・F4 側が正(P20)になるのが教科書)
    my %bval; $bval{ lc $labels[$file_i[$_]] } = $bc->at($_) for 0 .. $ne-1;
    print "\n--- N20 topo b (CAR, at latency) ----------------------------\n";
    for my $L (qw(Fz Cz Pz F3 F4 C3 C4 P3 P4 FC3 FC4 CP3 CP4 CPz)) {
        printf "  %-4s % .3g\n", $L, $bval{lc $L} if exists $bval{lc $L};
    }

    # (2) 順投影相関(決定的): b が「期待 S1 の単一源トポ」と「inv ピークの単一源
    #     トポ」のどちらに似ているか。S1 相関が高いのに inv が別を指す=inverse/正則化
    #     の問題。S1 相関が低い=b 自体が S1 らしくない(潜時/符号/電極対応/データ)。
    my ($ex,$ey,$ez) = (34,-27,51);                     # 右 S1(hand knob 近傍)
    ($ex = -34) if defined $side && uc($side) eq 'R';   # 右手刺激なら左 S1
    ($ex,$ey,$ez) = split /,/, $expect if defined $expect;
    my $tgt = pdl($ex,$ey,$ez);
    my $s1v = ( (($vc - $tgt->dummy(0))**2)->mv(1,0)->sumover )->minimum_ind;
    my $bcn = $bc - $bc->avg;
    my $corr = sub { my ($a,$b)=@_; my $an=$a-$a->avg; my $bn=$b-$b->avg;
        my $d = sqrt( ($an*$an)->sum * ($bn*$bn)->sum ); $d==0 ? 0 : ($an*$bn)->sum/$d };
    my $fwd = sub { my $v=shift; my $j=zeroes($Ns); $j->set($v,1); forward_project($Kc,$j) };
    printf "\n--- forward check: corr(b, single-source topo) --------------\n";
    printf "  expect S1  v%-6d (%.0f,%.0f,%.0f) %-30s corr=% .3f\n",
        $s1v, _mni($s1v), ($ny->area_of_vertex($s1v)//'?'), $corr->($bcn,$fwd->($s1v));
    printf "  inv peak   v%-6d (%.0f,%.0f,%.0f) %-30s corr=% .3f\n",
        $peak, _mni($peak), ($ny->area_of_vertex($peak)//'?'), $corr->($bcn,$fwd->($peak));

    # 全皮質を単一ダイポール適合でスキャン: b と各源トポ(K の行)の |corr| 最大=最良単一源。
    # これが postcentral なら「データは正常・inverse が前頭に引っ張られている」、前頭なら
    # 「b 自体が S1 らしくない(潜時/符号/電極対応/データ)」と決定的に切り分く。
    my $Kr  = $Kc - $Kc->avg(1)->dummy(1);              # 各源行を電極軸で中心化
    my $dot = ($Kr * $bcn->dummy(0))->mv(1,0)->sumover; # (Ns)
    my $nk  = sqrt( ($Kr*$Kr)->mv(1,0)->sumover );
    my $nb  = sqrt( ($bcn*$bcn)->sum );
    my $cc  = $dot / ($nk * $nb + 1e-30);
    my $bf  = $cc->abs->maximum_ind;
    printf "  best 1-dip  v%-6d (%.0f,%.0f,%.0f) %-30s |corr|=% .3f\n",
        $bf, _mni($bf), ($ny->area_of_vertex($bf)//'?'), $cc->at($bf);

    # (3) 上位頂点
    my $order = $pv->qsorti;                            # 昇順 index
    printf "\n--- top %d source vertices (pow / MNI / HO) ------------------\n", $topn;
    for my $r (0 .. $topn-1) {
        my $v = $order->at($pv->nelem-1-$r);
        printf "  #%2d v%-6d %.3g  (%.0f,%.0f,%.0f)  %s\n",
            $r+1, $v, $pv->at($v), _mni($v), ($ny->area_of_vertex($v)//'?');
    }

    # (4) postcentral / precentral の最良頂点(順位と最大比)
    my $ho = eval { $ny->ho_labels };  my $ai = eval { $ny->atlas_index };
    if (ref($ho) eq 'ARRAY' && defined $ai) {
        for my $nm ('Postcentral', 'Precentral') {
            my @lv = grep { defined $ho->[$_] && $ho->[$_] =~ /\Q$nm\E/i } 0 .. $#$ho;
            next unless @lv;
            my $mask = zeroes(long, $ai->nelem);
            $mask = $mask | ($ai == ($_ + 1)) for @lv;   # atlas は 1-based
            my $mp = $pv->where($mask);
            next unless $mp->nelem;
            my $mx = $mp->max;
            my $vi = which( ($pv >= $mx) & $mask )->at(0);
            my $rank = ($pv > $mx)->sum + 1;
            printf "  best %-12s v%-6d %.3g = %.1f%% of max, rank #%d  (%.0f,%.0f,%.0f)\n",
                $nm, $vi, $mx, 100*$mx/$gmax, $rank, _mni($vi);
        }
    }
}

# ----------------------------------------------------- write & report ---
my $xyz  = $vc->slice("($peak),:");
my $area = $ny->area_of_vertex($peak) // '(no HO label)';
my $hemi = ($xyz->at(0) < 0) ? 'Left' : 'Right';

open my $of, '>', $out or die "open $out: $!";
printf $of "%.8g\n", $pv->at($_) for 0 .. $pv->nelem-1;
close $of;

printf "\n--- result ---------------------------------------------------\n";
printf "peak vertex   : %d\n", $peak;
printf "peak MNI (mm) : x=%.1f y=%.1f z=%.1f  (%s hemisphere)\n",
    $xyz->at(0), $xyz->at(1), $xyz->at(2), $hemi;
printf "peak HO area  : %s\n", $area;
printf "stim side     : %s -> expect peak in contralateral Postcentral Gyrus\n",
    ($side // '?') if defined $side;
printf "wrote %s (%d source values, cortex75K order)\n", $out, $pv->nelem;
printf "\nnext: demo_gs3d_nyhead.pl --mat %s --surf /sa/cortex10K "
     . "--no-electrodes --source-power %s --cmap inferno\n", $mat, $out;
