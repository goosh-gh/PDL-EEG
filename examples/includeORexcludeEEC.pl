#!/usr/bin/env perl
# erp_proxy_eloreta.pl — 「代理電極ありなしで eLORETA 解がどう変わるか」の A/B。
#
# 実データに実際の外耳道電極 EEC がある。EEC の実測電位を、NYHead リードフィールドの
# 近傍既存電極(既定 Ex19 / LPA)の行に載せて(=代理)、その電極を「入れた解」と
# 「入れない解」を eLORETA で作り、定量比較する。
#
#   perl erp_proxy_eloreta.pl --mat ~/src/NYHead/sa_nyhead.mat \
#        --eeg 58tgtBN0.63_LP30Hz.renamedr2_EEC --latency 130 --proxy Ex19,LPA
#   perl erp_proxy_eloreta.pl --self-test        # 合成リードフィールドでコア検証(.mat/eeg.pm不要)
#
# 代理の意味(結果に添えること): 動くのは (a)耳領域センサ増の真の効き + (b)EEC真位置
# (外耳道内、Ex19 より ~12mm 内側)と代理位置ズレの順モデル誤差、が混在。(b) を消して
# (a) だけ取り出すのはモデル再構築(パターン3, pending)。今回の A/B は「耳領域チャネルが
# 解を動かすか」の一次判定。ERP は源が広く低 dipolarity ゆえ差は微妙 → 定量で出す。
#
# 依存: PDL, PDL::EEG::Inverse::MinimumNorm。実データ時は PDL::IO::NYHead と eeg.pm。
use strict; use warnings;
use Getopt::Long; use PDL; # use PDL::NiceSlice;
use lib '/Users/goosh/src/PDL_IO_NYHead/lib'; use PDL::IO::NYHead;
use lib '/Users/goosh/src/PDL_EEG/lib'; use PDL::EEG::Inverse::MinimumNorm qw(inverse_operator source_power forward_project avg_reference);
use eeg; use mat;

my %o = (
    mat      => "$ENV{HOME}/src/NYHead/sa_nyhead.mat",
    eeg      => undef,
    trig     => 'trgCH',      # NK_READ_ONEAVGFILE2 の第2引数
    arg3     => 0,            # 第3引数
    ear_label=> 'EEC',        # 外耳道電極のラベル
    proxy    => 'Ex19,LPA',   # 代理に使う既存電極(カンマ区切りで複数=各々別々に A/B)
    latency  => undef,        # ms。未指定なら post-stim GFP ピークを自動採用
    method   => 'eloreta',
    reg_frac => 0.05,
    out      => 'proxy',      # 出力接頭辞: <out>_without.dat / <out>_with_<proxy>.dat
    eeg_lib  => undef,        # eeg.pm のあるディレクトリ(必要なら)
    self_test=> 0,
);
GetOptions(
    'mat=s'=>\$o{mat}, 'eeg=s'=>\$o{eeg}, 'trig=s'=>\$o{trig}, 'arg3=i'=>\$o{arg3},
    'ear-label=s'=>\$o{ear_label}, 'proxy=s'=>\$o{proxy}, 'latency=f'=>\$o{latency},
    'method=s'=>\$o{method}, 'reg-frac=f'=>\$o{reg_frac}, 'out=s'=>\$o{out},
    'eeg-lib=s'=>\$o{eeg_lib}, 'self-test'=>\$o{self_test},
) or die "bad args\n";


# ============================================================================
# エンジン・コア: リードフィールド $lf=(Ne_all, Ns) と 電極 index piddle $idx、
# 電位 $b(選択電極, $idx と同順) から標準化源パワー (Ns) を返す。実データ・合成 共通。
# ============================================================================
sub run_inverse {
    my ($lf, $idx, $b, %opt) = @_;
    my $K = $lf->dice_axis(0, $idx)->transpose->sever;   # (Ns, Ne_sub)
    my ($Ns,$Ne) = $K->dims;
    die "need Ns>Ne, got ($Ns,$Ne)\n" unless $Ns > $Ne;
    $K = avg_reference($K);                                # subset 上で再 CAR
    my $op = inverse_operator($K, method=>$opt{method}//'eloreta', reg_frac=>$opt{reg_frac}//0.05);
    my $bb = $b - $b->avg;                                 # データも平均参照(centering でなく mean 減算)
    return source_power($op, $bb);                         # (Ns)
}

# 2つの源パワーマップを比較(定量)。$vc=(Ns,3) があればピークの MNI/距離も出す。
sub compare {
    my ($pw0, $pw1, $vc) = @_;
    my $p0 = $pw0->maximum_ind; my $p1 = $pw1->maximum_ind;   # ピーク頂点 index
    my $r  = _corr($pw0, $pw1);                              # 2解の相関(全頂点)
    # 正規化差(各々 max=1 に規格化してから)
    my $n0 = $pw0/($pw0->max||1); my $n1 = $pw1/($pw1->max||1);
    my $dmax = ($n1-$n0)->abs->max;
    my %R = (peak0=>$p0->sclr, peak1=>$p1->sclr, corr=>$r, dmax_norm=>$dmax->sclr);
    if (defined $vc) {
        my $m0 = $vc->dice_axis(0,$p0)->flat; my $m1 = $vc->dice_axis(0,$p1)->flat;
        $R{mni0} = [map {sprintf'%.1f',$_} $m0->list];
        $R{mni1} = [map {sprintf'%.1f',$_} $m1->list];
        $R{peak_shift_mm} = sqrt((($m1-$m0)**2)->sum)->sclr;
    }
    return \%R;
}
sub _corr {
    my ($a,$b)=@_; my $am=$a-$a->avg; my $bm=$b-$b->avg;
    my $den = sqrt(($am*$am)->sum * ($bm*$bm)->sum); return 0 if $den==0;
    (($am*$bm)->sum/$den)->sclr;
}

# ============================================================================
# 実データ経路
# ============================================================================
sub real_run {
    require PDL::IO::NYHead;
    unshift @INC, $o{eeg_lib} if $o{eeg_lib};
    require eeg;

    my $ny = PDL::IO::NYHead->new($o{mat});
    my $lf = $ny->leadfield;                    # (Ne=231, Ns=74382)
    my $elabs = $ny->electrode_labels;          # 231 ラベル
    my %eidx; for my $i (0..$#$elabs){ $eidx{uc $elabs->[$i]} = $i; }
    my %alias = (T3=>'T7',T4=>'T8',T5=>'P7',T6=>'P8',A1=>'', A2=>'');  # 記録→NYHead 別名
    my $nyfind = sub { my $l=uc shift; $l=$alias{$l} if exists $alias{$l}; return undef unless $l; $eidx{$l} };
    my $vc = eval { $ny->surface('cortex75K')->{vc} };   # (Ns,3) MNI (peak 位置報告用)

    # --- eeg.pm 読込 ---
    defined $o{eeg} or die "--eeg <NKファイル基底名> が要ります\n";
    my ($rd,$rh) = eeg::NK_READ_ONEAVGFILE2($o{eeg}, $o{trig}, $o{arg3});
    my $NCH  = $rh->{general}{ChannelN};
    my $sf   = $rh->{average}{sampling_freq};
    my $pre  = $rh->{average}{pretrigger};
    my $data = pdl([ map { $rd->[$_] } 0..$NCH-1 ]);      # (Nt, NCH)? -> 実は (行=時間? ) 要注意
    # $rd->[電極][時間] なので pdl(\@rows) は (Nt, NCH): dim0=時間, dim1=電極。
    # 各電極行 = $data(:,($i)) = 時間系列。潜時 sample での値を取る。
    printf STDERR "eeg: %d ch, sfreq=%g, pretrigger=%d, samples=%d\n", $NCH,$sf,$pre,$data->dim(0);

    # --- チャンネル→NYHead 対応(standard名でNYHeadに在るものだけ採用、EECは代理) ---
    my (@std_idx, @std_row, @std_lab, $eec_row);
    my @skipped;
    for my $i (0..$NCH-1) {
        my $lab = $rh->{ChannelName}{$i};
        next unless defined $lab && $lab =~ /\S/;
        if (uc $lab eq uc $o{ear_label}) { $eec_row = $i; next; }
        my $ni = $nyfind->($lab);
        if (defined $ni) { push @std_idx,$ni; push @std_row,$i; push @std_lab,$lab; }
        else { push @skipped, $lab; }
    }
    printf STDERR "matched %d scalp ch to NYHead; skipped(非対応): %s\n",
        scalar(@std_idx), (@skipped?join(",",@skipped):'(なし)');
    die "EEC ラベル '$o{ear_label}' が見つからない\n" unless defined $eec_row;
    printf STDERR "ear electrode '%s' = data row %d\n", $o{ear_label}, $eec_row;

    # --- 潜時 sample ---
    my $lat = $o{latency};
    if (!defined $lat) {                          # GFP ピーク(post-stim)を自動
        my $sel = pdl(@std_row);
        my $D   = $data->dice_axis(1,$sel);       # (Nt, Nsel)
        my $gfp = (($D - $D->xchg(0,1)->avgover->dummy(1))**2)->xchg(0,1)->avgover->sqrt; # (Nt)
        $gfp->slice("0:$pre") .= 0;                         # 刺激前は無視
        $lat = (($gfp->maximum_ind->sclr) - $pre)*1000.0/$sf;
        printf STDERR "auto latency (GFP peak) = %.1f ms\n", $lat;
    }
    my $samp = $pre + sprintf("%.0f", $lat*$sf/1000.0);
    printf STDERR "latency %.1f ms -> sample %d\n", $lat, $samp;

    my $b_std = $data->slice("($samp),:")->dice_axis(0, pdl(@std_row))->sever;  # 標準ch電位
    my $b_eec = $data->at($samp, $eec_row);                                     # EEC電位(スカラ)

    # --- WITHOUT ---
    my $idx0 = pdl(@std_idx);
    my $pw0  = run_inverse($lf, $idx0, $b_std, method=>$o{method}, reg_frac=>$o{reg_frac});
    _write("$o{out}_without.dat", $pw0);

    # --- WITH: 代理電極ごと ---
    for my $pl (split /,/, $o{proxy}) {
        my $pi = $nyfind->($pl) // $eidx{uc $pl};
        unless (defined $pi) { warn "proxy '$pl' が NYHead に無い(skip)\n"; next; }
        my $idx1 = pdl(@std_idx, $pi);
        my $b1   = $b_std->append(pdl($b_eec));
        my $pw1  = run_inverse($lf, $idx1, $b1, method=>$o{method}, reg_frac=>$o{reg_frac});
        _write("$o{out}_with_$pl.dat", $pw1);
        my $R = compare($pw0, $pw1, $vc);
        _report($pl, $R, $ny);
    }
}

sub _report {
    my ($pl,$R,$ny)=@_;
    my $area0 = eval { $ny->area_of_vertex($R->{peak0}) } // '?';
    my $area1 = eval { $ny->area_of_vertex($R->{peak1}) } // '?';
    print  "\n==== 代理 = $pl ====\n";
    printf "without: peak v%d  %s  [%s]\n", $R->{peak0}, $area0, join(",",@{$R->{mni0}//[]});
    printf "with   : peak v%d  %s  [%s]\n", $R->{peak1}, $area1, join(",",@{$R->{mni1}//[]});
    printf "peak shift = %.1f mm\n", $R->{peak_shift_mm} if defined $R->{peak_shift_mm};
    printf "corr(without,with) = %.4f    max|Δ|(norm) = %.3f\n", $R->{corr}, $R->{dmax_norm};
    print  "  → shift 小&corr≈1: 耳ch は解をほぼ動かさない / shift 大 or corr 低下: 動かす\n";
}
sub _write { my ($f,$pw)=@_; open my $fh,'>',$f or die "$f: $!"; print $fh join("\n",$pw->list),"\n"; close $fh;
    printf STDERR "wrote %s (%d rows)\n",$f,$pw->nelem; }

# ============================================================================
# 合成セルフテスト: 実 MinimumNorm で run_inverse/compare を検証(.mat/eeg.pm 不要)
# ============================================================================
sub self_test {
    my ($Ns,$Ne) = (200, 20);
    my $lf = grandom($Ne,$Ns);                    # (Ne, Ns) 合成リードフィールド
    # 源 s0 を1点立て、全電極へ順投影して合成電位を作る
    my $s0 = 111;
    my $j  = zeroes($Ns); $j->slice("($s0)") .= 1.0;
    my $Kall = $lf->transpose->sever;             # (Ns,Ne)
    my $b_all = forward_project(avg_reference($Kall), $j);  # (Ne) 平均参照系
    # WITHOUT = 電極 0..Ne-2、WITH = さらに最後の電極(=代理相当)を足す
    my $idx0 = sequence($Ne-1);                   # 0..Ne-2
    my $idx1 = sequence($Ne);                     # 0..Ne-1
    my $b0 = $b_all->dice_axis(0,$idx0)->sever;
    my $b1 = $b_all->sever;
    my $pw0 = run_inverse($lf,$idx0,$b0);
    my $pw1 = run_inverse($lf,$idx1,$b1);
    my $R = compare($pw0,$pw1);
    printf "self-test: Ns=%d Ne=%d, planted source s0=%d\n",$Ns,$Ne,$s0;
    printf "  WITHOUT peak v%d (=s0? %s)\n", $R->{peak0}, ($R->{peak0}==$s0?'YES':'no');
    printf "  WITH    peak v%d (=s0? %s)\n", $R->{peak1}, ($R->{peak1}==$s0?'YES':'no');
    printf "  corr=%.4f  max|Δ|norm=%.3f\n", $R->{corr}, $R->{dmax_norm};
    print  $R->{peak0}==$s0 && $R->{peak1}==$s0 ? "CORE OK: 両解が s0 に局在、比較指標も算出\n"
                                                : "NOTE: 合成ランダム LF では厳密局在しないことがある(コア経路が通れば可)\n";
}

$o{self_test} ? self_test() : real_run();



