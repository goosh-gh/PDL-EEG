#!/usr/bin/env perl
# ---------------------------------------------------------------------------
# sep_n20_sweep.pl  (P:EEG22)
#
#   潜時スイープ: 加算 SEP を読み、lat-min..lat-max を lat-step おきに走査して
#   各潜時のトポ b から源パワーを求め、<outdir>/powers.<lat>.dat に書き出す。
#   逆作用素は K だけで決まり b に依存しないので **1回だけ**構築し、潜時ごとには
#   source_power(matrix×vector) を回すだけ=451 潜時でも高速。
#
#   出力:
#     <outdir>/powers.0NN.N.dat        各潜時の源パワー(74382 行, cortex75K 順)
#     <outdir>/peak_by_latency.tsv     潜時→ピーク頂点/MNI/HO野/Postcentral比 の時系列
#
#   例:
#     perl -I<PDL_EEG/lib> -I<PDL_IO_NYHead/lib> sep_n20_sweep.pl \
#       --mat sa_nyhead.mat --data sep_ave.txt --labels sep_labels.txt \
#       --sfreq 10000 --tmin -0.05 --method eloreta \
#       --lat-min 15 --lat-max 60 --lat-step 0.1 --outdir ./LORETA_output
#
#   ある潜時を見る:
#     demo_gs3d_nyhead.pl --mat sa_nyhead.mat --surf /sa/cortex10K \
#       --no-electrodes --source-power ./LORETA_output/powers.020.1.dat \
#       --cmap inferno --render zbuffer
# ---------------------------------------------------------------------------
use strict;
use warnings;
use PDL;
use PDL::IO::NYHead;
use File::Path qw(make_path);
use Getopt::Long;

my @INV_FN = qw(avg_reference inverse_operator source_power);
my $HAVE_INV = eval { require PDL::EEG::Inverse::MinimumNorm;
                      PDL::EEG::Inverse::MinimumNorm->import(@INV_FN); 1 };

my ($data_file, $labels_file, $times_file, $mat);
my $sfreq   = 0;
my $tmin    = 0;
my $avg_ms  = 0;
my $method  = "eloreta";
my $reg     = 0.05;
my $side;
my $lat_min = 15;
my $lat_max = 60;
my $lat_step= 0.1;
my $outdir  = "./LORETA_output";
my $summary_only = 0;     # 個別 powers.<lat>.dat を書かず peak_by_latency.tsv だけ
my $dump_list;            # "20,20.1,21" 等: この潜時だけ powers.<lat>.dat を書く
GetOptions(
    "data=s"     => \$data_file,
    "labels=s"   => \$labels_file,
    "times=s"    => \$times_file,
    "sfreq=f"    => \$sfreq,
    "tmin=f"     => \$tmin,
    "mat=s"      => \$mat,
    "avg-ms=f"   => \$avg_ms,
    "method=s"   => \$method,
    "reg-frac=f" => \$reg,
    "side=s"     => \$side,
    "lat-min=f"  => \$lat_min,
    "lat-max=f"  => \$lat_max,
    "lat-step=f" => \$lat_step,
    "outdir=s"   => \$outdir,
    "summary-only" => \$summary_only,
    "dump=s"       => \$dump_list,
) or die "bad options\n";
$data_file   or die "need --data ch-by-time text matrix\n";
$labels_file or die "need --labels one-label-per-line\n";
$mat         or die "need --mat sa_nyhead.mat\n";
$HAVE_INV    or die "PDL::EEG::Inverse::MinimumNorm が見つかりません(\@INC に PDL-EEG)。\n";
$lat_step > 0 or die "--lat-step must be > 0\n";

# ---------------------------------------------------------- read waveform ---
my @rows;
open my $df, '<', $data_file or die "open $data_file: $!";
while (<$df>) { next unless /\S/; push @rows, [ split ' ' ]; }
close $df;
my $data = pdl(\@rows);                 # (n_time, n_ch)
my ($ntime, $nch) = ($data->dim(0), $data->dim(1));

open my $lf, '<', $labels_file or die "open $labels_file: $!";
my @labels = map { s/\s+$//r } grep { /\S/ } <$lf>;
close $lf;
@labels == $nch or die "labels(${\ scalar @labels}) != data rows($nch)\n";

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
my $ms = $times * 1000;
printf "loaded %s: %d ch x %d samp,  %.1f..%.1f ms\n",
    $data_file, $nch, $ntime, $ms->at(0), $ms->at($ntime-1);

# -------------------------------------------- map electrodes -> NYHead ---
my $ny  = PDL::IO::NYHead->new($mat);
my $ecl = $ny->electrode_labels;                 # 231
my %ny_idx; $ny_idx{ lc $ecl->[$_] } = $_ for 0 .. $#$ecl;
my %alias = (t3=>'t7', t4=>'t8', t5=>'p7', t6=>'p8', a1=>'', a2=>'');
my (@file_i, @ny_i, @unmatched);
for my $c (0 .. $nch-1) {
    my $key = lc $labels[$c]; $key =~ s/\s+//g;
    $key = $alias{$key} if exists $alias{$key};
    if ($key ne '' && exists $ny_idx{$key}) { push @file_i, $c; push @ny_i, $ny_idx{$key}; }
    else { push @unmatched, $labels[$c]; }
}
my $ne = @ny_i;
$ne >= 4 or die "only $ne electrodes matched NYHead — check labels\n";
printf "matched %d/%d scalp electrodes (dropped %d: %s)\n",
    $ne, $nch, scalar(@unmatched),
    (@unmatched > 12 ? join(',',@unmatched[0..11]).",..." : join(',',@unmatched));

my $sel = pdl(long, \@ny_i);
my $fidx= pdl(long, \@file_i);
my $K   = $ny->leadfield->dice_axis(0, $sel)->transpose->sever;  # (Ns, Ne)
my $Ns  = $K->dim(0);

# ------------------------------------------------ build operator ONCE ---
my $Kc = avg_reference($K);                                      # subset CAR
my $op = inverse_operator($Kc, method => $method, ref => 'car', reg_frac => $reg);
printf "inverse operator: method=%s, Ns=%d, Ne=%d (built once, reused per latency)\n",
    $method, $Ns, $ne;

# Postcentral マスク(S1 比の指標。無ければ比は N/A)
my $ho = eval { $ny->ho_labels };  my $ai = eval { $ny->atlas_index };
my $pc_mask;
if (ref($ho) eq 'ARRAY' && defined $ai) {
    my @lv = grep { defined $ho->[$_] && $ho->[$_] =~ /Postcentral/i } 0 .. $#$ho;
    if (@lv) { $pc_mask = zeroes(long, $ai->nelem); $pc_mask = $pc_mask | ($ai == ($_+1)) for @lv; }
}
my $vc = $ny->cortex('cortex75K')->{vc};                        # (Ns,3)

# ------------------------------------------------------ sweep latencies ---
make_path($outdir) unless -d $outdir;
my $sumf = "$outdir/peak_by_latency.tsv";
open my $sf, '>', $sumf or die "open $sumf: $!";
print $sf "latency_ms\tpeak_v\tx\ty\tz\tHO_area\tpostcentral_frac_pct\n";

my $wsamp = ($avg_ms > 0 && $sfreq > 0) ? int($avg_ms/1000*$sfreq + 0.5) : 0;
my $lo10 = int($lat_min*10 + 0.5);
my $hi10 = int($lat_max*10 + 0.5);
my $st10 = int($lat_step*10 + 0.5); $st10 = 1 if $st10 < 1;
my %dump = map { sprintf("%.1f",$_) => 1 } grep { length } split /,/, ($dump_list // '');

# 潜時リスト + 全潜時のトポを (Ne, Nlat) に一括構築
my @lats; for (my $l = $lo10; $l <= $hi10; $l += $st10) { push @lats, $l/10; }
my $Nlat = @lats;
my $B = zeroes($ne, $Nlat);
for my $k (0 .. $Nlat-1) {
    my $lat = $lats[$k];
    my $ti  = int(($lat/1000 - $tmin) * $sfreq + 0.5);
    $ti = 0 if $ti < 0;  $ti = $ntime-1 if $ti > $ntime-1;
    my $topo;
    if ($wsamp > 0) {
        my $a = $ti-$wsamp; $a = 0 if $a < 0;
        my $bb= $ti+$wsamp; $bb= $ntime-1 if $bb > $ntime-1;
        $topo = $data->slice("$a:$bb,:")->mv(0,0)->average;
    } else {
        $topo = $data->slice("($ti),:");
    }
    $B->slice(":,($k)") .= $topo->index($fidx);
}
my $Bc = $B - $B->average->dummy(0);        # 各列(潜時)を電極軸で CAR 中心化

# 一括解: 作用素は b 非依存なので全潜時を 1 回の行列積で解く。source_power が (Ne,Nlat)
# を受ければ (Ns,Nlat) が一撃で返る(=BLAS gemm)。受けなければ潜時ループにフォールバック。
my $POW;
my $batched = eval {
    my $p = source_power($op, $Bc);
    ($p->ndims == 2 && $p->dim(0) == $Ns && $p->dim(1) == $Nlat) ? do { $POW = $p; 1 } : 0;
};
if (!$batched) {
    warn "source_power is not vectorized over columns -> per-latency loop (operator still reused)\n";
    $POW = zeroes($Ns, $Nlat);
    for my $k (0 .. $Nlat-1) {
        $POW->slice(":,($k)") .= source_power($op, $Bc->slice(":,($k)")->sever)->flat;
    }
}
printf "computed %d latencies %s (Ns=%d, Ne=%d)\n",
    $Nlat, ($batched ? "in one batched solve" : "per-latency"), $Ns, $ne;

# サマリ量は潜時方向の縮約でベクトル化
my $peakv = $POW->maximum_ind;                       # (Nlat) 各潜時のピーク頂点
my $tot   = $POW->sumover;                           # (Nlat)
my $pcsum = defined $pc_mask ? ($POW * $pc_mask->dummy(1))->sumover : undef;   # (Nlat)

my $nfile = 0;
for my $k (0 .. $Nlat-1) {
    my $lat = $lats[$k];
    my $v   = $peakv->at($k);
    my ($px,$py,$pz) = ($vc->at($v,0), $vc->at($v,1), $vc->at($v,2));
    my $area = $ny->area_of_vertex($v) // '?';
    my $frac = defined $pcsum
        ? sprintf("%.1f", $tot->at($k) > 0 ? 100*$pcsum->at($k)/$tot->at($k) : 0) : 'NA';
    printf $sf "%.1f\t%d\t%.1f\t%.1f\t%.1f\t%s\t%s\n", $lat,$v,$px,$py,$pz,$area,$frac;

    next if $summary_only;                                    # 個別 .dat を書かない
    next if %dump && !$dump{ sprintf("%.1f",$lat) };          # --dump: 指定潜時だけ
    $POW->slice(":,($k)")->sever->wcols(sprintf("%s/powers.%05.1f.dat", $outdir, $lat));
    $nfile++;
}
close $sf;

printf "\nwrote %d powers file(s) to %s/%s\n", $nfile, $outdir,
    ($summary_only ? "  (--summary-only: individual maps skipped)" : "");
printf "summary: %s  (latency -> peak vertex / MNI / HO / Postcentral%%)\n", $sumf;
