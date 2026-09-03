#!/usr/bin/env perl
# ---------------------------------------------------------------------------
# sep_gof_sweep.pl
#
# 単一ダイポール適合度(best-1-dip goodness of fit)を潜時掃引して 2 列ファイル
# "latency_ms<TAB>gof" を書く。sep_latency_movie.pl の --gof にそのまま渡せる。
#
# GoF の定義は sep_n20_inverse.pl の --diag と同じ = 各潜時のトポ b と、NYHead
# リードフィールドの全皮質源トポ K_v との |相関| の最大(=最良単一固定向きダイ
# ポール)。逆問題は解かない（$Kc と b だけで出る）。
#   gof = max_v | corr( b , K_v ) |            (--metric corr, 既定)
#   gof = max_v   corr( b , K_v )^2            (--metric r2, 分散説明率)
#
# 電極対応・CAR は sep_n20_sweep.pl と同一（同じ結果になるように)。
#
# 使い方:
#   perl sep_gof_sweep.pl --mat sa_nyhead.mat \
#     --data sep_ave.txt --labels sep_labels.txt --sfreq 10000 --tmin -0.05 \
#     --lat-min 10 --lat-max 32 --lat-step 0.2 --out gof.tsv
#
# それから:
#   perl sep_latency_movie.pl --ave sep_ave.txt --names sep_labels.txt \
#     --sfreq 10000 --tmin -0.05 --montage standard_1020_eog_nose.elc \
#     --gof gof.tsv --anim-min-ms 10 --anim-max-ms 32 --step-ms 0.2 \
#     --out sep_latency.mp4
# ---------------------------------------------------------------------------
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::NYHead;
use PDL::EEG::Inverse::MinimumNorm qw(avg_reference);
use Getopt::Long;

my ($mat, $data_file, $labels_file, $times_file);
my $sfreq   = 10000;
my $tmin    = -0.05;
my $avg_ms  = 0;          # >0 なら ±avg_ms を平均してトポを作る（既定=点）
my $lat_min = 10;
my $lat_max = 32;
my $lat_step= 0.2;
my $metric  = 'corr';     # corr | r2
my $out     = 'gof.tsv';
GetOptions(
    'mat=s'     => \$mat,
    'data=s'    => \$data_file,
    'labels=s'  => \$labels_file,
    'times=s'   => \$times_file,
    'sfreq=f'   => \$sfreq,
    'tmin=f'    => \$tmin,
    'avg-ms=f'  => \$avg_ms,
    'lat-min=f' => \$lat_min,
    'lat-max=f' => \$lat_max,
    'lat-step=f'=> \$lat_step,
    'metric=s'  => \$metric,
    'out=s'     => \$out,
) or die "bad options\n";
$mat && $data_file && $labels_file or die "--mat, --data, --labels required\n";

# -------------------------------------------------- load evoked (ch x time) ---
my @rows;
open my $df, '<', $data_file or die "open $data_file: $!";
while (<$df>) { next unless /\S/; push @rows, [split] }
close $df;
my $data = pdl(\@rows)->sever;           # (ntime, nch)  ← sep_n20_sweep.pl と同じ向き
my ($ntime, $nch) = $data->dims;

my @labels;
open my $lf, '<', $labels_file or die "open $labels_file: $!";
while (<$lf>) { s/\s+$//; next unless /\S/; push @labels, $_ }
close $lf;
@labels == $nch or die "labels (".scalar(@labels).") != channels ($nch)\n";

# 潜時[ms]（--times があれば読む。無ければ tmin/sfreq から）
my $ms;
if ($times_file) {
    my @t; open my $tf,'<',$times_file or die $!;
    while (<$tf>){ next unless /\S/; push @t,(split)[0] } close $tf;
    $ms = pdl(\@t) * 1000;
} else {
    $ms = ($tmin + sequence($ntime)/$sfreq) * 1000;
}
printf "loaded %s: %d ch x %d samp, %.1f..%.1f ms\n",
    $data_file, $nch, $ntime, $ms->at(0), $ms->at($ntime-1);

# -------------------------------------------- map electrodes -> NYHead ---
# （sep_n20_sweep.pl と同一のマッピング）
my $ny  = PDL::IO::NYHead->new($mat);
my $ecl = $ny->electrode_labels;                 # 231
my %ny_idx; $ny_idx{ lc $ecl->[$_] } = $_ for 0 .. $#$ecl;
my %alias = (t3=>'t7', t4=>'t8', t5=>'p7', t6=>'p8', a1=>'', a2=>'');
my (@file_i, @ny_i, @unmatched);
for my $c (0 .. $nch-1) {
    my $key = lc $labels[$c]; $key =~ s/\s+//g;
    $key = $alias{$key} if exists $alias{$key};
    if ($key ne '' && exists $ny_idx{$key}) { push @file_i,$c; push @ny_i,$ny_idx{$key}; }
    else { push @unmatched, $labels[$c]; }
}
my $ne = @ny_i;
$ne >= 4 or die "only $ne electrodes matched NYHead — check labels\n";
printf "matched %d/%d scalp electrodes (dropped %d: %s)\n",
    $ne, $nch, scalar(@unmatched),
    (@unmatched > 12 ? join(',',@unmatched[0..11]).",..." : join(',',@unmatched));

my $sel  = pdl(long, \@ny_i);
my $fidx = pdl(long, \@file_i);
my $K    = $ny->leadfield->dice_axis(0, $sel)->transpose->sever;  # (Ns, Ne)
my $Ns   = $K->dim(0);
my $Kc   = avg_reference($K);                                     # subset CAR (源ごと電極平均を除く)

# best-1-dip 用に源ノルムを 1 回だけ用意（Kc は既に電極中心化済み＝corr の分母）
my $Kn   = sqrt( ($Kc*$Kc)->mv(1,0)->sumover );                   # (Ns)

# ------------------------------------------------------ sweep latencies ---
my $wsamp = ($avg_ms > 0 && $sfreq > 0) ? int($avg_ms/1000*$sfreq + 0.5) : 0;
my $lo10  = int($lat_min*10 + 0.5);
my $hi10  = int($lat_max*10 + 0.5);
my $st10  = int($lat_step*10 + 0.5); $st10 = 1 if $st10 < 1;
my @lats; for (my $l=$lo10; $l<=$hi10; $l+=$st10) { push @lats, $l/10 }

open my $of, '>', $out or die "open $out: $!";
print $of "#latency_ms\tgof\n";
for my $lat (@lats) {
    my $ti = int(($lat/1000 - $tmin)*$sfreq + 0.5);
    $ti = 0 if $ti < 0; $ti = $ntime-1 if $ti > $ntime-1;
    my $topo;
    if ($wsamp > 0) {
        my $a=$ti-$wsamp; $a=0 if $a<0;
        my $b=$ti+$wsamp; $b=$ntime-1 if $b>$ntime-1;
        $topo = $data->slice("$a:$b,:")->mv(0,0)->average;   # 窓平均 (nch)
    } else {
        $topo = $data->slice("($ti),:");                     # 点 (nch)
    }
    my $b  = $topo->index($fidx)->sever;      # 対応電極だけ (Ne)
    my $bc = $b - $b->avg;                     # CAR（電極平均を除く）
    my $bn = sqrt( ($bc*$bc)->sum );
    my $dot= ($Kc * $bc->dummy(0))->mv(1,0)->sumover;   # (Ns)  dot[v]=Kc_v·bc
    my $cc = $dot / ($Kn*$bn + 1e-30);                  # (Ns)  = corr(b,K_v)
    my $bf = $cc->abs->maximum_ind;
    my $gof= $metric eq 'r2' ? $cc->at($bf)**2 : $cc->abs->at($bf);
    printf $of "%.1f\t%.4f\n", $lat, $gof;
}
close $of;
printf "wrote %d latencies -> %s  (metric=%s)\n", scalar(@lats), $out, $metric;
