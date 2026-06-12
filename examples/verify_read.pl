#!/usr/bin/env perl
# examples/verify_read.pl
#
# 実データでの動作確認スクリプト。make test とは独立して使う。
#
# Usage:
#   perl -Ilib examples/verify_read.pl patient.eeg
#   perl -Ilib examples/verify_read.pl patient.eeg --plot   # Cairo でプロット
#   perl -Ilib examples/verify_read.pl                      # 合成データで確認

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PDL::EEG::IO::NihonKohden qw(read_nk);

my $file  = shift;
my $plot  = grep { $_ eq '--plot' } @ARGV;
my ($fs)  = map { /^--fs=?(\d+)$/ ? $1 : () } @ARGV;
$fs     //= do { my ($i) = grep { ($ARGV[$_]//'') eq '--fs' } 0..$#ARGV;
                 defined $i ? $ARGV[$i+1] : undef };

# 引数なし → 合成データを使う
unless ($file && -f $file) {
    my $synthetic = "$Bin/../t/data/test01.eeg";
    if (-f $synthetic) {
        warn "No file given; using synthetic test data: $synthetic\n";
        $file = $synthetic;
        $fs //= 1000;
    } else {
        die "Usage: $0 patient.eeg [--fs 1000] [--plot]\n"
          . "  (Or run 'perl t/mk_synthetic_nk.pl' first to generate synthetic data)\n";
    }
}

# ---------------------------------------------------------------
# Step 1: ファイル読み込み
# ---------------------------------------------------------------
print "=== Step 1: read_nk ===\n";
my $rec = eval { read_nk($file, $fs ? (fs => $fs) : ()) };
if ($@) {
    # PDL missing → try metadata-only mode
    if ($@ =~ /PDL is required/) {
        die "FAILED: PDL not installed. Install PDL first:\n"
          . "  cpanm PDL\n  (or: port install p5-pdl)\n";
    }
    die "FAILED: $@\n";
}
print "OK\n\n";

# ---------------------------------------------------------------
# Step 2: メタデータ表示
# ---------------------------------------------------------------
print "=== Step 2: Metadata ===\n";
printf "  Device        : %s\n",   $rec->{device} // 'unknown';
printf "  Sampling rate : %g Hz\n",   $rec->{fs};
printf "  Valid channels: %d\n",      $rec->{n_ch_valid} // scalar(@{$rec->{labels}})-1;
printf "  Total ch (data): %d\n",     scalar @{$rec->{labels}};
printf "  Wfm blocks    : %d\n",      $rec->{n_blocks};
printf "  Start time    : %s\n",      $rec->{t_start} // 'unknown';
printf "  Labels        : %s\n",      join(', ', @{$rec->{labels}});
if (@{$rec->{events}}) {
    printf "  Events (%d)   :\n", scalar @{$rec->{events}};
    printf "    t=%gs  %s\n", $_->{t}, $_->{label} for @{$rec->{events}};
} else {
    print  "  Events        : none\n";
}
print "\n";

# ---------------------------------------------------------------
# Step 3: PDL データ確認
# ---------------------------------------------------------------
print "=== Step 3: PDL data ===\n";
eval { require PDL } or do { print "PDL not available, skipping\n\n"; goto PLOT };

my $data = $rec->{data};
printf "  Shape         : [%s]  (n_ch × n_samples)\n", join(', ', $data->dims);
printf "  Type          : %s\n",   $data->type;
printf "  Min / Max     : %.2f / %.2f µV\n", $data->min, $data->max;

# チャネルごとの RMS
for my $i (0 .. $data->dim(0) - 1) {
    my $ch   = $data->slice("($i),:");
    my $rms  = sqrt(($ch**2)->avg);
    my $pk   = $ch->abs->max;
    printf "  ch%02d %-12s  RMS=%6.1f µV  peak=%6.1f µV\n",
        $i, $rec->{labels}[$i], $rms, $pk;
}
print "\n";

# ---------------------------------------------------------------
# Step 4: プロット (--plot 時のみ)
# ---------------------------------------------------------------
PLOT:
if ($plot) {
    print "=== Step 4: Plot (PDL::Graphics::Cairo) ===\n";
    eval {
        require PDL::Graphics::Cairo;
        PDL::Graphics::Cairo->import(qw(figure subplots));
    } or do { print "PDL::Graphics::Cairo not available\n"; exit };

    my $n_ch   = $data->dim(0);
    my $n_plot = $n_ch < 8 ? $n_ch : 8;
    my $n_sec  = 10;
    my $n_samp = $n_sec * $rec->{fs};
    $n_samp    = $data->dim(1) if $n_samp > $data->dim(1);

    my ($fig, @ax) = subplots($n_plot, 1, figsize => [14, $n_plot * 1.4]);
    my $t = PDL::sequence($n_samp) / $rec->{fs};

    for my $i (0 .. $n_plot - 1) {
        my $ch = $data->slice("($i),0:@{[$n_samp-1]}");
        $ax[$i]->line($t, $ch, color => '#2196F3', lw => 0.7);
        $ax[$i]->ylabel($rec->{labels}[$i], fontsize => 7);
        $ax[$i]->axis('off') if $i < $n_plot - 1;
    }
    $ax[-1]->xlabel('Time (s)');

    # イベントマーカー
    for my $evt (@{$rec->{events}}) {
        next if $evt->{t} > $n_sec;
        $_->axvline($evt->{t}, color => 'red', lw => 0.8, alpha => 0.5) for @ax;
    }

    $fig->suptitle("PDL::EEG — $file", fontsize => 10);
    $fig->tight_layout;
    $fig->show;
    print "Plot displayed.\n";
}
