#!/usr/bin/env perl
# Convert a Nihon Kohden .EEG recording to BESA ASCII multiplexed (.mul).
#
#   perl -Ilib examples/nk_to_mul.pl subject.EEG
#   perl -Ilib examples/nk_to_mul.pl subject.EEG --out out.mul --suffix -BN
#   perl -Ilib examples/nk_to_mul.pl subject.EEG --block 2
#   perl -Ilib examples/nk_to_mul.pl subject.EEG --allblocks --trig-width 4
#
#   --cut SPEC         write ONE .mul per range in DATA-coordinate seconds:
#                      "a-b[:name],c-d[:name]"  e.g. "21-376:b0b1_21_376".
#   --cut-clock SPEC   same, but ranges are WALL-CLOCK "HH:MM:SS" (as shown in
#                      the vendor viewer), e.g. "14:06:14-14:07:15:task2".
#                      A single range is fine; repeat the command for more.
#                      Wall-clock mapping assumes a continuous recording; if
#                      --allblocks spans a discontinuity, use data-coordinate
#                      --cut instead.
#   --bne[=PROP]       re-reference to the balanced non-cephalic (BNE) reference
#                      before writing: y = x - (PROP*BN1 + (1-PROP)*BN2). PROP
#                      defaults to 0.5 (the BN balance is set in analog hardware,
#                      so 0.5 = no extra digital re-balance; BN1=V, BN2=S). The
#                      acquisition reference (SystemReference, e.g. C3,C4) cancels
#                      because the weights sum to 1. Re-referenced channels are
#                      tagged -BN (or --suffix); DC/Trigger pass through; BN1/BN2
#                      are dropped. WITHOUT --bne the data are written as recorded
#                      (no re-reference) -- this is the default.
#
# Output name in cut mode: <base>_<name>.mul (<base>_cutNN.mul if no :name).

use strict;
use warnings;
use lib 'lib';
use PDL;
use PDL::EEG::IO::NihonKohden qw(read_nk);
use PDL::EEG::IO::BESA::ASCII qw(write_mul);
use PDL::EEG::Derivation      qw(bne);
use Getopt::Long;

my ($out, $block, $allblocks, $suffix, $trigw, $decimals, $cut, $cutclock, $bne);
GetOptions(
    'out=s'                => \$out,
    'block=i'              => \$block,
    'allblocks'            => \$allblocks,
    'suffix=s'             => \$suffix,       # e.g. -BN ; scalp channels only
    'trig-width=i'         => \$trigw,        # 4 reproduces the compact-trigger look
    'decimals=i'           => \$decimals,
    'cut=s'                => \$cut,
    'cutclock|cut-clock=s' => \$cutclock,
    'bne:s'                => \$bne,          # optional value: --bne or --bne=0.6
) or die "bad options\n";

my $in = shift @ARGV
    or die "usage: $0 [--out FILE] [--block N|--allblocks] [--suffix -BN] "
         . "[--trig-width 4] [--decimals 2] [--bne[=0.5]]\n"
         . "          [--cut a-b[:name]] [--cut-clock HH:MM:SS-HH:MM:SS[:name]] file.EEG\n";

# --- local, reader-independent range helpers (shared design with edf_to_mul) -
sub select_range {                 # end-exclusive: keep samples [lo, hi)
    my ($rec, $lo, $hi) = @_;
    my $n = $rec->{data}->dim(1);
    $lo = 0  if $lo < 0;
    $hi = $n if $hi > $n;
    my $data = $rec->{data}->slice(":," . $lo . ":" . ($hi - 1))->sever;
    return { %$rec, data => $data };
}
sub clock_to_samp { my ($rec, $sec) = @_; int($sec * $rec->{fs} + 0.5) }

# apply montage suffix + trig/decimal opts, then write one .mul
sub write_out {
    my ($r, $file) = @_;
    my @labels = @{ $r->{labels} };
    if (defined $suffix) {
        # montage suffix on scalp/aux channels; skip trigger/DC and the
        # reference/common channels (X1 IS suffixed, per the vendor .mul).
        @labels = map { /^(?:DC\d+|Trigger|A[12]_ref|COM|E)$/ ? $_ : "$_$suffix" }
                  @labels;
    }
    my %wopt = (labels => \@labels);
    $wopt{trig_width} = $trigw    if defined $trigw;
    $wopt{decimals}   = $decimals if defined $decimals;
    $wopt{segment}    = sprintf('BNE_prop%g', $r->{bne_prop})
        if defined $r->{bne_prop};    # provenance -> SegmentName (BESA-standard field)
    write_mul($r, $file, %wopt);
    printf "  wrote %s : %d ch x %d samp\n",
        $file, $r->{data}->dim(0), $r->{data}->dim(1);
}

my %ropt;
$ropt{block}      = $block if defined $block;
$ropt{all_blocks} = 1      if $allblocks;
my $rec = read_nk($in, %ropt);

# ----- optional BNE re-reference (default: OFF -> written as recorded) --------
if (defined $bne) {
    my $prop = length($bne) ? $bne + 0 : 0.5;
    my $suf  = defined $suffix ? $suffix : '-BN';
    $rec = eval { bne($rec, prop => $prop, suffix => $suf) };
    die "$0: --bne failed (need BN1/BN2 electrodes): $@" unless $rec;
    printf "BNE re-reference: y = x - (%.3g*BN1 + %.3g*BN2)"
         . "  [SystemReference=%s cancels]\n",
        $prop, 1 - $prop, ($rec->{system_reference} // 'unknown');
    $suffix = undef;   # bne() already applied the montage suffix; don't double it
}

# ----- cut mode: arbitrary ranges (data-coordinate or wall-clock) ------------
if (defined $cut || defined $cutclock) {
    my $fs = $rec->{fs};
    my $n  = $rec->{data}->dim(1);
    printf "read %s : %d ch @ %g Hz, %.1f s data\n",
        $in, $rec->{data}->dim(0), $fs, $n / $fs;

    my $is_clock = defined $cutclock;
    my $spec     = $is_clock ? $cutclock : $cut;

    my $start_clk = 0;
    if ($is_clock && ($rec->{t_start} // '') =~ /(\d{2}):(\d{2}):(\d{2})/) {
        $start_clk = $1 * 3600 + $2 * 60 + $3;
    }
    my $to_samp = sub {
        my $tok = shift;
        if ($tok =~ /^(\d+):(\d+):(\d+)$/) {            # HH:MM:SS (wall-clock)
            my $wall = ($1 * 3600 + $2 * 60 + $3) - $start_clk;
            return clock_to_samp($rec, $wall);
        }
        return $is_clock ? clock_to_samp($rec, 0 + $tok)   # bare sec = wall-clock
                         : int((0 + $tok) * $fs + 0.5);    # bare sec = data coord
    };

    my $default_base = ($in =~ s/\.[Ee][Ee][Gg]$//r) . '.mul';
    (my $base = ($out // $default_base)) =~ s/\.mul$//i;

    my $k = 0;
    my $T = qr/\d{1,2}:\d{2}:\d{2}|\d+(?:\.\d+)?/;   # HH:MM:SS or (fractional) seconds
    for my $part (split /\s*,\s*/, $spec) {
        next unless length $part;
        my ($a, $b, $name);
        if ($part =~ /^\s*($T)\s*-\s*($T)(?::(.*))?$/) {
            ($a, $b, $name) = ($1, $2, $3);
        } else {
            warn "skip bad range '$part'\n"; next;
        }
        my ($lo, $hi) = ($to_samp->($a), $to_samp->($b));
        ($lo, $hi) = ($hi, $lo) if $lo > $hi;
        if ($lo >= $n) {
            warn "skip '$part': start sample $lo >= data length $n\n"; next;
        }
        my $sub  = select_range($rec, $lo, $hi);
        my $file = (defined $name && length $name)
                 ? sprintf('%s_%s.mul', $base, $name)
                 : sprintf('%s_cut%02d.mul', $base, $k);
        printf "cut %s: samples %d..%d  (%.3f-%.3f s)\n",
            ($name // "#$k"), $lo, $hi, $lo / $fs, $hi / $fs;
        write_out($sub, $file);
        $k++;
    }
    exit 0;
}

# ----- whole-recording mode --------------------------------------------------
$out //= ($in =~ s/\.[Ee][Ee][Gg]$//r) . '.mul';
write_out($rec, $out);
printf "wrote %s : %d ch x %d samp @ %g Hz\n",
    $out, $rec->{data}->dim(0), $rec->{data}->dim(1), $rec->{fs};
