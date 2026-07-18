#!/usr/bin/env perl
# t/mk_synthetic_nk.pl — synthetic EEG-1100C test file
#
# Confirmed layout from real hardware (subject.EEG):
#   +0x00        : 0x01
#   +0x01..+0x10 : "TIME164330000000"
#   +0x11        : 0x00
#   +0x12..+0x13 : 0x02 0x02 (version)
#   +0x14..+0x19 : BCD timestamp
#   +0x26        : n_ch_entries
#   +0x2F..      : n_ch_entries × 10-byte channel entries
#   +0x171       : data (uint16 LE offset binary, center=0x8000)
# n_ch_data = n_ch_entries + 1 (trailing zero-pad ch)

use strict; use warnings;
use File::Path qw(make_path);
use POSIX qw(floor);
use Getopt::Long;

# --long[=SEC] also writes t/data/*_long.eeg: same layouts, but long enough that
# the viewer's time slider actually has somewhere to go.
#
# These are NOT committed. 38 channels x 2 bytes x 1000 Hz is 76 KB per second of
# recording, so anything long enough to scroll through is megabytes -- t/data is
# git-tracked and the test suite does not need it. Generate it when you want to
# drive the viewer; add t/data/*_long.eeg to .gitignore.
my $LONG = 0;
GetOptions('long:i' => \$LONG) or die "usage: $0 [--long[=SECONDS]]\n";
$LONG = 300 if defined $LONG && $LONG == 0 && grep { /^--long/ } @ARGV;   # bare --long

make_path('t/data');

my $N_CH_ENTRIES = 4;    # valid channels (FP1,FP2,CZ,Fz equiv)
my $N_CH_DATA    = $N_CH_ENTRIES + 1;   # +1 zero-pad
my $N_SAMP       = 1000;
my $FS           = 1000;

# Channel indices (EEG-1100 lookup table)
my @CH_IDX = (0x01, 0x02, 0x10, 0x11);  # Fp1,Fp2,Fz,Cz (1-indexed)

my $DATA_REL  = 0x171;
my $WFM_ADDR  = 0x1800;  # must be > 0x17FE
my $DATA_ABS  = $WFM_ADDR + $DATA_REL;
my $TOTAL     = $DATA_ABS + $N_CH_DATA * $N_SAMP * 2;

my $buf = "\x00" x $TOTAL;
sub put { substr($buf, $_[0], length($_[1])) = $_[1] }
sub bcd { my $v=$_[0]; chr((($v/10)<<4)|($v%10)) }

# File header
put(0x0000, "EEG-1100C V01.00");
put(0x0081, "EEG-1100C V01.00");
put(0x0091, pack('C', 1));
put(0x0092, pack('V', 0x0200));

# Control block @0x0200
put(0x0200+17, pack('C', 1));
put(0x0200+18, pack('V', $WFM_ADDR));

# Waveform block @0x1800
put($WFM_ADDR+0x00, pack('C', 0x01));                  # block type
put($WFM_ADDR+0x01, "TIME164330000000");               # TIME string
put($WFM_ADDR+0x11, "\x00");                            # NUL
put($WFM_ADDR+0x12, pack('CC', 0x02, 0x02));           # version
put($WFM_ADDR+0x14, bcd(25).bcd(12).bcd(21));          # YY MM DD
put($WFM_ADDR+0x17, bcd(16).bcd(43).bcd(30));          # HH MM SS
put($WFM_ADDR+0x1A, pack('v', 1000 | (3<<14)));        # fs: lower14=1000Hz, upper2=flags
put($WFM_ADDR+0x1C, pack('v', 340));                   # DMA chunk hint
put($WFM_ADDR+0x26, pack('C', $N_CH_ENTRIES));         # n_ch_entries

# Channel table: 10-byte entries at +0x2F
for my $i (0 .. $N_CH_ENTRIES-1) {
    my $eoff = $WFM_ADDR + 0x2F + $i*10;
    put($eoff, pack('CC', 0x10, 0x05));
    put($eoff+2, pack('C', $CH_IDX[$i]));
}

# Sample data: uint16 LE offset binary (0x8000 = 0µV)
# Synthetic: sine waves
my @samples;
for my $s (0..$N_SAMP-1) {
    for my $c (0..$N_CH_ENTRIES-1) {
        my $hz  = ($c+1)*5;
        my $amp = 100;   # 100µV amplitude
        my $val = int($amp * sin(2*3.14159265*$hz*$s/$FS)) + 0x8000;
        push @samples, $val;
    }
    push @samples, 0x8000;   # zero-pad ch: 0µV = 0x8000
}
put($DATA_ABS, pack('v*', @samples));

# Verify
die "BUG: sig check" unless unpack('C', substr($buf, $WFM_ADDR,1)) == 0x01;
die "BUG: 0x17FE not in wfmblock" unless $WFM_ADDR == 0x1800;

open my $fh, '>:raw', 't/data/test01.eeg' or die $!;
print $fh $buf; close $fh;

printf "Written: t/data/test01.eeg (%d bytes)\n", $TOTAL;
printf "  wfmblock=0x%04X  data=0x%04X\n", $WFM_ADDR, $DATA_ABS;
printf "  n_ch_entries=%d  n_ch_data=%d  N_SAMP=%d\n",
    $N_CH_ENTRIES, $N_CH_DATA, $N_SAMP;

# LOG file
my @events = ({t=>10, label=>'EYES OPEN'},{t=>30, label=>'HYPERVENT'});
my $LOG_LB   = 0x0200;
my $LOG_SIZE = $LOG_LB + 0x0014 + @events*45 + 64;
my $lbuf = "\x00" x $LOG_SIZE;
sub lput { substr($lbuf,$_[0],length($_[1]))=$_[1] }
lput(0x0000,"EEG-1100C V01.00");
lput(0x0091,pack('C',1));
lput(0x0092,pack('V',$LOG_LB));
lput($LOG_LB+0x0012,pack('C',scalar@events));
for my $i (0..$#events) {
    my $off = $LOG_LB + 0x0014 + $i * 45;
    lput($off,      sprintf("%-20s", $events[$i]{label}));
    lput($off+20, sprintf("%06d",  $events[$i]{t}));   # ASCII 6-digit seconds at byte 20
}
open $fh,'>:raw','t/data/test01.LOG' or die $!;
print $fh $lbuf; close $fh;
printf "Written: t/data/test01.LOG (%d bytes)\n", $LOG_SIZE;

# ===========================================================================
# t/data/test02.eeg — MULTI-BLOCK EEG-1100C file (for block_extents / --cut).
#
# test01 has a single waveform block, so it cannot exercise block boundaries,
# the gap_samples padding read_nk inserts between blocks, or a partial read that
# straddles a boundary. test02 has THREE blocks of DELIBERATELY UNEQUAL length,
# so uniform-length assumptions and off-by-one errors cannot hide.
#
# n_samp is NOT stored on disk: the reader derives it from the gap to the next
# block address (EOF for the last). The addresses below are laid out so that
# (next_addr - waddr - 0x171) is exactly n_ch_data * 2 * n_samp.
#
# Samples encode their own ABSOLUTE (concatenated) index, so a test can prove a
# glued partial read landed on the right samples:
#     raw(abs, ch) = 0x8000 + ((abs*7 + ch*13) % 2000) - 1000
# ===========================================================================
{
    my @BLK_LEN  = (1000, 2500, 700);                  # samples per block
    my @BLK_TIME = ([25,12,21,16,43,30],               # YY MM DD HH MM SS
                    [25,12,21,16,44,30],
                    [25,12,21,16,50,00]);
    my $CTL      = 0x0200;
    my $W0       = 0x1800;                             # must be > 0x17FE

    my @waddr;
    my $a = $W0;
    for my $b (0 .. $#BLK_LEN) {
        push @waddr, $a;
        $a += $DATA_REL + $N_CH_DATA * 2 * $BLK_LEN[$b];
    }
    my $total = $a;                                    # EOF == end of last block

    my $b2 = "\x00" x $total;
    my $put2 = sub { substr($b2, $_[0], length $_[1]) = $_[1] };

    $put2->(0x0000, "EEG-1100C V01.00");
    $put2->(0x0081, "EEG-1100C V01.00");
    $put2->(0x0091, pack('C', 1));                     # one control block
    $put2->(0x0092, pack('V', $CTL));
    $put2->(0x03EE, pack('V', 0));                     # ext_address = 0 -> wfmblock

    $put2->($CTL + 17, pack('C', scalar @waddr));      # data block count
    $put2->($CTL + 18 + $_ * 20, pack('V', $waddr[$_])) for 0 .. $#waddr;

    my $abs = 0;                                       # concatenated sample index
    for my $b (0 .. $#BLK_LEN) {
        my $w = $waddr[$b];
        $put2->($w + 0x00, pack('C', 0x01));
        $put2->($w + 0x01, "TIME164330000000");
        $put2->($w + 0x12, pack('CC', 0x02, 0x02));
        $put2->($w + 0x14, join '', map { bcd($_) } @{ $BLK_TIME[$b] });
        $put2->($w + 0x1A, pack('v', $FS | (3 << 14)));
        $put2->($w + 0x1C, pack('v', 340));
        $put2->($w + 0x26, pack('C', $N_CH_ENTRIES));
        for my $i (0 .. $N_CH_ENTRIES - 1) {
            my $eoff = $w + 0x2F + $i * 10;
            $put2->($eoff,     pack('CC', 0x10, 0x05));
            $put2->($eoff + 2, pack('C', $CH_IDX[$i]));
        }

        my @s;
        for my $t (0 .. $BLK_LEN[$b] - 1) {
            my $x = $abs + $t;
            push @s, 0x8000 + ((($x * 7 + $_ * 13) % 2000) - 1000)
                for 0 .. $N_CH_ENTRIES - 1;
            push @s, 0x8000;                           # zero-pad channel
        }
        $put2->($w + $DATA_REL, pack('v*', @s));
        $abs += $BLK_LEN[$b];
    }

    open my $fh2, '>:raw', 't/data/test02.eeg' or die $!;
    print {$fh2} $b2; close $fh2;

    printf "Written: t/data/test02.eeg (%d bytes)\n", $total;
    printf "  %d blocks, n_samp = %s (unequal on purpose)\n",
        scalar @BLK_LEN, join('/', @BLK_LEN);
    printf "  wfmblock addrs = %s\n", join(', ', map { sprintf('0x%04X', $_) } @waddr);
}

# ===========================================================================
# t/data/test03.eeg — MULTI-SEGMENT EXTBLOCK (EEG-1200A) file.
#
# This is the fixture that would have caught the bug that ate a week.
#
# An extblock file is NOT one contiguous data block. At every recording break
# the recorder writes a fresh copy of the channel-info block -- 72 + (n_ch-1)*10
# bytes, identical to the one at eb3 except for its timestamp -- straight into
# the sample stream and carries on. read_nk() used to read those headers as if
# they were EEG samples, so the channel phase slipped by hdr_len % stride bytes
# (442 % 76 = 62 = 31 channels) at every break, and from the first break onward
# every channel label sat on another channel's data.
#
#   [ segment 0 samples ][ 442-byte header ][ segment 1 samples ][ header ] ...
#
# So the fixture needs:
#   * MORE THAN ONE segment (one segment cannot catch a per-break bug)
#   * segments of DIFFERENT lengths (equal lengths hide off-by-one errors)
#   * DC channels carrying SQUARE-WAVE TTL, and EEG channels carrying sine
#     -- if the phase slips, a trigger appears on an EEG label and vice versa,
#     which is exactly what the real recording showed
#   * a distinct wall-clock t_start per segment, with real dead time between
#     them, so event placement can be checked against the true anchors
#
# Layout (confirmed byte-for-byte against a real EEG-1200A recording):
#   channel-info block, 72 + (n_ch-1)*10 bytes:
#     +0x00        : 0x01  block type
#     +0x01..+0x10 : ASCII "TIME" + HHMMSS + zero padding
#     +0x12..+0x13 : 0x02 0x02  version
#     +0x14..+0x27 : ASCII "YYYYMMDDHHMMSS" + zero padding  <- segment start time
#     +0x28..+0x2B : u32 sampling rate
#     +0x44..+0x45 : u16 n_ch - 1
#     +0x48..      : channel table, (n_ch-1) x 10 bytes: [u16 hw-1][6 x 00][00 05]
# ===========================================================================
{
    my $E_NCH   = 38;                      # incl. STIM
    my $E_FS    = 1000;
    my $E_HDR   = 72 + ($E_NCH - 1) * 10;  # 442
    my ($EXT, $EB2, $EB3, $ECTL, $EDA) = (0x27CF, 0x2BFF, 0x43FB, 0x400, 0x17FE);
    my $E_REC   = $EB3 + $E_HDR;           # 0x45B5 -- first sample of segment 0
    my @E_HW    = (1..20, 23..30, 45, 46, 47, 48, 75, 76, 77, 78, 100);
    die "extblock fixture: hw list must have n_ch-1 entries" unless @E_HW == $E_NCH - 1;

    # segments: [ YYYYMMDDHHMMSS, n_samp ] -- unequal on purpose.
    # 14:03:03 + 5.000 s of data, then a 284 s break -> 14:07:52
    # 14:07:52 + 3.000 s of data, then a 219 s break -> 14:11:34
    my @E_SEG = ( [ '20260702140303', 5000 ],
                  [ '20260702140752', 3000 ],
                  [ '20260702141134', 2000 ] );

    my $mk_block = sub {                   # the 442-byte channel-info block
        my $ts = shift;                    # YYYYMMDDHHMMSS
        my $b  = "\x00" x $E_HDR;
        substr($b, 0x00, 1)    = pack('C', 0x01);
        substr($b, 0x01, 16)   = sprintf('%-16.16s', 'TIME' . substr($ts, 8, 6) . '000000');
        substr($b, 0x12, 2)    = pack('CC', 0x02, 0x02);
        substr($b, 0x14, 20)   = sprintf('%-20.20s', $ts . '000000');
        substr($b, 0x28, 4)    = pack('V', $E_FS);
        substr($b, 0x44, 2)    = pack('v', $E_NCH - 1);
        for my $i (0 .. $#E_HW) {
            my $o = 72 + $i * 10;
            substr($b, $o,     2) = pack('v', $E_HW[$i] - 1);   # 0-based hw index
            substr($b, $o + 8, 2) = pack('v', 0x0500);
        }
        return $b;
    };

    my $ebuf = "\x00" x $E_REC;
    my $eput = sub { substr($ebuf, $_[0], length $_[1]) = $_[1] };
    $eput->(0x0000, 'EEG-1200A V01.00');
    $eput->(0x03EE, pack('V', $EXT));
    $eput->($EXT + 18,  pack('V', $EB2));
    $eput->($EB2 + 20,  pack('V', $EB3));
    $eput->(0x0092,     pack('V', $ECTL));
    $eput->($ECTL + 18, pack('V', $EDA));
    $eput->($EDA + 0x1A, pack('v', $E_FS | (3 << 14)));
    $eput->($EDA + 0x14, join '', map { bcd($_) } (26, 7, 2, 14, 3, 3));   # BCD fallback
    $eput->($EB3, $mk_block->($E_SEG[0][0]));

    # DC channels are the 10th..7th from the end of the analog list: hw 45..48.
    # Find their DATA-COLUMN indices so the TTL really lands on the DC labels.
    my %col;  $col{ $E_HW[$_] } = $_ for 0 .. $#E_HW;
    my @DC = map { $col{$_} } (45, 46, 47, 48);

    my $body = '';
    for my $s (0 .. $#E_SEG) {
        $body .= $mk_block->($E_SEG[$s][0]) if $s;        # header at every break
        my ($ts, $n) = @{ $E_SEG[$s] };
        my @v;
        for my $t (0 .. $n - 1) {
            my @row;
            for my $c (0 .. $E_NCH - 2) {
                push @row, 0x8000 + int(200 * sin(2 * 3.14159265358979 * ($c + 1) * $t / $E_FS));
            }
            # square-wave TTL on the DC channels: 250-sample high/low
            $row[$_] = 0x8000 + ((int($t / 250) % 2) ? 9000 : 0) for @DC;
            push @row, 0;                                  # STIM
            push @v, @row;
        }
        $body .= pack('v*', @v);
    }

    open my $efh, '>:raw', 't/data/test03.eeg' or die $!;
    print {$efh} $ebuf . $body;
    close $efh;

    my $tot = 0; $tot += $_->[1] for @E_SEG;
    printf "Written: t/data/test03.eeg (%d bytes)\n", length($ebuf) + length($body);
    printf "  extblock, %d segments, n_samp = %s (unequal on purpose), total %d\n",
        scalar @E_SEG, join('/', map { $_->[1] } @E_SEG), $tot;
    printf "  channel-info block = 72 + (%d-1)*10 = %d bytes; %d %% %d = %d (the phase slip)\n",
        $E_NCH, $E_HDR, $E_HDR, $E_NCH * 2, $E_HDR % ($E_NCH * 2);
    printf "  square-wave TTL on data columns %s (hw 45-48 = DC01-DC04)\n", join(',', @DC);
}

# ===========================================================================
# --long : the same two layouts, scaled up so the viewer has something to scroll.
# Not part of the test suite; not committed.
# ===========================================================================
if ($LONG) {
    my $sec = $LONG;
    printf "\n--- --long: generating %d s versions (NOT for git) ---\n", $sec;

    # ---- wfmblock (EEG-1100C), 4 blocks of unequal length -------------------
    {
        my $FS_L    = 500;                     # halve the rate: same duration, half the bytes
        my @FRAC    = (0.30, 0.35, 0.20, 0.15);
        my @LEN     = map { int($sec * $FS_L * $_) } @FRAC;
        my $HDR_L   = 0x171;
        my $NCH_D   = $N_CH_DATA;
        my @waddr;
        my $a = 0x1800;
        for my $b (0 .. $#LEN) { push @waddr, $a; $a += $HDR_L + $NCH_D * 2 * $LEN[$b] }
        my $buf = "\x00" x $a;
        my $p = sub { substr($buf, $_[0], length $_[1]) = $_[1] };
        $p->(0x0000, 'EEG-1100C V01.00');
        $p->(0x0081, 'EEG-1100C V01.00');
        $p->(0x0091, pack('C', 1));
        $p->(0x0092, pack('V', 0x0200));
        $p->(0x03EE, pack('V', 0));
        $p->(0x0200 + 17, pack('C', scalar @waddr));
        $p->(0x0200 + 18 + $_ * 20, pack('V', $waddr[$_])) for 0 .. $#waddr;
        my @T = ([25,12,21,16,43,30], [25,12,21,16,50,10],
                 [25,12,21,17,05,00], [25,12,21,17,20,45]);
        my $abs = 0;
        for my $b (0 .. $#LEN) {
            my $w = $waddr[$b];
            $p->($w, pack('C', 0x01));
            $p->($w + 0x01, 'TIME164330000000');
            $p->($w + 0x12, pack('CC', 0x02, 0x02));
            $p->($w + 0x14, join '', map { bcd($_) } @{ $T[$b] });
            $p->($w + 0x1A, pack('v', $FS_L | (3 << 14)));
            $p->($w + 0x26, pack('C', $N_CH_ENTRIES));
            for my $i (0 .. $N_CH_ENTRIES - 1) {
                my $e = $w + 0x2F + $i * 10;
                $p->($e, pack('CC', 0x10, 0x05));
                $p->($e + 2, pack('C', $CH_IDX[$i]));
            }
            my @s;
            for my $t (0 .. $LEN[$b] - 1) {
                my $x = $abs + $t;
                for my $c (0 .. $N_CH_ENTRIES - 1) {
                    push @s, 0x8000 + int(300 * sin(2 * 3.14159265358979 * ($c + 1) * 3 * $x / $FS_L));
                }
                push @s, 0x8000;
            }
            $p->($w + $HDR_L, pack('v*', @s));
            $abs += $LEN[$b];
        }
        open my $fh2, '>:raw', 't/data/test02_long.eeg' or die $!;
        print {$fh2} $buf; close $fh2;
        printf "Written: t/data/test02_long.eeg (%.1f MB)  wfmblock, %d blocks, %s s @ %d Hz\n",
            $a / 1048576, scalar @LEN,
            join('/', map { sprintf('%.0f', $_ / $FS_L) } @LEN), $FS_L;
    }

    # ---- extblock (EEG-1200A), 4 segments, TTL on the DC channels ------------
    {
        my $E_NCH = 38;
        my $E_FS  = 500;
        my $E_HDR = 72 + ($E_NCH - 1) * 10;
        my ($EXT, $EB2, $EB3, $ECTL, $EDA) = (0x27CF, 0x2BFF, 0x43FB, 0x400, 0x17FE);
        my $E_REC = $EB3 + $E_HDR;
        my @E_HW  = (1..20, 23..30, 45, 46, 47, 48, 75, 76, 77, 78, 100);
        my @FRAC  = (0.30, 0.30, 0.25, 0.15);
        my @LEN   = map { int($sec * $E_FS * $_) } @FRAC;
        my @TS    = ('20260702140303', '20260702140752',
                     '20260702141134', '20260702141641');

        my $mk = sub {
            my $ts = shift;
            my $b  = "\x00" x $E_HDR;
            substr($b, 0x00, 1)  = pack('C', 0x01);
            substr($b, 0x01, 16) = sprintf('%-16.16s', 'TIME' . substr($ts, 8, 6) . '000000');
            substr($b, 0x12, 2)  = pack('CC', 0x02, 0x02);
            substr($b, 0x14, 20) = sprintf('%-20.20s', $ts . '000000');
            substr($b, 0x28, 4)  = pack('V', $E_FS);
            substr($b, 0x44, 2)  = pack('v', $E_NCH - 1);
            for my $i (0 .. $#E_HW) {
                my $o = 72 + $i * 10;
                substr($b, $o, 2)     = pack('v', $E_HW[$i] - 1);
                substr($b, $o + 8, 2) = pack('v', 0x0500);
            }
            return $b;
        };

        my $eb = "\x00" x $E_REC;
        my $ep = sub { substr($eb, $_[0], length $_[1]) = $_[1] };
        $ep->(0x0000, 'EEG-1200A V01.00');
        $ep->(0x03EE, pack('V', $EXT));
        $ep->($EXT + 18,  pack('V', $EB2));
        $ep->($EB2 + 20,  pack('V', $EB3));
        $ep->(0x0092,     pack('V', $ECTL));
        $ep->($ECTL + 18, pack('V', $EDA));
        $ep->($EDA + 0x1A, pack('v', $E_FS | (3 << 14)));
        $ep->($EDA + 0x14, join '', map { bcd($_) } (26, 7, 2, 14, 3, 3));
        $ep->($EB3, $mk->($TS[0]));

        my %col; $col{ $E_HW[$_] } = $_ for 0 .. $#E_HW;
        my @DC = map { $col{$_} } (45, 46, 47, 48);

        my $body = '';
        for my $s (0 .. $#LEN) {
            $body .= $mk->($TS[$s]) if $s;
            my @v;
            for my $t (0 .. $LEN[$s] - 1) {
                my @row;
                for my $c (0 .. $E_NCH - 2) {
                    push @row, 0x8000 + int(200 * sin(2 * 3.14159265358979 * ($c + 1) * $t / $E_FS));
                }
                # 9000 counts x 366.3 uV/bit = 3.3 V -- a real TTL, not a 3.3 mV one
                $row[$_] = 0x8000 + ((int($t / ($E_FS / 4)) % 2) ? 9000 : 0) for @DC;
                push @row, 0;
                push @v, @row;
            }
            $body .= pack('v*', @v);
        }
        open my $fh3, '>:raw', 't/data/test03_long.eeg' or die $!;
        print {$fh3} $eb . $body; close $fh3;
        my $tot = 0; $tot += $_ for @LEN;
        printf "Written: t/data/test03_long.eeg (%.1f MB)  extblock, %d segments, %s s @ %d Hz\n",
            (length($eb) + length($body)) / 1048576, scalar @LEN,
            join('/', map { sprintf('%.0f', $_ / $E_FS) } @LEN), $E_FS;
    }
    print "  (add t/data/*_long.eeg to .gitignore -- these are not test fixtures)\n";
}
