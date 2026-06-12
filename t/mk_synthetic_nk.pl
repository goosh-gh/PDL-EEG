#!/usr/bin/env perl
# t/mk_synthetic_nk.pl — synthetic EEG-1100C test file
#
# Confirmed layout from real hardware (YJ0394VB.EEG):
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
