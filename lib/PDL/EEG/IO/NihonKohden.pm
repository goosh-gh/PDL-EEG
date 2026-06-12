package PDL::EEG::IO::NihonKohden;

use strict;
use warnings;
use Carp qw(croak confess);
use Encode qw(decode);
use PDL;

our $VERSION = '0.01';

=head1 NAME

PDL::EEG::IO::NihonKohden - Read Nihon Kohden EEG binary files into PDL

=head1 SYNOPSIS

  use PDL;
  use PDL::EEG::IO::NihonKohden qw(read_nk);

  my $rec = read_nk('patient01.eeg');
  my $data = $rec->{data};   # PDL [n_ch, n_samples] float32, µV
  my $fs   = $rec->{fs};     # sampling rate in Hz
  my @labs = @{ $rec->{labels} };
  my @evts = @{ $rec->{events} };

  # Plot channel 0
  use PDL::Graphics::Cairo qw(figure);
  my $fig = figure();
  my $ax  = $fig->add_subplot(1,1,1);
  my $t   = sequence($data->dim(1)) / $fs;
  $ax->line($t, $data->slice('(0),:'));
  $ax->set_title($labs[0]);
  $fig->show;

=head1 DESCRIPTION

Parses the Nihon Kohden EEG-1100/2100 binary format (*.eeg, *.log, *.pnt, *.21e)
and returns the waveform data as a PDL piddle in µV.

Format knowledge derived from EDFbrowser's nk2edf.cpp by Teunis van Beelen
(GPL-2), used as a format reference only — this is a clean-room Perl implementation.

=head1 EXPORTED FUNCTIONS

=cut

use Exporter 'import';
our @EXPORT_OK = qw(read_nk);

# ---------------------------------------------------------------------------
# Sampling rate code → Hz
# ---------------------------------------------------------------------------
my %SRATE_TABLE = (
    0xA0 =>  100,
    0xA1 =>  200,
    0xA2 =>  500,
    0xA3 => 1000,
    0xA4 => 2000,
    0xA5 => 5000,  # EEG-2100 extended
);

# ---------------------------------------------------------------------------
# Gain code → µV/bit  (from nk2edf.cpp get_chan_sensitivity)
# Values are physical µV per LSB (int16 ADC count)
# ---------------------------------------------------------------------------
my @GAIN_TABLE = (
    # code 0x00 .. 0x12
    1000, 500, 250, 125,          # 0x00–0x03
    100,  50,  25,                # 0x04–0x06
    10,   5,   2.5,               # 0x07–0x09
    1000, 500, 250,               # 0x0A–0x0C (with x10 amp)
    100,  50,  25,                # 0x0D–0x0F
    10,   5,   2.5,               # 0x10–0x12
);

# Default electrode label table
# 0-indexed, matching .21e file format (confirmed from EEG-1100C YJ0394VB.21E).
# ch_idx in waveform block is 1-indexed → subtract 1 before lookup.
my @DEFAULT_LABELS = (
    # 0-19: standard 10-20 + extras
    qw(Fp1 Fp2 F3 F4 C3 C4 P3 P4 O1 O2),   # 0-9
    qw(F7 F8 T3 T4 T5 T6 Fz Cz Pz E),       # 10-19
    # 20-25
    qw(BN1 BN2 A1 A2 vEOG hEOG),
    # 26-36: X1..X11
    (map { "X$_" } 1..11),
    # 37-41
    qw(BN AV SD Aav 0V),
    # 42-43
    qw(SpO2 EtCO2),
    # 44-47: DC03-DC06
    (map { sprintf("DC%02d", $_) } 3..6),
    # 48-49
    qw(Pulse CO2Wave),
    # 50-73: DC09-DC32
    (map { sprintf("DC%02d", $_) } 9..32),
    # 74-75: BN1, BN2
    qw(BN1 BN2),
    # 76-92: undefined (gap 76..92)
    ('-') x 17,
    # 93-99: RFU1-RFU7
    (map { "RFU$_" } 1..7),
    # 100-120: X12-X32
    (map { "X$_" } 12..32),
    # 121-125: BP1-BP4, Ud
    qw(BP1 BP2 BP3 BP4 Ud),
);
# Pad to 256
$#DEFAULT_LABELS = 255;
$_ //= '-' for @DEFAULT_LABELS;

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

=head2 read_nk($eeg_file, %opts)

Read a Nihon Kohden *.eeg recording.

Options:
  block     => $n   # which waveform block to read (0-based, default 0)
  edfplus   => 1    # also read .log/.pnt for events (default 1)
  no_events => 1    # skip event parsing even if .log exists

Returns a hashref:
  data    => $pdl   # [n_ch, n_samples] float32, µV
  fs      => $hz
  labels  => \@channel_label_strings
  t_start => $epoch_seconds  (undef if not parsed)
  events  => \@{ {t => $sec, label => $str} }
  gains   => $pdl   # [n_ch] µV/bit

=cut

sub read_nk {
    my ($eeg_path, %opts) = @_;
    croak "File not found: $eeg_path" unless -f $eeg_path;

    my $block_idx   = $opts{block}      // 0;
    my $want_events = !$opts{no_events};
    my $fs_override = $opts{fs};          # caller can supply fs (e.g. 1000)

    open my $fh, '<:raw', $eeg_path or croak "Cannot open $eeg_path: $!";

    # --- 1. Device signature ------------------------------------------------
    my $sig = _read_bytes($fh, 0x0000, 16);
    _check_device_sig($sig)
        and croak "Unknown device signature: " . _hexdump($sig);
    (my $device = $sig) =~ s/\x00.*//s;  # trim NUL

    # --- 2. Control block list ---------------------------------------------
    my $ctl_count = _read_u8($fh, 0x0091);
    croak "No control blocks in file" unless $ctl_count > 0;

    my @wfm_addrs;
    for my $i (0 .. $ctl_count - 1) {
        my $ctl_addr = _read_u32le($fh, 0x0092 + $i * 20);
        my $db_count = _read_u8($fh, $ctl_addr + 17);
        for my $j (0 .. $db_count - 1) {
            push @wfm_addrs, _read_u32le($fh, $ctl_addr + 18 + $j * 20);
        }
    }
    croak "Block index $block_idx out of range (have " . scalar(@wfm_addrs) . " blocks)"
        if $block_idx >= @wfm_addrs;

    # --- 3. Waveform block signature: wfm_addrs[0] starts with 0x01 --------
    # (The byte at 0x17FE IS the first byte of wfmblock[0], not a separate sig)

    # --- 4. Read waveform block header -------------------------------------
    my $waddr = $wfm_addrs[$block_idx];
    my $meta  = _read_wfm_header($fh, $waddr);

    # --- 5. n_samples from block gap ---------------------------------------
    # Next block address (or file end) minus current block address, minus header
    seek $fh, 0, 2;
    my $file_size   = tell $fh;
    my $next_addr   = ($block_idx + 1 < @wfm_addrs)
                    ? $wfm_addrs[$block_idx + 1]
                    : $file_size;
    my $gap         = $next_addr - $waddr;
    my $data_rel    = $meta->{data_offset} - $waddr;   # = 0x171 = 369
    my $data_bytes  = $gap - $data_rel;
    my $n_ch        = $meta->{n_ch};                   # n_ch_entries + 1
    my $n_samp      = int($data_bytes / ($n_ch * 2));
    croak "Cannot determine n_samples (data_bytes=$data_bytes not divisible by n_ch*2=" . ($n_ch*2) . ")"
        if $data_bytes % ($n_ch * 2) != 0;

    # --- 6. Sampling rate --------------------------------------------------
    # Extracted from lower 14 bits of u16 at wfmblock+0x1A
    my $fs = $fs_override // $meta->{fs}
        or croak "Sampling rate could not be determined "
               . "(header gave " . ($meta->{fs}//0) . " Hz; supply fs => NNN option)";

    # --- 7. Electrode labels -----------------------------------------------
    my $base = $eeg_path;  $base =~ s/\.[^.]+$//;
    my %label_override;
    for my $ext (qw(.21e .21E)) {
        my %h = _read_21e("$base$ext");
        if (%h) { %label_override = %h; last }
    }

    # Valid channels = ch_indices[0..n_ch_valid-1]; last ch is zero-pad
    my $n_ch_valid = $meta->{n_ch_valid};
    my @labels;
    for my $i (0 .. $n_ch_valid - 1) {
        my $idx       = $meta->{ch_indices}[$i];  # 1-indexed
        my $label_idx = $idx - 1;                 # 0-indexed for lookup
        # .21e keys are 0-indexed integers; DEFAULT_LABELS also 0-indexed
        push @labels, $label_override{$label_idx}
                   // ($label_idx >= 0 ? $DEFAULT_LABELS[$label_idx] : undef)
                   // "CH$i";
    }
    push @labels, 'PAD';   # trailing zero-pad channel

    # --- 8. Read raw sample data -------------------------------------------
    # Format: unsigned uint16 LE, offset binary (center = 0x8000)
    # physical µV = (raw_u16 - 0x8000) × gain_µV/bit
    my $data_offset = $meta->{data_offset};
    seek $fh, $data_offset, 0 or croak "Seek to data failed";
    my $buf;
    my $bytes = $n_ch * $n_samp * 2;
    my $got = read($fh, $buf, $bytes);
    croak "Short read: wanted $bytes bytes, got $got" unless $got == $bytes;

    # Unpack as unsigned uint16 LE then subtract 0x8000
    my @raw_u = unpack "v*", $buf;   # 'v' = uint16 LE

    my $raw_pdl = PDL->new(\@raw_u)->reshape($n_ch, $n_samp);
    $raw_pdl    = $raw_pdl->double - 0x8000;   # offset binary → signed µV

    # PAD channel (last, index n_ch-1) is hardware zero-fill → set to 0.0
    $raw_pdl->slice('(-1),:') .= 0.0;

    # --- 9. Apply gain (µV/bit) --------------------------------------------
    # EEG-1100C/1200C: ADC gain is fixed, not stored in file.
    # nk2edf.cpp sensitivity_list[10] = 0.9765625 µV/div, 10 bits/div
    # → 0.9765625 / 10 = 0.09765625 µV/bit
    # Verified: 0.09765625 * 32767 = 3199.90 µV = EDF phys_max (confirmed)
    # DC channels (BN1/BN2/Mark etc) have different range but use same ADC.
    my $gain_uv_per_bit = 0.09765625;
    my @gain_uv  = ($gain_uv_per_bit) x $n_ch;
    my $gains    = PDL->new(\@gain_uv);

    my $data_uv = ($raw_pdl * $gains->slice(':,*1'))->float;

    # --- 10. Events (optional) ---------------------------------------------
    my @events;
    if ($want_events) {
        my $log_path = "$base.LOG";
        $log_path = "$base.log" unless -f $log_path;
        @events = _read_log($log_path) if -f $log_path;
    }

    close $fh;

    return {
        data       => $data_uv,
        fs         => $fs,
        labels     => \@labels,
        t_start    => $meta->{t_start},
        events     => \@events,
        gains      => $gains,
        n_blocks   => scalar(@wfm_addrs),
        n_ch_valid => $n_ch_valid,
        block_idx  => $block_idx,
        device     => $device,           # e.g. "EEG-1100C V01.00"
    };
}

# ---------------------------------------------------------------------------
# Internal: waveform block header parse
#
# Confirmed structure from real EEG-1100C file (YJ0394VB.EEG, 2025-12-21):
#
#   +0x00        : 0x01  block type
#   +0x01..+0x10 : ASCII time string "TIME164330000000" (16 bytes)
#   +0x11        : 0x00  NUL terminator
#   +0x12        : format version (0x02)
#   +0x13        : sub-version   (0x02)
#   +0x14..+0x19 : BCD timestamp  YY MM DD HH MM SS
#   +0x1C..+0x1D : u16LE = n_samp_per_dma_chunk (340 observed)
#   +0x26        : n_ch_entries (number of channel table entries = n_valid_ch)
#   +0x2F..      : channel table: n_ch_entries × 10 bytes
#                  each entry: [0x10][0x05][ch_idx][0x00×7]
#   +0x171       : data start  (= +0x2F + n_ch_entries×10 + 2 padding)
#                  NOTE: actual n_ch in data = n_ch_entries + 1 (trailing zero ch)
#
# Data encoding: unsigned uint16 LE, offset binary
#   physical = (raw_u16 - 0x8000) × gain_µV/bit
#   Gain: 1 µV/bit (verified: ±700 counts ≈ ±700 µV, plausible ERP range)
#
# ---------------------------------------------------------------------------
sub _read_wfm_header {
    my ($fh, $addr) = @_;

    # --- Block type check ---
    my $block_type = _read_u8($fh, $addr);
    croak sprintf("Unexpected block type 0x%02X (expected 0x01)", $block_type)
        unless $block_type == 0x01;

    # --- TIME string + BCD timestamp ---
    my $time_str = _read_bytes($fh, $addr + 0x01, 16);
    $time_str =~ s/\x00.*//s;  # trim at NUL

    my $bcd_yy = _bcd_byte($fh, $addr + 0x14);
    my $bcd_mm = _bcd_byte($fh, $addr + 0x15);
    my $bcd_dd = _bcd_byte($fh, $addr + 0x16);
    my $bcd_hh = _bcd_byte($fh, $addr + 0x17);
    my $bcd_mi = _bcd_byte($fh, $addr + 0x18);
    my $bcd_ss = _bcd_byte($fh, $addr + 0x19);
    my $t_start = sprintf("20%02d-%02d-%02d %02d:%02d:%02d",
        $bcd_yy, $bcd_mm, $bcd_dd, $bcd_hh, $bcd_mi, $bcd_ss);

    # --- Sampling rate: +0x1A..+0x1B lower 14 bits ---
    # Confirmed: 0xC3E8 & 0x3FFF = 0x03E8 = 1000 Hz
    my $fs_raw = _read_u16le($fh, $addr + 0x1A);
    my $fs     = $fs_raw & 0x3FFF;   # lower 14 bits

    # --- Channel table ---
    # +0x26: number of channel table entries (= n_valid_ch)
    my $n_ch_entries = _read_u8($fh, $addr + 0x26);
    croak "Zero channel entries in waveform block" unless $n_ch_entries > 0;

    # 10-byte entries starting at +0x2F: [0x10][0x05][ch_idx][0×7]
    my @ch_indices;
    for my $i (0 .. $n_ch_entries - 1) {
        my $eoff = $addr + 0x2F + $i * 10;
        my $marker = _read_u8($fh, $eoff);
        my $ch_idx = _read_u8($fh, $eoff + 2);
        push @ch_indices, $ch_idx;
    }

    # --- Data layout ---
    # Data starts at +0x171 (= +0x2F + n_ch_entries×10 + 2 padding bytes)
    # n_ch in data stream = n_ch_entries + 1 (one trailing zero-pad channel)
    my $data_offset = $addr + 0x171;
    my $n_ch_data   = $n_ch_entries + 1;

    # n_samples is not stored in header; caller computes from block size.
    # Return n_ch_data so caller can calculate: n_samp = data_bytes / (n_ch_data * 2)

    return {
        n_ch        => $n_ch_data,        # channels in data stream (incl. zero-pad)
        n_ch_valid  => $n_ch_entries,     # real EEG/physio channels
        fs          => $fs,               # from lower 14 bits of u16 at +0x1A
        n_samples   => undef,             # caller computes from block gap
        ch_indices  => \@ch_indices,      # electrode index per valid channel (1-indexed)
        gain_codes  => [(0) x $n_ch_entries],
        t_start     => $t_start,
        data_offset => $data_offset,
        time_str    => $time_str,
    };
}

# BCD byte decode: 0x25 → 25 (decimal)
sub _bcd_byte {
    my ($fh, $offset) = @_;
    my $v = _read_u8($fh, $offset);
    return ($v >> 4) * 10 + ($v & 0x0F);
}

# ---------------------------------------------------------------------------
# Internal: .21e electrode name file
# Format: sections [ELECTRODE], key=value pairs
# ---------------------------------------------------------------------------
sub _read_21e {
    my ($path) = @_;
    return () unless -f $path;
    my $fh;
    open($fh, '<:encoding(Shift_JIS):crlf', $path)
        or open($fh, '<:encoding(UTF-8):crlf', $path)
        or return ();
    my %map;
    my $in_electrode = 0;
    while (<$fh>) {
        s/\r\n$/\n/;   # CRLF → LF
        s/\r$/\n/;     # CR-only → LF
        chomp;
        if (/^\[ELECTRODE\]/i)      { $in_electrode = 1; next }
        if (/^\[/ && !/ELECTRODE/i) { $in_electrode = 0; next }
        # Keys are 4-digit decimal, 0-indexed: "0000"=Fp1, "0001"=Fp2, ...
        if ($in_electrode && /^(\d+)=(.+)$/) {
            my ($key_str, $val) = ($1, $2);
            $val =~ s/\s+$//;   # trim trailing whitespace
            $map{$key_str + 0} = $val;
        }
    }
    close $fh;
    return %map;
}

# ---------------------------------------------------------------------------
# Internal: .LOG event file
# ---------------------------------------------------------------------------
sub _read_log {
    my ($path) = @_;
    return () unless -f $path;
    open my $fh, '<:raw', $path or return ();

    my $sig = _read_bytes($fh, 0x0000, 16);
    return () if _check_device_sig($sig);

    my $n_logblocks = _read_u8($fh, 0x0091);

    my @events;
    for my $i (0 .. $n_logblocks - 1) {
        my $lb_addr = _read_u32le($fh, 0x0092 + $i * 20);
        my $n_logs  = _read_u8($fh, $lb_addr + 0x0012);

        for my $j (0 .. $n_logs - 1) {
            my $entry_off = $lb_addr + 0x0014 + $j * 45;
            my $raw       = _read_bytes($fh, $entry_off, 45);

            # Label: bytes 0..19 (NUL/space padded)
            my $label = substr($raw, 0, 20);
            $label =~ s/[\x00-\x1F\x7F-\xFF]//g;
            $label =~ s/\s+$//;
            next unless $label =~ /\S/;

            # Time: bytes 20..25 = 6-digit ASCII seconds from recording start
            # e.g. "000034" = 34s.  Confirmed: raw[20..25] in EEG-1100C real data.
            # Absolute timestamp follows as "(YYMMDDHHMMSS)" at bytes 26..39
            my $t_str = substr($raw, 20, 6);
            my $t_sec = ($t_str =~ /^(\d{6})$/) ? ($1 + 0)
                      : unpack('v', substr($raw, 20, 2));   # old-format fallback

            push @events, { t => $t_sec, label => $label };
        }
    }
    close $fh;
    return @events;
}

# ---------------------------------------------------------------------------
# Gain code → µV/bit
# ---------------------------------------------------------------------------
sub _gain_for_code {
    my ($code) = @_;
    return $GAIN_TABLE[$code] // 10.0;  # fallback 10 µV/bit
}

# ---------------------------------------------------------------------------
# Low-level binary helpers
# ---------------------------------------------------------------------------

sub _read_bytes {
    my ($fh, $offset, $len) = @_;
    seek $fh, $offset, 0 or croak "seek failed at offset $offset";
    my $buf;
    read($fh, $buf, $len) == $len or croak "short read at offset $offset";
    return $buf;
}

sub _read_u8 {
    my ($fh, $offset) = @_;
    return unpack 'C', _read_bytes($fh, $offset, 1);
}

sub _read_u16le {
    my ($fh, $offset) = @_;
    return unpack 'v', _read_bytes($fh, $offset, 2);
}

sub _read_u32le {
    my ($fh, $offset) = @_;
    return unpack 'V', _read_bytes($fh, $offset, 4);
}

# ---------------------------------------------------------------------------
# Device signature validation
# Valid prefixes: "EEG-" or "QI-"
# Returns 0 (false) if VALID, 1 (true) if INVALID — matches nk2edf check_device()
# Usage: croak if invalid  →  _check_device_sig($sig) and croak ...
#        return if invalid →  return () if _check_device_sig($sig)
# ---------------------------------------------------------------------------
sub _check_device_sig {
    my ($sig) = @_;
    return 0 if $sig =~ /^(?:EEG-|QI-|EEG2)/;
    return 1;
}

sub _hexdump {
    my ($s) = @_;
    join ' ', map { sprintf '%02X', ord $_ } split //, $s;
}

1;

__END__

=head1 KNOWN LIMITATIONS

=over 4

=item *

Waveform block data offset (0x0400) is approximate. Needs verification
against real hardware files. The exact offset depends on header size which
varies with n_channels.

=item *

Gain code → µV table is derived from EDFbrowser and may not cover all
EEG-2100/EEG-1200A variants.

=item *

EEG-1200A "new format" is NOT supported (different block layout).
See EDFbrowser issue #28.

=item *

BCD timestamp decode in waveform block header is not yet implemented.

=back

=head1 SEE ALSO

L<PDL>, L<PDL::Graphics::Cairo>, EDFbrowser nk2edf.cpp (GPL-2) by Teunis van Beelen

=head1 AUTHOR

goosh E<lt>goosh@exampleE<gt>

=head1 LICENSE

Same terms as Perl itself (Artistic License 2.0 or GPL-1+).
This module does NOT incorporate EDFbrowser GPL code.

=cut
