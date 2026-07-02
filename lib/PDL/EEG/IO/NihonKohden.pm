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

  my $rec = read_nk('subject.eeg');
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
our @EXPORT_OK = qw(read_nk nk_layout nk_format_hint);

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
    # 44-47: DC channels — NOTE the jack numbering differs by FORMAT/generation:
    #   EEG-1100C  -> DC03..DC06   EEG-1200A (JE-92NX) -> DC01..DC04
    # These defaults (DC03..DC06) are the 1100C convention and are only used
    # when a .21e is absent. The real .21e always wins, so live files get the
    # correct names (e.g. DC01..DC04 on EEG-1200A). Stimulus TTL triggers live
    # on these DC channels.
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
# On-disk layout is determined by the FORMAT SIGNATURE string at 0x0000,
# NOT by the physical recorder model (which is not stored in the file):
#   'wfmblock' : legacy layout, channel table inside the wfmblock
#                (EEG-1100x, EEG-2100, QI-403A, DAE-2100D, ...)
#   'extblock' : newer layout, channel info via ext_address(@0x3EE) chain
#                (EEG-1200A ...)
# ---------------------------------------------------------------------------
our %FORMAT_LAYOUT = (
    'EEG-1100A V01.00' => 'wfmblock', 'EEG-1100B V01.00' => 'wfmblock',
    'EEG-1100C V01.00' => 'wfmblock',
    'EEG-1100A V02.00' => 'wfmblock', 'EEG-1100B V02.00' => 'wfmblock',
    'EEG-1100C V02.00' => 'wfmblock',
    'QI-403A V01.00'   => 'wfmblock', 'QI-403A V02.00'   => 'wfmblock',
    'EEG-2100 V01.00'  => 'wfmblock', 'EEG-2100 V02.00'  => 'wfmblock',
    'DAE-2100D V01.30' => 'wfmblock', 'DAE-2100D V02.00' => 'wfmblock',
    'EEG-1200A V01.00' => 'extblock',
);

# Physical recorder model -> likely format signature (NON-authoritative HINT).
# The model is not in the file; a firmware update can change the emitted format
# (e.g. EEG-1290 -> 1200B/C, or EEG-1214 migrating 1100C -> 1200A). The real
# signature in the file always wins; this only helps a human guess beforehand.
our %DEVICE_FORMAT_HINT = (
    'EEG-1290' => 'EEG-1200A V01.00',   # extblock; may become 1200B/C after update
    'EEG-1214' => 'EEG-1100C V01.00',   # wfmblock; may move to 1200A after update
    'EEG-1200' => 'EEG-1200A V01.00',
    'EEG-1100' => 'EEG-1100C V01.00',
    'EEG-2100' => 'EEG-2100 V01.00',
);

# extblock ADC gains (µV/bit), hardware-fixed, not stored in file.
# EEG/micro channels use the ±3200 µV range; DC/other use the ±12002 range.
my $EXT_GAIN_UV = (3199.902 + 3200.0)  / (32767 + 32768);   # ~0.09765624
my $EXT_GAIN_MV = (12002.56 + 12002.9) / (32767 + 32768);   # ~0.36629984
# 1-based hardware channel index -> "micro" (µV) vs DC/other range
sub _ext_is_micro {
    my $c = shift;
    return 1 if $c >= 1 && $c <= 42;
    return 1 if $c == 75 || $c == 76;
    return 1 if $c >= 79 && $c <= 1096;
    return 0;
}

=head2 nk_layout($eeg_file)

Resolve the on-disk layout from the file's own format signature, with a safe
fallback so future/unseen signatures still read when possible.

  1. exact match in %FORMAT_LAYOUT            -> authoritative
  2. structural: ext_address(@0x3EE) != 0     -> extblock, else wfmblock
  3. name-family prefix used only to cross-check / annotate the decision
  4. signature not NK-like at all             -> undef (caller should stop)

Returns C<($sig, $layout, $how)>. C<$layout> is undef for unknown signatures.

=cut

sub nk_layout {
    my ($path) = @_;
    open my $fh, '<:raw', $path or croak "Cannot open $path: $!";
    (my $sig = _read_bytes($fh, 0x0000, 16)) =~ s/\x00.*//s;

    if (my $l = $FORMAT_LAYOUT{$sig}) { close $fh; return ($sig, $l, 'table'); }

    unless ($sig =~ /^(EEG-|QI-|DAE-)/) {
        close $fh; return ($sig, undef, 'unknown-signature');
    }

    my $ext = _read_u32le($fh, 0x03EE);
    close $fh;
    my $struct = $ext ? 'extblock' : 'wfmblock';
    my $name = ($sig =~ /^EEG-12/) ? 'extblock'
             : ($sig =~ /^EEG-11/) ? 'wfmblock' : undef;
    my $how = sprintf('fallback:ext_address(0x%X)', $ext);
    if    (defined $name && $name ne $struct) {
        $how .= " WARNING name-family suggests $name but structure says $struct";
    } elsif (defined $name) {
        $how .= " (name-family agrees: $name)";
    }
    return ($sig, $struct, $how);
}

=head2 nk_format_hint($model)

Non-authoritative guess of the format signature/layout for a physical recorder
model name (e.g. 'EEG-1290'). Returns C<($guessed_sig, $guessed_layout, $note)>.
Always confirm against the actual file via C<nk_layout()>.

=cut

sub nk_format_hint {
    my ($model) = @_;
    $model =~ s/^\s+|\s+$//g;
    my $sig = $DEVICE_FORMAT_HINT{$model};
    if (!$sig) {
        for my $k (sort { length($b) <=> length($a) } keys %DEVICE_FORMAT_HINT) {
            if (index($model, $k) == 0) { $sig = $DEVICE_FORMAT_HINT{$k}; last; }
        }
    }
    return (undef, undef, "no hint for model '$model'") unless $sig;
    return ($sig, $FORMAT_LAYOUT{$sig},
            "HINT ONLY for '$model' -> expect $sig ($FORMAT_LAYOUT{$sig}); "
          . "confirm via the file signature");
}

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

    # --- Layout dispatch (signature-based, with fallback) ------------------
    my (undef, $layout, $how) = nk_layout($eeg_path);
    croak "Unknown Nihon Kohden signature in $eeg_path" unless defined $layout;
    warn "read_nk: layout resolved via $how\n" if $how ne 'table' && $ENV{NK_DEBUG};
    return _read_extblock($eeg_path, %opts) if $layout eq 'extblock';
    # else: fall through to the legacy 'wfmblock' path below

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
# Internal: extblock layout reader (EEG-1200A and family)
#
# Confirmed from real data JJ0090J6.EEG (recorder EEG-1290, MMN, 38ch, 1000Hz;
# format signature "EEG-1200A V01.00"). Channel info lives in the extended
# block chain, not the wfmblock:
#
#   ext  = u32(@0x03EE)                 # != 0 for extblock
#   eb2  = u32(ext + 18)
#   eb3  = u32(eb2 + 20)                 # wfmblock-like subheader
#   n_ch = u16(eb3 + 68) + 1             # +1 appended STIM/marker channel
#   hw[i]= u16(eb3 + 72 + i*10) + 1      # 1-based hardware index (i=0..n_ch-2)
#   rec  = eb3 + 72 + (n_ch-1)*10        # waveform start
#
# Data: sample-interleaved uint16 LE, offset binary (center 0x8000). n_samples
# is not stored; computed from file size. Gain fixed (µV for EEG/micro indices,
# mV-range for DC). Last channel (STIM/marker) is raw (no offset).
# Same return contract as read_nk().  Ref: Brainstorm in_fopen_nk.m (independent
# re-implementation, no code copied).
# ---------------------------------------------------------------------------
sub _read_extblock {
    my ($eeg_path, %opts) = @_;
    my $fs_override = $opts{fs};
    my $want_events = !$opts{no_events};

    open my $fh, '<:raw', $eeg_path or croak "Cannot open $eeg_path: $!";
    (my $device = _read_bytes($fh, 0x0000, 16)) =~ s/\x00.*//s;

    my $ext = _read_u32le($fh, 0x03EE)
        or croak "extblock: ext_address is 0 in $eeg_path";
    my $ctl0      = _read_u32le($fh, 0x0092);
    my $data_addr = _read_u32le($fh, $ctl0 + 18);          # 0x17FE

    # sampling rate: lower 14 bits of u16 at data_addr+0x1A
    my $fs = $fs_override // (_read_u16le($fh, $data_addr + 0x1A) & 0x3FFF)
        or croak "extblock: sampling rate not determined (supply fs => NNN)";

    # extended block chain -> channels
    my $eb2 = _read_u32le($fh, $ext + 18);
    my $eb3 = _read_u32le($fh, $eb2 + 20);
    my $n_ch = _read_u16le($fh, $eb3 + 68) + 1;            # +1 STIM
    my @hw;                                                # 1-based hw indices
    push @hw, _read_u16le($fh, $eb3 + 72 + $_ * 10) + 1 for 0 .. $n_ch - 2;
    my $rec = $eb3 + 72 + ($n_ch - 1) * 10;

    # n_samples from file size
    seek $fh, 0, 2; my $file_size = tell $fh;
    my $n_samp = int(($file_size - $rec) / $n_ch / 2);
    croak "extblock: no samples (rec=$rec, size=$file_size)" if $n_samp <= 0;

    # --- timestamp (BCD at data_addr+0x14) ---
    my $t_start = sprintf("20%02d-%02d-%02d %02d:%02d:%02d",
        map { _bcd_byte($fh, $data_addr + 0x14 + $_) } 0 .. 5);

    # --- read interleaved uint16 straight into a ushort piddle [n_ch, n_samp] ---
    seek $fh, $rec, 0;
    my $buf; my $want = $n_ch * $n_samp * 2;
    my $got = read($fh, $buf, $want);
    $n_samp = int($got / $n_ch / 2) if $got < $want;       # tolerate short file
    my $u16 = zeroes(ushort, $n_ch, $n_samp);              # (ch, t)
    ${ $u16->get_dataref } = substr($buf, 0, $n_ch * $n_samp * 2);
    $u16->upd_data;

    # --- per-channel gain (µV/bit) + digital offset ---
    my (@gain_uv, @offset);
    for my $c (@hw) {
        push @gain_uv, _ext_is_micro($c) ? $EXT_GAIN_UV : $EXT_GAIN_MV;
        push @offset,  32768;
    }
    push @gain_uv, 1.0; push @offset, 0;                   # STIM: raw code
    my $gains = PDL->new(\@gain_uv);
    my $offs  = PDL->new(\@offset);

    my $data_uv = (($u16->double - $offs->slice(':,*1')) * $gains->slice(':,*1'))->float;

    # --- labels from .21e (fall back to DEFAULT_LABELS), STIM appended ---
    my $base = $eeg_path; $base =~ s/\.[^.]+$//;
    my %ov;
    for my $ext21 (qw(.21e .21E)) {
        my %h = _read_21e("$base$ext21");
        if (%h) { %ov = %h; last }
    }
    my @labels;
    for my $c (@hw) {
        my $li = $c - 1;                                   # 0-indexed lookup
        push @labels, $ov{$li} // ($li >= 0 ? $DEFAULT_LABELS[$li] : undef) // "CH$li";
    }
    push @labels, 'STIM';                                  # appended marker channel

    # --- events (.LOG) via the shared helper ---
    my @events;
    if ($want_events) {
        my $lp = "$base.LOG"; $lp = "$base.log" unless -f $lp;
        @events = _read_log($lp) if -f $lp;
    }

    close $fh;

    return {
        data             => $data_uv,          # [n_ch, n_samp] float32 µV
        fs               => $fs,
        labels           => \@labels,
        t_start          => $t_start,
        events           => \@events,
        gains            => $gains,            # [n_ch] µV/bit (STIM=1)
        n_blocks         => 1,
        n_ch_valid       => $n_ch - 1,         # analog channels (excl. STIM)
        block_idx        => 0,
        device           => $device,
        layout           => 'extblock',
        ch_hw_idx        => \@hw,              # 1-based hardware indices
        stim_index       => $n_ch,            # 1-based; last row
        t_block_starts   => [0],
        n_samp_per_block => [$n_samp],
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
# Format: sections [ELECTRODE], [REFERENCE], ... key=value pairs.
# Keys are 4-digit decimal electrode indices, 0-based ("0000"=Fp1, ...).
# [ELECTRODE] holds normal electrode names; [REFERENCE] holds reference-input
# channels with a '$' prefix (e.g. 0076=$A1, 0077=$A2, 0080=$Cz). Some montages
# (e.g. EEG-1290/JE-92NX) record reference channels whose indices are absent
# from [ELECTRODE] but present in [REFERENCE]; we fill those in as a fallback.
# [ELECTRODE] always wins on key collision.
# ---------------------------------------------------------------------------
sub _read_21e {
    my ($path) = @_;
    return () unless -f $path;
    my $fh;
    open($fh, '<:encoding(Shift_JIS):crlf', $path)
        or open($fh, '<:encoding(UTF-8):crlf', $path)
        or return ();
    my (%elec, %ref);
    my $section = '';
    while (<$fh>) {
        s/\r\n$/\n/;   # CRLF → LF
        s/\r$/\n/;     # CR-only → LF
        chomp;
        if (/^\[([^\]]+)\]/) { $section = uc $1; next }
        next unless $section eq 'ELECTRODE' || $section eq 'REFERENCE';
        # Keys are 4-digit decimal, 0-indexed: "0000"=Fp1, "0001"=Fp2, ...
        if (/^(\d+)=(.+)$/) {
            my ($key, $val) = ($1 + 0, $2);
            $val =~ s/\s+$//;   # trim trailing whitespace
            if ($section eq 'ELECTRODE') { $elec{$key} = $val }
            else {
                # NK marks reference derivations with a leading '$' (e.g. $A1,
                # $AV, $Cz). '$' is hazardous in Perl/filenames and, stripped
                # bare, would collide with the plain electrode (A1). Normalize to
                # a safe, self-documenting, collision-free "<name>_ref" suffix.
                $val =~ s/^\$(.+)/${1}_ref/;
                $ref{$key} = $val;
            }
        }
    }
    close $fh;
    # [REFERENCE] fills only holes that are ALSO blank ('-') in the built-in
    # defaults, so a meaningful DEFAULT name is never hidden by a $-reference
    # name (preserves pre-existing wfmblock/1100C behavior). [ELECTRODE] always
    # wins. The 1200A reference channels (indices 76/77, whose defaults are '-')
    # still pick up $A1/$A2.
    for my $k (keys %ref) {
        next if defined $elec{$k};
        my $def = $DEFAULT_LABELS[$k];
        $elec{$k} = $ref{$k} if !defined $def || $def eq '-';
    }
    return %elec;
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

Two on-disk layouts are supported, dispatched by the file's format signature
via C<nk_layout()>: the legacy 'wfmblock' layout (EEG-1100x / EEG-2100 / QI-403A
signatures) and the newer 'extblock' layout (EEG-1200A signature; used e.g. by
the EEG-1290 recorder). Unseen signatures fall back to a structural check
(ext_address at 0x03EE) so future variants (EEG-1200B/C, updated firmware) read
without code changes; truly non-NK signatures are rejected.

=item *

In the extblock layout, the appended last channel is the STIM/marker channel
(returned raw). Experiment triggers on EEG-1290 recordings are TTL levels on the
DC channels (per the montage's .21e labels), not on this STIM channel.

=back

=head1 SEE ALSO

L<PDL>, L<PDL::Graphics::Cairo>, EDFbrowser nk2edf.cpp (GPL-2) by Teunis van Beelen

=head1 AUTHOR

goosh E<lt>goosh@exampleE<gt>

=head1 LICENSE

Same terms as Perl itself (Artistic License 2.0 or GPL-1+).
This module does NOT incorporate EDFbrowser GPL code.

=cut
