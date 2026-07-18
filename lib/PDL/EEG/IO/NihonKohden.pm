package PDL::EEG::IO::NihonKohden;

use strict;
use warnings;
use Carp qw(croak confess carp);
use Encode qw(decode);
use POSIX ();
use PDL;

our $VERSION = '0.02';

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

  # Whole recording (all waveform blocks concatenated, gaps marked)
  my $all = read_nk('subject.eeg', all_blocks => 1);
  # $all->{gap_bounds}  -> [[start_samp, end_samp], ...] (recording breaks)

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
our @EXPORT_OK = qw(read_nk nk_layout nk_format_hint block_extents block_ranges select_block select_range clock_to_samp);

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
# 0-indexed, matching .21e file format (confirmed from EEG-1100C subject.21E).
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
    'EEG-1200C V01.00' => 'extblock',   # seen in the wild (4 kHz, 18 segments);
                                        # was already reading via the ext_address
                                        # fallback, but make it explicit
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

# extblock ADC gains, hardware-fixed, not stored in the file.
#
# The two ranges are quoted by the vendor in DIFFERENT units:
#   EEG / "micro" channels : +/- 3200    uV
#   DC / other channels    : +/- 12002.9 mV   -- i.e. +/- 12 V, not +/- 12 mV
#
# An EEG-1200A DC input is a +/-12 V line. A saturated DC channel reads exactly
# -12002.9 on the vendor's scale, which is the rail; a 3.3-5 V trigger pulse sits
# around 3.6-4.2 on it. Earlier versions of this reader used the mV figure as
# though it were uV/bit, so every DC channel came out 1000x too small and a 4 V
# TTL looked like a 4 mV wobble.
#
# {data} and {gains} are now uniformly MICROVOLTS, for every channel. The DC gain
# is the vendor's mV/bit x 1000.
#
# {units} is a separate thing: it is the dimension each channel should be WRITTEN
# in when exported. A DC channel cannot go into EDF in uV -- +/-12002900 uV needs
# nine characters and EDF's physical_min field is eight -- so it is exported in
# mV. write_edf() reads {units} and does that conversion; nothing else should.
my $EXT_GAIN_UV = (3199.902 + 3200.0)  / (32767 + 32768);          # ~0.0977 uV/bit
my $EXT_GAIN_MV = (12002.56 + 12002.9) / (32767 + 32768) * 1000;   # ~366.3  uV/bit
# ---------------------------------------------------------------------------
# DC TRIGGER CHANNEL NUMBERING
#
# The four trigger DC inputs are hardware electrode codes 45-48. What they are
# CALLED on the front panel is not the same on every recorder:
#
#   EEG-1100A/B/C : DC03 DC04 DC05 DC06     (what @DEFAULT_LABELS says)
#   EEG-1200A     : DC01 DC02 DC03 DC04
#
# and the physical recorder model is NOT IN THE FILE -- see %FORMAT_LAYOUT above:
# the 16-byte string at 0x0000 is a FORMAT signature, and the same recorder can
# emit different signatures after a firmware update. So the best we have is the
# format signature, and the mapping below is keyed on that.
#
# Signatures whose DC numbering has NOT been confirmed against a real recording
# are deliberately ABSENT. We do not guess: a wrong DC label silently attaches
# the stimulus triggers to the wrong channel, and that error survives all the way
# into an analysis. read_nk() croaks and asks for the .21e, or for an explicit
# dc_base.
#
# A .21e names the channels outright and always wins. When one is present AND it
# names the DC codes, we cross-check it against this table and carp on
# disagreement -- which is how an unconfirmed signature gets confirmed, or the
# table gets corrected, the first time a real file of that kind is read.
#
# The rule is by FAMILY, and each half of it has a source:
#
#   EEG-1100*  -> 3   This is just @DEFAULT_LABELS, the vendor electrode-code
#                     table out of nk2edf.cpp, which names code 45 "DC03". Not a
#                     guess -- the documented table. Confirmed in practice on
#                     EEG-1100C.
#
#   EEG-1200*  -> 1   The 1200 line DEVIATES from that table and starts at DC01.
#                     Confirmed on real recordings from two different signatures,
#                     EEG-1200A and EEG-1200C, whose own .21e files name codes
#                     45-48 DC01..DC04.
#
# Anything else -- EEG-2100, QI-403A, DAE-2100D, or a signature that does not
# exist yet -- follows NEITHER convention as far as we know, so read_nk() croaks
# rather than pick one. A wrong DC label silently attaches the stimulus triggers
# to the wrong channel, and that error survives all the way into an analysis.
#
# Beware "EEG2100": it is ALSO the directory the vendor software exports into
# (NKT/EEG2100/...), and a file in that folder is routinely signed
# "EEG-1200A V01.00" -- the recording this was all debugged on is exactly that.
# The directory names the product line; the 16 bytes at 0x0000 name the FORMAT.
# They are different things and they disagree.
my @DC_FAMILY = (
    [ qr/^EEG-1200/ => 1 ],
    [ qr/^EEG-1100/ => 3 ],
);
my @DC_CODES = (45, 46, 47, 48);

sub _dc_labels {
    my ($sig, $dc_base) = @_;
    my ($fmt) = ($sig // '') =~ /^(\S+)/;               # "EEG-1200C V01.00" -> "EEG-1200C"
    my $base = $dc_base;
    if (!defined $base && defined $fmt) {
        for my $f (@DC_FAMILY) {
            if ($fmt =~ $f->[0]) { $base = $f->[1]; last }
        }
    }
    return undef unless defined $base;
    return { map { $DC_CODES[$_] => sprintf('DC%02d', $base + $_) } 0 .. $#DC_CODES };
}

# Refuse to guess when it matters: this file HAS trigger DC channels, we have no
# numbering for its signature, and nothing else names them.
sub _check_dc_known {
    my ($sig, $dc, $hw, $ov) = @_;
    return if $dc;
    my %code = map { $_ => 1 } @DC_CODES;
    my @unnamed = grep { $code{$_} && !defined $ov->{ $_ - 1 } } @$hw;
    return unless @unnamed;
    croak
        "DC channel numbering is unknown for format signature '" . ($sig // '?') . "'.\n"
      . "  Hardware codes 45-48 are the four trigger DC inputs, but their panel\n"
      . "  names differ by family:  EEG-1100* -> DC03..DC06\n"
      . "                           EEG-1200* -> DC01..DC04\n"
      . "  and this signature is in neither family.\n"
      . "  and this file has no .21e naming them. Refusing to guess: a wrong DC\n"
      . "  label puts the stimulus triggers on the wrong channel.\n"
      . "  Either put the .21e next to the .EEG, or say which panel it is:\n"
      . "      read_nk(\$file, dc_base => 1)   # DC01..DC04  (EEG-1200A)\n"
      . "      read_nk(\$file, dc_base => 3)   # DC03..DC06  (EEG-1100x)\n";
}

# When the .21e DOES name the DC codes, check our table against it. This is the
# only way %DC_BASE ever gets validated, so it is worth the four comparisons.
sub _verify_dc_base {
    my ($sig, $dc, $hw, $ov) = @_;
    return unless $dc;
    my %code = map { $_ => 1 } @DC_CODES;
    my @bad;
    for my $c (grep { $code{$_} } @$hw) {
        my $named = $ov->{ $c - 1 } or next;             # .21e did not name it
        next unless $named =~ /^DC\d+$/i;                # a custom name proves nothing
        push @bad, "code $c: .21e says $named, table says $dc->{$c}"
            if uc($named) ne uc($dc->{$c});
    }
    carp "DC numbering for '" . ($sig // '?') . "' disagrees with the .21e:\n"
       . join('', map { "    $_\n" } @bad)
       . "  The .21e wins (it is the recording's own montage), but %DC_BASE in "
       . __PACKAGE__ . " is wrong for this signature -- please fix it."
        if @bad;
}

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
  block       => $n   # which waveform block to read (0-based, default 0)
  all_blocks  => 1    # concatenate ALL waveform blocks into one recording
                      #   (wfmblock only; extblock is already a single block).
                      #   Blocks are butt-joined: NO synthetic samples are
                      #   inserted at a recording break. Break positions are in
                      #   t_block_starts / block_meta, which also carry each
                      #   block's wall-clock t_start, so the real elapsed gap is
                      #   recoverable. The data are NOT padded to elapsed time.
  dc_base     => 1|3  # what hardware code 45 is CALLED on the front panel:
                      #   1 -> DC01..DC04  (EEG-1200A)
                      #   3 -> DC03..DC06  (EEG-1100A/B/C)
                      # Only needed when the file's format signature is not in
                      # %DC_BASE and there is no .21e naming the DC channels; in
                      # that case read_nk croaks rather than guess, because a
                      # wrong DC label puts the triggers on the wrong channel.
                      # A .21e always wins over this.
  gap_samples => $n   # DEPRECATED. Insert $n zero samples at each recording
                      #   break (and report them in gap_bounds). Default 0.
                      #   Zeros are not data: they ring through filters, corrupt
                      #   spectra and distort waveform metrics. Only set this if
                      #   you have code that still keys off gap_bounds.
  no_events   => 1    # skip event parsing even if .log exists
  label_map   => \%h  # override channel labels, keyed by 1-based ch_idx
                      #   (electrode index in the file). Highest priority, above
                      #   the .21e and the built-in DEFAULT_LABELS. Use this to
                      #   name trigger/DC channels correctly per headbox/montage,
                      #   e.g. label_map => { 45=>'DC03', 74=>'DC06' }. The DC/
                      #   trigger channel<->name mapping is headbox-specific and
                      #   NOT reliably derivable from the format; supply it here
                      #   (see PDL::EEG::Trigger for signal+.PTN based resolution).

Returns a hashref:
  data           => $pdl   # [n_ch, n_samples] float32, µV
  fs             => $hz
  labels         => \@channel_label_strings
  ch_indices     => \@idx  # 1-based ch_idx (electrode index) per label (both layouts)
  t_start        => "YYYY-MM-DD HH:MM:SS"  (first block)
  events         => \@{ {t => $sec, label => $str} }
  gains          => $pdl   # [n_ch] µV/bit
  n_blocks       => total waveform blocks in file
  gap_bounds     => \@[ [start_samp, end_samp], ... ]  # concat coords (empty if none)
  t_block_starts => \@start_samp                       # concat start of each block
  block_meta     => \@{ start_samp, n_samp, t_start }  # per block

=cut

sub read_nk {
    my ($eeg_path, %opts) = @_;
    croak "File not found: $eeg_path" unless -f $eeg_path;

    my $block_idx   = $opts{block}      // 0;
    my $want_events = !$opts{no_events};
    my $fs_override = $opts{fs};          # caller can supply fs (e.g. 1000)
    my $all_blocks  = $opts{all_blocks} ? 1 : 0;
    my $gap_samples = $opts{gap_samples} // 0;

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
        if !$all_blocks && $block_idx >= @wfm_addrs;

    seek $fh, 0, 2;
    my $file_size = tell $fh;

    # --- 3-9. Read waveform sample data ------------------------------------
    # Data encoding (per block): uint16 LE, offset binary (center 0x8000);
    # µV = (raw - 0x8000) * 0.09765625.  n_samples is derived from the gap to
    # the next block address (or EOF for the last block).
    my ($data_uv, $n_ch, $n_ch_valid, $fs, $meta);
    my (@gap_bounds, @t_block_starts, @block_meta);

    if ($all_blocks) {
        # Read every block, then concatenate with a short zero gap marker so
        # recording breaks stay visible without padding to real elapsed time.
        my @blk;
        for my $b (0 .. $#wfm_addrs) {
            my $next = ($b + 1 <= $#wfm_addrs) ? $wfm_addrs[$b+1] : $file_size;
            my $r = _read_wfm_block($fh, $wfm_addrs[$b], $next, $fs_override);
            $meta       //= $r->{meta};
            $n_ch       //= $r->{n_ch};
            $n_ch_valid //= $r->{n_ch_valid};
            $fs         //= $r->{fs};
            croak "all_blocks: channel count differs at block $b ($r->{n_ch} vs $n_ch)"
                if $r->{n_ch} != $n_ch;
            push @blk, $r;
        }
        my $total = 0; $total += $_->{n_samp} for @blk;
        $total += $gap_samples * $#blk if @blk > 1;
        $data_uv = zeroes(float, $n_ch, $total);       # gap regions stay 0.0
        my $pos = 0;
        for my $b (0 .. $#blk) {
            my $ns = $blk[$b]{n_samp};
            $data_uv->slice(",${pos}:${\($pos+$ns-1)}") .= $blk[$b]{data} if $ns > 0;
            push @t_block_starts, $pos;
            push @block_meta,
                { start_samp => $pos, n_samp => $ns, t_start => $blk[$b]{meta}{t_start} };
            $pos += $ns;
            if ($b < $#blk && $gap_samples > 0) {      # optional gap marker
                push @gap_bounds, [ $pos, $pos + $gap_samples - 1 ];
                $pos += $gap_samples;
            }
        }
    } else {
        my $next = ($block_idx + 1 < @wfm_addrs) ? $wfm_addrs[$block_idx+1] : $file_size;
        my $r = _read_wfm_block($fh, $wfm_addrs[$block_idx], $next, $fs_override);
        ($data_uv, $n_ch, $n_ch_valid, $fs, $meta) =
            ($r->{data}, $r->{n_ch}, $r->{n_ch_valid}, $r->{fs}, $r->{meta});
        push @t_block_starts, 0;
        push @block_meta,
            { start_samp => 0, n_samp => $r->{n_samp}, t_start => $meta->{t_start} };
    }

    # --- Electrode labels (block-independent) ------------------------------
    my $base = $eeg_path;  $base =~ s/\.[^.]+$//;
    my %label_override;
    for my $ext (qw(.21e .21E)) {
        my %h = _read_21e("$base$ext");
        if (%h) { %label_override = %h; last }
    }

    # Valid channels = ch_indices[0..n_ch_valid-1]; last ch is zero-pad
    my @hw_idx = @{ $meta->{ch_indices} }[ 0 .. $n_ch_valid - 1 ];
    my $dc = _dc_labels($device, $opts{dc_base});
    _check_dc_known($device, $dc, \@hw_idx, \%label_override);
    _verify_dc_base($device, $dc, \@hw_idx, \%label_override);

    my @labels;
    for my $i (0 .. $n_ch_valid - 1) {
        my $idx       = $meta->{ch_indices}[$i];  # 1-indexed
        my $label_idx = $idx - 1;                 # 0-indexed for lookup
        my $user      = $opts{label_map} ? $opts{label_map}{$idx} : undef;  # ch_idx key
        push @labels, $user
                   // $label_override{$label_idx}
                   // ($dc ? $dc->{$idx} : undef)   # this signature's DC panel
                   // ($label_idx >= 0 ? $DEFAULT_LABELS[$label_idx] : undef)
                   // "CH$i";
    }
    push @labels, 'PAD';   # trailing zero-pad channel

    # --- Gain report (data already in µV) ----------------------------------
    # EEG-1100C: ADC gain fixed, not stored in file. 0.9765625 µV/div / 10 =
    # 0.09765625 µV/bit (0.09765625 * 32767 = 3199.90 µV = EDF phys_max).
    my $gains = PDL->new([ (0.09765625) x $n_ch ]);

    # --- Events (optional) -------------------------------------------------
    my @events;
    if ($want_events) {
        my $log_path = "$base.LOG";
        $log_path = "$base.log" unless -f $log_path;
        @events = _read_log($log_path) if -f $log_path;
    }
    # Place events at their data-sample position, but only for a concatenated
    # multi-block session (all_blocks): there the .LOG's REC START markers
    # delimit the blocks and block_meta carries each block's exact start_samp,
    # so _attach_recstart_samp is exact (same as the extblock path; within a
    # block wall-clock and data advance 1:1). It self-guards on REC START count
    # == block count and is a no-op otherwise. For a single-block read the .LOG
    # is session-wide while the data is one block, so we do NOT place events
    # (they keep {t,label,epoch} only, as before) rather than guess.
    _attach_recstart_samp(\@events, \@block_meta, $fs) if @block_meta > 1;

    close $fh;

    return {
        data           => $data_uv,
        fs             => $fs,
        labels         => \@labels,
        t_start        => $meta->{t_start},
        events         => \@events,
        gains          => $gains,
        n_blocks       => scalar(@wfm_addrs),
        n_ch_valid     => $n_ch_valid,
        ch_indices     => [ @{ $meta->{ch_indices} }[0 .. $n_ch_valid - 1] ], # 1-based ch_idx per label

        block_idx      => ($all_blocks ? -1 : $block_idx),
        all_blocks     => $all_blocks,
        device         => $device,           # e.g. "EEG-1100C V01.00"
        layout         => 'wfmblock',
        system_reference => $label_override{system_reference},  # e.g. "C3,C4" (.21e [SYSTEM_SETUP])
        last_pattern     => $label_override{last_pattern},      # recording montage no. (.21e [LASTPATTERN])
        gap_bounds     => \@gap_bounds,       # [[start_samp,end_samp],...] concat coords
        t_block_starts => \@t_block_starts,   # concat start sample of each block
        block_meta     => \@block_meta,       # per-block {start_samp,n_samp,t_start}
        n_samp_per_block => [ map { $_->{n_samp} } @block_meta ],
        units          => [ ('uV') x $n_ch ],            # wfmblock is uV throughout
    };
}

# ---------------------------------------------------------------------------
# Internal: read one wfmblock's samples -> µV piddle [n_ch, n_samp] + metadata.
# Used by both the single-block (block => N) and all_blocks paths.
# ---------------------------------------------------------------------------
sub _read_wfm_block {
    my ($fh, $waddr, $next_addr, $fs_override) = @_;
    my $meta       = _read_wfm_header($fh, $waddr);
    my $data_rel   = $meta->{data_offset} - $waddr;      # = 0x171
    my $data_bytes = ($next_addr - $waddr) - $data_rel;
    my $n_ch       = $meta->{n_ch};                      # n_ch_entries + 1
    croak "Cannot determine n_samples (data_bytes=$data_bytes not divisible by "
        . "n_ch*2=" . ($n_ch*2) . ")" if $data_bytes % ($n_ch * 2) != 0;
    my $n_samp = int($data_bytes / ($n_ch * 2));

    my $fs = $fs_override // $meta->{fs}
        or croak "Sampling rate could not be determined "
               . "(header gave " . ($meta->{fs}//0) . " Hz; supply fs => NNN option)";

    seek $fh, $meta->{data_offset}, 0 or croak "Seek to data failed";
    my $buf;
    my $bytes = $n_ch * $n_samp * 2;
    my $got = read($fh, $buf, $bytes);
    croak "Short read: wanted $bytes bytes, got $got" unless $got == $bytes;

    # Interpret the little-endian uint16 buffer directly into a PDL, the same
    # way _read_extblock does: no multi-million-element Perl list, no double
    # blow-up. ~4.5x faster and far lighter on memory at full recording scale.
    # (On little-endian hosts ushort matches "v"; both reader paths assume LE.)
    my $u16 = zeroes(ushort, $n_ch, $n_samp);
    ${ $u16->get_dataref } = substr($buf, 0, $n_ch * $n_samp * 2);
    $u16->upd_data;
    my $raw = $u16->double - 0x8000;                     # offset binary → signed
    $raw->slice('(-1),:') .= 0.0;                        # PAD channel zero-fill
    my $data = ($raw * 0.09765625)->float;               # fixed µV/bit

    return {
        data       => $data,
        n_samp     => $n_samp,
        n_ch       => $n_ch,
        n_ch_valid => $meta->{n_ch_valid},
        meta       => $meta,
        fs         => $fs,
    };
}

# ---------------------------------------------------------------------------
# Internal: extblock layout reader (EEG-1200A and family)
#
# Confirmed from real data subject.EEG (recorder EEG-1290, MMN, 38ch, 1000Hz;
# format signature "EEG-1200A V01.00"). Channel info lives in the extended
# block chain, not the wfmblock. Data: sample-interleaved uint16 LE, offset
# binary (center 0x8000). n_samples computed from file size. Gain fixed (µV for
# EEG/micro indices, mV-range for DC). Last channel (STIM/marker) is raw.
# One data block per file, so all_blocks is a no-op here.
# Same return contract as read_nk().  Ref: Brainstorm in_fopen_nk.m.
# ---------------------------------------------------------------------------
sub _read_extblock {
    my ($eeg_path, %opts) = @_;
    my $fs_override = $opts{fs};
    my $want_events = !$opts{no_events};
    my $block_idx   = $opts{block};
    my $all_blocks  = $opts{all_blocks} ? 1 : 0;

    open my $fh, '<:raw', $eeg_path or croak "Cannot open $eeg_path: $!";
    (my $device = _read_bytes($fh, 0x0000, 16)) =~ s/\x00.*//s;

    my $ext = _read_u32le($fh, 0x03EE)
        or croak "extblock: ext_address is 0 in $eeg_path";
    my $ctl0      = _read_u32le($fh, 0x0092);
    my $data_addr = _read_u32le($fh, $ctl0 + 18);          # 0x17FE

    my $fs = $fs_override // (_read_u16le($fh, $data_addr + 0x1A) & 0x3FFF)
        or croak "extblock: sampling rate not determined (supply fs => NNN)";

    my $eb2  = _read_u32le($fh, $ext + 18);
    my $eb3  = _read_u32le($fh, $eb2 + 20);
    my $n_ch = _read_u16le($fh, $eb3 + 68) + 1;            # +1 STIM
    my @hw;
    push @hw, _read_u16le($fh, $eb3 + 72 + $_ * 10) + 1 for 0 .. $n_ch - 2;

    my $hdr_len = 72 + ($n_ch - 1) * 10;   # size of the channel-info block (e.g. 442)
    my $stride  = $n_ch * 2;               # bytes per sample (e.g. 76)
    my $rec     = $eb3 + $hdr_len;         # first sample of the FIRST segment

    seek $fh, 0, 2; my $file_size = tell $fh;
    croak "extblock: no samples (rec=$rec, size=$file_size)" if $file_size <= $rec;

    # -----------------------------------------------------------------------
    # SEGMENTS.
    #
    # An extblock file is NOT one contiguous data block. At every recording
    # break the recorder writes a fresh copy of the channel-info block --
    # 72 + (n_ch-1)*10 bytes, byte-identical to the one at eb3 except for its
    # timestamp -- straight into the sample stream, and then carries on.
    #
    # Earlier versions of this reader assumed a single block running to EOF, so
    # they read each embedded header as if it were 5.8 samples of EEG. The
    # channel phase then slipped by hdr_len % stride bytes (442 % 76 = 62 bytes
    # = 31 channels) at every break, and from the first break onward every
    # channel label sat on another channel's data. It also over-reported the
    # sample count by hdr_len/stride per break.
    #
    # The header is unmistakable and is always sample-aligned:
    #     +0x00        0x01                       block type
    #     +0x01..+0x04 "TIME"
    #     +0x12,+0x13  0x02 0x02                  version
    #     +0x14..      "YYYYMMDDHHMMSS" (ASCII)   this segment's start time
    #     +0x44        u16 = n_ch - 1             channel count, must agree
    #     +0x48..      the channel table again
    # -----------------------------------------------------------------------
    my @seg = _ext_segments($fh, $file_size, $rec, $hdr_len, $stride, $n_ch,
                           _ext_hdr_time($fh, $eb3) // _bcd_time($fh, $data_addr));
    croak "extblock: no data segments found" unless @seg;

    my $n_blocks = scalar @seg;
    croak "extblock: block index $block_idx out of range (have $n_blocks)"
        if defined $block_idx && ($block_idx < 0 || $block_idx >= $n_blocks);

    # which segments to read
    my @want = (defined $block_idx) ? ($seg[$block_idx])
             : ($n_blocks == 1 || $all_blocks) ? @seg
             : ($seg[0]);                          # default: first segment only

    # -----------------------------------------------------------------------
    # Read the wanted segments and concatenate (butt-joined; gap_samples is
    # honoured for callers that still want the old zero marker).
    # -----------------------------------------------------------------------
    my $gap_samples = $opts{gap_samples} // 0;
    my $n_samp = 0; $n_samp += $_->{n_samp} for @want;
    $n_samp   += $gap_samples * $#want if @want > 1;
    croak "extblock: no samples" if $n_samp <= 0;

    my $u16 = zeroes(ushort, $n_ch, $n_samp);
    my $dst = $u16->get_dataref;
    my $pos = 0;
    my (@t_block_starts, @block_meta, @gap_bounds, @n_per);
    for my $i (0 .. $#want) {
        my $s = $want[$i];
        seek $fh, $s->{off}, 0;
        my $buf;
        my $want_b = $s->{n_samp} * $stride;
        my $got = read($fh, $buf, $want_b);
        croak "extblock: short read in segment $i" unless $got == $want_b;
        substr($$dst, $pos * $stride, $want_b) = $buf;
        push @t_block_starts, $pos;
        push @n_per, $s->{n_samp};
        push @block_meta, { start_samp => $pos, n_samp => $s->{n_samp},
                            t_start => $s->{t_start} };
        $pos += $s->{n_samp};
        if ($i < $#want && $gap_samples > 0) {
            push @gap_bounds, [ $pos, $pos + $gap_samples - 1 ];
            $pos += $gap_samples;
        }
    }
    $u16->upd_data;

    my $t_start = $want[0]{t_start};

    my (@gain_uv, @offset, @units);
    for my $c (@hw) {
        my $micro = _ext_is_micro($c);
        push @gain_uv, $micro ? $EXT_GAIN_UV : $EXT_GAIN_MV;   # both uV/bit
        push @units,   $micro ? 'uV' : 'mV';                   # EXPORT dimension
        push @offset,  32768;
    }
    push @gain_uv, 1.0; push @offset, 0; push @units, 'code';   # STIM: raw code
    my $gains = PDL->new(\@gain_uv);
    my $offs  = PDL->new(\@offset);

    my $data_uv = (($u16->double - $offs->slice(':,*1')) * $gains->slice(':,*1'))->float;

    my $base = $eeg_path; $base =~ s/\.[^.]+$//;
    my %ov;
    for my $ext21 (qw(.21e .21E)) {
        my %h = _read_21e("$base$ext21");
        if (%h) { %ov = %h; last }
    }
    # DC numbering keys on the FORMAT SIGNATURE, not on the layout. extblock does
    # not imply EEG-1200A, and the recorder model is not in the file at all.
    my $dc = _dc_labels($device, $opts{dc_base});
    _check_dc_known($device, $dc, \@hw, \%ov);
    _verify_dc_base($device, $dc, \@hw, \%ov);

    my @labels;
    for my $c (@hw) {
        my $li = $c - 1;
        my $user = $opts{label_map} ? $opts{label_map}{$c} : undef;   # ch_idx (hw index) key
        push @labels, $user
                   // $ov{$li}                                  # .21e override
                   // ($dc ? $dc->{$c} : undef)                 # model's DC panel
                   // ($li >= 0 ? $DEFAULT_LABELS[$li] : undef)
                   // "CH$li";
    }
    push @labels, 'STIM';

    my @events;
    if ($want_events) {
        my $lp = "$base.LOG"; $lp = "$base.log" unless -f $lp;
        @events = _read_log($lp) if -f $lp;
    }
    # Segment boundaries and their wall-clock starts are now known exactly, from
    # the embedded block headers. Prefer _attach_recstart_samp: it delimits
    # segments by the .LOG REC START markers, anchors each to the exact header
    # boundary (block_meta start_samp), and places within a segment by the .LOG
    # time offset from that segment's REC START -- exact, because within a
    # segment wall-clock and data advance 1:1. Fall back to the uniform-epoch
    # page map, then to the wall-clock seg map. (_attach_seg_samp assumes .LOG
    # time equals the block-header wall-clock; on real recordings the .LOG clock
    # also counts paused setup time between blocks, so that assumption drifts by
    # tens to hundreds of seconds and misfiles late-segment events -- see
    # t/08_epoch.t.)
    _attach_recstart_samp(\@events, \@block_meta, $fs)
        or _attach_epoch_samp(\@events, $n_samp, $fs)
        or _attach_seg_samp(\@events, \@block_meta, $fs);

    close $fh;

    return {
        data             => $data_uv,          # [n_ch, n_samp] float32 µV
        fs               => $fs,
        labels           => \@labels,
        t_start          => $t_start,
        events           => \@events,
        gains            => $gains,            # [n_ch] uV/bit, ALL channels
        units            => \@units,           # [n_ch] EXPORT dimension:
                                               #   'uV' | 'mV' (DC) | 'code' (STIM).
                                               #   {data} is uV regardless.
        n_blocks         => $n_blocks,
        n_ch_valid       => $n_ch - 1,         # analog channels (excl. STIM)
        block_idx        => (defined $block_idx ? $block_idx : ($all_blocks ? -1 : 0)),
        all_blocks       => $all_blocks,
        device           => $device,
        layout           => 'extblock',
        system_reference => $ov{system_reference},   # e.g. "C3,C4" (.21e [SYSTEM_SETUP])
        last_pattern     => $ov{last_pattern},        # recording montage no. (.21e [LASTPATTERN])
        ch_hw_idx        => \@hw,              # 1-based hardware indices
        ch_indices       => \@hw,             # alias: 1-based ch_idx per label (matches wfmblock)
        stim_index       => $n_ch,            # 1-based; last row
        gap_bounds       => \@gap_bounds,
        t_block_starts   => \@t_block_starts,
        block_meta       => \@block_meta,
        n_samp_per_block => \@n_per,
        hdr_len          => $hdr_len,         # embedded channel-info block size
    };
}

# ---------------------------------------------------------------------------
# _ext_segments($fh, $file_size, $rec, $hdr_len, $stride, $n_ch, $t0)
#
# Walk an extblock's data region and return its segments.
#
#   [ { off, n_samp, t_start }, ... ]
#
# The sample grid RESTARTS after every embedded header, because the header is
# hdr_len bytes and hdr_len % stride != 0 (442 % 76 = 62). So a header can only
# be tested for alignment against the START OF THE CURRENT SEGMENT -- testing it
# against the first sample of the file (as a first cut of this code did) rejects
# every header after the first, which is exactly the kind of off-by-one that
# produced the original bug.
# ---------------------------------------------------------------------------
sub _ext_segments {
    my ($fh, $file_size, $rec, $hdr_len, $stride, $n_ch, $t0) = @_;

    # candidate offsets: every "\x01TIME" in the data region
    my @cand;
    my $CHUNK = 1 << 22;
    for (my $b = $rec; $b < $file_size; $b += $CHUNK - $hdr_len) {
        my $len = $CHUNK;
        $len = $file_size - $b if $b + $len > $file_size;
        last if $len <= 0;
        seek $fh, $b, 0;
        my $buf;
        read($fh, $buf, $len) or last;
        my $p = -1;
        # CORE::index -- PDL::Slices exports an index() that shadows the builtin
        # ("Usage: PDL::index(a, ind)"). Same trap as t/09_i18n.t.
        while (($p = CORE::index($buf, "\x01TIME", $p + 1)) >= 0) {
            my $off = $b + $p;
            next if $off + $hdr_len > $file_size;
            push @cand, $off unless @cand && $cand[-1] == $off;
        }
    }

    my (@seg, @hdr);
    my $seg_off = $rec;
    my $seg_t   = $t0;
    for my $off (@cand) {
        next if $off < $seg_off;                          # inside a header we took
        next if ($off - $seg_off) % $stride;              # not on THIS segment's grid
        next unless _ext_is_hdr($fh, $off, $n_ch);
        my $n = int(($off - $seg_off) / $stride);
        push @seg, { off => $seg_off, n_samp => $n, t_start => $seg_t } if $n > 0;
        $seg_t   = _ext_hdr_time($fh, $off) // $seg_t;
        $seg_off = $off + $hdr_len;
    }
    my $n = int(($file_size - $seg_off) / $stride);
    push @seg, { off => $seg_off, n_samp => $n, t_start => $seg_t } if $n > 0;
    return @seg;
}

# Is there an embedded channel-info block at $off? The signature is strong:
# 0x01, "TIME", version 0x02 0x02, and a channel count that agrees with eb3.
sub _ext_is_hdr {
    my ($fh, $off, $n_ch) = @_;
    my $b = _read_bytes($fh, $off, 0x46);
    return 0 unless length($b) >= 0x46;
    return 0 unless substr($b, 0, 5) eq "\x01TIME";
    return 0 unless substr($b, 0x12, 2) eq "\x02\x02";
    return 0 unless substr($b, 0x14, 14) =~ /^\d{14}$/;      # YYYYMMDDHHMMSS
    return 0 unless unpack('v', substr($b, 0x44, 2)) == $n_ch - 1;
    return 1;
}

# The block header carries its segment's start time as ASCII YYYYMMDDHHMMSS at
# +0x14 -- no BCD, no epoch guessing.
sub _ext_hdr_time {
    my ($fh, $off) = @_;
    my $b = _read_bytes($fh, $off + 0x14, 14);
    return undef unless defined $b && $b =~ /^(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})$/;
    return "$1-$2-$3 $4:$5:$6";
}

sub _bcd_time {
    my ($fh, $addr) = @_;
    return sprintf("20%02d-%02d-%02d %02d:%02d:%02d",
        map { _bcd_byte($fh, $addr + 0x14 + $_) } 0 .. 5);
}

# Place .LOG events using the REC START markers as segment delimiters. Each
# REC START opens a segment; its exact data-sample start is block_meta's
# start_samp (from the embedded block header), and within a segment wall-clock
# and data advance 1:1 with no gaps, so an event's data offset is simply its
# .LOG time minus the opening REC START's .LOG time. This is exact and does not
# depend on what the .LOG's absolute clock means (it counts paused setup time
# between blocks, which is why _attach_seg_samp's wall-clock assumption drifts).
#
# Requires the REC START count to equal the segment count (1:1 delimiters);
# otherwise returns 0 so the caller falls back to the epoch-page map. Events
# before the first REC START anchor to segment 0.
sub _attach_recstart_samp {
    my ($events, $meta, $fs) = @_;
    return 0 unless $events && @$events && $meta && @$meta && $fs;
    my @rec = grep { ($_->{label} // '') =~ /REC\s*START/i } @$events;
    return 0 unless @rec == @$meta;                    # 1:1 delimiters required

    my $seg = -1;
    my $anchor_t    = 0;
    my $anchor_samp = $meta->[0]{start_samp};
    my $seg_end     = $meta->[0]{start_samp} + $meta->[0]{n_samp} - 1;
    for my $e (@$events) {
        if (($e->{label} // '') =~ /REC\s*START/i) {
            $seg++;
            $anchor_t    = $e->{t} // 0;
            $anchor_samp = $meta->[$seg]{start_samp};
            $seg_end     = $anchor_samp + $meta->[$seg]{n_samp} - 1;
        }
        my $off = ($e->{t} // 0) - $anchor_t;
        $off = 0 if $off < 0;
        my $s = $anchor_samp + int($off * $fs + 0.5);
        $s = $anchor_samp if $s < $anchor_samp;
        $s = $seg_end     if $s > $seg_end;            # clamp into this segment
        $e->{samp}   = $s;
        $e->{t_data} = $s / $fs;
    }
    return 1;
}

# Place .LOG events using the real segment anchors: within a segment wall-clock
# and data advance 1:1, so a wall-clock time t inside segment b maps to
#     samp = start_samp[b] + (t - t_start[b]) * fs
# This replaces _attach_epoch_samp()'s assumption that every segment is the same
# length -- an assumption that put events seconds away from where they belong.
sub _attach_seg_samp {
    my ($events, $meta, $fs) = @_;
    return 0 unless $events && @$events && $meta && @$meta && $fs;
    my @ep = map { _epoch_of($_->{t_start}) } @$meta;
    return 0 if grep { !defined } @ep;
    for my $e (@$events) {
        my $wall = $e->{epoch_sec};
        if (!defined $wall) {
            # .LOG t is seconds from the recording start
            next unless defined $e->{t};
            $wall = $ep[0] + $e->{t};
        }
        my $b = 0;
        for my $i (0 .. $#ep) { $b = $i if $wall >= $ep[$i] }
        my $s = $meta->[$b]{start_samp} + int(($wall - $ep[$b]) * $fs + 0.5);
        $s = $meta->[$b]{start_samp} if $s < $meta->[$b]{start_samp};
        my $end = $meta->[$b]{start_samp} + $meta->[$b]{n_samp} - 1;
        $s = $end if $s > $end;
        $e->{samp}   = $s;
        $e->{t_data} = $s / $fs;
    }
    return 1;
}

sub _epoch_of {
    my $ts = shift // '';
    my ($Y,$M,$D,$h,$m,$s) = $ts =~ /^(\d{4})-(\d{2})-(\d{2})[ T](\d{2}):(\d{2}):(\d{2})/
        or return undef;
    my $e = POSIX::mktime($s, $m, $h, $D, $M - 1, $Y - 1900, 0, 0, -1);
    return (defined $e && $e >= 0) ? $e : undef;
}

# ---------------------------------------------------------------------------
# Internal: waveform block header parse
#
# Confirmed structure from real EEG-1100C file (subject.EEG, 2025-12-21):
#   +0x00        : 0x01  block type
#   +0x01..+0x10 : ASCII time string "TIME164330000000" (16 bytes)
#   +0x14..+0x19 : BCD timestamp  YY MM DD HH MM SS
#   +0x1A..+0x1B : u16LE, lower 14 bits = sampling rate
#   +0x26        : n_ch_entries (channel table entries = n_valid_ch)
#   +0x2F..      : channel table: n_ch_entries × 10 bytes ([0x10][0x05][ch_idx][0×7])
#   +0x171       : data start (n_ch in stream = n_ch_entries + 1, trailing zero ch)
# ---------------------------------------------------------------------------
sub _read_wfm_header {
    my ($fh, $addr) = @_;

    my $block_type = _read_u8($fh, $addr);
    croak sprintf("Unexpected block type 0x%02X (expected 0x01)", $block_type)
        unless $block_type == 0x01;

    my $time_str = _read_bytes($fh, $addr + 0x01, 16);
    $time_str =~ s/\x00.*//s;  # trim at NUL

    my $t_start = sprintf("20%02d-%02d-%02d %02d:%02d:%02d",
        map { _bcd_byte($fh, $addr + 0x14 + $_) } 0 .. 5);

    # Sampling rate: +0x1A..+0x1B lower 14 bits (0xC3E8 & 0x3FFF = 1000)
    my $fs = _read_u16le($fh, $addr + 0x1A) & 0x3FFF;

    # +0x26: number of channel table entries (= n_valid_ch)
    my $n_ch_entries = _read_u8($fh, $addr + 0x26);
    croak "Zero channel entries in waveform block" unless $n_ch_entries > 0;

    # 10-byte entries starting at +0x2F: [0x10][0x05][ch_idx][0×7]
    my @ch_indices;
    for my $i (0 .. $n_ch_entries - 1) {
        my $eoff = $addr + 0x2F + $i * 10;
        push @ch_indices, _read_u8($fh, $eoff + 2);
    }

    return {
        n_ch        => $n_ch_entries + 1,     # data-stream channels (incl. zero-pad)
        n_ch_valid  => $n_ch_entries,         # real EEG/physio channels
        fs          => $fs,
        ch_indices  => \@ch_indices,          # electrode index per valid ch (1-indexed)
        t_start     => $t_start,
        data_offset => $addr + 0x171,
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
# Sections [ELECTRODE], [REFERENCE], ... key=value; keys are 0-based electrode
# indices. [REFERENCE] holds '$'-prefixed reference derivations; those are
# normalized to a Perl/filename-safe "<name>_ref" suffix and only fill indices
# that are blank in both [ELECTRODE] and the built-in defaults. [ELECTRODE] wins.
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
        # System-level fields. Keys here are names (not electrode indices), so
        # they ride along in %elec under string keys and never collide with the
        # integer electrode lookups; read_nk pulls them out into the record.
        if ($section eq 'SYSTEM_SETUP' && /^SystemReference=(.+)$/) {
            (my $v = $1) =~ s/\s+$//;
            $elec{system_reference} = $v;          # e.g. "C3,C4"
            next;
        }
        if ($section eq 'LASTPATTERN' && /^PATTERN=(-?\d+)/) {
            $elec{last_pattern} = $1 + 0;          # recording-time display montage no.
            next;
        }
        next unless $section eq 'ELECTRODE' || $section eq 'REFERENCE';
        if (/^(\d+)=(.+)$/) {
            my ($key, $val) = ($1 + 0, $2);
            $val =~ s/\s+$//;   # trim trailing whitespace
            if ($section eq 'ELECTRODE') { $elec{$key} = $val }
            else {
                $val =~ s/^\$(.+)/${1}_ref/;   # $A1 -> A1_ref (safe, no collision)
                $ref{$key} = $val;
            }
        }
    }
    close $fh;
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

            my $raw_label = substr($raw, 0, 20);
            $raw_label =~ s/\x00+$//;                       # trim NUL padding
            # Nihon Kohden labels are Shift_JIS (CP932). Decode to Unicode so
            # Japanese task names (e.g. 安静開眼) survive; fall back to the raw
            # bytes if decoding fails, then strip only control characters.
            my $label = eval { decode('cp932', $raw_label, Encode::FB_DEFAULT) };
            $label = $raw_label unless defined $label;
            $label =~ s/[\x00-\x08\x0B\x0C\x0E-\x1F\x7F]//g; # control chars only
            $label =~ s/\s+$//;
            next unless $label =~ /\S/;

            # bytes 20..25 = 6-digit ASCII seconds from recording start
            my $t_str = substr($raw, 20, 6);
            my $t_sec = ($t_str =~ /^(\d{6})$/) ? ($1 + 0)
                      : unpack('v', substr($raw, 20, 2));   # old-format fallback

            # byte 42 (u16 LE) = epoch number * 256 (0-based); +1 -> 1-based epoch.
            # Used to place events at their true data-sample position (extblock/
            # EEG-1200A stores gap-removed continuous data, so wall-clock t does
            # not map to samples; epoch does).
            my $ep_raw = unpack('v', substr($raw, 42, 2));
            my $epoch  = int($ep_raw / 256) + 1;

            push @events, { t => $t_sec, label => $label, epoch => $epoch };
        }
    }
    close $fh;
    return @events;
}

# Attach data-sample position to events from their epoch number.
# Epoch is a fixed-size acquisition block; a REC START marks a segment boundary
# and sits at an epoch boundary, so its sample = (epoch-1) * (n_total/max_epoch).
# Within a segment there are no gaps (wall-clock advances 1:1 with data), so
# other events (task markers, etc.) are placed at the preceding REC START's
# sample plus their wall-clock offset from that REC START. This keeps segment
# boundaries exact while giving markers second-accuracy (matching the vendor's
# "elapsed/積算時間"), instead of rounding everything to the epoch boundary.
# Sets $e->{samp} and $e->{t_data}=samp/fs.
sub _attach_epoch_samp {
    my ($events, $n_total, $fs) = @_;
    return unless $events && @$events && $n_total && $fs;
    my $max_ep = 0;
    for my $e (@$events) {
        my $ep = $e->{epoch} || 0;
        $max_ep = $ep if $ep > $max_ep;
    }
    return unless $max_ep > 1;                 # no usable epoch information
    my $L = $n_total / $max_ep;                # samples per epoch (float)

    my $clamp = sub {
        my $s = shift;
        $s = $n_total - 1 if $s > $n_total - 1;
        $s = 0            if $s < 0;
        return $s;
    };

    my $anchor;                                # last REC START { samp, t }
    for my $e (@$events) {
        next unless defined $e->{epoch} && $e->{epoch} >= 1;
        my $is_rec = ($e->{label} // '') =~ /REC\s*START/i;
        my $s;
        if ($is_rec) {
            $s = $clamp->(int(($e->{epoch} - 1) * $L + 0.5));   # epoch boundary
            $anchor = { samp => $s, t => ($e->{t} // 0) };
        } elsif ($anchor) {
            # same segment (events between this and the next REC START): add the
            # within-segment wall-clock offset (no gaps inside a segment).
            my $off = (($e->{t} // 0) - $anchor->{t}) * $fs;
            $off = 0 if $off < 0;
            $s = $clamp->(int($anchor->{samp} + $off + 0.5));
        } else {
            $s = $clamp->(int(($e->{epoch} - 1) * $L + 0.5));   # no anchor yet
        }
        $e->{samp}   = $s;
        $e->{t_data} = $s / $fs;
    }
    return;
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
sub _read_u8    { my ($fh,$o)=@_; unpack 'C', _read_bytes($fh,$o,1) }
sub _read_u16le { my ($fh,$o)=@_; unpack 'v', _read_bytes($fh,$o,2) }
sub _read_u32le { my ($fh,$o)=@_; unpack 'V', _read_bytes($fh,$o,4) }

# ---------------------------------------------------------------------------
# Device signature validation. Valid prefixes: "EEG-" / "QI-".
# Returns 0 (false) if VALID, 1 (true) if INVALID — matches nk2edf check_device()
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

# ---------------------------------------------------------------------------
# Block / segment selection helpers
# ---------------------------------------------------------------------------

=head2 block_extents($eeg_file, %opt)

Return the per-block extents of a file B<without reading any sample data>: only
the control-block address table and each waveform block's header are read.
C<n_samp> is derived from the address gap, exactly as the reader does.

  [ { index => 0, addr => 0x..., n_samp => 34000,
      start_samp => 0,     end_samp => 34000,
      t_start => "...", fs => 1000, n_ch => 34, n_ch_valid => 33 },
    { index => 1, ..., start_samp => 34100, end_samp => 407100, ... }, ... ]

C<start_samp>/C<end_samp> are in the same coordinate system that
C<read_nk($f, all_blocks =E<gt> 1)> produces. Both default to C<gap_samples =E<gt> 0>
(blocks butt-joined); if you pass a non-zero C<gap_samples> to one, pass the same
to the other or the coordinates will not line up.

C<t_start> is each block's own wall-clock start, so the elapsed gap at a break is
C<epoch(t_start[b+1]) - epoch(t_start[b]) - n_samp[b]/fs> -- the information the
old zero-padding was standing in for, without putting fake samples in the data.

Use this to plan a partial read -- "which blocks does 100-500 s touch?" -- and
then C<read_nk($f, block =E<gt> $i)> only those. C<block_ranges> cannot do this:
it operates on a record that has already been read.

extblock (EEG-1200A) files hold a single data block, so one entry is returned.

=head2 block_ranges($rec, %opt)

Return the block/segment boundaries of a record as an arrayref of hashrefs
(0-based, end-exclusive samples):

  [ { index => 0, start => 0,     end => 34000, t_start => "..." },
    { index => 1, start => 34000, end => 68000, t_start => "..." }, ... ]

Source of the boundaries (option C<source>, default C<'auto'>):

  physical  use the reader's block_meta (real waveform blocks; wfmblock)
  log       split by ".LOG" "REC START" markers in $rec->{events}
            (for extblock/EEG-1200A, which is one physical block per file)
  auto      physical if >1 block, else log, else a single segment

=head2 select_block($rec, $i, %opt)

Return a new record (same contract as read_nk) containing only block/segment
C<$i> (0-based). Data is sliced, events are filtered to the segment and rebased
to segment-local time, and per-segment C<t_start> is set. Useful for writing one
EDF per recording segment.

=cut

sub _add_seconds {
    my ($tstr, $sec) = @_;
    return $tstr unless defined $tstr
        && $tstr =~ /^(\d{4})-(\d{2})-(\d{2})[ T](\d{2}):(\d{2}):(\d{2})/;
    my $ep = POSIX::mktime($6, $5, $4, $3, $2 - 1, $1 - 1900, 0, 0, -1);
    return $tstr if !defined $ep || $ep < 0;
    my @lt = localtime($ep + $sec);
    return POSIX::strftime("%Y-%m-%d %H:%M:%S", @lt);
}

# ---------------------------------------------------------------------------
# block_extents($eeg_path, %opt) -> arrayref of per-block extents, HEADER ONLY.
#
# Reads the control-block address table and each waveform block's header, but NO
# sample data. n_samp is derived from the address gap exactly as _read_wfm_block
# does, so start_samp/end_samp are in the SAME coordinate system that
# read_nk(all_blocks => 1) produces. Both default to gap_samples => 0 (blocks
# butt-joined); a non-zero value must be passed to BOTH or they will not line up.
#
# This lets a caller plan a partial read ("which blocks does 100-500 s touch?")
# without loading the whole recording. block_ranges() cannot do this: it needs a
# record that has already been read.
#
#   [ { index, addr, n_samp, start_samp, end_samp, t_start, fs, n_ch, n_ch_valid }, ... ]
#
# extblock (EEG-1200A) files hold a single data block, so this returns one entry.
# ---------------------------------------------------------------------------
sub block_extents {
    my ($eeg_path, %opt) = @_;
    croak "File not found: $eeg_path" unless -f $eeg_path;
    my $gap_samples = $opt{gap_samples} // 0;

    my (undef, $layout) = nk_layout($eeg_path);
    croak "Unknown Nihon Kohden signature in $eeg_path" unless defined $layout;

    open my $fh, '<:raw', $eeg_path or croak "Cannot open $eeg_path: $!";
    seek $fh, 0, 2;
    my $file_size = tell $fh;

    if ($layout eq 'extblock') {
        # extblock is NOT a single data block: the recorder re-emits the
        # channel-info block (72 + (n_ch-1)*10 bytes) into the sample stream at
        # every recording break. Walk them the same way _read_extblock does.
        my $ext  = _read_u32le($fh, 0x03EE)
            or croak "extblock: ext_address is 0 in $eeg_path";
        my $ctl0      = _read_u32le($fh, 0x0092);
        my $data_addr = _read_u32le($fh, $ctl0 + 18);
        my $fs        = $opt{fs} // (_read_u16le($fh, $data_addr + 0x1A) & 0x3FFF);
        my $eb2  = _read_u32le($fh, $ext + 18);
        my $eb3  = _read_u32le($fh, $eb2 + 20);
        my $n_ch = _read_u16le($fh, $eb3 + 68) + 1;
        my $hdr_len = 72 + ($n_ch - 1) * 10;
        my $stride  = $n_ch * 2;
        my $rec     = $eb3 + $hdr_len;

        my @seg = _ext_segments($fh, $file_size, $rec, $hdr_len, $stride, $n_ch,
                               _ext_hdr_time($fh, $eb3) // _bcd_time($fh, $data_addr));
        my @out;
        my $pos = 0;
        for my $i (0 .. $#seg) {
            push @out, { index      => $i,
                         addr       => $seg[$i]{off},
                         n_samp     => $seg[$i]{n_samp},
                         start_samp => $pos,
                         end_samp   => $pos + $seg[$i]{n_samp},
                         t_start    => $seg[$i]{t_start},
                         fs         => $fs,
                         n_ch       => $n_ch,
                         n_ch_valid => $n_ch - 1 };
            $pos += $seg[$i]{n_samp};
            $pos += $gap_samples if $i < $#seg;
        }
        close $fh;
        return \@out;
    }

    # --- wfmblock: control block list -> waveform block addresses -------------
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
    croak "No waveform blocks in $eeg_path" unless @wfm_addrs;

    my @out;
    my $pos = 0;
    for my $b (0 .. $#wfm_addrs) {
        my $waddr = $wfm_addrs[$b];
        my $next  = ($b + 1 <= $#wfm_addrs) ? $wfm_addrs[$b + 1] : $file_size;
        my $h     = _read_wfm_header($fh, $waddr);            # header only
        my $data_bytes = ($next - $waddr) - ($h->{data_offset} - $waddr);
        croak "block_extents: block $b data_bytes=$data_bytes not divisible by "
            . "n_ch*2=" . ($h->{n_ch} * 2)
            if $data_bytes % ($h->{n_ch} * 2) != 0;
        my $n_samp = int($data_bytes / ($h->{n_ch} * 2));
        push @out, { index      => $b,
                     addr       => $waddr,
                     n_samp     => $n_samp,
                     start_samp => $pos,                       # all_blocks coords
                     end_samp   => $pos + $n_samp,
                     t_start    => $h->{t_start},
                     fs         => $h->{fs},
                     n_ch       => $h->{n_ch},
                     n_ch_valid => $h->{n_ch_valid} };
        $pos += $n_samp;
        $pos += $gap_samples if $b < $#wfm_addrs;              # matches read_nk
    }
    close $fh;
    return \@out;
}

sub block_ranges {
    my ($rec, %opt) = @_;
    croak "block_ranges: need record hashref" unless ref $rec eq 'HASH';
    my $n_total = eval { $rec->{data}->dim(1) } // 0;
    my $fs      = $rec->{fs} || 1;
    my $src     = $opt{source} || 'auto';

    my $meta = $rec->{block_meta} || [];
    my $use_physical = ($src eq 'physical')
                    || ($src eq 'auto' && @$meta > 1);

    my @out;
    if ($use_physical && @$meta) {
        for my $b (0 .. $#$meta) {
            my $lo = $meta->[$b]{start_samp};
            my $ns = $meta->[$b]{n_samp};
            push @out, { index => $b, start => $lo, end => $lo + $ns,
                         t_start => $meta->[$b]{t_start} };
        }
    } else {
        # Segment by .LOG "REC START" markers. Prefer the epoch-derived data
        # sample position ($e->{samp}); fall back to wall-clock t*fs only if no
        # epoch info (which for gap-removed extblock data would be wrong -> we
        # then refuse to split; see below).
        my @starts;
        my $have_samp = 0;
        for my $e (@{ $rec->{events} || [] }) {
            my $lab = ref $e eq 'HASH' ? ($e->{label} // '') : "$e";
            next unless $lab =~ /REC\s*START/i;
            if (defined $e->{samp}) { push @starts, $e->{samp}; $have_samp = 1; }
            else                    { push @starts, int(($e->{t} || 0) * $fs + 0.5); }
        }
        my %seen; @starts = grep { !$seen{$_}++ } sort { $a <=> $b } @starts;
        unshift @starts, 0 unless @starts && $starts[0] == 0;
        @starts = (0) unless @starts;

        # Without epoch (samp) info, wall-clock starts may exceed the gap-removed
        # data length; we cannot split reliably, so return a single segment.
        if (!$have_samp && grep { $_ >= $n_total } @starts) {
            carp "block_ranges: .LOG times exceed data length and no epoch info "
               . "(gap-removed data). Returning a single segment.";
            return [ { index => 0, start => 0, end => $n_total,
                       t_start => $rec->{t_start} } ];
        }

        for my $i (0 .. $#starts) {
            my $lo = $starts[$i];
            next if $lo >= $n_total;
            my $hi = $i < $#starts ? $starts[$i + 1] : $n_total;
            $hi = $n_total if $hi > $n_total;
            next if $hi <= $lo;
            push @out, { index => scalar @out, start => $lo, end => $hi,
                         t_start => _add_seconds($rec->{t_start}, $lo / $fs) };
        }
    }
    return \@out;
}

sub select_block {
    my ($rec, $i, %opt) = @_;
    my $ranges = block_ranges($rec, %opt);
    croak "select_block: index $i out of range (have " . scalar(@$ranges) . " block(s))"
        if $i < 0 || $i > $#$ranges;
    my $r = $ranges->[$i];
    return select_range($rec, $r->{start}, $r->{end}, t_start => $r->{t_start});
}

# Extract an arbitrary [$lo, $hi) sample range as a new record (rebased events).
sub select_range {
    my ($rec, $lo, $hi, %opt) = @_;
    my $fs = $rec->{fs} || 1;
    my $n  = eval { $rec->{data}->dim(1) } // 0;
    $lo = 0  if $lo < 0;
    $hi = $n if $hi > $n;
    croak "select_range: empty/invalid range [$lo,$hi)" if $hi <= $lo;
    my $len = $hi - $lo;

    my $data = $rec->{data}->slice(":," . $lo . ":" . ($hi - 1))->sever;

    my @events;
    for my $e (@{ $rec->{events} || [] }) {
        my $s = defined $e->{samp} ? $e->{samp} : ($e->{t} || 0) * $fs;
        next unless $s >= $lo && $s < $hi;
        my %ne = %$e;
        if (defined $e->{samp}) {
            $ne{samp}   = $e->{samp} - $lo;
            $ne{t_data} = $ne{samp} / $fs;
            $ne{t}      = $ne{samp} / $fs;
        } else {
            $ne{t} = ($s - $lo) / $fs;
        }
        push @events, \%ne;
    }

    my $tstart = defined $opt{t_start} ? $opt{t_start}
               : _add_seconds($rec->{t_start}, $lo / $fs);

    return {
        %$rec,
        data             => $data,
        events           => \@events,
        t_start          => $tstart,
        n_blocks         => 1,
        all_blocks       => 0,
        gap_bounds       => [],
        t_block_starts   => [0],
        n_samp_per_block => [$len],
        block_meta       => [ { start_samp => 0, n_samp => $len, t_start => $tstart } ],
    };
}

# Anchors for the wall-clock -> data-sample map, as { t, samp[, end] } with t in
# seconds from the recording start. The map is 1:1 within a segment with a jump
# at each recording break, so a single anchor per segment (its wall-clock start
# and its data-sample start) is enough.
#
# Source, in preference order:
#   1. block_meta when it resolves >1 segment. Each block's t_start is a real
#      wall-clock start (wfmblock: per block; extblock: from the embedded block
#      header), so this is authoritative and needs no .LOG. This is the same
#      map _attach_seg_samp() uses to place events. It is what fixes wall-clock
#      queries on multi-block recordings, where there are no REC START events
#      and the old code silently fell back to a 1:1 (wall==data) map.
#   2. REC START markers, for a hand-built rec or a file whose segment headers
#      did not resolve but whose .LOG carries REC START events with {samp}.
# With neither, the caller gets a 1:1 map (correct for a single continuous
# segment).
sub _clock_anchors {
    my ($rec) = @_;
    my $meta = $rec->{block_meta} || [];
    if (@$meta > 1) {
        my @ep = map { _epoch_of($_->{t_start}) } @$meta;
        unless (grep { !defined } @ep) {
            my $t0 = $ep[0];
            return map {
                { t    => $ep[$_] - $t0,
                  samp => $meta->[$_]{start_samp},
                  end  => $meta->[$_]{start_samp} + $meta->[$_]{n_samp} - 1 }
            } 0 .. $#$meta;
        }
    }
    my @ev = grep { defined $_->{samp} && defined $_->{t}
                    && ($_->{label} // '') =~ /REC\s*START/i }
             @{ $rec->{events} || [] };
    return sort { $a->{t} <=> $b->{t} }
           map { { t => $_->{t}, samp => $_->{samp} } } @ev;
}

# Convert a wall-clock time (seconds from recording start) to a data-sample
# index. Uses the block_meta piecewise-linear map (see _clock_anchors); within
# a segment wall-clock and data advance 1:1, and a query that lands in a
# recording gap clamps to the last real sample of the preceding segment rather
# than inventing samples that are not in the data.
sub clock_to_samp {
    my ($rec, $wall_sec) = @_;
    my $fs = $rec->{fs} || 1;
    my $n  = eval { $rec->{data}->dim(1) } // 0;

    my @anch = _clock_anchors($rec);
    return int($wall_sec * $fs + 0.5) unless @anch;   # no anchors: assume 1:1

    my $best;
    for my $a (@anch) { $best = $a if $a->{t} <= $wall_sec; }
    $best ||= $anch[0];

    my $s   = $best->{samp} + int(($wall_sec - $best->{t}) * $fs + 0.5);
    my $end = defined $best->{end} ? $best->{end} : ($n ? $n - 1 : $s);
    $s = $best->{samp} if $s < $best->{samp};   # clamp into this segment
    $s = $end          if $s > $end;
    $s = 0             if $s < 0;
    $s = $n - 1        if $n && $s > $n - 1;     # never past the data
    return $s;
}

1;

__END__

=head1 KNOWN LIMITATIONS

=over 4

=item *

Two on-disk layouts are supported, dispatched by the file's format signature
via C<nk_layout()>: the legacy 'wfmblock' layout (EEG-1100x / EEG-2100 / QI-403A
signatures) and the newer 'extblock' layout (EEG-1200A signature; used e.g. by
the EEG-1290 recorder). Unseen signatures fall back to a structural check
(ext_address at 0x03EE) so future variants (EEG-1200B/C, updated firmware) read
without code changes; truly non-NK signatures are rejected.

=item *

C<all_blocks =E<gt> 1> concatenates all wfmblock waveform blocks. Blocks are
butt-joined and the data are NOT padded to real elapsed time: a recording break
is a discontinuity, not a stretch of samples. Break positions are in
C<t_block_starts> / C<block_meta> (concatenated-sample coordinates), and each
block's wall-clock C<t_start> is there too, so the elapsed gap at break I<b> is

  epoch(t_start[b+1]) - epoch(t_start[b]) - n_samp[b] / fs

Epoching must reject segments that straddle a C<t_block_starts> boundary.

B<Changed in 0.2:> earlier versions inserted C<gap_samples> (default 100) zero
samples at each break and reported them in C<gap_bounds>, so that a break was
visible as a flat stretch. Those zeros are not data -- they ring through
filters, corrupt spectra, and skew waveform-similarity metrics -- so the default
is now 0 and C<gap_bounds> is empty unless C<gap_samples> is set explicitly.
Concatenated-sample coordinates therefore shifted by 100 samples per preceding
break. Channel layout is assumed constant across blocks; a mismatch is a fatal
error. extblock files are a single data block, so C<all_blocks> is a no-op
there.

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
