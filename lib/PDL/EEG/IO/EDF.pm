package PDL::EEG::IO::EDF;

use strict;
use warnings;
use Carp qw(croak carp);
use PDL;
use POSIX qw(floor);
use Encode qw(encode decode);
use Exporter 'import';

our @EXPORT_OK = qw(write_edf read_edf clean_edf_label);
our $VERSION   = '0.02';

=head1 NAME

PDL::EEG::IO::EDF - Write EEG records to EDF / EDF+ (European Data Format)

=head1 SYNOPSIS

  use PDL::EEG::IO::NihonKohden qw(read_nk);
  use PDL::EEG::IO::EDF         qw(write_edf);

  my $rec = read_nk('subject.EEG');       # { data,fs,labels,events,... }
  write_edf($rec, 'out.edf');              # EDF+C, events -> annotations

  # Plain EDF (no annotation channel), byte-compatible physical scaling:
  write_edf($rec, 'out.edf', plus => 0, phys => 'gain');

=head1 DESCRIPTION

Writes the hash returned by C<read_nk> as a 16-bit EDF file.  The record is
expected to contain:

  data    PDL [n_ch, n_samples] float32, physical units (uV) -- ALWAYS uV,
          for every channel, DC included. A DC channel cannot be WRITTEN in uV
          (EDF's physical_min is 8 ASCII chars and "-12002900" is 9), so it is
          written in mV and declared as such; {units} says which dimension each
          signal goes out in / came in as. read_edf normalises back to uV.
  fs      sampling rate in Hz (integer)
  labels  arrayref of n_ch channel names
  events  (optional) arrayref of events -> EDF+ annotations
  start   (optional) [Y,M,D,h,m,s] or epoch seconds; else $rec->{start_datetime}
  t_block_starts (optional, from all_blocks=>1) sample indices of block starts

By default an B<EDF+C> (continuous) file is written with an "EDF Annotations"
channel so that C<events> survive the round trip. Pass C<< plus => 0 >> for a
plain EDF file (no annotation channel).

=head2 Options

  plus       => 1        1 = EDF+C (default), 0 = plain EDF
  phys       => 'auto'   per-channel min/max -> full 16-bit range (default)
             => 'gain'   fixed +/- 32768*gain uV (exact NK round-trip, may clip)
             => <number> fixed symmetric physical range +/- <number> uV
  gain       => <num>    force one uV/bit for phys eq 'gain'; default uses
                         $rec->{gains} per channel (DC/STIM differ), else 0.09765625
  record_dur => 1.0      seconds per data record
  unit       => 'uV'     force ONE physical dimension on every signal. Without
                         it, $rec->{units} (from read_nk) is used per signal, so
                         EEG goes out as uV and DC as mV -- which is what EDF+
                         wants and what the 8-char physical_min field requires.
  subject    => ''       subject CODE (written to the EDF+ patient-id field);
                         'patient' accepted as a deprecated alias
  sex        => 'X'      EDF+: 'M' | 'F' | 'X'
  birthdate  => 'X'      EDF+: 'dd-MMM-yyyy' (e.g. 02-MAY-1951) or 'X'
  name       => 'X'      EDF+: subject name (spaces -> '_')
  recording  => ''       EDF+: investigation/admin CODE subfield (else free text)
  technician => 'X'      EDF+: technician code
  equipment  => 'X'      EDF+: equipment/device (default: vendor + $rec->{device},
                         e.g. Nihon_Kohden_EEG-1200A_V01.00)
                         (EDF+ builds "code sex bd name" / "Startdate date admin tech equip";
                          plain EDF, plus=>0, leaves subject/recording as free text)
  start      => undef     recording start; default parses $rec->{t_start}
                          ("YYYY-MM-DD HH:MM:SS"). Also accepts [Y,M,D,h,m,s] or epoch.
  labels     => undef     override channel labels (arrayref)
  annotate_blocks => 1   mark multi-block boundaries as annotations
  annot_encoding => 'utf8'  annotation text encoding: 'utf8' (default, keeps
                         Japanese labels e.g. 安静開眼) or 'ascii' (non-ASCII
                         bytes replaced with '_')

=head1 CAVEATS

Discontinuous multi-block recordings (C<all_blocks=>1> with time gaps) are
written as B<continuous> EDF+C; block starts are added as annotations so the
gaps are at least visible. True EDF+D (per-record start times) is not yet
implemented.

=cut

# ---------------------------------------------------------------------------
# byte-exact field formatting
# ---------------------------------------------------------------------------

sub _ascii {                       # left-justified, space-padded, hard-truncated
    my ($s, $len) = @_;
    $s = '' unless defined $s;
    $s = substr($s, 0, $len);
    return $s . (' ' x ($len - length $s));
}

sub _num {                         # numeric field that MUST fit in $len (<=8) chars
    my ($x, $len) = @_;
    $len ||= 8;
    if ($x == int($x) && abs($x) < 10 ** ($len - ($x < 0 ? 1 : 0))) {
        return _ascii(sprintf('%d', $x), $len);
    }
    for my $dp (reverse 0 .. $len - 1) {
        my $s = sprintf('%.*f', $dp, $x);
        if ($s =~ /\./) { $s =~ s/0+$//; $s =~ s/\.$//; }
        return _ascii($s, $len) if length($s) <= $len;
    }
    return _ascii(substr(sprintf('%.*g', $len - 2, $x), 0, $len), $len);
}

sub _onset {                       # EDF+ onset: always signed, compact
    my $t = shift;
    my $s = sprintf('%.6f', $t);
    $s =~ s/0+$//; $s =~ s/\.$//;
    return ($t < 0 ? '' : '+') . $s;
}

sub _tal {                         # one Time-stamped Annotation List
    my ($onset, $text, $dur) = @_;
    my $s = _onset($onset);
    if (defined $dur) {
        my $d = sprintf('%.6f', $dur); $d =~ s/0+$//; $d =~ s/\.$//;
        $s .= "\x15$d";
    }
    $s .= defined $text ? "\x14$text\x14\x00" : "\x14\x14\x00";
    return $s;
}

# ---------------------------------------------------------------------------
# start date/time
# ---------------------------------------------------------------------------

sub _startdatetime {
    my ($start) = @_;
    my @lt;
    if (!defined $start || $start eq '') {
        return ('01.01.85', '00.00.00', 0);   # placeholder; caller warns
    } elsif (ref $start eq 'ARRAY') {
        @lt = @$start;                          # [Y,M,D,h,m,s]
    } elsif ($start =~ /^\s*(\d{4})-(\d{2})-(\d{2})[ T](\d{2}):(\d{2}):(\d{2})/) {
        @lt = ($1, $2, $3, $4, $5, $6);         # read_nk t_start string
    } elsif ($start =~ /^\d+$/) {
        my @t = localtime($start);              # epoch seconds
        @lt = ($t[5] + 1900, $t[4] + 1, $t[3], $t[2], $t[1], $t[0]);
    } else {
        return ('01.01.85', '00.00.00', 0);     # unrecognised
    }
    my ($Y, $M, $D, $h, $m, $s) = @lt;
    my $yy = $Y % 100;
    return (sprintf('%02d.%02d.%02d', $D, $M, $yy),
            sprintf('%02d.%02d.%02d', $h, $m, $s), 1, $Y, $M, $D);
}

# EDF+ subfield sanitiser: spaces -> underscore, empty -> 'X'
sub _sub {
    my $s = shift;
    return 'X' unless defined $s && length $s;
    $s =~ s/\s+/_/g;
    return $s;
}

my @MON3 = qw(XXX JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC);

# dd-MMM-yyyy (EDF+), or 'X' if unknown
sub _ddmmmyyyy {
    my ($have, $Y, $M, $D) = @_;
    return 'X' unless $have && $M >= 1 && $M <= 12;
    return sprintf('%02d-%s-%04d', $D, $MON3[$M], $Y);
}

# EDF+ local patient identification field (called "subject" in this API):
# "code sex birthdate name"
sub _edfplus_subject {
    my %o = @_;
    my $sex = (defined $o{sex} && $o{sex} =~ /^[MF]$/i) ? uc $o{sex} : 'X';
    my $bd  = (defined $o{birthdate} && length $o{birthdate}) ? _sub($o{birthdate}) : 'X';
    return join ' ', _sub($o{code}), $sex, $bd, _sub($o{name});
}

# EDF+ local recording identification: "Startdate dd-MMM-yyyy admin tech equip"
sub _edfplus_recording {
    my ($startdate, %o) = @_;
    return join ' ', 'Startdate', ($startdate // 'X'),
                     _sub($o{admin}), _sub($o{tech}), _sub($o{equip});
}

# --- equipment / vendor identification -------------------------------------
# The EDF+ recording-id "equipment" subfield should name the acquisition
# system (e.g. Nihon_Kohden_EEG-1200A_V01.00). The device signature stored in
# $rec->{device} carries only the model+version (e.g. "EEG-1200A V01.00"); the
# vendor name is not in the file, so it is inferred from the model prefix via
# the table below. First matching pattern wins. Extend as new readers appear.
my @VENDOR_TABLE = (
    [ qr/^EEG-1[0-9]{3}/i => 'Nihon Kohden' ],   # EEG-1100C, EEG-1200A, ...
    [ qr/^Neurofax/i      => 'Nihon Kohden' ],
    [ qr/^QP-/i           => 'Nihon Kohden' ],    # Neurofax QP series
);

# default vendor/equipment when the device signature is absent or unrecognised.
# At present no other acquisition system is supported end-to-end, so unknown
# sources are labelled as Brain Products / BrainAmp rather than left blank.
my $DEFAULT_EQUIP = 'Brain Products BrainAmp';

sub _vendor_for {
    my $dev = shift;
    return undef unless defined $dev;
    for my $r (@VENDOR_TABLE) { return $r->[1] if $dev =~ $r->[0] }
    return undef;
}

# build the equipment string. Explicit override (--equipment) always wins;
# else prefix the recognised vendor to the device signature; else fall back
# to $DEFAULT_EQUIP. _sub() later turns spaces into '_' for the EDF+ field.
sub _equip_string {
    my ($rec, $override) = @_;
    return $override if defined $override && length $override;

    my $dev = defined $rec->{device} ? $rec->{device} : '';
    $dev =~ s/^\s+//; $dev =~ s/\s+$//;
    return $DEFAULT_EQUIP unless length $dev;

    return $dev if $dev =~ /nihon|kohden|brain/i;      # already vendor-prefixed
    my $vendor = _vendor_for($dev);
    return defined $vendor ? "$vendor $dev" : $DEFAULT_EQUIP;
}

# ---------------------------------------------------------------------------
# normalise events -> [ [onset_sec, text], ... ]
# ---------------------------------------------------------------------------

sub _norm_events {
    my ($events, $fs, $enc) = @_;
    $enc = 'utf8' unless defined $enc;
    my @out;
    return \@out unless $events && @$events;
    for my $e (@$events) {
        my ($t, $txt);
        if (ref $e eq 'HASH') {
            for my $k (qw(t_data time t onset sec sample)) {
                if (exists $e->{$k}) { $t = $e->{$k}; $t /= $fs if $k eq 'sample'; last }
            }
            for my $k (qw(label name text desc annotation)) {
                if (defined $e->{$k}) { $txt = $e->{$k}; last }
            }
        } elsif (ref $e eq 'ARRAY') {
            ($t, $txt) = @$e;
        } else {
            $txt = $e;
        }
        next unless defined $t;
        $txt = 'event' unless defined $txt && length $txt;
        # Encode to octets so downstream TAL byte-length math is correct.
        # EDF+ is nominally Latin-1, but UTF-8 annotations are widely read
        # (MNE, recent EDFbrowser); 'ascii' mode keeps only 7-bit bytes.
        if ($enc eq 'ascii') {
            $txt = encode('UTF-8', $txt) if utf8::is_utf8($txt);
            $txt =~ s/[^\x20-\x7e]/_/g;
        } else {                                  # utf8 (default)
            $txt = encode('UTF-8', $txt);
        }
        $txt =~ s/[\x00\x14\x15]/ /g;            # strip TAL delimiters
        push @out, [ 0 + $t, "$txt" ];
    }
    return [ sort { $a->[0] <=> $b->[0] } @out ];
}

# ---------------------------------------------------------------------------
# main entry
# ---------------------------------------------------------------------------

sub write_edf {
    my ($rec, $path, %o) = @_;
    croak "write_edf: need a record hashref"   unless ref $rec eq 'HASH';
    croak "write_edf: need an output path"      unless defined $path && length $path;

    my $data = $rec->{data};
    croak "write_edf: \$rec->{data} is not a PDL" unless eval { $data->isa('PDL') };
    my $fs   = $rec->{fs} or croak "write_edf: \$rec->{fs} missing";
    my $n_ch = $data->dim(0);
    my $n_s  = $data->dim(1);
    croak "write_edf: empty data" if $n_ch < 1 || $n_s < 1;

    my $plus  = exists $o{plus} ? $o{plus} : 1;
    my $phys  = defined $o{phys} ? $o{phys} : 'auto';
    my $gain  = defined $o{gain} ? $o{gain} : 0.09765625;
    my $rdur  = defined $o{record_dur} ? $o{record_dur} : 1.0;
    my $unit  = defined $o{unit} ? $o{unit} : 'uV';
    my $ab    = exists $o{annotate_blocks} ? $o{annotate_blocks} : 1;
    my $aenc  = defined $o{annot_encoding} ? $o{annot_encoding} : 'utf8';

    my @labels = $o{labels} ? @{ $o{labels} }
               : $rec->{labels} ? @{ $rec->{labels} }
               : map { "ch$_" } 1 .. $n_ch;
    croak "write_edf: label count ($#labels+1) != n_ch ($n_ch)"
        if @labels != $n_ch;

    my $spr = int($fs * $rdur + 0.5);            # signal samples per data record
    croak "write_edf: fs*record_dur must be a positive integer (got $spr)"
        if $spr < 1;
    my $n_rec = int(($n_s + $spr - 1) / $spr);   # ceil -> pad last record

    # per-channel gains (µV/bit) for phys=>'gain': explicit opt overrides,
    # else use $rec->{gains} from read_nk (DC/STIM channels differ), else scalar.
    my @gain_ch;
    if (defined $o{gain}) {
        @gain_ch = ($o{gain}) x $n_ch;
    } elsif (eval { $rec->{gains}->isa('PDL') } && $rec->{gains}->nelem == $n_ch) {
        @gain_ch = $rec->{gains}->list;
    } else {
        @gain_ch = ($gain) x $n_ch;
    }

    # --- PER-SIGNAL physical dimension --------------------------------------
    #
    # EDF gives every signal its own physical dimension, and here it is not a
    # nicety -- it is forced. {data} is uniformly µV, but an EEG-1200A DC input is
    # a ±12 V line, i.e. ±12002900 µV. EDF's physical_min field is EIGHT ASCII
    # characters, and "-12002900" is nine. A DC channel simply CANNOT be written
    # in µV. In mV it is "-12002.9" -- eight characters, exactly.
    #
    # So: {units} (from read_nk) says which dimension each channel should be
    # written in; the samples are divided by $UV_PER{dim} on the way out, and the
    # declared physical range is expressed in the same dimension. An explicit
    # unit => '...' still forces one dimension on everything, as before.
    my %UV_PER = ('uV' => 1, 'uv' => 1, 'µV' => 1, 'mV' => 1_000, 'mv' => 1_000);
    my @dim_ch;
    if (defined $o{unit}) {
        @dim_ch = ($o{unit}) x $n_ch;                      # caller forced one unit
    } elsif (ref $rec->{units} eq 'ARRAY' && @{ $rec->{units} } == $n_ch) {
        @dim_ch = @{ $rec->{units} };
    } else {
        @dim_ch = ($unit) x $n_ch;
    }
    my @uvper = map { $UV_PER{ $dim_ch[$_] } // 1 } 0 .. $n_ch - 1;  # 'code' -> 1

    # --- per-channel physical/digital ranges -------------------------------
    my $DMIN = -32768;
    my $DMAX =  32767;
    my (@pmin, @pmax, @dig);                      # @dig: padded short piddles
    for my $c (0 .. $n_ch - 1) {
        # express this signal in ITS declared dimension (µV / uvper). For a µV
        # channel uvper is 1 and nothing changes; for a DC channel written in mV
        # it is 1000, and both the samples and the declared range shrink by 1000.
        my $uvp = $uvper[$c] || 1;
        my $ch  = $uvp == 1 ? $data->slice("($c),:")
                            : $data->slice("($c),:") / $uvp;
        my ($pmn, $pmx);
        if ($phys eq 'auto') {
            ($pmn, $pmx) = $ch->minmax;
            if ($pmn == $pmx) { $pmn -= 1; $pmx += 1; }    # avoid /0 on flat ch
        } elsif ($phys eq 'gain') {
            my $g = ($gain_ch[$c] || $gain) / $uvp;        # µV/bit -> dim/bit
            $pmn = $DMIN * $g;                             # e.g. -3200 µV, -12002.9 mV
            $pmx = $DMAX * $g;
        } else {                                          # explicit +/- range
            my $r = abs($phys) || 1;
            $pmn = -$r; $pmx = $r;
        }
        my $scale = ($DMAX - $DMIN) / ($pmx - $pmn);
        my $d = ((($ch - $pmn) * $scale) + $DMIN)->rint;
        my $nclip = ($d < $DMIN)->sum + ($d > $DMAX)->sum;
        carp sprintf("write_edf: channel %d (%s): %d sample(s) clipped to +/-32k",
                     $c, $labels[$c], $nclip) if $nclip;
        $d = $d->clip($DMIN, $DMAX);
        # pad to n_rec*spr with the digital code nearest physical 0
        my $zero = sprintf('%.0f', (0 - $pmn) * $scale + $DMIN);
        $zero = $DMIN if $zero < $DMIN; $zero = $DMAX if $zero > $DMAX;
        my $full = zeroes(short, $n_rec * $spr) + $zero;
        $full->slice("0:" . ($n_s - 1)) .= $d->short;
        push @pmin, $pmn; push @pmax, $pmx; push @dig, $full;
    }

    # --- annotations per record (EDF+ only) --------------------------------
    my $annot_spr = 0;
    my @rec_tal;                                  # $rec_tal[$r] = concatenated TAL bytes
    if ($plus) {
        @rec_tal = ('') x $n_rec;
        # timekeeping TAL first in every record
        $rec_tal[$_] = _tal($_ * $rdur, undef) for 0 .. $n_rec - 1;

        my $ev = _norm_events($rec->{events}, $fs, $aenc);
        for my $e (@$ev) {
            my ($t, $txt) = @$e;
            next if $t < 0 || $t >= $n_rec * $rdur;
            my $r = int($t / $rdur); $r = $n_rec - 1 if $r > $n_rec - 1;
            $rec_tal[$r] .= _tal($t, $txt);
        }
        if ($ab && $rec->{t_block_starts} && @{ $rec->{t_block_starts} } > 1) {
            my @bs = @{ $rec->{t_block_starts} };
            for my $i (1 .. $#bs) {               # skip block 0 (== recording start)
                my $t = $bs[$i] / $fs;
                next if $t >= $n_rec * $rdur;
                my $r = int($t / $rdur); $r = $n_rec - 1 if $r > $n_rec - 1;
                $rec_tal[$r] .= _tal($t, sprintf('Block %d', $i + 1));
            }
        }
        my $maxlen = 0;
        for (@rec_tal) { $maxlen = length($_) if length($_) > $maxlen; }
        my $bytes = $maxlen < 2 ? 2 : $maxlen;
        $bytes++ if $bytes % 2;                   # 2-byte samples
        $annot_spr = $bytes / 2;
    }

    # --- header ------------------------------------------------------------
    my $start_src = defined $o{start} ? $o{start} : $rec->{t_start};
    my ($sdate, $stime, $have_start, $sY, $sM, $sD) = _startdatetime($start_src);
    carp "write_edf: no usable start date/time; using placeholder 01.01.85 00.00.00"
        unless $have_start;

    my $ns_total = $n_ch + ($plus ? 1 : 0);
    my $hbytes   = 256 * (1 + $ns_total);
    my $reserved = $plus ? 'EDF+C' : '';

    # Subject / recording identification. Note: the EDF+ spec names the first
    # 80-byte field "local patient identification"; we use 'subject' throughout
    # the API (matching MRI/BIDS convention) and write it into that field.
    # 'patient' is still accepted as a deprecated alias for backward compat.
    my $subject_opt = defined $o{subject}   ? $o{subject}   : $o{patient};
    my $subject_rec = defined $rec->{subject} ? $rec->{subject} : $rec->{patient};

    my ($subject_fld, $recording_fld);
    if ($plus) {
        $subject_fld = _edfplus_subject(
            code => (defined $subject_opt ? $subject_opt : $subject_rec),
            sex  => $o{sex}, birthdate => $o{birthdate}, name => $o{name},
        );
        $recording_fld = _edfplus_recording(
            _ddmmmyyyy($have_start, $sY, $sM, $sD),
            admin => (defined $o{recording} ? $o{recording} : $rec->{recording}),
            tech  => $o{technician},
            equip => _equip_string($rec, $o{equipment}),
        );
    } else {
        $subject_fld   = defined $subject_opt  ? $subject_opt  : ($subject_rec      // '');
        $recording_fld = defined $o{recording} ? $o{recording} : ($rec->{recording} // '');
    }

    my $hdr = '';
    $hdr .= _ascii('0', 8);                                   # version
    $hdr .= _ascii($subject_fld,   80);
    $hdr .= _ascii($recording_fld, 80);
    $hdr .= _ascii($sdate, 8);
    $hdr .= _ascii($stime, 8);
    $hdr .= _ascii($hbytes, 8);
    $hdr .= _ascii($reserved, 44);
    $hdr .= _ascii($n_rec, 8);
    $hdr .= _num($rdur, 8);
    $hdr .= _ascii($ns_total, 4);

    # per-signal fields (all signals for a given field are contiguous)
    my @sig_label = (map { _ascii($labels[$_], 16) } 0 .. $n_ch - 1);
    my @sig_trans = (_ascii('', 80)) x $n_ch;
    my @sig_dim   = (map { _ascii($dim_ch[$_], 8) } 0 .. $n_ch - 1);
    my @sig_pmin  = (map { _num($pmin[$_], 8) } 0 .. $n_ch - 1);
    my @sig_pmax  = (map { _num($pmax[$_], 8) } 0 .. $n_ch - 1);
    my @sig_dmin  = (($DMIN eq -32768 ? _ascii('-32768',8) : _num($DMIN,8))) x $n_ch;
    my @sig_dmax  = (_num($DMAX, 8)) x $n_ch;
    my @sig_pref  = (_ascii('', 80)) x $n_ch;
    my @sig_nspr  = (map { _ascii($spr, 8) } 0 .. $n_ch - 1);
    my @sig_rsv   = (_ascii('', 32)) x $n_ch;

    if ($plus) {
        push @sig_label, _ascii('EDF Annotations', 16);
        push @sig_trans, _ascii('', 80);
        push @sig_dim,   _ascii('', 8);
        push @sig_pmin,  _num(-1, 8);
        push @sig_pmax,  _num( 1, 8);
        push @sig_dmin,  _ascii('-32768', 8);
        push @sig_dmax,  _num(32767, 8);
        push @sig_pref,  _ascii('', 80);
        push @sig_nspr,  _ascii($annot_spr, 8);
        push @sig_rsv,   _ascii('', 32);
    }

    $hdr .= join('', @sig_label, @sig_trans, @sig_dim, @sig_pmin, @sig_pmax,
                     @sig_dmin,  @sig_dmax,  @sig_pref, @sig_nspr, @sig_rsv);

    croak sprintf("write_edf: internal header size %d != %d", length($hdr), $hbytes)
        if length($hdr) != $hbytes;

    # --- write -------------------------------------------------------------
    open my $fh, '>:raw', $path or croak "write_edf: cannot open $path: $!";
    print {$fh} $hdr;

    for my $r (0 .. $n_rec - 1) {
        my $lo = $r * $spr;
        my $hi = $lo + $spr - 1;
        for my $c (0 .. $n_ch - 1) {
            print {$fh} ${ $dig[$c]->slice("$lo:$hi")->short->sever->get_dataref };
        }
        if ($plus) {
            my $a = $rec_tal[$r];
            $a .= "\x00" x ($annot_spr * 2 - length $a);
            print {$fh} $a;
        }
    }
    close $fh or croak "write_edf: close failed: $!";
    return $path;
}

# ---------------------------------------------------------------------------
# read_edf($path, %opt) -> record hashref matching the read_nk contract:
#   { data => PDL[n_ch,n_samp] float (uV), fs, labels=>[...],
#     t_start => "YYYY-MM-DD HH:MM:SS", events=>[ {onset=>sec, label=>str}, ... ],
#     device, recording, edf_type }
# Annotation ("EDF Annotations") signals are parsed into {events} and kept out
# of {data}. Assumes all non-annotation signals share one sample rate; if not,
# it carps and uses the first signal's rate.
# ---------------------------------------------------------------------------
my %_MON = (JAN=>1,FEB=>2,MAR=>3,APR=>4,MAY=>5,JUN=>6,
            JUL=>7,AUG=>8,SEP=>9,OCT=>10,NOV=>11,DEC=>12);

sub read_edf {
    my ($path, %opt) = @_;
    open my $fh, '<:raw', $path or croak "read_edf: $path: $!";

    my $hdr; read($fh, $hdr, 256) == 256 or croak "read_edf: short header";
    my $trim = sub { my $s = shift; $s =~ s/\s+$//; $s };
    my $recording_fld = $trim->(substr($hdr,  88, 80));
    my $startdate     = $trim->(substr($hdr, 168,  8));   # dd.mm.yy
    my $starttime     = $trim->(substr($hdr, 176,  8));   # hh.mm.ss
    my $reserved      = substr($hdr, 192, 44);
    my $nrec          = 0 + $trim->(substr($hdr, 236, 8));
    my $recdur        = 0 + $trim->(substr($hdr, 244, 8));
    my $ns            = 0 + $trim->(substr($hdr, 252, 4));
    croak "read_edf: bad signal count" if $ns <= 0;
    my $edf_type = ($reserved =~ /EDF\+([CD])/) ? "EDF+$1" : 'EDF';

    my $sh; read($fh, $sh, $ns*256) == $ns*256 or croak "read_edf: short signal header";
    my $col = sub {                       # extract field i of width w from block base
        my ($base, $w, $i) = @_; $trim->(substr($sh, $base + $i*$w, $w));
    };
    my @label  = map { $col->(0,          16, $_) } 0..$ns-1;
    my @pdim   = map { $trim->($col->(96*$ns, 8, $_)) } 0..$ns-1;   # 16+80
    my @pmin   = map { $col->(104*$ns,     8, $_) } 0..$ns-1;   # 16+80+8
    my @pmax   = map { $col->(112*$ns,     8, $_) } 0..$ns-1;
    my @dmin   = map { $col->(120*$ns,     8, $_) } 0..$ns-1;
    my @dmax   = map { $col->(128*$ns,     8, $_) } 0..$ns-1;
    my @spr    = map { 0 + $col->(216*$ns, 8, $_) } 0..$ns-1;   # samples per record
    # field offsets: label16 transducer80 physdim8 pmin8 pmax8 dmin8 dmax8
    #                prefilter80 spr8 reserved32  -> cumulative *$ns blocks

    # classify annotation vs data signals
    my (@sig, @annot);
    for my $s (0..$ns-1) {
        if ($label[$s] eq 'EDF Annotations') { push @annot, $s } else { push @sig, $s }
    }
    croak "read_edf: no data signals" unless @sig;

    # sample rate: from first data signal; warn if others differ
    my $fs = $recdur ? $spr[$sig[0]] / $recdur : $spr[$sig[0]];
    for my $s (@sig) {
        my $f = $recdur ? $spr[$s]/$recdur : $spr[$s];
        if ($f != $fs) {
            carp "read_edf: signal '$label[$s]' rate ${f}Hz != ${fs}Hz; ".
                 "read_nk contract needs uniform rate (using ${fs}Hz)";
            last;
        }
    }

    # read all data records; concatenate each signal's raw int16 bytes (no
    # per-sample Perl lists), annotation bytes separately.
    my $rec_samples = 0; $rec_samples += $_ for @spr;
    my @sigbytes = ('') x $ns;
    my @annbytes;
    for my $r (0 .. $nrec-1) {
        my $rb; my $got = read($fh, $rb, $rec_samples*2);
        last unless defined $got && $got >= $rec_samples*2;
        my $p = 0;
        for my $s (0..$ns-1) {
            my $cnt = $spr[$s];
            if ($label[$s] eq 'EDF Annotations') {
                push @annbytes, substr($rb, $p*2, $cnt*2);
            } else {
                $sigbytes[$s] .= substr($rb, $p*2, $cnt*2);
            }
            $p += $cnt;
        }
    }
    close $fh;

    # build [n_ch, n_samp] float piddle in uV; interpret each signal's bytes
    # directly as int16 LE via get_dataref (~90x faster than unpack->pdl), scale.
    #
    # The physical dimension is PER SIGNAL and must be honoured. A DC channel is
    # written in mV (it cannot be written in µV: EDF's physical_min is 8 ASCII
    # chars and "-12002900" is 9), so its physical values come out 1000x smaller
    # than everything else. Reading them as µV -- which this used to do -- turns a
    # 4 V trigger into a 4 mV one. Normalise everything to µV here, and report the
    # dimension each signal was stored in via {units}.
    my %UV_PER = ('uV' => 1, 'uv' => 1, 'µV' => 1, 'mV' => 1_000, 'mv' => 1_000);
    my $n_ch   = scalar @sig;
    my $n_samp = length($sigbytes[$sig[0]]) / 2;
    my $data   = zeroes(float, $n_ch, $n_samp);
    my (@labels, @units);
    for my $i (0 .. $#sig) {
        my $s = $sig[$i];
        push @labels, $label[$s];
        push @units,  ($pdim[$s] // 'uV');
        my $uvp  = $UV_PER{ $pdim[$s] // 'uV' } // 1;      # unknown dim -> as-is
        my $dd   = ($dmax[$s] - $dmin[$s]);
        my $gain = $dd ? ($pmax[$s] - $pmin[$s]) / $dd : 1;
        my $col  = zeroes(short, length($sigbytes[$s]) / 2);
        ${ $col->get_dataref } = $sigbytes[$s];
        $col->upd_data;
        $data->slice("($i),:") .=
            (($col->float - $dmin[$s]) * $gain + $pmin[$s]) * $uvp;
    }

    # parse EDF+ annotations -> events
    my @events;
    if (@annbytes) {
        my $blob = join('', @annbytes);
        for my $tal (split /\x00/, $blob) {
            next unless length $tal;
            my @parts = split /\x14/, $tal, -1;
            my $timing = shift @parts;
            next unless defined $timing && $timing =~ /^([+-][\d.]+)/;
            my $onset = 0 + $1;
            for my $txt (@parts) {
                next unless defined $txt && length $txt;
                my $label = eval { decode('UTF-8', $txt, Encode::FB_CROAK) };
                $label = $txt if $@;                 # fall back to raw bytes
                push @events, { onset => $onset, label => $label };
            }
        }
    }

    # t_start: prefer full year from "Startdate DD-MMM-YYYY" in recording field
    my ($Y,$Mo,$D);
    if ($recording_fld =~ /Startdate\s+(\d{2})-([A-Z]{3})-(\d{4})/i) {
        ($D,$Mo,$Y) = ($1, $_MON{uc $2} // 1, $3);
    } elsif ($startdate =~ /(\d{2})\.(\d{2})\.(\d{2})/) {
        ($D,$Mo,my $yy) = ($1,$2,$3);
        $Y = $yy >= 85 ? 1900+$yy : 2000+$yy;        # EDF clipping rule
    }
    my ($H,$Mi,$Se) = $starttime =~ /(\d{2})\.(\d{2})\.(\d{2})/ ? ($1,$2,$3) : (0,0,0);
    my $t_start = (defined $Y)
        ? sprintf('%04d-%02d-%02d %02d:%02d:%02d', $Y,$Mo,$D,$H,$Mi,$Se)
        : undef;

    return {
        data      => $data,          # ALWAYS uV -- mV signals are normalised
        fs        => $fs,
        labels    => \@labels,
        units     => \@units,        # the dimension each signal was STORED in
        t_start   => $t_start,
        events    => \@events,
        recording => $recording_fld,
        edf_type  => $edf_type,
        n_records => $nrec,
        rec_dur   => $recdur,
    };
}

# ---------------------------------------------------------------------------
# clean_edf_label($raw) -> montage-ready electrode name
#
# EDF+ signal labels are "TYPE electrode" (e.g. "EEG Fp1-Ref", "POL DC01"), and
# referential montages append a "-Ref" marker. Downstream text formats such as
# BESA .mul use a whitespace-delimited label row, so any space inside a label
# splits into an extra token and desyncs the reader from the channel count.
# This normalises an EDF+ label to a single clean token:
#   * strip a leading recognised EDF+ type code   (EEG Fp1-Ref -> Fp1-Ref)
#   * drop a trailing referential marker           (Fp1-Ref    -> Fp1)
#   * map Nihon Kohden "$A1" reference channels    ($A1        -> A1_ref)
#   * collapse any residual whitespace to '_'
# Unknown prefixes are left untouched; undef passes through unchanged.
# ---------------------------------------------------------------------------
my %_EDF_TYPE = map { $_ => 1 }
    qw(EEG EOG ECG EKG EMG ERG MEG POL DC SaO2 SpO2 Status Trig Temp Resp);

sub clean_edf_label {
    my ($lab) = @_;
    return $lab unless defined $lab;
    $lab =~ s/^\s+//; $lab =~ s/\s+$//;
    $lab = $2 if $lab =~ /^(\S+)\s+(.+)$/ && $_EDF_TYPE{$1};   # strip EDF+ type code
    $lab =~ s/-Ref$//i;                                        # drop referential marker
    $lab =~ s/^\$(.+)$/${1}_ref/;                              # NK ref chan: $A1 -> A1_ref
    $lab =~ s/\s+/_/g;                                         # any residual whitespace
    return $lab;
}

1;

=head1 SEE ALSO

L<PDL::EEG::IO::NihonKohden>, EDF/EDF+ spec (Kemp & Olivan 2003).

=cut
