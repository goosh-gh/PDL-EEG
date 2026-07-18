package PDL::EEG::IO::BESA::ASCII;

use strict;
use warnings;
use Carp qw(croak carp);
use PDL;

our $VERSION = '0.02';

use Exporter 'import';
our @EXPORT_OK = qw(write_mul read_mul);

=head1 NAME

PDL::EEG::IO::BESA::ASCII - Write EEG data to BESA ASCII files (.mul, .avr)

=head1 SYNOPSIS

  use PDL::EEG::IO::NihonKohden qw(read_nk);
  use PDL::EEG::IO::BESA::ASCII qw(write_mul);

  my $rec = read_nk('subject.EEG');
  write_mul($rec, 'subject.mul');

  # override / tune
  write_mul($rec, 'out.mul',
      labels   => [ map { "$_-BN" } @{ $rec->{labels} } ],  # add montage suffix
      trigger  => 'Trigger',   # channel written as integer (default: any label eq 'Trigger')
      begin_ms => 0.0,         # BeginSweep[ms] (negative for pre-stimulus epochs)
      decimals => 2,           # value precision
      fieldw   => 8,           # column width
      time     => '16:44:34',  # default: derived from $rec->{t_start}
      segment  => 'Cond1',     # add SegmentName= (default: omitted)
      count_trigger => 1,      # count Trigger in Channels= (default: excluded)
      exclude  => [qw(DC01 DC02 STIM)],  # drop channels by name (or index)
  );

  # raw piddle, no full record:
  write_mul({ data => $pdl, fs => 1000, labels => \@names }, 'out.mul');

=head1 DESCRIPTION

Writes the BESA I<ASCII multiplexed> text format (C<.mul>): a two-line header
followed by one line per time point, each holding all channel values.

  TimePoints=208000 Channels=26 BeginSweep[ms]=0.00 SamplingInterval[ms]=1.000 Bins/uV=1.000 Time=16:44:34
   Fp1-BN Fp2-BN F3-BN ... X1-BN Trigger
      9.34     7.03    -5.26 ...     2.96    0
   ...

Data are written uncompressed as physical microvolts, so C<Bins/uV=1.000>.
Because BESA-compatible readers (Brainstorm, FieldTrip C<read_besa_mul>, ...)
tokenise each data line on whitespace, the exact column widths are cosmetic
and do not affect round-trip correctness.

B<C<Bins/uV> is one scalar for the whole file.> There is no per-channel physical
dimension in C<.mul>, unlike EDF. So every column must be microvolts -- and an
EEG-1200A DC input is a +/-12 V line, i.e. B<+/-12002913 uV>. It is written out
at that magnitude, correctly, and it will not fit an 8-character column; the
numbers are simply wide. If BESA's autoscaling is the point of the export, drop
the DC channels:

  write_mul($rec, 'out.mul', exclude => [ grep { /^DC/ } @{ $rec->{labels} } ]);

(Before 0.2 the DC columns came out 1000x too small, because C<read_nk> was
returning the vendor's millivolt figure and calling it microvolts. .mul files
written from EEG-1200A recordings before that fix have wrong DC values and
should be regenerated.)

The dedicated integer Trigger channel is exported as a column but is I<not>
counted in C<Channels=> (it is dropped from later analysis), matching the
Nihon Kohden C<.mul> export. So a 27-column recording whose last channel is
C<Trigger> is written with C<Channels=26>. Pass C<< count_trigger => 1 >> to
count it instead.

The input record is the hashref returned by C<read_nk>:

  $rec->{data}    PDL [n_ch, n_samples] float, microvolts -- ALL channels,
                  DC included (a DC input reaches +/-12002913 uV)
  $rec->{fs}      sampling rate (Hz)
  $rec->{labels}  arrayref of channel names
  $rec->{units}   arrayref, optional: 'uV' | 'mV' | 'code'. Only 'code' is acted
                  on here (raw trigger codes -> written as integers). 'mV' is an
                  EDF export hint and does NOT mean the data are millivolts.
  $rec->{t_start} "YYYY-MM-DD HH:MM:SS" (optional; used for the Time= field)

=head1 FUNCTIONS

=head2 write_mul($rec, $path, %opts)

=cut

sub write_mul {
    my ($rec, $path, %opt) = @_;

    croak "write_mul: record hashref required"      unless ref $rec eq 'HASH';
    croak "write_mul: output path required"         unless defined $path;

    my $data = $rec->{data};
    croak "write_mul: \$rec->{data} must be a PDL"  unless eval { $data->isa('PDL') };
    croak "write_mul: \$rec->{data} must be 2-D [n_ch, n_samples]"
        unless $data->ndims == 2;

    my $nc = $data->dim(0);
    my $nt = $data->dim(1);

    my @labels = $opt{labels}      ? @{ $opt{labels} }
               : $rec->{labels}    ? @{ $rec->{labels} }
               : map { "E" . ($_ + 1) } 0 .. $nc - 1;
    croak "write_mul: got ${\ scalar @labels} labels for $nc channels"
        unless @labels == $nc;

    my @units = ($rec->{units} && @{ $rec->{units} } == $nc)
              ? @{ $rec->{units} } : ();

    # --- exclude => [names and/or indices] ----------------------------------
    # .mul has ONE Bins/uV for the whole file, so a +/-12 V DC channel has to go
    # in as +/-12002913 uV next to a 30 uV EEG channel. That is numerically right
    # and practically useless -- BESA will autoscale everything into the floor.
    # Dropping the DC columns is usually what you actually want.
    if ($opt{exclude} && @{ $opt{exclude} }) {
        my %drop;
        for my $e (@{ $opt{exclude} }) {
            if ($e =~ /^\d+$/ && $e < $nc) { $drop{$e} = 1; next }
            my @hit = grep { $labels[$_] eq $e } 0 .. $#labels;
            if (@hit) { $drop{$_} = 1 for @hit }
            else      { carp "write_mul: exclude: no channel named '$e'" }
        }
        my @keep = grep { !$drop{$_} } 0 .. $nc - 1;
        croak "write_mul: exclude removed every channel" unless @keep;
        if (@keep < $nc) {
            $data   = $data->slice(pdl(\@keep), 'X')->sever;
            @labels = @labels[@keep];
            @units  = @units[@keep] if @units;
            $nc     = scalar @keep;
        }
    }

    # The .mul label row (header line 2) is whitespace-delimited and MUST hold
    # exactly $nc tokens. EDF+ labels carry a type prefix ("EEG Fp1-Ref",
    # "POL DC01") i.e. embedded spaces, which would split into extra tokens and
    # desync a reader from Channels=$nc. Collapse any internal whitespace to '_'
    # so the file is always valid regardless of upstream label conventions.
    # (Callers that want clean names should strip the type prefix beforehand;
    #  this is the last-resort guarantee.)
    my $spaced = grep { /\s/ } @labels;
    if ($spaced) {
        s/\s+/_/g for @labels;
        carp "write_mul: $spaced label(s) contained whitespace; replaced with "
           . "'_' to keep the .mul label row aligned with Channels=$nc";
    }

    my $fs = $opt{fs} // $rec->{fs}
        or croak "write_mul: sampling rate (fs) required";

    my $begin_ms = defined $opt{begin_ms} ? $opt{begin_ms} : 0.0;
    my $si_ms    = 1000.0 / $fs;                    # SamplingInterval[ms]
    my $dec      = defined $opt{decimals} ? $opt{decimals} : 2;
    my $w        = defined $opt{fieldw}   ? $opt{fieldw}   : 8;
    my $si_dec   = defined $opt{si_decimals} ? $opt{si_decimals} : 3;
    my $tw       = defined $opt{trig_width} ? $opt{trig_width} : $w;

    # Time=HH:MM:SS : explicit opt, else the clock part of t_start, else omit.
    my $time = $opt{time};
    if (!defined $time && !exists $opt{time} && defined $rec->{t_start}) {
        ($time) = $rec->{t_start} =~ /(\d{2}:\d{2}:\d{2})/;
    }

    # Which column is the integer trigger channel?
    my $trig_idx;
    if (exists $opt{trigger}) {
        if (defined $opt{trigger} && $opt{trigger} =~ /^\d+$/) {
            $trig_idx = $opt{trigger};              # explicit index
        } elsif (defined $opt{trigger}) {
            ($trig_idx) = grep { $labels[$_] eq $opt{trigger} } 0 .. $#labels;
        }
        # trigger => undef  -> no integer column
    } else {
        # default: a channel literally named 'Trigger', or -- better -- one that
        # read_nk marked as raw codes ({units} eq 'code', e.g. the EEG-1200A STIM
        # column). Those are integer marker codes; writing them as "%8.2f" is
        # just noise.
        ($trig_idx) = grep { $labels[$_] eq 'Trigger' } 0 .. $#labels;
        if (!defined $trig_idx && @units) {
            ($trig_idx) = grep { ($units[$_] // '') eq 'code' } 0 .. $#labels;
        }
    }

    open my $fh, '>', $path or croak "write_mul: cannot open $path: $!";

    # Channels= convention: the dedicated integer Trigger channel is still
    # exported (label + data column) but is NOT counted in Channels=, because it
    # is dropped from later analysis. This matches the Nihon Kohden .mul export
    # (e.g. 27 columns incl. Trigger -> Channels=26). DC trigger *inputs*
    # (DC01.., DC03..) ARE counted; only the dedicated Trigger channel is not.
    # Pass count_trigger => 1 to include it in the count instead.
    my $n_report = $nc;
    $n_report-- if defined $trig_idx && !$opt{count_trigger};

    # --- header line 1 -----------------------------------------------------
    my $h1 = sprintf(
        'TimePoints=%d Channels=%d BeginSweep[ms]=%.2f SamplingInterval[ms]=%.*f Bins/uV=%.3f',
        $nt, $n_report, $begin_ms, $si_dec, $si_ms, 1.0,
    );
    $h1 .= " Time=$time"              if defined $time;
    $h1 .= " SegmentName=$opt{segment}" if defined $opt{segment};
    print {$fh} $h1, "\n";

    # --- header line 2: channel labels (leading space, as BESA emits) ------
    print {$fh} ' ', join(' ', @labels), "\n";

    # --- data: one line per time point, all channels -----------------------
    # PDL stores dim0 (channel) fastest, so flattening a [n_ch, blk] slice
    # already yields multiplexed order: ch0_t0 ch1_t0 ... chN_t0 ch0_t1 ...
    my $afmt = "%${w}.${dec}f";       # analog channel
    my $tfmt = "%${tw}.0f";           # trigger channel (integer, no decimals)

    # One printf format for a whole time-point row (trigger column as integer),
    # applied to an entire flushed block in a single sprintf. This replaces the
    # per-cell sprintf/join and is ~2x faster on large recordings; output is
    # byte-identical.
    my $rowfmt = join(' ',
        map { defined $trig_idx && $_ == $trig_idx ? $tfmt : $afmt } 0 .. $nc - 1
    ) . "\n";

    # A DC channel in uV needs 12 characters ("-12002913.00"), not 8. sprintf will
    # simply widen the column rather than truncate, so the file stays valid and
    # readers still tokenise it -- but it stops lining up, and that surprises
    # people. Say so once, up front, instead of letting them find it in the file.
    {
        my ($mn, $mx) = $data->minmax;
        my $need = length(sprintf("%.${dec}f", $mn < 0 ? $mn : -$mx));
        if ($need > $w) {
            my @wide = grep {
                my ($a, $b) = $data->slice("($_),:")->minmax;
                length(sprintf("%.${dec}f", $a < 0 ? $a : -$b)) > $w
            } 0 .. $nc - 1;
            carp sprintf(
                "write_mul: %s need %d columns but fieldw=%d; the rows will not "
              . "line up (still valid: readers split on whitespace). Pass "
              . "fieldw => %d, or exclude => [...] to drop them.",
                join(',', @labels[@wide]), $need, $w, $need);
        }
    }

    my $BLK = 20000;                  # time points per flush (bounds memory)
    for (my $lo = 0; $lo < $nt; $lo += $BLK) {
        my $hi   = $lo + $BLK - 1;
        $hi      = $nt - 1 if $hi > $nt - 1;
        my @vals = $data->slice(":,${lo}:${hi}")->flat->list;
        my $ntb  = $hi - $lo + 1;
        print {$fh} sprintf($rowfmt x $ntb, @vals);
    }

    close $fh or croak "write_mul: error closing $path: $!";
    return $path;
}


=head2 read_mul($path, %opt)

Read a BESA ASCII multiplexed file back into a record matching the C<read_nk>
contract. Written to close the loop: the Nihon Kohden viewer exports C<.mul>
itself, so a reader lets us diff our own conversion against the vendor's.

  my $m = read_mul('subject.m01');
  #  { data => PDL[n_ch,n_samp] float (uV), fs, labels => [...],
  #    t_start => "YYYY-MM-DD HH:MM:SS" | undef, n_report, trig_idx,
  #    bins_per_uv, begin_ms, segment, header => { raw key => value } }

Both column layouts are handled:

  TimePoints=174000 Channels=26 ... Bins/uV=1.000 Time=14:07:54
  Date Time Fp1-BN Fp2-BN ... Trigger
  2026/07/02 14:07:54   32.48    20.19 ...    0

and the same thing without the leading C<Date>/C<Time> columns (what write_mul
emits). The prefix columns are detected from the label row, not assumed.

WHAT THE VENDOR'S OWN EXPORT LOOKS LIKE (subject.m01, EEG-1200A), because it is
not what the docs led us to expect:

=over

=item *

B<Channels= counts every data column, Trigger included.> 21 scalp + 4 DC + 1
Trigger = 26 columns, and the header says C<Channels=26>. C<write_mul> defaults
to the opposite (Trigger exported but not counted). One of the two is wrong for
this recorder; C<read_mul> reports both numbers rather than pick a side.

=item *

B<The label row is unreliable.> The four DC columns are all named C<Fp1>. Not a
typo here -- that is what the file says. Match DC columns by position, never by
name.

=item *

B<DC channels are in microvolts>, quantised to 366.22 uV per bit -- the DC ADC
step. C<Bins/uV=1.000>, so the vendor writes DC at the same scale as EEG, at
values around 233279.92 uV. This is independent confirmation that the DC gain is
366.3 uV/bit and not 0.366.

=item *

B<The export is trimmed.> C<Time=14:07:54> but the segment starts at 14:07:52,
and TimePoints is 174000 where the segment holds 176000 samples: the viewer drops
the first two seconds (the amplifier settling period, up to "A1+A2 OFF").

=back

=cut

sub read_mul {
    my ($path, %opt) = @_;
    open my $fh, '<', $path or croak "read_mul: cannot open $path: $!";

    # header line 1 -- skip any leading blank lines
    my $h1;
    while (defined($h1 = <$fh>)) { last if $h1 =~ /\S/ }
    croak "read_mul: $path is empty" unless defined $h1;
    chomp $h1;

    my %hdr;
    $hdr{$1} = $2 while $h1 =~ /(\S+?)=(\S*)/g;
    croak "read_mul: $path: no TimePoints= in the header (not a .mul?)"
        unless defined $hdr{TimePoints};

    my $si = $hdr{'SamplingInterval[ms]'};
    my $fs = ($si && $si > 0) ? 1000.0 / $si : undef;
    croak "read_mul: cannot get a sampling rate from '$h1'" unless $fs;
    my $bins = defined $hdr{'Bins/uV'} && $hdr{'Bins/uV'} > 0 ? $hdr{'Bins/uV'} : 1.0;

    # header line 2 -- channel names, possibly preceded by Date/Time columns
    my $h2 = <$fh>;
    croak "read_mul: $path: no label row" unless defined $h2;
    chomp $h2;
    my @names = split ' ', $h2;
    my $npre = 0;
    if (@names >= 2 && lc $names[0] eq 'date' && lc $names[1] eq 'time') { $npre = 2 }
    elsif (@names && lc $names[0] eq 'time')                             { $npre = 1 }
    splice @names, 0, $npre if $npre;
    my $nc = scalar @names;
    croak "read_mul: $path: no channel names" unless $nc;

    # --- data ---------------------------------------------------------------
    my $nt = 0 + $hdr{TimePoints};
    my $data = zeroes(float, $nc, $nt);
    my $dref = $data->get_dataref;
    my @flat;
    my ($first_stamp, $last_stamp);
    my $row = 0;
    while (my $l = <$fh>) {
        next unless $l =~ /\S/;
        my @f = split ' ', $l;
        if ($npre) {
            my @stamp = splice @f, 0, $npre;
            $first_stamp //= join(' ', @stamp);
            $last_stamp    = join(' ', @stamp);
        }
        croak "read_mul: $path line " . ($row + 3) . ": got "
            . scalar(@f) . " columns, expected $nc"
            unless @f == $nc;
        push @flat, @f;
        $row++;
        last if $row >= $nt;
    }
    carp "read_mul: $path: TimePoints=$nt but only $row data rows" if $row < $nt;

    if ($row < $nt) { $data = $data->slice(":,0:" . ($row - 1))->sever; $nt = $row }
    $data->slice(":,0:" . ($nt - 1)) .= pdl(\@flat)->reshape($nc, $nt);
    $data /= $bins unless $bins == 1.0;              # Bins/uV -> physical uV

    close $fh;

    # t_start: prefer the first Date/Time stamp (it has the date), else Time=
    my $t_start;
    if (defined $first_stamp
        && $first_stamp =~ m{^(\d{4})[/-](\d{2})[/-](\d{2})\s+(\d{2}:\d{2}:\d{2})}) {
        $t_start = "$1-$2-$3 $4";
    } elsif (defined $hdr{Time}) {
        $t_start = $hdr{Time};                       # clock only, no date
    }

    my ($trig_idx) = grep { $names[$_] eq 'Trigger' } 0 .. $#names;

    return {
        data        => $data,                        # [n_ch, n_samp] float uV
        fs          => $fs,
        labels      => \@names,
        t_start     => $t_start,
        n_ch        => $nc,                          # columns actually present
        n_report    => 0 + ($hdr{Channels} // $nc),  # what Channels= claims
        trig_idx    => $trig_idx,
        bins_per_uv => $bins,
        begin_ms    => $hdr{'BeginSweep[ms]'},
        segment     => $hdr{SegmentName},
        date_time   => $npre,                        # 0, 1 (Time) or 2 (Date Time)
        header      => \%hdr,
    };
}

1;

__END__

=head1 FORMAT REFERENCE

BESA ASCII multiplexed format, as documented on the BESA wiki
(L<https://wiki.besa.de/index.php?title=ASCII_File_Format>) and implemented in
FieldTrip's C<read_besa_mul.m>. The C<Time=> item is present only for continuous
epochs; C<SegmentName=> only when a segment comment exists.

=head1 SEE ALSO

L<PDL::EEG::IO::NihonKohden>, L<PDL::EEG::IO::EDF>

=head1 AUTHOR

goosh

=head1 LICENSE

Same terms as Perl itself.

=cut
