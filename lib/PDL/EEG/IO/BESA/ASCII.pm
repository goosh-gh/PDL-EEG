package PDL::EEG::IO::BESA::ASCII;

use strict;
use warnings;
use Carp qw(croak carp);
use PDL;

our $VERSION = '0.01';

use Exporter 'import';
our @EXPORT_OK = qw(write_mul);

=head1 NAME

PDL::EEG::IO::BESA::ASCII - Write EEG data to BESA ASCII files (.mul, .avr)

=head1 SYNOPSIS

  use PDL::EEG::IO::NihonKohden qw(read_nk);
  use PDL::EEG::IO::BESA::ASCII qw(write_mul);

  my $rec = read_nk('patient.EEG');
  write_mul($rec, 'patient.mul');

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

The dedicated integer Trigger channel is exported as a column but is I<not>
counted in C<Channels=> (it is dropped from later analysis), matching the
Nihon Kohden C<.mul> export. So a 27-column recording whose last channel is
C<Trigger> is written with C<Channels=26>. Pass C<< count_trigger => 1 >> to
count it instead.

The input record is the hashref returned by C<read_nk>:

  $rec->{data}    PDL [n_ch, n_samples] float, microvolts
  $rec->{fs}      sampling rate (Hz)
  $rec->{labels}  arrayref of channel names
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
        # default: auto-detect a channel literally named 'Trigger'
        ($trig_idx) = grep { $labels[$_] eq 'Trigger' } 0 .. $#labels;
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
