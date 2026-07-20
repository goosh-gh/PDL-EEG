package PDL::EEG::IO::ASA;

use strict;
use warnings;
use Carp qw(croak carp);
use PDL;

our $VERSION = '0.01';

use Exporter 'import';
our @EXPORT_OK = qw(read_elc parse_ELEC_POS3D_ASA_4AdventCalendar);

=head1 NAME

PDL::EEG::IO::ASA - Read ASA electrode position files (.elc)

=head1 SYNOPSIS

  use PDL::EEG::IO::ASA qw(read_elc);

  my $mon = read_elc('standard_1020.elc');

  $mon->{coords}      # PDL [3, N] float, native unit, order = $mon->{labels}
  $mon->{labels}      # arrayref of N electrode names (trimmed)
  $mon->{pos}{Cz}     # { x=>.., y=>.., z=>.., index=>.. }  name lookup
  $mon->{unit}        # 'mm'   (from UnitPosition)
  $mon->{reference}   # 'avg'  (from ReferenceLabel)
  $mon->{comment}     # first line (FileComment)
  $mon->{n}           # N (actual electrode count)
  $mon->{fiducials}   # { LPA=>[x,y,z], RPA=>[x,y,z], Nz=>[x,y,z] } if present

  # 3D-ready (3,N) piddle straight into a viewer:
  my ($x, $y, $z) = $mon->{coords}->dog;   # three 1-D piddles

=head1 DESCRIPTION

Reads the ASA C<.elc> electrode-coordinate text format (as shipped with
mne-python, e.g. C<standard_1020.elc>). An ASA file has an C<xyz> coordinate
block in its upper half and a matching label block below, with
C<Positions[k]> corresponding to C<Labels[k]>:

  # ASA electrode file
  ReferenceLabel   avg
  UnitPosition     mm
  NumberPositions= 97
  Positions
  -86.0761 -19.9897 -47.9860
  ...
  Labels
  LPA
  ...

This reader is the library form of the parser described in the PDL Advent
Calendar 2024 (Day 12), rebuilt around three points:

=over 4

=item * B<Robust whitespace.> Lines are chomped (CR stripped) and leading
whitespace removed before splitting, so indented coordinate blocks parse
correctly. The Advent C<parse_ASCII> split on C</\s+/> without either, which
on an indented C<.elc> shifts C<x> into an empty field and silently drops
C<z>; the C<< map $_||0 >> guard only masked it by turning C<x> into 0.

=item * B<Vectorised coordinates.> The coordinate block is parsed into one
flat list and reshaped to a C<(3,N)> piddle in a single step, rather than
looping row-by-row into a pre-zeroed piddle.

=item * B<Block location by keyword.> C<Positions>/C<Labels> are found by
their marker lines rather than a fixed header count, so extra metadata lines
do not throw off the offsets.

=back

The declared C<NumberPositions=> is cross-checked against the actual count
and a mismatch is C<carp>ed (not fatal), matching the Advent behaviour.

=head1 FUNCTIONS

=head2 read_elc($path)

Returns a hashref describing the montage (keys listed in the SYNOPSIS).
C<coords> is a C<(3,N)> float piddle in file order; C<labels> is the parallel
name list; C<pos> is a name-keyed lookup. Croaks on a malformed file.

=cut

# Canonical fiducial names. ASA files lead with the three head landmarks.
my %FID_CANON = (
    LPA     => 'LPA', RPA     => 'RPA',
    NZ      => 'Nz',  NAS     => 'Nz', NASION => 'Nz',
);

sub _fid_canon {
    my ($name) = @_;
    return $FID_CANON{ uc $name };
}

sub read_elc {
    my ($path) = @_;
    croak "read_elc: no file given" unless defined $path;
    croak "read_elc: '$path' is empty or missing" unless -s $path;

    open my $fh, '<', $path or croak "read_elc: cannot open '$path': $!";
    my @lines = <$fh>;
    close $fh;
    for (@lines) { s/\r?\n?$//; }          # chomp + strip CR (CRLF-safe)

    my %meta;
    $meta{comment} = @lines ? $lines[0] : '';
    for my $ln (@lines) {
        $meta{reference} = $1 if $ln =~ /^ReferenceLabel\s+(\S+)/;
        $meta{unit}      = $1 if $ln =~ /^UnitPosition\s+(\S+)/;
        $meta{declared}  = $1 if $ln =~ /^NumberPositions=?\s+(\d+)/;
    }

    my ($pos_i) = grep { $lines[$_] =~ /^Positions\b/ } 0 .. $#lines;
    my ($lab_i) = grep { $lines[$_] =~ /^Labels\b/    } 0 .. $#lines;
    croak "read_elc: no 'Positions' block in '$path'" unless defined $pos_i;
    croak "read_elc: no 'Labels' block in '$path'"    unless defined $lab_i;
    croak "read_elc: 'Labels' precedes 'Positions' in '$path'" if $lab_i <= $pos_i;

    # Coordinate rows: between Positions and Labels markers.
    my @nums;
    my $ncoord = 0;
    for my $i ($pos_i + 1 .. $lab_i - 1) {
        my $ln = $lines[$i];
        $ln =~ s/^\s+//; $ln =~ s/\s+$//;
        next if $ln eq '';
        my @f = split /\s+/, $ln;
        croak "read_elc: coordinate line ", $i + 1, " has ", scalar @f,
              " field(s), need >= 3: '$lines[$i]'" if @f < 3;
        push @nums, @f[0, 1, 2];
        $ncoord++;
    }
    croak "read_elc: no coordinates found in '$path'" unless $ncoord;

    # Labels: as many rows as coordinates, first token of each.
    my @labels;
    for my $i ($lab_i + 1 .. $#lines) {
        last if @labels >= $ncoord;
        my $ln = $lines[$i];
        $ln =~ s/^\s+//; $ln =~ s/\s+$//;
        next if $ln eq '';
        my @f = split /\s+/, $ln;
        push @labels, $f[0];
    }
    croak "read_elc: found $ncoord coordinates but only ", scalar @labels,
          " label(s) in '$path'" if @labels != $ncoord;

    if (defined $meta{declared} && $meta{declared} != $ncoord) {
        carp "read_elc: NumberPositions=$meta{declared} but $ncoord "
           . "coordinate rows in '$path'";
    }

    # Vectorised (3,N): flat [x0,y0,z0,x1,...] -> reshape.
    my $coords = pdl(float, \@nums)->reshape(3, $ncoord);

    my (%pos, %fid);
    for my $i (0 .. $ncoord - 1) {
        my ($x, $y, $z) = map { $coords->at($_, $i) } 0 .. 2;
        $pos{ $labels[$i] } = { x => $x, y => $y, z => $z, index => $i };
        if (my $c = _fid_canon($labels[$i])) { $fid{$c} = [ $x, $y, $z ]; }
    }

    return {
        coords           => $coords,
        labels           => \@labels,
        pos              => \%pos,
        fiducials        => \%fid,
        unit             => $meta{unit}      // '',
        reference        => $meta{reference} // '',
        comment          => $meta{comment},
        n                => $ncoord,
        number_positions => $meta{declared},
        file             => $path,
    };
}

=head2 parse_ELEC_POS3D_ASA_4AdventCalendar($path)

Backward-compatible shim reproducing the four-value return of the PDL Advent
Calendar 2024 (Day 12) parser, for scripts written against it:

  my ($r_h, $r_epos, $labels, $coords) = parse_ELEC_POS3D_ASA_4AdventCalendar($file);

  $r_h              header hash: FileComment, ReferenceLabel, UnitPosition,
                    NumberPositions, plus per-electrode $r_h->{Cz}{x|y|z|DeviceCh}
  $r_epos           arrayref of { name, x, y, z } in file order
  $labels           arrayref of names, each prefixed with two spaces
                    (the Advent script's TriD label padding)
  $coords           (3,N) float piddle

Built on C<read_elc>, so it inherits the robust parse. Note: the Advent
original returned C<< $coords->using(0,1,2) >>; C<using> is not a core PDL
method, so this returns the plain C<(3,N)> piddle instead.

=cut

sub parse_ELEC_POS3D_ASA_4AdventCalendar {
    my ($path) = @_;
    my $m = read_elc($path);

    my @epos;
    my %h = (
        filename        => $path,
        FileComment     => $m->{comment},
        ReferenceLabel  => $m->{reference},
        UnitPosition    => $m->{unit},
        NumberPositions => (defined $m->{number_positions}
                              ? $m->{number_positions} : $m->{n}),
        N_Coords        => $m->{n},
    );
    for my $i (0 .. $m->{n} - 1) {
        my $name = $m->{labels}[$i];
        my ($x, $y, $z) = map { $m->{coords}->at($_, $i) } 0 .. 2;
        $epos[$i] = { name => $name, x => $x, y => $y, z => $z };
        $h{$name} = { DeviceCh => $i, x => $x, y => $y, z => $z };
    }
    my @labels = map { '  ' . $_->{name} } @epos;   # Advent 2-space padding

    return (\%h, \@epos, \@labels, $m->{coords});
}

=head1 SEE ALSO

L<PDL::EEG>, L<PDL::EEG::IO::NihonKohden>,
L<https://pdl.perl.org/advent/blog/2024/12/12/eeg/>

=head1 AUTHOR

goosh

=head1 LICENSE

Same terms as Perl itself.

=cut

1;
