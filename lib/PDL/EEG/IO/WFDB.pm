package PDL::EEG::IO::WFDB;

# Read PhysioNet WFDB records with PDL: the header (.hea), format-16 signals in a
# .dat file or in a Matlab v4 .mat file (as in the PhysioNet / Computing in
# Cardiology Challenge 2018), and MIT-format annotation files (.atr, .arousal, ...).

use strict;
use warnings;
use Carp qw(croak);
use File::Basename qw(dirname);
use File::Spec;
use PDL;
use Exporter 'import';

our @EXPORT_OK = qw(read_wfdb read_wfdb_annotations wfdb_stages_events);

# ---------------------------------------------------------------------------
# $rec = read_wfdb($stem, %opt)
#
#   $stem : record path without extension (/path/tr03-0005); the .hea is
#           "$stem.hea" unless hea => FILE; the signal file is the one the
#           header names, in the header's directory, unless dat => FILE
#   data  => 1   also build ->{data} (nch, nt) double in physical units
#
# Returns a hash reference:
#   raw (nch, nt) short   the samples as stored
#   fs, nt, name, labels [..], units [..], gain [..], base [..]
#   chan => sub ($label, from => s0, to => s1)   one channel in physical
#           units (double), samples s0 .. s1-1 (default the whole record)
#
# Only format 16 (16-bit little-endian, channels interleaved) in one signal
# file is read. physical = (raw - baseline) / gain; a gain of 0 means 200
# ADC units per unit (WFDB default); the baseline defaults to the ADC zero.
# A signal file cut short (e.g. a partial download) is read up to its last
# whole sample, with a warning.
# ---------------------------------------------------------------------------
sub read_wfdb {
    my ($stem, %o) = @_;
    my $hea = $o{hea} // "$stem.hea";
    open my $h, '<', $hea or croak "read_wfdb: $hea: $!";
    my @l = grep { /\S/ && !/^#/ } <$h>;
    close $h;
    chomp @l;
    my ($name, $nsig, $fsf, $nt) = split ' ', shift @l;
    croak "read_wfdb: $hea: multi-segment records are not read" if $name =~ m{/};
    my ($fs) = ($fsf // 250) =~ /^([\d.eE+]+)/;
    croak "read_wfdb: $hea: $nsig signals but " . scalar(@l) . " signal lines" if @l < $nsig;
    my (@lab, @unit, @gain, @base, %file, $off);
    for my $i (0 .. $nsig - 1) {
        my @f = split ' ', $l[$i];
        $file{ $f[0] } = 1;
        my ($fmt) = $f[1] =~ /^(\d+)/;
        croak "read_wfdb: $hea: signal format $f[1] (only 16 is read)" unless $fmt == 16;
        ($off) = $f[1] =~ /\+(\d+)/;
        $off //= 0;
        my ($g, $bl, $u) = ($f[2] // '0') =~ /^([-\d.eE+]+)(?:\(([-\d]+)\))?(?:\/(\S+))?$/
            or croak "read_wfdb: $hea: gain field '$f[2]'";
        $g = 200 if $g == 0;
        push @gain, $g;
        push @base, defined $bl ? $bl : ($f[4] // 0);
        push @unit, $u // 'mV';
        push @lab,  join(' ', @f[8 .. $#f]);
    }
    croak "read_wfdb: $hea: signals in more than one file (" . join(', ', sort keys %file) . ")" if keys %file > 1;
    my $sig = $o{dat} // File::Spec->catfile(dirname($hea), (keys %file)[0]);

    open my $m, '<:raw', $sig or croak "read_wfdb: $sig: $!";
    my $size = -s $sig;
    if ($sig =~ /\.mat$/i) {                            # Matlab v4: type, mrows, ncols, imagf, namlen, name
        my $hdr;
        read($m, $hdr, 20) == 20 or croak "read_wfdb: $sig: short header";
        my ($type, $mrows, $ncols) = unpack('l<3', $hdr);
        croak "read_wfdb: $sig: $mrows rows, the header says $nsig signals" unless $mrows == $nsig;
        $nt = $ncols if !defined $nt || $nt > $ncols;
    }
    my $n = int(($size - $off) / (2 * $nsig));
    $nt //= $n;
    if ($n < $nt) { warn "read_wfdb: $sig holds $n of $nt samples\n"; $nt = $n }
    seek($m, $off, 0);
    my $raw = zeroes(short, $nsig, $nt);
    my $buf;
    read($m, $buf, 2 * $nsig * $nt) == 2 * $nsig * $nt or croak "read_wfdb: $sig: short read";
    close $m;
    ${ $raw->get_dataref } = $buf;
    undef $buf;
    $raw->upd_data;
    $raw->bswap2 if unpack('S', pack('s<', 1)) != 1;     # the file is little-endian

    my %ix;
    @ix{@lab} = 0 .. $#lab;
    my $chan = sub {
        my ($lb, %c) = @_;
        my $i = $ix{$lb};
        croak "read_wfdb: no channel '$lb' (have: @lab)" unless defined $i;
        my ($s0, $s1) = ($c{from} // 0, ($c{to} // $nt) - 1);
        return ($raw->slice("($i),$s0:$s1")->double - $base[$i]) / $gain[$i];
    };
    my $rec = { raw => $raw, fs => $fs, labels => \@lab, units => \@unit, gain => \@gain, base => \@base,
                nt => $nt, name => $name, chan => $chan };
    $rec->{data} = ($raw->double - pdl(@base)->dummy(1)) / pdl(@gain)->dummy(1) if $o{data};
    return $rec;
}

# ---------------------------------------------------------------------------
# $ann = read_wfdb_annotations($file)   [ [sample, code, aux text], ... ]
#
# MIT format: 16-bit little-endian words, code A = word >> 10, interval
# I = word & 1023.  A = 59 (SKIP): the next two words are a signed 32-bit
# interval, high word first; A = 63 (AUX): I bytes of text follow (padded to
# even) and belong to the annotation before; A = 60, 61, 62 (NUM, SUB, CHN)
# carry I as a value; A = 0 with I = 0 ends the file.  Any other A is an
# annotation at time += I.  Annotations with code 0 are dropped.
# ---------------------------------------------------------------------------
sub read_wfdb_annotations {
    my ($file) = @_;
    open my $f, '<:raw', $file or croak "read_wfdb_annotations: $file: $!";
    local $/;
    my $buf = <$f>;
    close $f;
    my @w = unpack('v*', $buf);
    my ($t, $i, @ann) = (0, 0);
    while ($i < @w) {
        my $A = $w[$i] >> 10;
        my $I = $w[$i] & 1023;
        $i++;
        if ($A == 0 && $I == 0) { last }
        elsif ($A == 59) {
            my $v = ($w[$i] << 16) + $w[$i + 1];
            $v -= 4294967296 if $v >= 2147483648;       # signed: the Challenge 2018 files use -1
            $t += $v;
            $i += 2;
        }
        elsif ($A == 63) {
            my $s = substr($buf, 2 * $i, $I);
            $s =~ s/\0+$//;
            $ann[-1][2] = $s if @ann;
            $i += int(($I + 1) / 2);
        }
        elsif ($A >= 60 && $A <= 62) { }
        else { $t += $I; push @ann, [ $t, $A, '' ] }
    }
    return [ grep { $_->[1] != 0 } @ann ];
}

# ---------------------------------------------------------------------------
# ($stage, $names, $events, $other) = wfdb_stages_events($ann, $nt)
#
#   $ann : read_wfdb_annotations output; sleep stages are aux texts at stage
#          changes, events are "(name" at their start and "name)" at their
#          end (as in the Challenge 2018 .arousal files)
#   $stage  (nt) byte: 0 unscored, 1 W, 2 N1, 3 N2, 4 N3, 5 R
#   $names  { code => name }
#   $events [ [name, start sample, end sample], ... ]
#   $other  { text => count } texts that are neither a stage nor an event
#           (lines starting with ## are header notes and are skipped)
#
# Stage texts (any case): W wake S0 -> W; N1 nonrem1 S1 -> N1; N2 nonrem2
# S2 -> N2; N3 N4 nonrem3 nonrem4 S3 S4 -> N3; R REM -> R; undefined ? U
# unscored -> unscored.
# ---------------------------------------------------------------------------
sub _stage_code {
    my $x = lc shift;
    return 1 if $x =~ /^(w|wake|s0)$/;
    return 2 if $x =~ /^(n1|nonrem1|s1)$/;
    return 3 if $x =~ /^(n2|nonrem2|s2)$/;
    return 4 if $x =~ /^(n3|n4|nonrem3|nonrem4|s3|s4)$/;
    return 5 if $x =~ /^(r|rem)$/;
    return 0 if $x =~ /^(undefined|\?|unscored|u)$/;
    return undef;
}

sub wfdb_stages_events {
    my ($ann, $nt) = @_;
    my %name = (0 => 'unscored', 1 => 'W', 2 => 'N1', 3 => 'N2', 4 => 'N3', 5 => 'R');
    my $st = zeroes(byte, $nt);
    my (%other, @chg, %open, @ev);
    for my $an (@$ann) {
        my $c = _stage_code($an->[2]);
        if (defined $c) { push @chg, [ $an->[0], $c ] }
        elsif ($an->[2] !~ /^\(|\)$|^##/) { $other{ $an->[2] eq '' ? "(code $an->[1], no text)" : $an->[2] }++ }
    }
    for my $k (0 .. $#chg) {
        my $s0 = $chg[$k][0];
        my $s1 = $k < $#chg ? $chg[$k + 1][0] - 1 : $nt - 1;
        next if $s0 >= $nt;
        $s1 = $nt - 1 if $s1 >= $nt;
        $st->slice("$s0:$s1") .= $chg[$k][1];
    }
    for my $an (@$ann) {
        my $x = $an->[2];
        if ($x =~ /^\((.+)$/) { $open{$1} = $an->[0] }
        elsif ($x =~ /^(.+)\)$/ && defined $open{$1}) { push @ev, [ $1, delete $open{$1}, $an->[0] ] }
    }
    return ($st, \%name, \@ev, \%other);
}

1;

__END__

=head1 NAME

PDL::EEG::IO::WFDB - read PhysioNet WFDB records (format 16, .dat or Matlab v4 .mat) and MIT annotation files

=head1 SYNOPSIS

    use PDL::EEG::IO::WFDB qw(read_wfdb read_wfdb_annotations wfdb_stages_events);

    my $rec = read_wfdb('/data/tr03-0005');            # tr03-0005.hea + the signal file it names
    my $ecg = $rec->{chan}->('ECG');                   # physical units
    my $abd = $rec->{chan}->('ABD', from => 0, to => 200 * 3600);

    my $ann = read_wfdb_annotations('/data/tr03-0005.arousal');
    my ($stage, $names, $events, $other) = wfdb_stages_events($ann, $rec->{nt});

=head1 DESCRIPTION

C<read_wfdb> reads a single-segment record whose signals are all in format 16
in one file: a WFDB C<.dat> file (offset from the format field, e.g. C<16+24>)
or a Matlab v4 C<.mat> file holding one int16 matrix (the layout of the
PhysioNet / Computing in Cardiology Challenge 2018). The samples stay int16 in
memory (2 bytes per value); C<chan> converts one channel to physical units
when asked. C<data =E<gt> 1> converts all channels at once.

C<read_wfdb_annotations> reads an MIT-format annotation file; the interval of
a SKIP is signed. C<wfdb_stages_events> turns the aux texts into a stage per
sample and a list of events.

Matlab v5 / v7.3 files and other WFDB signal formats are not read.

=cut
