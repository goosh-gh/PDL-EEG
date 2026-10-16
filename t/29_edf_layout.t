use strict;
use warnings;
use Test::More;
use File::Temp qw(tempfile);
use PDL;
use PDL::EEG::IO::EDF qw(write_edf read_edf);

# Byte layout of the data records, checked for every sample.
# The EDF specification: header 256 + 256*ns bytes; then nrec data records,
# each holding, signal after signal, spr[s] little-endian int16 samples.
# Plain Perl (pack / unpack) on both sides of write_edf and read_edf.

sub slurp { my ($p) = @_; open my $fh, '<:raw', $p or die $!; local $/; my $b = <$fh>; return $b }
sub spew  { my ($p, $b) = @_; open my $fh, '>:raw', $p or die $!; print {$fh} $b; close $fh or die $! }
sub fld   { my ($s) = @_; $s =~ s/\s+$//; return $s }

# parse a file: header fields, and per signal the concatenated samples of all records
sub parse_edf {
    my ($buf) = @_;
    my $ns   = 0 + fld(substr($buf, 252, 4));
    my $nrec = 0 + fld(substr($buf, 236, 8));
    my $base = 256;
    my @label = map { fld(substr($buf, $base + 16 * $_, 16)) } 0 .. $ns - 1;
    my @spr   = map { 0 + fld(substr($buf, $base + 216 * $ns + 8 * $_, 8)) } 0 .. $ns - 1;
    my $p = 256 * (1 + $ns);
    my @samples = map { [] } 0 .. $ns - 1;
    my @bytes   = ('') x $ns;
    for my $r (0 .. $nrec - 1) {
        for my $s (0 .. $ns - 1) {
            my $chunk = substr($buf, $p, 2 * $spr[$s]); $p += 2 * $spr[$s];
            push @{ $samples[$s] }, unpack('s<*', $chunk);
            $bytes[$s] .= $chunk;
        }
    }
    return { ns => $ns, nrec => $nrec, label => \@label, spr => \@spr,
             samples => \@samples, bytes => \@bytes, end => $p, len => length $buf };
}

# integer data in [-30000, 30000]; with phys => 'gain', gain => 1 the digital
# code of a sample is the sample itself, and the padding code is 0.
sub value { my ($c, $i) = @_; return ((($c + 1) * 1000 + 7 * $i + 13 * $c * $i) % 60001) - 30000 }

for my $case ([3, 250, 100, 1], [5, 1000, 100, 0], [1, 7, 4, 0], [4, 333, 50, 1], [2, 1, 10, 1]) {
    my ($nch, $ns_samp, $fs, $plus) = @$case;
    my $tag = "write_edf: $nch ch x $ns_samp samples, fs $fs, plus $plus";
    my @rows = map { my $c = $_; [ map { value($c, $_) } 0 .. $ns_samp - 1 ] } 0 .. $nch - 1;
    my $data = pdl(float, \@rows)->transpose;                          # (n_ch, n_samp)
    is_deeply([$data->dims], [$nch, $ns_samp], "$tag: input dims") or next;
    my $rec = { data => $data, fs => $fs, labels => [ map { "E$_" } 1 .. $nch ],
                events => [ { time => 0.0, label => 'start' } ], t_start => '2026-10-16 10:00:00' };
    my (undef, $path) = tempfile(SUFFIX => '.edf', UNLINK => 1);
    write_edf($rec, $path, phys => 'gain', gain => 1, plus => $plus);
    my $f = parse_edf(slurp($path));

    my $nrec = int(($ns_samp + $fs - 1) / $fs);
    is($f->{nrec}, $nrec, "$tag: record count");
    is($f->{ns}, $nch + ($plus ? 1 : 0), "$tag: signal count");
    is($f->{end}, $f->{len}, "$tag: records fill the file exactly");
    for my $c (0 .. $nch - 1) {
        my @want = (@{ $rows[$c] }, (0) x ($nrec * $fs - $ns_samp));
        is_deeply($f->{samples}[$c], \@want, "$tag: channel $c, every sample in its record");
    }
    if ($plus) {
        is($f->{label}[-1], 'EDF Annotations', "$tag: annotation signal last");
        like($f->{bytes}[-1], qr/^\+0\x14\x14/, "$tag: first record starts with its time-keeping TAL");
    }

    my $back = read_edf($path);
    is_deeply([$back->{data}->dims], [$nch, $nrec * $fs], "$tag: read_edf dims");
    my @got = map { [ $back->{data}->slice("($_),:")->list ] } 0 .. $nch - 1;
    my @want = map { [ @{ $rows[$_] }, (0) x ($nrec * $fs - $ns_samp) ] } 0 .. $nch - 1;
    is_deeply(\@got, \@want, "$tag: read_edf returns every sample");
}

# ---- samples outside the digital range: clipped, counted in a warning ------------
{
    my @v = (0, 40000, -40000, 32767, -32768, 32767.6, -32768.6, 5);
    my $data = pdl(float, [\@v, [ (1) x @v ]])->transpose;   # ch 0 clips, ch 1 does not
    my (undef, $path) = tempfile(SUFFIX => '.edf', UNLINK => 1);
    my @w; local $SIG{__WARN__} = sub { push @w, @_ };
    write_edf({ data => $data, fs => 4, labels => [qw(A B)], t_start => '2026-10-16 10:00:00' }, $path,
              phys => 'gain', gain => 1, plus => 0);
    my $f = parse_edf(slurp($path));
    is_deeply($f->{samples}[0], [0, 32767, -32768, 32767, -32768, 32767, -32768, 5],
              'clipping: out-of-range samples become the digital limits');
    is(scalar @w, 1, 'clipping: one warning, for the clipped channel only');
    like($w[0] // '', qr/channel 0 \(A\): 4 sample\(s\) clipped/, 'clipping: the warning counts 4 samples');
}

# ---- read_edf on a file written here: annotation signal FIRST, other spr ----------
{
    my ($nrec, $spr, $sprA) = (4, 6, 9);
    my @lab  = ('EDF Annotations', 'Fz', 'Cz');
    my @sprs = ($sprA, $spr, $spr);
    my $ns = @lab;
    my $hdr = sprintf('%-8s%-80s%-80s%-8s%-8s%-8d%-44s%-8d%-8s%-4d',
        '0', 'X X X X', 'Startdate 16-OCT-2026 X X X', '16.10.26', '10.00.00',
        256 * (1 + $ns), 'EDF+C', $nrec, '1', $ns);
    $hdr .= join '', map { sprintf '%-16s', $_ } @lab;
    $hdr .= ' ' x (80 * $ns);                                          # transducer
    $hdr .= join '', map { sprintf '%-8s', $_ } ('', 'uV', 'uV');      # physical dim
    $hdr .= join '', map { sprintf '%-8s', $_ } (-1, -32768, -32768);  # physical min
    $hdr .= join '', map { sprintf '%-8s', $_ } (1, 32767, 32767);     # physical max
    $hdr .= join '', map { sprintf '%-8s', $_ } (-32768) x 3;          # digital min
    $hdr .= join '', map { sprintf '%-8s', $_ } (32767) x 3;           # digital max
    $hdr .= ' ' x (80 * $ns);                                          # prefilter
    $hdr .= join '', map { sprintf '%-8s', $_ } @sprs;                 # samples per record
    $hdr .= ' ' x (32 * $ns);
    is(length $hdr, 256 * (1 + $ns), 'hand-made header length');

    my (@fz, @cz, $body);
    for my $r (0 .. $nrec - 1) {
        my $tal = "+$r\x14\x14\x00" . ($r == 2 ? "+2.5\x14blink\x14\x00" : '');
        $tal .= "\x00" x (2 * $sprA - length $tal);
        my @a = map { 100 * $r + $_ } 0 .. $spr - 1;
        my @b = map { -1000 * $r - 3 * $_ } 0 .. $spr - 1;
        push @fz, @a; push @cz, @b;
        $body .= $tal . pack('s<*', @a) . pack('s<*', @b);
    }
    my (undef, $path) = tempfile(SUFFIX => '.edf', UNLINK => 1);
    spew($path, $hdr . $body);
    my $back = read_edf($path);
    is_deeply($back->{labels}, ['Fz', 'Cz'], 'annotation first: labels of the data signals');
    is_deeply([ $back->{data}->slice('(0),:')->list ], \@fz, 'annotation first: every Fz sample');
    is_deeply([ $back->{data}->slice('(1),:')->list ], \@cz, 'annotation first: every Cz sample');
    is_deeply([ map { [ $_->{onset}, $_->{label} ] } @{ $back->{events} } ], [[2.5, 'blink']],
              'annotation first: the event of record 2');

    # a short last record is dropped
    spew($path, $hdr . substr($body, 0, length($body) - 5));
    my $cut = read_edf($path);
    is($cut->{data}->dim(1), 3 * $spr, 'short last record: 3 complete records read');
    is_deeply([ $cut->{data}->slice('(1),:')->list ], [ @cz[0 .. 3 * $spr - 1] ],
              'short last record: every sample of the complete ones');
}

done_testing;
