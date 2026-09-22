package PDL::EEG::Spans;
# Read/write time-span annotations (as produced by plot_raw_nihonkohden.pl) and
# turn them into per-sample boolean masks for covariance construction (GED S/R,
# ICA segment selection).
#
# A span is a hashref { kind => STR, t0 => SEC, t1 => SEC } where t0/t1 are in
# the recording's data-coordinate seconds (the same axis the viewer prints, i.e.
# absolute sample index / fs). kind is a free string; the viewer uses
# 'blink', 'hsaccade', 'bad', 'clean'.
#
# Sidecar file format: tab-separated, one span per line, sorted by t0:
#   # kind    t0        t1
#   blink     12.340    12.610
#   hsaccade  45.100    45.480
# Lines beginning with '#' and blank lines are ignored on read.

use strict;
use warnings;
use PDL;
use Exporter 'import';
our @EXPORT_OK = qw(read_spans write_spans spans_to_mask kinds_present);

sub write_spans {
    my ($path, $spans) = @_;
    open my $fh, '>', $path or die "write_spans: cannot write $path: $!\n";
    print {$fh} "# kind\tt0\tt1\n";
    for my $s (sort { $a->{t0} <=> $b->{t0} } @$spans) {
        printf {$fh} "%s\t%.4f\t%.4f\n", $s->{kind}, $s->{t0}, $s->{t1};
    }
    close $fh;
    return scalar @$spans;
}

sub read_spans {
    my ($path) = @_;
    my @spans;
    open my $fh, '<', $path or die "read_spans: cannot read $path: $!\n";
    while (my $line = <$fh>) {
        next if $line =~ /^\s*#/ || $line =~ /^\s*$/;
        chomp $line;
        my ($kind, $t0, $t1) = split /\t/, $line;
        next unless defined $t1;
        push @spans, { kind => $kind, t0 => $t0 + 0, t1 => $t1 + 0 };
    }
    close $fh;
    return \@spans;
}

# spans_to_mask(\@spans, fs => Hz, nt => n_samples, abs0 => 0, kind => undef|STR)
# Returns a byte PDL of length nt with 1 inside any matching span, else 0.
# Time t (data-coord sec) maps to local sample: round(t*fs) - abs0.
sub spans_to_mask {
    my ($spans, %o) = @_;
    my $fs   = $o{fs}   or die "spans_to_mask: fs required\n";
    my $nt   = $o{nt}   or die "spans_to_mask: nt required\n";
    my $abs0 = $o{abs0} // 0;
    my $want = $o{kind};                      # undef = any kind
    my $mask = zeroes(byte, $nt);
    for my $s (@$spans) {
        next if defined $want && $s->{kind} ne $want;
        my $s0 = int($s->{t0} * $fs + 0.5) - $abs0;
        my $s1 = int($s->{t1} * $fs + 0.5) - $abs0;
        ($s0, $s1) = ($s1, $s0) if $s0 > $s1;
        $s0 = 0       if $s0 < 0;
        $s1 = $nt - 1 if $s1 > $nt - 1;
        next if $s1 < 0 || $s0 > $nt - 1;
        $mask->slice("$s0:$s1") .= 1;
    }
    return $mask;
}

sub kinds_present {
    my ($spans) = @_;
    my %seen;
    $seen{ $_->{kind} }++ for @$spans;
    return sort keys %seen;
}

1;
