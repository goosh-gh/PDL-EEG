#!/usr/bin/env perl
# examples/find_bn_balance.pl
#
#   perl examples/find_bn_balance.pl JJ.EEG=0.7093 IJ.EEG=0.6394 [more...]
#
# We now KNOW the BN balance for several recordings, measured against the vendor's
# own .mul export (mul_to_nk.pl --solve-bne):
#
#     subject  (EEG-1200A)  p = 0.7093
#     IJ0200S9  (EEG-1260Next / 1200-family) p = 0.6394, logged as 0.65
#
# The balance is a hardware pot setting, so it MUST be stored somewhere for the
# viewer to re-reference with it -- either in the .EEG header, the .21e, or the
# .PNT. This looks for it, using the fact that we know the target number in each
# file: scan every plausible fixed-point / float encoding at every header offset,
# and keep the offsets where EVERY file holds ITS OWN known value.
#
# That "every file agrees at the same offset" test is what makes this reliable:
# a single file has thousands of coincidental byte matches, but the probability
# that the same offset encodes 0.7093 in one file and 0.6394 in another by chance
# is negligible.

use strict;
use warnings;

my @spec = @ARGV;
@spec >= 1 or die
    "usage: $0 file1.EEG=0.7093 [file2.EEG=0.6394 ...]\n" .
    "  give each .EEG with its measured BN balance (from mul_to_nk --solve-bne)\n";

# The value we MEASURED (e.g. 0.6394) is the effective mixing ratio -- pot setting
# plus contact-impedance asymmetry. What the operator DIALLED IN, and what the
# machine is likely to store, is a round number (0.65). So for each file we build
# a set of candidate set-values: the measured value, and every round value near
# it. The scan then accepts an offset if each file holds ANY of its candidates.
sub candidates {
    my $p = shift;
    my %c = ($p => 1);
    # round to 2 and 1 decimals, and the neighbours at 0.05 / 0.01 spacing
    for my $q (sprintf('%.2f', $p) + 0, sprintf('%.1f', $p) + 0) { $c{$q} = 1 }
    for my $step (0.05, 0.01) {
        my $base = int($p / $step + 0.5) * $step;
        $c{ sprintf('%.4f', $base + $step * $_) + 0 } = 1 for -1 .. 1;
    }
    return [ sort { $a <=> $b } grep { $_ > 0 && $_ < 1 } keys %c ];
}

my @F;
for (@spec) {
    my ($path, $p) = split /=/, $_, 2;
    defined $p or die "need file=value, got '$_'\n";
    open my $fh, '<:raw', $path or die "$path: $!";
    read $fh, my $buf, 1 << 16;
    close $fh;
    my $cand = candidates($p + 0);
    push @F, { path => $path, p => $p + 0, cand => $cand, buf => $buf };
    printf "%-40s measured %.4f   candidates {%s}   (%d bytes)\n",
        ($path =~ m{([^/]+)$})[0], $p + 0,
        join(',', map { sprintf('%.2f', $_) } @$cand), length $buf;
}
my $minlen = $F[0]{buf};
$minlen = length($_->{buf}) < length($minlen) ? $_->{buf} : $minlen for @F;
my $L = length $F[0]{buf};
$L = length($_->{buf}) < $L ? length($_->{buf}) : $L for @F;

# tolerance: the pot is set by hand to ~2 decimals, and our measurement itself
# has ~0.005 of scatter, so accept a value within 0.01 of the target.
my $TOL = 0.008;

# ---- the encodings a balance could be stored as ----------------------------
# Each: a name, the number of bytes, and a decoder from a substr to a number.
my @dec = (
    [ 'u8/255'      , 1, sub { unpack('C', $_[0]) / 255 } ],
    [ 'u8/256'      , 1, sub { unpack('C', $_[0]) / 256 } ],
    [ 'u8/100'      , 1, sub { unpack('C', $_[0]) / 100 } ],
    [ 'u16le/65535' , 2, sub { unpack('v', $_[0]) / 65535 } ],
    [ 'u16le/65536' , 2, sub { unpack('v', $_[0]) / 65536 } ],
    [ 'u16le/32768' , 2, sub { unpack('v', $_[0]) / 32768 } ],
    [ 'u16le/1000'  , 2, sub { unpack('v', $_[0]) / 1000 } ],
    [ 'u16le/10000' , 2, sub { unpack('v', $_[0]) / 10000 } ],
    [ 'u16be/65535' , 2, sub { unpack('n', $_[0]) / 65535 } ],
    [ 's16le/1000'  , 2, sub { my $v = unpack('v', $_[0]); ($v > 32767 ? $v - 65536 : $v) / 1000 } ],
    [ 'u32le/1e6'   , 4, sub { unpack('V', $_[0]) / 1e6 } ],
    [ 'f32le'       , 4, sub { unpack('f<', $_[0]) } ],
    [ 'f32be'       , 4, sub { unpack('f>', $_[0]) } ],
    [ 'f64le'       , 8, sub { unpack('d<', $_[0]) } ],
    [ 'ascii4'      , 4, sub { $_[0] =~ /^[\d.]+$/ ? $_[0] + 0 : undef } ],
    [ 'ascii5'      , 5, sub { $_[0] =~ /^[\d.]+$/ ? $_[0] + 0 : undef } ],
    [ 'ascii6'      , 6, sub { $_[0] =~ /^[\d.]+$/ ? $_[0] + 0 : undef } ],
);

# ---- scan ------------------------------------------------------------------
# For each (encoding, offset), does every file decode to its own p (within TOL)?
my @hit;
for my $d (@dec) {
    my ($name, $w, $fn) = @$d;
    OFF: for my $o (0 .. $L - $w) {
        my @matched;                        # the candidate each file matched
        for my $f (@F) {
            my $val = $fn->(substr($f->{buf}, $o, $w));
            next OFF unless defined $val;
            my ($best) = sort { abs($val - $a) <=> abs($val - $b) } @{ $f->{cand} };
            next OFF unless abs($val - $best) <= $TOL;
            push @matched, $best;
        }
        # the files must actually hold DIFFERENT values here, or it is a constant
        if (@F > 1) {
            my ($lo, $hi) = ($matched[0], $matched[0]);
            for (@matched) { $lo = $_ if $_ < $lo; $hi = $_ if $_ > $hi }
            next OFF if $hi - $lo < 0.02;
        }
        push @hit, [ $o, $name, $w, [@matched] ];
    }
}

# ---- report ----------------------------------------------------------------
if (!@hit) {
    print "\nNo offset encodes each file's own balance in any tried format.\n";
    print "The balance may be:\n";
    print "  - in the .21e or .PNT, not the .EEG (run this on those too)\n";
    print "  - stored only as a montage/derivation the viewer re-applies\n";
    print "  - in an encoding not in the list above\n";
    print "  - genuinely not persisted (the .mul carries the RESULT, not the setting)\n";
    exit 0;
}

printf "\n%d candidate offset(s) encode every file's own balance:\n\n", scalar @hit;
for my $h (sort { $a->[0] <=> $b->[0] } @hit) {
    my ($o, $name, $w, $matched) = @$h;
    printf "  offset 0x%04X (%5d)  as %-12s :", $o, $o, $name;
    my ($n2, $ww, $fn) = @{ (grep { $_->[0] eq $name } @dec)[0] };
    for my $i (0 .. $#F) {
        my $val = $fn->(substr($F[$i]{buf}, $o, $ww));
        printf "  %s: %.4f~%.2f", ($F[$i]{path} =~ m{([^/]+)$})[0], $val, $matched->[$i];
    }
    print "\n";
}

print "\nNext: read the value at that offset in a THIRD file whose balance you\n";
print "know, to confirm. If it holds, that offset is the BN balance, and\n";
print "read_nk() can pull the correct prop straight out of the header.\n";
