#!/usr/bin/env perl
# examples/find_bn_diff.pl
#
#   perl find_bn_diff.pl JJ.21e=0.71 IJ.21e=0.65 [--width N] [--tol 0.02]
#
# The value-scanning finder (find_bn_balance.pl) assumes an encoding. This one does
# not: it walks the two files byte for byte and, at every offset, asks a simpler
# question -- "could THESE bytes be file 1's balance while the SAME bytes in file 2
# are file 2's balance?" -- across a big menu of interpretations including ones the
# other tool lacks (BCD especially: NK stores timestamps in BCD, so 0.71 -> 0x71,
# 0.65 -> 0x65 is very much in play, and a plain u8 read never sees it).
#
# We know the two targets differ (0.71 vs 0.65), so the offset that encodes each
# file's own value is highly unlikely to arise by chance.

use strict;
use warnings;

my ($tol, $wmax, $cap) = (0.02, 4, 1 << 16);
my @spec;
for (@ARGV) {
    if    (/^--tol=(.+)/)   { $tol  = $1 }
    elsif (/^--width=(.+)/) { $wmax = $1 }
    elsif (/^--cap=(.+)/)   { $cap  = $1 }        # bytes to scan (0 = whole file)
    elsif (/^--all/)        { $cap  = 0 }
    elsif (/^--/)           { die "unknown option $_\n" }
    else                    { push @spec, $_ }
}
@spec == 2 or die "usage: $0 fileA=0.71 fileB=0.65 [--tol=0.02] [--width=4]\n";

# Only the HEADER matters. A balance is a setting; the rest of a .EEG is a
# gigabyte of waveform, and scanning it just turns up BCD-shaped noise (a run of
# sample bytes will read as "0.72" somewhere). Cap the read -- pass --all to
# override, but you almost never want to.
my @F;
for (@spec) {
    my ($p, $v) = split /=/, $_, 2;
    open my $fh, '<:raw', $p or die "$p: $!";
    my $sz = -s $p;
    my $take = $cap && $cap < $sz ? $cap : $sz;
    read $fh, my $buf, $take;
    close $fh;
    push @F, { path => $p, p => $v + 0, buf => $buf };
    printf "%-40s target %.4f   scanning %d of %d bytes%s\n",
        ($p =~ m{([^/]+)$})[0], $v + 0, $take, $sz,
        ($take < $sz ? "  (header only; --all for the lot)" : '');
}
die "the two targets are identical; need different balances\n"
    if abs($F[0]{p} - $F[1]{p}) < 0.02;

my $L = length($F[0]{buf}) < length($F[1]{buf})
      ? length($F[0]{buf}) : length($F[1]{buf});

# ---- interpretations of the bytes at an offset -----------------------------
# each returns a number (or undef) for a byte string
my @enc = (
    [ 'u8/100'      , 1, sub { unpack('C', $_[0]) / 100 } ],
    [ 'u8/255'      , 1, sub { unpack('C', $_[0]) / 255 } ],
    [ 'u8-pct'      , 1, sub { unpack('C', $_[0]) / 100 } ],           # 71 -> 0.71
    [ 'u8-bcd'      , 1, sub { my $b = unpack('C', $_[0]);             # 0x71 -> 0.71
                               my ($hi, $lo) = ($b >> 4, $b & 0xF);
                               ($hi <= 9 && $lo <= 9) ? ($hi * 10 + $lo) / 100 : undef } ],
    [ 'u16le/65535' , 2, sub { unpack('v', $_[0]) / 65535 } ],
    [ 'u16le/65536' , 2, sub { unpack('v', $_[0]) / 65536 } ],
    [ 'u16le/32768' , 2, sub { unpack('v', $_[0]) / 32768 } ],
    [ 'u16le/1000'  , 2, sub { unpack('v', $_[0]) / 1000 } ],          # 710 -> 0.71
    [ 'u16le/100'   , 2, sub { unpack('v', $_[0]) / 100 } ],           # 71  -> 0.71
    [ 'u16be/65535' , 2, sub { unpack('n', $_[0]) / 65535 } ],
    [ 'u16be/1000'  , 2, sub { unpack('n', $_[0]) / 1000 } ],
    [ 'u16le-bcd'   , 2, sub { my @b = unpack('C2', $_[0]);            # 0x0071 or 0x7100
                               my $n = 0;
                               for my $b (reverse @b) {
                                   my ($h, $l) = ($b >> 4, $b & 0xF);
                                   return undef if $h > 9 || $l > 9;
                                   $n = $n * 100 + $h * 10 + $l;
                               }
                               $n / 1000 } ],
    [ 'f32le'       , 4, sub { unpack('f<', $_[0]) } ],
    [ 'f32be'       , 4, sub { unpack('f>', $_[0]) } ],
    [ 'f64le'       , 8, sub { unpack('d<', $_[0]) } ],
    [ 'ascii3'      , 3, \&_an ], [ 'ascii4', 4, \&_an ],
    [ 'ascii5'      , 5, \&_an ], [ 'ascii6', 6, \&_an ],
);
sub _an {
    my $t = $_[0];
    return undef unless $t =~ /^\d?\.\d+$/ || $t =~ /^\d\.\d+$/;
    my $v = $t + 0;
    ($v > 0 && $v < 1) ? $v : undef;
}

# ---- scan ------------------------------------------------------------------
my @hit;
for my $e (@enc) {
    my ($name, $w, $fn) = @$e;
    next if $w > 8;
    for my $o (0 .. $L - $w) {
        my $va = $fn->(substr($F[0]{buf}, $o, $w)); next unless defined $va;
        next unless abs($va - $F[0]{p}) <= $tol;
        my $vb = $fn->(substr($F[1]{buf}, $o, $w)); next unless defined $vb;
        next unless abs($vb - $F[1]{p}) <= $tol;
        push @hit, [ $o, $name, $w, $va, $vb ];
    }
}

if (!@hit) {
    print "\nNo offset holds both targets in any encoding tried",
          " (incl. BCD, percent, per-mil).\n",
          "Widen --tol, or the value is stored somewhere these two files do not share.\n";
    # Still useful: show where each file INDIVIDUALLY has a byte that could be its
    # balance in BCD or percent, so the eye can look for structure.
    for my $fi (0, 1) {
        my $tgt = $F[$fi]{p};
        my @pct;
        for my $o (0 .. $L - 1) {
            my $b = unpack('C', substr($F[$fi]{buf}, $o, 1));
            my ($h, $l) = ($b >> 4, $b & 0xF);
            my $bcd = ($h <= 9 && $l <= 9) ? ($h * 10 + $l) / 100 : -1;
            push @pct, sprintf('0x%04X', $o)
                if abs($b / 100 - $tgt) <= $tol || abs($bcd - $tgt) <= $tol;
        }
        printf "\n  %s: bytes that alone could be %.2f (percent or BCD): %s%s\n",
            ($F[$fi]{path} =~ m{([^/]+)$})[0], $tgt,
            join(' ', @pct[0 .. ($#pct < 30 ? $#pct : 30)]),
            (@pct > 31 ? sprintf(' ... (%d total)', scalar @pct) : '');
    }
    exit 0;
}

printf "\n%d offset(s) hold EACH file's own balance:\n\n", scalar @hit;
my $exact = 0;
for my $h (sort { $a->[0] <=> $b->[0] || $a->[1] cmp $b->[1] } @hit) {
    my ($o, $name, $w, $va, $vb) = @$h;
    my $ha = uc unpack('H*', substr($F[0]{buf}, $o, $w));
    my $hb = uc unpack('H*', substr($F[1]{buf}, $o, $w));
    my $tight = (abs($va - $F[0]{p}) <= 0.005 && abs($vb - $F[1]{p}) <= 0.005);
    $exact++ if $tight;
    printf "  0x%05X (%6d)  %-12s  A=%.4f [%s]   B=%.4f [%s]%s\n",
        $o, $o, $name, $va, $ha, $vb, $hb, ($tight ? '   <== exact' : '');
}
print "\n";
if ($exact) {
    print "The '<== exact' line(s) hit both targets dead on -- those are the ones to\n"
        . "confirm with a third file of known balance.\n";
} else {
    print "NB: every hit is off-target by ~0.01 and they are scattered. That is what\n"
        . "BCD-shaped noise in a data region looks like, not a real field. Treat as\n"
        . "negative unless a third file reproduces the SAME offset exactly.\n";
}
