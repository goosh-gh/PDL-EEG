use strict;
use warnings;
use Test::More;
use PDL::EEG::IO::NihonKohden::PTN qw(parse_ptn find_montage_file);
use File::Temp qw(tempfile);
use FindBin qw($Bin);

# ---------------------------------------------------------------------------
# Build a minimal but format-valid .PTN in memory and round-trip it, so the
# parser is tested with no external binary. Then, if the real 21A fixture is
# present (t/data/Pattern_032.PTN), verify against it too.
# ---------------------------------------------------------------------------

sub rec80 {
    my (%a) = @_;                                 # g1,g2,sens,pos,name
    my $r = "\0" x 80;
    substr($r, 0, 1) = chr($a{g1}   // 0);
    substr($r, 1, 1) = chr($a{g2}   // 0);
    substr($r, 2, 1) = chr($a{sens} // 0);
    substr($r, 3, 7) = "\x0a\x0e\x09\x0a\x01\x04\x03";  # channel-type signature
    substr($r, 12, 2) = pack('v', $a{pos} // 1);        # display x-position
    if (defined $a{name}) {
        my $nm = substr($a{name}, 0, 30);
        substr($r, 14, length($nm) + 1) = $nm . "\0";
    }
    return $r;
}

sub build_ptn {
    my $b = "\0" x 1040;
    substr($b, 0, 31)    = "EEG-1000/9000 Pattern Info File";
    substr($b, 0x80, 3)  = "TST";
    my @recs = (
        rec80(g1 => 0, g2 => 0x25, sens => 0x0f, pos => 0x20, name => 'Fp1'),  # inline
        rec80(g1 => 1, g2 => 0x25, sens => 0x0f, pos => 0x40, name => 'Fp2'),  # inline
        rec80(g1 => 2, g2 => 0x25, sens => 0x0f, pos => 0x60),                 # no inline
        rec80(g1 => 0, g2 => 0,    sens => 0x05, pos => 0x80, name => 'TrigA'),# trigger
        rec80(g1 => 0, g2 => 0,    sens => 0x05, pos => 0xA0, name => 'TrigB'),# trigger
    );
    # end-of-list sentinel: display-position word == 0xFFFF
    my $end = "\0" x 80; substr($end, 12, 2) = pack('v', 0xFFFF);
    return $b . join('', @recs) . $end . ("\0" x 160);
}

my ($fh, $path) = tempfile(SUFFIX => '.PTN', UNLINK => 1);
binmode $fh; print {$fh} build_ptn(); close $fh;

my $m = parse_ptn($path);
is($m->{name}, 'TST', 'synthetic NAME parsed');
is($m->{n}, 5, 'stops at 0xFFFF sentinel (5 channels, padding ignored)');
is(scalar @{ $m->{triggers} }, 2, 'two trigger channels found');
is($m->{channels}[0]{inline}, 'Fp1', 'inline name #0');
is($m->{channels}[0]{ch_idx}, 1, 'EEG ch_idx = G1+1');
is($m->{channels}[2]{inline}, undef, 'no inline name where absent');
is($m->{channels}[2]{ch_idx}, 3, 'ch_idx still derived from G1');
ok($m->{channels}[3]{trigger}, 'trigger flagged (G1=0,G2=0,SENS=0x05)');
is($m->{channels}[3]{ch_idx}, undef, 'trigger has no electrode index');
is($m->{channels}[4]{inline}, 'TrigB', 'trigger inline name');

# --- optional: real 21A fixture --------------------------------------------
my $fix = "$Bin/data/Pattern_032.PTN";
SKIP: {
    skip "real fixture $fix not present", 4 unless -f $fix;
    my $r = parse_ptn($fix);
    is($r->{name}, '21A', 'real fixture NAME = 21A');
    is($r->{n}, 25, 'real fixture has 25 display channels');
    is(scalar @{ $r->{triggers} }, 4, 'real fixture has 4 trigger channels');
    is($r->{channels}[$r->{triggers}[0]]{inline}, 'TrigBit0', 'first trigger = TrigBit0');
}

done_testing();
