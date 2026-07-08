use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::IO::NihonKohden::Montage qw(montage_from_log resolve_labels);
use File::Temp qw(tempdir);

# --- montage_from_log -------------------------------------------------------
is(montage_from_log([
    { t => 0,  label => 'REC START IIA EEG' },
    { t => 1,  label => 'A1+A2 OFF' },
]), 'IIA', 'montage name from REC START ... EEG');

is(montage_from_log([
    { t => 0, label => 'REC START 21A CAL' },
]), '21A', 'montage name from a CAL marker too');

is(montage_from_log([ { t => 0, label => 'A1+A2 OFF' } ]), undef,
   'undef when no REC START marker');

# --- build a matching .PTN (NAME "TST", 4 triggers T0..T3) ------------------
sub rec80 {
    my (%a) = @_;
    my $r = "\0" x 80;
    substr($r, 0, 1)  = chr($a{g1}   // 0);
    substr($r, 1, 1)  = chr($a{g2}   // 0);
    substr($r, 2, 1)  = chr($a{sens} // 0);
    substr($r, 3, 7)  = "\x0a\x0e\x09\x0a\x01\x04\x03";
    substr($r, 12, 2) = pack('v', $a{pos} // 1);
    if (defined $a{name}) { substr($r, 14, length($a{name}) + 1) = $a{name} . "\0" }
    return $r;
}
sub build_ptn {
    my $b = "\0" x 1040;
    substr($b, 0, 31)   = "EEG-1000/9000 Pattern Info File";
    substr($b, 0x80, 3) = "TST";
    my @recs = (
        rec80(g1 => 0, g2 => 0x25, sens => 0x0f, pos => 0x20, name => 'Fp1'),
        rec80(g1 => 0, g2 => 0, sens => 0x05, pos => 0x40, name => 'T0'),
        rec80(g1 => 0, g2 => 0, sens => 0x05, pos => 0x50, name => 'T1'),
        rec80(g1 => 0, g2 => 0, sens => 0x05, pos => 0x60, name => 'T2'),
        rec80(g1 => 0, g2 => 0, sens => 0x05, pos => 0x70, name => 'T3'),
    );
    my $end = "\0" x 80; substr($end, 12, 2) = pack('v', 0xFFFF);
    return $b . join('', @recs) . $end;
}
my $dir = tempdir(CLEANUP => 1);
open my $fh, '>:raw', "$dir/Pattern_001.PTN" or die $!;
print {$fh} build_ptn(); close $fh;

# --- synthetic record: 6 EEG + 4 trigger channels, known ch_indices ---------
my $n = 2000; srandom(7);
my @rows;
push @rows, grandom($n) * (15 + 8 * $_) for 0 .. 5;      # 6 EEG (small ranges)
for my $k (0 .. 3) {                                       # 4 triggers
    my $ch = zeroes(double, $n);
    my $s = 150 + $k * 350;
    $ch->slice("$s:" . ($s + 50)) .= 3188;
    push @rows, $ch;
}
my $data = cat(@rows)->transpose;                          # [10, n]

# recorded ch_idx: EEG at 1..6, triggers at 45,46,47,74 (deliberately non-contiguous)
my @ch_indices = (1, 2, 3, 4, 5, 6, 45, 46, 47, 74);

my $rec = {
    data       => $data->float,
    fs         => 1000,
    labels     => [ (map { "E$_" } 1 .. 6), qw(x x x x) ],
    ch_indices => \@ch_indices,
    events     => [ { t => 0, label => 'REC START TST EEG' } ],
};

my $r = resolve_labels($rec, ptn_dir => $dir, skip_sec => 0, apply => 1);

is($r->{montage}, 'TST', 'resolved montage name from events');
ok($r->{ptn}, 'matching .PTN located by NAME');
is(scalar @{ $r->{triggers} }, 4, 'four triggers resolved');

# names (montage order T0..T3) zipped onto triggers sorted by ch_idx (45,46,47,74)
is_deeply($r->{label_map},
    { 45 => 'T0', 46 => 'T1', 47 => 'T2', 74 => 'T3' },
    'label_map binds montage names to detected ch_idx in order');

# apply => 1 rewrote the in-place labels via ch_indices
is($rec->{labels}[6], 'T0', 'labels updated for ch_idx 45');
is($rec->{labels}[9], 'T3', 'labels updated for ch_idx 74');

done_testing();
