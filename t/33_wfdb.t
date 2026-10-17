use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use PDL;
use PDL::EEG::IO::WFDB qw(read_wfdb read_wfdb_annotations wfdb_stages_events);

# Synthetic records and annotation files, written here (no data files shipped).
my $dir = tempdir(CLEANUP => 1);
my ($nsig, $nt) = (3, 500);
my $raw = (sequence(long, $nsig, $nt) * 7 % 2001 - 1000)->short;       # (nsig, nt), values -1000 .. 1000
$raw->slice('(2),:') .= (sequence($nt) % 100 - 50)->short;
my $bytes = pack('s<*', $raw->list);                                     # channels interleaved, little-endian
my @gain = (1, 200, 0);                                                  # 0 -> WFDB default 200
my @base = (-20, 5, 0);
my $sigline = sub {
    my ($file, $fmt) = @_;
    return (
        "$file $fmt 1(-20)/uV 16 0 -20 0 0 F3-M2",
        "$file $fmt 200(5)/mV 16 0 5 0 0 ECG",
        "$file $fmt 0/mV 16 0 0 0 0 ABD",
    );
};
my $write = sub { my ($f, $s) = @_; open my $h, '>:raw', "$dir/$f" or die; print $h $s; close $h };

# ---- Matlab v4 .mat (Challenge 2018 layout) ------------------------------------
my $mat = pack('l<5', 20, $nsig, $nt, 0, 4) . "val\0" . $bytes;          # 24-byte header
$write->('rec1.mat', $mat);
$write->('rec1.hea', join("\n", "rec1 $nsig 200 $nt", $sigline->('rec1.mat', '16+24')) . "\n");
my $r1 = read_wfdb("$dir/rec1");
is_deeply([$r1->{raw}->dims], [$nsig, $nt], 'read_wfdb .mat: (nsig, nt)');
ok(all($r1->{raw} == $raw), 'read_wfdb .mat: samples as written');
is($r1->{fs}, 200, 'read_wfdb: fs');
is_deeply($r1->{labels}, [qw(F3-M2 ECG ABD)], 'read_wfdb: labels');
is_deeply($r1->{units}, [qw(uV mV mV)], 'read_wfdb: units');
is_deeply($r1->{gain}, [1, 200, 200], 'read_wfdb: gain 0 -> 200');
my $ecg = $r1->{chan}->('ECG');
cmp_ok(($ecg - ($raw->slice('(1),:')->double - 5) / 200)->abs->max, '<', 1e-12, 'chan: (raw - baseline) / gain');
my $part = $r1->{chan}->('F3-M2', from => 100, to => 150);
is($part->nelem, 50, 'chan: from .. to-1');
cmp_ok(($part - ($raw->slice('(0),100:149')->double + 20))->abs->max, '<', 1e-12, 'chan: the right samples');
eval { $r1->{chan}->('EOG') };
like($@, qr/no channel 'EOG'/, 'chan: unknown label croaks');
my $rd = read_wfdb("$dir/rec1", data => 1);
cmp_ok(($rd->{data}->slice('(1),:') - $ecg)->abs->max, '<', 1e-12, 'data => 1: all channels in physical units');

# ---- .dat with the same samples, offset 0 ----------------------------------------
$write->('rec2.dat', $bytes);
$write->('rec2.hea', join("\n", "# a comment line", "rec2 $nsig 250 $nt", $sigline->('rec2.dat', '16')) . "\n");
my $r2 = read_wfdb("$dir/rec2");
ok(all($r2->{raw} == $raw), 'read_wfdb .dat: samples as written');
is($r2->{fs}, 250, 'read_wfdb .dat: fs');

# ---- a .mat cut short is read up to its last whole sample ------------------------
$write->('rec3.mat', substr($mat, 0, 24 + 2 * $nsig * 300 + 3));
$write->('rec3.hea', join("\n", "rec3 $nsig 200 $nt", $sigline->('rec3.mat', '16+24')) . "\n");
my $warn = '';
my $r3 = do { local $SIG{__WARN__} = sub { $warn .= shift }; read_wfdb("$dir/rec3") };
is($r3->{nt}, 300, 'cut-short file: 300 whole samples');
like($warn, qr/holds 300 of 500 samples/, 'cut-short file: warning');
ok(all($r3->{raw} == $raw->slice(':,0:299')), 'cut-short file: samples as written');

# ---- refusals ----------------------------------------------------------------------
$write->('rec4.hea', "rec4 1 200 10\nrec4.dat 212 200/mV 12 0 0 0 0 X\n");
eval { read_wfdb("$dir/rec4") };
like($@, qr/format 212/, 'read_wfdb: other formats refused');
$write->('rec5.hea', "rec5 2 200 10\na.dat 16 200/mV 16 0 0 0 0 X\nb.dat 16 200/mV 16 0 0 0 0 Y\n");
eval { read_wfdb("$dir/rec5") };
like($@, qr/more than one file/, 'read_wfdb: signals in several files refused');

# ---- annotations -------------------------------------------------------------------
# word = A << 10 | I.  Events: W at 0, "(resp_hypopnea" at 300, a SKIP of -1 then
# "resp_hypopnea)" at 300 + 900 - 1, N2 at 2000 (via a SKIP of 1500), a NUM word,
# a code-0 annotation, the end word.
my $aux = sub { my $s = shift; my $p = length($s) % 2 ? "$s\0" : $s; return (63 << 10 | length($s), unpack('v*', $p)) };
my @w = (
    22 << 10 | 0, $aux->('W'),
    22 << 10 | 300, $aux->('(resp_hypopnea'),
    59 << 10, 0xFFFF, 0xFFFF,                    # SKIP -1
    22 << 10 | 900, $aux->('resp_hypopnea)'),
    59 << 10, 0, 1500 - 900 + 1,                 # SKIP to 300 + 899 + 601 = 1800
    22 << 10 | 200, $aux->('N2'),                # 2000
    60 << 10 | 3,                                # NUM, no annotation
    0 << 10 | 10,                                # code 0 at 2010: dropped
    22 << 10 | 100, $aux->('## a header note'),  # 2110
    22 << 10 | 50, $aux->('Lights on'),          # 2160
    0,
);
$write->('rec1.arousal', pack('v*', @w));
my $ann = read_wfdb_annotations("$dir/rec1.arousal");
is_deeply([map { [ @$_[0, 2] ] } @$ann],
          [[0, 'W'], [300, '(resp_hypopnea'], [1199, 'resp_hypopnea)'], [2000, 'N2'], [2110, '## a header note'], [2160, 'Lights on']],
          'read_wfdb_annotations: times (signed SKIP) and aux texts; code 0 dropped');
my ($st, $names, $ev, $other) = wfdb_stages_events($ann, 2500);
is($st->nelem, 2500, 'stages: one per sample');
is_deeply([$st->at(0), $st->at(1999), $st->at(2000), $st->at(2499)], [1, 1, 3, 3], 'stages: W until 1999, N2 from 2000');
is($names->{3}, 'N2', 'stages: code names');
is_deeply($ev, [['resp_hypopnea', 300, 1199]], 'events: name, start, end');
is_deeply($other, { 'Lights on' => 1 }, 'other texts counted, ## notes skipped');

done_testing();
