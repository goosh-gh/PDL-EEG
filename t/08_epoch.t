use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::IO::NihonKohden qw(block_ranges select_range clock_to_samp);

# data[c,s] = s  over 900 samples @100 Hz (9 s of gap-removed data)
my $fs   = 100;
my $n    = 900;
my $data = sequence($n)->dummy(0, 3)->sever;

# Simulate what read_nk attaches for extblock: each event carries an epoch-based
# data-sample position ({samp}/{t_data}) that is DIFFERENT from wall-clock {t}.
# Three recording segments of 3 s data each, separated by 10 s wall-clock gaps.
my $rec = {
    data   => $data->float,
    fs     => $fs,
    labels => [qw(a b c)],
    events => [
        { label => 'REC START MMN EEG', t => 0,  samp => 0,   t_data => 0 },
        { label => 'task1',             t => 1,  samp => 100, t_data => 1 },
        { label => 'REC START MMN EEG', t => 13, samp => 300, t_data => 3 },
        { label => 'task2',             t => 14, samp => 400, t_data => 4 },
        { label => 'REC START MMN EEG', t => 26, samp => 600, t_data => 6 },
        { label => 'task3',             t => 27, samp => 700, t_data => 7 },
    ],
    t_start => '2026-07-02 14:03:03',
    block_meta => [ { start_samp => 0, n_samp => $n, t_start => '2026-07-02 14:03:03' } ],
    n_blocks   => 1,
};

# --- block_ranges must use epoch samp (not wall-clock t) --------------------
my $r = block_ranges($rec);
is scalar @$r, 3, 'three segments from REC START samp';
is $r->[1]{start}, 300, 'segment 1 starts at data sample 300 (not t*fs=1300)';
is $r->[1]{end},   600, 'segment 1 ends at 600';
is $r->[2]{end},   900, 'segment 2 ends at data length';

# --- select_range on a segment ---------------------------------------------
my $seg = select_range($rec, 300, 600);
is $seg->{data}->dim(1), 300, 'segment length 300';
is $seg->{data}->at(0, 0), 300, 'segment starts at original sample 300';
# events inside [300,600): the REC START@300 and task2@400 -> rebased
my @lab = map { $_->{label} } @{ $seg->{events} };
ok( (grep { /task2/ } @lab), 'task2 present in segment' );
my ($t2) = grep { $_->{label} eq 'task2' } @{ $seg->{events} };
ok abs($t2->{t_data} - 1.0) < 1e-9, 'task2 rebased to 1.0 s within segment';
ok !(grep { /task3/ } @lab), 'task3 (sample 700) excluded from [300,600)';

# --- clock_to_samp: wall-clock -> data sample via anchors -------------------
# query wall 14 s -> anchor REC START{t=13,samp=300} -> 300 + (14-13)*100 = 400
is clock_to_samp($rec, 14), 400, 'wall 14 s maps to data sample 400';
# query wall 27 s -> anchor {t=26,samp=600} -> 600 + 100 = 700
is clock_to_samp($rec, 27), 700, 'wall 27 s maps to data sample 700';
# before first anchor -> 0
is clock_to_samp($rec, 0), 0, 'wall 0 maps to 0';

# --- _attach_epoch_samp: epoch boundary for REC START, wall-clock offset within
#     segment for other markers -------------------------------------------------
{
    my $nt = 970000;          # pretend total; max_epoch below is 97 -> L~10000
    my @ev = (
        { label => 'REC START MMN CAL', epoch => 1,  t => 0   },
        { label => 'task1',             epoch => 2,  t => 16  },  # segment 0, +16s
        { label => 'REC START MMN EEG', epoch => 21, t => 325 },
        { label => 'task2',             epoch => 21, t => 327 },  # +2s into seg
        { label => 'REC START MMN CAL', epoch => 97, t => 1602 }, # sets max_epoch=97
    );
    PDL::EEG::IO::NihonKohden::_attach_epoch_samp(\@ev, $nt, 1000);
    my $L = $nt / 97;

    is $ev[0]{samp}, 0, 'REC START epoch1 -> sample 0';
    # task1 is in segment 0 (anchor REC START epoch1 @0), +16 s -> 16000
    is $ev[1]{samp}, 16000, 'task1 = anchor(0) + 16 s wall-clock offset';
    # REC START epoch21 -> epoch boundary
    is $ev[2]{samp}, int((21 - 1) * $L + 0.5), 'REC START epoch21 at epoch boundary';
    # task2 = that boundary + (327-325) s
    is $ev[3]{samp}, $ev[2]{samp} + 2000, 'task2 = anchor + 2 s within segment';
    ok $ev[4]{t_data} < $nt/1000, 'final CAL within data length';
}

# --- clock_to_samp via block_meta (wfmblock multi-block, NO events) -----------
# P:EEG13: with per-block wall-clock t_start in block_meta and no REC START
# events, wall-clock must follow the piecewise map, not a 1:1 map. Three 3 s
# segments @100 Hz, real 10 s wall-clock gaps between them; data butt-joined to
# 0..299 / 300..599 / 600..899. (clock_to_samp uses only epoch DIFFERENCES, so
# these are timezone-independent.)
{
    my $n2 = 900;
    my $d2 = sequence($n2)->dummy(0, 3)->sever->float;
    my $recw = {
        data => $d2, fs => 100, labels => [qw(a b c)], events => [],
        t_start    => '2026-07-02 14:03:03',
        block_meta => [
            { start_samp => 0,   n_samp => 300, t_start => '2026-07-02 14:03:03' },
            { start_samp => 300, n_samp => 300, t_start => '2026-07-02 14:03:16' },
            { start_samp => 600, n_samp => 300, t_start => '2026-07-02 14:03:29' },
        ],
        n_blocks => 3,
    };
    is clock_to_samp($recw, 1),  100, 'block_meta: wall 1s  -> 100 (seg0, 1:1 inside)';
    is clock_to_samp($recw, 14), 400, 'block_meta: wall 14s -> 400 (seg1, NOT 1400)';
    is clock_to_samp($recw, 26), 600, 'block_meta: wall 26s -> 600 (seg2 start)';
    is clock_to_samp($recw, 8),  299, 'block_meta: wall 8s in a gap -> 299 (clamp seg0 end)';
    is clock_to_samp($recw, 40), 899, 'block_meta: past end -> 899 (data end)';
}

# --- _attach_recstart_samp: REC START-delimited exact placement --------------
# Real-data regression (subject.EEG): the .LOG clock counts paused setup time
# between blocks, so wall-clock placement drifts and late-segment events get
# misfiled / clamped to the data end. REC START-delimited placement anchors each
# segment to its exact header boundary and is drift-free.
{
    my $fs2  = 1000;
    my $meta = [
        { start_samp=>0,      n_samp=>205000, t_start=>'2026-07-02 14:03:03' },
        { start_samp=>205000, n_samp=>176000, t_start=>'2026-07-02 14:07:52' },
        { start_samp=>381000, n_samp=>30000,  t_start=>'2026-07-02 14:11:34' },
        { start_samp=>411000, n_samp=>214000, t_start=>'2026-07-02 14:12:21' },
        { start_samp=>625000, n_samp=>206000, t_start=>'2026-07-02 14:16:41' },
        { start_samp=>831000, n_samp=>62000,  t_start=>'2026-07-02 14:20:34' },
        { start_samp=>893000, n_samp=>69000,  t_start=>'2026-07-02 14:21:48' },
        { start_samp=>962000, n_samp=>19000,  t_start=>'2026-07-02 14:23:02' },
    ];
    my @ev = map { { t=>$_->[0], label=>$_->[1] } } (
        [16,'task1'], [325,'REC START MMN EEG'], [327,'task2'],
        [621,'REC START MMN EEG'], [625,'task3 practice'],
        [651,'REC START MMN EEG'], [654,'task4'],
        [1025,'REC START MMN EEG'], [1033,'task5'],
        [1351,'REC START MMN EEG'], [1355,'安静開眼'],
        [1453,'REC START MMN EEG'], [1456,'安静閉眼'],
        [1602,'REC START MMN CAL'],
    );
    # prepend seg0's REC START so REC START count (8) == segment count (8)
    unshift @ev, { t=>0, label=>'REC START MMN CAL' };

    my $ok = PDL::EEG::IO::NihonKohden::_attach_recstart_samp(\@ev, $meta, $fs2);
    ok $ok, '_attach_recstart_samp fired (8 REC STARTs == 8 segments)';

    my %s = map { $_->{label} => $_->{samp} } @ev;
    is $s{task1},           16000,  'task1  -> 16000 (seg0)';
    is $s{task2},           207000, 'task2  -> 207000 (seg1, +2s; NOT 243000)';
    is $s{'task3 practice'},385000, 'task3  -> 385000 (seg2, correct segment)';
    is $s{task4},           414000, 'task4  -> 414000 (seg3)';
    is $s{task5},           633000, 'task5  -> 633000 (seg4; not clamped to end)';
    is $s{'安静開眼'},       835000, '安静開眼 -> 835000 (seg5)';
    is $s{'安静閉眼'},       896000, '安静閉眼 -> 896000 (seg6)';

    # every event lands inside its own segment (no cross-segment misfiling)
    my $inside = 1;
    for my $e (@ev) {
        my $in = 0;
        for my $m (@$meta) {
            $in = 1 if $e->{samp} >= $m->{start_samp}
                    && $e->{samp} <  $m->{start_samp} + $m->{n_samp};
        }
        $inside = 0 unless $in;
    }
    ok $inside, 'all events land inside a real segment (none clamped to data end)';

    # count-mismatch -> returns 0 so the caller can fall back
    my @ev2 = ({ t=>0, label=>'task only' });
    my $r2  = PDL::EEG::IO::NihonKohden::_attach_recstart_samp(\@ev2, $meta, $fs2);
    is $r2, 0, 'REC START count != segment count -> 0 (fall back)';
}

done_testing();
