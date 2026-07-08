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

done_testing();
