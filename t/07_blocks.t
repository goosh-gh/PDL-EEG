use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::IO::NihonKohden qw(block_ranges select_block);

# data[c,s] = s  (so we can verify which columns a slice picked)
my $fs   = 100;
my $n    = 300;
my $data = sequence($n)->dummy(0, 3)->sever;      # [3, 300], each [.,s] == s

# ---------------------------------------------------------------------------
# (A) physical blocks: two waveform blocks via block_meta
# ---------------------------------------------------------------------------
my $rec = {
    data       => $data->float,
    fs         => $fs,
    labels     => [qw(a b c)],
    ch_indices => [1, 2, 3],
    events     => [ { t => 0.5, label => 'X' },   # s=50  -> block 0
                    { t => 1.5, label => 'Y' },   # s=150 -> block 1
                    { t => 2.5, label => 'Z' } ], # s=250 -> block 1
    t_start    => '2025-12-21 16:43:30',
    block_meta => [ { start_samp => 0,   n_samp => 150, t_start => '2025-12-21 16:43:30' },
                    { start_samp => 150, n_samp => 150, t_start => '2025-12-21 16:45:10' } ],
    n_blocks   => 2,
};

my $ranges = block_ranges($rec);
is scalar @$ranges, 2, 'two physical blocks';
is $ranges->[1]{start}, 150, 'block 1 start sample';
is $ranges->[1]{end},   300, 'block 1 end sample';
is $ranges->[1]{t_start}, '2025-12-21 16:45:10', 'block 1 per-block t_start';

my $b0 = select_block($rec, 0);
is $b0->{data}->dim(1), 150, 'block 0 length';
is $b0->{data}->at(0, 0),   0,   'block 0 starts at sample 0';
is $b0->{data}->at(0, 149), 149, 'block 0 ends at sample 149';
is scalar @{ $b0->{events} }, 1, 'block 0 has 1 event';
is $b0->{events}[0]{label}, 'X', 'block 0 event = X';
is $b0->{t_start}, '2025-12-21 16:43:30', 'block 0 t_start';
is $b0->{n_blocks}, 1, 'sub-record reports single block';

my $b1 = select_block($rec, 1);
is $b1->{data}->at(0, 0), 150, 'block 1 starts at original sample 150';
is scalar @{ $b1->{events} }, 2, 'block 1 has 2 events';
is $b1->{events}[0]{label}, 'Y', 'first is Y';
ok abs($b1->{events}[0]{t} - 0.0) < 1e-9, 'Y rebased to t=0';
ok abs($b1->{events}[1]{t} - 1.0) < 1e-9, 'Z rebased to t=1.0';
is $b1->{ch_indices}[0], 1, 'ch_indices carried through';

eval { select_block($rec, 5) };
like $@, qr/out of range/, 'select_block croaks on bad index';

# ---------------------------------------------------------------------------
# (B) single physical block -> segment by .LOG "REC START" markers
# ---------------------------------------------------------------------------
my $rec2 = {
    data       => $data->float,
    fs         => $fs,
    labels     => [qw(a b c)],
    t_start    => '2025-12-21 16:43:30',
    events     => [ { t => 0,   label => 'REC START IIA EEG' },
                    { t => 1.0, label => 'REC START IIA EEG' },
                    { t => 0.5, label => 'stim' } ],
    block_meta => [ { start_samp => 0, n_samp => 300, t_start => '2025-12-21 16:43:30' } ],
    n_blocks   => 1,
};

my $lr = block_ranges($rec2);
is scalar @$lr, 2, 'two .LOG segments';
is $lr->[0]{start}, 0,   'segment 0 start';
is $lr->[0]{end},   100, 'segment 0 end (next REC START at 1.0 s)';
is $lr->[1]{start}, 100, 'segment 1 start';
is $lr->[1]{end},   300, 'segment 1 end (to EOF)';
is $lr->[1]{t_start}, '2025-12-21 16:43:31', 'segment 1 t_start = base + 1 s';

done_testing();
