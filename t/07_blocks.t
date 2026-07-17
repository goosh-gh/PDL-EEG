use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::IO::NihonKohden qw(block_ranges select_block block_extents read_nk);

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

# ---------------------------------------------------------------------------
# (C) block_extents(): plan a partial read from the FILE, without reading samples
#
# block_ranges() above takes a RECORD -- the data must already be in memory.
# block_extents() takes a PATH and reads only the control-block address table
# and the per-block headers, so a caller can decide which blocks a time range
# touches before loading anything. Its coordinates must line up, sample for
# sample, with read_nk(all_blocks => 1). If they ever drift apart, a --cut range
# read off the viewer stops matching what nk_to_mul.pl writes out.
#
# Since 0.2 both default to gap_samples => 0: blocks are BUTT-JOINED and no
# synthetic samples enter the data. gap_samples > 0 is still honoured (and still
# tested here), but it is deprecated -- zeros are not data.
#
# Fixture: t/data/test02.eeg (3 blocks, 1000/2500/700 samples -- unequal on
# purpose, so uniform-length assumptions and off-by-one errors cannot hide).
# ---------------------------------------------------------------------------
my $F2 = 't/data/test02.eeg';

SKIP: {
    skip "no $F2 (run: perl t/mk_synthetic_nk.pl)", 45 unless -f $F2;

    my $GAP = 100;                          # deprecated, but must still work
    my $ext = block_extents($F2, gap_samples => $GAP);

    my @LEN  = (1000, 2500, 700);
    my @TIME = ('2025-12-21 16:43:30', '2025-12-21 16:44:30', '2025-12-21 16:50:00');

    is ref $ext,     'ARRAY',  'block_extents returns an arrayref';
    is scalar @$ext, 3,        'one entry per waveform block';
    is_deeply [ map { $_->{index}      } @$ext ], [0,1,2],       'index 0-based, ordered';
    is_deeply [ map { $_->{n_samp}     } @$ext ], \@LEN,         'n_samp derived from address gaps';
    is_deeply [ map { $_->{t_start}    } @$ext ], \@TIME,        'per-block t_start (BCD) decoded';
    is_deeply [ map { $_->{fs}         } @$ext ], [(1000) x 3],  'fs decoded (14-bit mask applied)';
    is_deeply [ map { $_->{n_ch}       } @$ext ], [(5) x 3],     'n_ch includes the pad channel';
    is_deeply [ map { $_->{n_ch_valid} } @$ext ], [(4) x 3],     'n_ch_valid excludes it';

    is_deeply [ map { $_->{start_samp} } @$ext ], [0, 1100, 3700],
        'start_samp includes the gap_samples padding between blocks';
    is_deeply [ map { $_->{end_samp} } @$ext ], [1000, 3600, 4400],
        'end_samp = start_samp + n_samp (the gap is NOT part of the block)';

    # --- the DEFAULT: butt-joined, no synthetic samples ----------------------
    my $extd = block_extents($F2);                          # no gap_samples
    is_deeply [ map { $_->{start_samp} } @$extd ], [0, 1000, 3500],
        'default gap_samples => 0: blocks are butt-joined';
    is_deeply [ map { $_->{end_samp} } @$extd ], [1000, 3500, 4200],
        'default: end_samp is contiguous with the next start_samp';

    my $alld = read_nk($F2, all_blocks => 1);               # no gap_samples
    is $alld->{data}->dim(1), 4200,
        'default: concatenated length == sum of n_samp (no padding)';
    is_deeply $alld->{gap_bounds}, [],
        'default: gap_bounds is empty -- there are no synthetic samples to report';
    is_deeply $alld->{t_block_starts}, [0, 1000, 3500],
        'default: t_block_starts still marks the breaks';
    # Do NOT test this with "!= 0": the fixture's values are
    #     raw = 0x8000 + ((abs*7 + ch*13) % 2000 - 1000)
    # so a sample is LEGITIMATELY 0 uV whenever that lands on 1000. Test the
    # actual thing instead -- that samples 999 and 1000 hold the values belonging
    # to concatenated indices 999 and 1000, i.e. that nothing was inserted.
    my $uv_at = sub {
        my ($abs, $ch) = @_;
        return ((($abs * 7 + $ch * 13) % 2000) - 1000) * 0.09765625;
    };
    for my $k (998, 999, 1000, 1001) {          # 999|1000 straddles the break
        ok abs($alld->{data}->at(0, $k) - $uv_at->($k, 0)) < 1e-3,
            "default: sample $k is concatenated index $k (no padding inserted)";
    }

    # The elapsed wall-clock gap is what the old zero padding stood in for. It is
    # recoverable from the per-block t_start, which is why block_extents returns it:
    #     dt = epoch(t_start[b+1]) - epoch(t_start[b]) - n_samp[b]/fs
    # test02: block 0 starts 16:43:30 and is 1000 samp @1000Hz = 1 s of data;
    #         block 1 starts 16:44:30 -> 60 s apart -> 59 s of real dead time.
    require POSIX;
    my $ep = sub {
        my ($Y,$M,$D,$h,$m,$sec) = $_[0] =~ /^(\d{4})-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})/;
        POSIX::mktime($sec, $m, $h, $D, $M - 1, $Y - 1900, 0, 0, -1);
    };
    my $dt01 = $ep->($extd->[1]{t_start}) - $ep->($extd->[0]{t_start})
             - $extd->[0]{n_samp} / $extd->[0]{fs};
    my $dt12 = $ep->($extd->[2]{t_start}) - $ep->($extd->[1]{t_start})
             - $extd->[1]{n_samp} / $extd->[1]{fs};
    ok abs($dt01 - 59)    < 1e-9, 'elapsed gap at break 0/1 from t_start (60s - 1.0s of data)';
    ok abs($dt12 - 327.5) < 1e-9, 'elapsed gap at break 1/2 from t_start (330s - 2.5s of data)';

    # --- agreement with read_nk(all_blocks => 1): the coordinate contract -----
    my $all = read_nk($F2, all_blocks => 1, gap_samples => $GAP);

    is $all->{n_blocks}, 3,    'read_nk agrees on the block count';
    is $all->{fs},       1000, 'read_nk agrees on fs';
    is $ext->[-1]{end_samp}, $all->{data}->dim(1),
        'last end_samp == total concatenated length (no trailing gap)';
    is_deeply [ map { $_->{start_samp} } @$ext ], $all->{t_block_starts},
        'start_samp == read_nk t_block_starts';
    is_deeply [ map { $_->{n_samp} } @$ext ], $all->{n_samp_per_block},
        'n_samp == read_nk n_samp_per_block';
    is_deeply [ map { { start_samp => $_->{start_samp}, n_samp => $_->{n_samp},
                        t_start => $_->{t_start} } } @$ext ],
              [ map { { start_samp => $_->{start_samp}, n_samp => $_->{n_samp},
                        t_start => $_->{t_start} } } @{ $all->{block_meta} } ],
        'start_samp/n_samp/t_start == read_nk block_meta';
    is_deeply $all->{gap_bounds},
              [ map { [ $ext->[$_]{end_samp}, $ext->[$_+1]{start_samp} - 1 ] } 0 .. 1 ],
        'read_nk gap_bounds fill exactly the space between consecutive extents';

    # --- the point of the whole thing ----------------------------------------
    # A lazy partial read (block_extents -> read_nk(block => N) -> glue -> trim)
    # must reproduce the all_blocks buffer bit for bit, for ANY [lo,hi) -- also
    # when lo or hi lands inside a gap, where all_blocks has zeros.
    # This mirrors what examples/read_nihonkohden.pl does for --cut.
    my $lazy_read = sub {
        my ($lo, $hi) = @_;
        my $e    = block_extents($F2, gap_samples => $GAP);
        my $nch  = $e->[0]{n_ch};
        my @need = grep { $_->{end_samp} > $lo && $_->{start_samp} < $hi } @$e;

        my (@piece, $base);
        if (!@need) {                             # range lies entirely in a gap
            $base  = $lo;
            @piece = (zeroes(float, $nch, $hi - $lo));
        }
        else {
            for my $i (0 .. $#need) {
                if ($i > 0) {                     # re-insert the inter-block gap
                    my $g = $need[$i]{start_samp} - $need[$i-1]{end_samp};
                    push @piece, zeroes(float, $nch, $g) if $g > 0;
                }
                push @piece, read_nk($F2, block => $need[$i]{index})->{data};
            }
            $base = $need[0]{start_samp};
            if ($lo < $base) {                    # starts inside the preceding gap
                unshift @piece, zeroes(float, $nch, $base - $lo);
                $base = $lo;
            }
            if ($hi > $need[-1]{end_samp}) {      # ends inside the following gap
                push @piece, zeroes(float, $nch, $hi - $need[-1]{end_samp});
            }
        }
        my $d = shift @piece;
        $d = $d->glue(1, @piece) if @piece;
        return $d->slice(":," . ($lo - $base) . ":" . ($hi - $base - 1))->sever;
    };

    my @cases = (
        [ 'inside the first block',       10,   900 ],
        [ 'straddling the 0/1 boundary',  950,  1150 ],
        [ 'starting inside a gap',        1020, 1130 ],
        [ 'ending inside a gap',          970,  1040 ],
        [ 'entirely inside a gap',        1010, 1060 ],
        [ 'spanning all three blocks',    5,    4400 ],
        [ 'wholly inside the last block', 3800, 4100 ],
    );
    for my $c (@cases) {
        my ($what, $lo, $hi) = @$c;
        my $got = $lazy_read->($lo, $hi);
        my $ref = $all->{data}->slice(":," . $lo . ":" . ($hi - 1))->sever;
        is $got->dim(1), $hi - $lo, "lazy read length ($what)";
        ok all($got == $ref), "lazy read == all_blocks slice, sample for sample ($what)"
            or diag sprintf('max abs diff = %g', ($got - $ref)->abs->max);
    }

    # planner: read no more blocks than the range actually touches
    my $touch = sub {
        my ($lo, $hi) = @_;
        [ map { $_->{index} }
          grep { $_->{end_samp} > $lo && $_->{start_samp} < $hi } @$ext ];
    };
    is_deeply $touch->(3800, 4100), [2], 'a range inside one block plans a 1-block read';
    is_deeply $touch->(950,  1150), [0,1], 'a boundary-straddling range plans 2 adjacent blocks';
    is_deeply $touch->(1010, 1060), [],  'a range inside a gap plans no block read at all';
}

# ---------------------------------------------------------------------------
# (D) extblock (EEG-1200A): the recorder RE-EMITS the channel-info block into
#     the sample stream at every recording break.
#
# This is the regression test for the worst bug this module has had. read_nk()
# assumed an extblock file was one contiguous data block running to EOF. It is
# not: at every break the recorder writes a fresh 72 + (n_ch-1)*10 = 442-byte
# copy of the channel-info block straight into the samples. The old reader read
# those as EEG, so the channel phase slipped by 442 % 76 = 62 bytes = 31
# channels at each break. From the first break onward EVERY channel label sat on
# another channel's data -- a DC trigger line showed brain signal, an EEG
# electrode showed square waves -- and the sample count was over-reported by
# 442/76 = 5.8 samples per break.
#
# It survived for so long because there was no synthetic extblock fixture: the
# whole layout had only ever been exercised on real recordings, where a
# scrambled montage is not obvious. Hence t/data/test03.eeg.
# ---------------------------------------------------------------------------
my $F3 = 't/data/test03.eeg';

SKIP: {
    skip "no $F3 (run: perl t/mk_synthetic_nk.pl)", 25 unless -f $F3;

    my @LEN = (5000, 3000, 2000);                   # unequal on purpose
    my @TS  = ('2026-07-02 14:03:03',               # + 5 s data, 284 s break
               '2026-07-02 14:07:52',               # + 3 s data, 219 s break
               '2026-07-02 14:11:34');
    my $NCH    = 38;
    my $STRIDE = $NCH * 2;                          # 76
    my $HDRLEN = 72 + ($NCH - 1) * 10;              # 442
    my $TOTAL  = 0; $TOTAL += $_ for @LEN;          # 10000

    is $HDRLEN % $STRIDE, 62,
        'the embedded header is not a whole number of samples (442 % 76 = 62)';

    my $ext = block_extents($F3);
    is scalar @$ext, 3, 'extblock: block_extents finds all three segments';
    is_deeply [ map { $_->{n_samp}  } @$ext ], \@LEN, 'extblock: segment lengths';
    is_deeply [ map { $_->{t_start} } @$ext ], \@TS,
        'extblock: per-segment t_start comes from the embedded header, exactly';
    is_deeply [ map { $_->{start_samp} } @$ext ], [0, 5000, 8000],
        'extblock: start_samp is contiguous -- the headers are NOT samples';
    is_deeply [ map { $_->{n_ch} } @$ext ], [(38) x 3], 'extblock: n_ch';

    my $r = read_nk($F3, all_blocks => 1);
    is $r->{layout},   'extblock', 'extblock: layout dispatch';
    is $r->{n_blocks}, 3,          'extblock: n_blocks counts the segments';
    is $r->{data}->dim(1), $TOTAL,
        'extblock: sample count excludes the embedded headers'
        . " (the old reader reported $TOTAL + " . int(2 * $HDRLEN / $STRIDE) . ')';
    is_deeply $r->{t_block_starts},   [0, 5000, 8000], 'extblock: t_block_starts';
    is_deeply $r->{n_samp_per_block}, \@LEN,          'extblock: n_samp_per_block';
    is_deeply [ map { $_->{t_start} } @{ $r->{block_meta} } ], \@TS,
        'extblock: block_meta carries each segment start time';

    # --- THE test: is each label still on its own channel after a break? -----
    # The fixture puts a square-wave TTL on the DC channels and a sine on every
    # EEG channel. If the phase slips, a DC label shows a sine and an EEG label
    # shows a square wave. Nothing else in this file can produce that.
    my @labels = @{ $r->{labels} };
    my %ix; $ix{ $labels[$_] } = $_ for 0 .. $#labels;

    # Hardware codes 45-48. The electrode-code table numbers those DC03-DC06 --
    # the EEG-1100C front panel. The EEG-1200A panel calls them DC01-DC04, so
    # read_nk applies the 1200A numbering on the extblock path when there is no
    # .21e to say otherwise. This fixture has no .21e, so: DC01-DC04.
    for my $dc (qw(DC01 DC02 DC03 DC04)) {
        my $c = $ix{$dc};
        ok defined $c, "extblock: $dc is present";
        next unless defined $c;
        # sample the LAST segment -- i.e. after TWO breaks, where the old reader
        # was 62 channels adrift
        my $x = $r->{data}->slice("($c),8000:9999");
        my %lv; $lv{ sprintf('%.0f', $x->at($_)) } = 1 for map { $_ * 3 } 0 .. 600;
        ok scalar(keys %lv) <= 4,
            "extblock: $dc is STILL a step function after two breaks "
            . '(' . scalar(keys %lv) . ' levels)';
    }

    for my $eeg (qw(Fp1 Cz)) {
        my $c = $ix{$eeg};
        next unless defined $c;
        my $x = $r->{data}->slice("($c),8000:9999");
        my %lv; $lv{ sprintf('%.0f', $x->at($_)) } = 1 for map { $_ * 3 } 0 .. 600;
        ok scalar(keys %lv) > 10,
            "extblock: $eeg is STILL continuous after two breaks "
            . '(' . scalar(keys %lv) . ' levels)';
    }

    # single-segment reads must work too
    my $b1 = read_nk($F3, block => 1);
    is $b1->{data}->dim(1), $LEN[1], 'extblock: block => 1 reads just that segment';
    is $b1->{t_start}, $TS[1],       'extblock: block => 1 has that segment\'s t_start';
    eval { read_nk($F3, block => 9) };
    like $@, qr/out of range/,       'extblock: block index is range-checked';
}

done_testing();
