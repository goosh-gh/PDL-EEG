#!/usr/bin/env perl
use strict; use warnings;
use Test::More;
use FindBin qw($Bin);
use lib "$Bin/../lib";

use_ok 'PDL::EEG::IO::NihonKohden', 'read_nk';

my $eeg = "$Bin/data/test01.eeg";
my $log = "$Bin/data/test01.LOG";

# Auto-generate synthetic data if missing
unless (-f $eeg) {
    my $gen = "$Bin/mk_synthetic_nk.pl";
    # Run from the t/ directory so t/data/ is created correctly
    my $orig_dir = do { require Cwd; Cwd::getcwd() };
    chdir $Bin;
    system($^X, $gen) == 0 or die "Failed to generate test data: $gen";
    chdir $orig_dir;
}

ok -f $eeg, "test01.eeg exists";
ok -f $log, "test01.LOG exists";

SKIP: {
    skip "test01.eeg missing", 20 unless -f $eeg;

    # EEG-1100C requires fs to be supplied
    my $rec = eval { read_nk($eeg) };   # fs now auto-detected from header
    is $@, '', "read_nk lives";
    ok defined $rec, "read_nk returns value";
    isa_ok $rec, 'HASH';

    # Metadata
    is  $rec->{fs},           1000, 'fs = 1000 Hz';
    is  $rec->{n_ch_valid},   4,    'n_ch_valid = 4';
    is  $rec->{t_start}, '2025-12-21 16:43:30', 'BCD timestamp decoded';
    like $rec->{labels}[0], qr/FP1/i, 'label[0] ~ FP1';
    like $rec->{labels}[1], qr/FP2/i, 'label[1] ~ FP2';
    is  $rec->{labels}[-1], 'PAD',    'last label = PAD';
    ok  $rec->{n_blocks} >= 1,        'n_blocks >= 1';

    # PDL data: n_ch = n_ch_valid+1 = 5, n_samp = 1000
    my $data = $rec->{data};
    isa_ok $data, 'PDL', 'data is PDL';
    is $data->ndims,  2,    'data 2-D';
    is $data->dim(0), 5,    'dim(0) = n_ch = 5';
    is $data->dim(1), 1000, 'dim(1) = n_samples = 1000';

    # Valid channels should have signal, PAD channel = 0
    my $pad_ch  = $data->slice('(-1),:');
    my $max_pad = $pad_ch->abs->max;
    is $max_pad, 0, 'PAD channel is all zeros';

    my $ch0_rms = sqrt(($data->slice('(0),:')->pow(2))->avg);
    ok $ch0_rms > 0, "ch0 has signal (rms=$ch0_rms µV)";

    # Data in plausible µV range (±200µV synthetic sine)
    my $absmax = $data->slice('0:-2,:')->abs->max;
    ok $absmax < 500 && $absmax > 0, "data range plausible (absmax=$absmax µV)";

    # Events
    my @evts = @{ $rec->{events} };
    is scalar @evts, 2,           '2 events';
    is $evts[0]{t},  10,          'event[0].t=10s';
    like $evts[0]{label}, qr/EYES/i, 'event[0] label';

    # Error handling
    eval { read_nk('/no/such/file.eeg') };
    like $@, qr/not found/i, 'croak on missing file';

    # fs auto-detected from header (no longer needs fs option)
    ok $rec->{fs} > 0, 'fs auto-detected from header';
}

done_testing;
