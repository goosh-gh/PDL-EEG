use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::IO::BESA::ASCII qw(write_mul);
use File::Temp qw(tempfile);
# [n_ch=3, n_samples=2]; dim0 = channel, dim1 = time. Last channel is Trigger.
# In pdl([...]) the inner list is dim0, so each row here is one time point
# (all channels), which is exactly dims (3,2) = [n_ch, n_samples]. No transpose.
my $data   = pdl([ [9.34,  7.03, 0],     # t=0 : ch0 ch1 trig
                   [11.65, 7.55, 0] ]);  # t=1
is($data->dim(0), 3, 'n_ch');
is($data->dim(1), 2, 'n_samples');
my @labels = ('Fp1', 'Fp2', 'Trigger');
my $rec = { data => $data, fs => 1000, labels => \@labels,
            t_start => '2026-07-05 16:44:34' };
my (undef, $path) = tempfile(SUFFIX => '.mul', UNLINK => 1);
write_mul($rec, $path, trig_width => 4);
open my $in, '<', $path or die $!;
my @L = <$in>;
close $in;
# The last channel is Trigger, which is written as a column but NOT counted in
# Channels= (it is dropped from later analysis), so 3 channels -> Channels=2.
like($L[0], qr/^TimePoints=2 Channels=2 /,           'header: counts (Trigger not counted)');
like($L[0], qr/BeginSweep\[ms\]=0\.00 /,              'header: begin sweep');
like($L[0], qr/SamplingInterval\[ms\]=1\.000 /,       'header: interval (1000 Hz)');
like($L[0], qr{Bins/uV=1\.000 },                      'header: bins/uV');
like($L[0], qr/Time=16:44:34/,                        'header: time from t_start');
is  ($L[1], " Fp1 Fp2 Trigger\n",                     'label line, leading space');
is  ($L[2], "    9.34     7.03    0\n",               'row 0: float chans + int trigger');
is  ($L[3], "   11.65     7.55    0\n",               'row 1');
# Trigger auto-detection off when trigger => undef: all columns float.
write_mul($rec, $path, trigger => undef);
open $in, '<', $path or die $!;
@L = <$in>;
close $in;
like($L[2], qr/    0\.00\n$/, 'trigger => undef writes float 0.00');

# Trigger is exported as a data column but excluded from Channels= (the label
# row and data still include it); count_trigger => 1 restores the full count.
# [P:EEG11 convention]
{
    my (undef, $p) = tempfile(SUFFIX => '.mul', UNLINK => 1);
    write_mul($rec, $p, trig_width => 4);
    open my $fh, '<', $p or die $!;
    my @H = <$fh>;
    close $fh;
    like($H[0], qr/^TimePoints=2 Channels=2 /, 'trigger excluded from Channels=');
    is  ($H[1], " Fp1 Fp2 Trigger\n",          'label row still lists Trigger');
    is  (scalar(split ' ', $H[1]), 3,          'label row keeps all 3 tokens');

    write_mul($rec, $p, trig_width => 4, count_trigger => 1);
    open $fh, '<', $p or die $!;
    my @H2 = <$fh>;
    close $fh;
    like($H2[0], qr/^TimePoints=2 Channels=3 /, 'count_trigger => 1 counts it');
}

done_testing;
