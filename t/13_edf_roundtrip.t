#!/usr/bin/env perl
# t/13_edf_roundtrip.t
#
# (Delivered as its own file: a fragment dropped into t/ is picked up by
# `make test` via t/*.t and has to compile on its own. If you would rather have
# it inside t/02_edf.t, paste everything below the `use` lines in there instead.)

use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);

use PDL;
use PDL::EEG::IO::NihonKohden qw(read_nk);
use PDL::EEG::IO::EDF         qw(write_edf read_edf);

# ---------------------------------------------------------------------------
# EDF round trip: uV in, uV out -- INCLUDING the DC channels.
#
# The bug this pins down: an EEG-1200A DC input is a +/-12 V line. The vendor
# quotes its range as +/-12002.9 *mV*, and read_nk used that figure as though it
# were uV/bit, so every DC channel came out 1000x too small -- a 3.3 V trigger
# read as a 3.3 mV wobble. {data} is now uV for every channel.
#
# But a DC channel CANNOT be written to EDF in uV: +/-12002913 needs nine
# characters and EDF's physical_min field is eight. So write_edf uses
# $rec->{units} to give each signal its own physical dimension -- uV for EEG, mV
# for DC -- and read_edf normalises back to uV on the way in. If either half of
# that is missing, this test catches it, because the numbers come back 1000x out.
#
# t/data/test03.eeg puts a square-wave TTL of 9000 ADC counts on the DC channels.
# 9000 * 366.3 uV/bit = 3.30 V. That is what has to survive the round trip.
# ---------------------------------------------------------------------------
SKIP: {
    my $F3 = 't/data/test03.eeg';
    skip "no $F3 (run: perl t/mk_synthetic_nk.pl)", 14 unless -f $F3;

    my $rec = read_nk($F3, all_blocks => 1);
    my @lab = @{ $rec->{labels} };
    my %ix; $ix{ $lab[$_] } = $_ for 0 .. $#lab;

    # hw 45-48. With no .21e, read_nk falls back to the EEG-1200A DC numbering
    # for extblock, so these are DC01-DC04 (the electrode-code table would say
    # DC03-DC06, which is the EEG-1100C panel).
    my $dc  = $ix{DC01};
    my $eeg = $ix{Fp1};
    ok defined $dc,  'fixture has DC01 (EEG-1200A DC numbering)';
    ok defined $eeg, 'fixture has an EEG channel';

    # --- read_nk side: is the DC channel in VOLTS-worth of microvolts? -------
    is $rec->{units}[$dc],  'mV', 'read_nk marks DC for mV export';
    is $rec->{units}[$eeg], 'uV', 'read_nk marks EEG for uV export';

    my $ttl = $rec->{data}->slice("($dc),:")->max - $rec->{data}->slice("($dc),:")->min;
    ok abs($ttl - 9000 * 366.29984) < 5000,
        sprintf('read_nk: the DC TTL is %.2f V, not %.2f mV', $ttl / 1e6, $ttl / 1e3);
    ok $ttl > 3_000_000, 'read_nk: DC data are uV (a 3.3 V pulse is ~3.3e6), not mV';

    my $amp = $rec->{data}->slice("($eeg),:")->max - $rec->{data}->slice("($eeg),:")->min;
    ok $amp > 10 && $amp < 100, "read_nk: the EEG channel is still tens of uV ($amp)";

    # --- write, then read back ----------------------------------------------
    my $dir = tempdir(CLEANUP => 1);
    my $edf = "$dir/rt.edf";
    write_edf($rec, $edf, phys => 'gain');
    ok -s $edf, 'write_edf produced a file';

    # the physical dimension must differ PER SIGNAL, or the DC range cannot fit
    open my $fh, '<:raw', $edf or die $!;
    read $fh, my $hdr, 256;
    my $ns = 0 + substr($hdr, 252, 4);
    read $fh, my $sig, 256 * $ns;
    my $dim = sub { my $i = shift; my $d = substr($sig, 96 * $ns + 8 * $i, 8); $d =~ s/\s+$//; $d };
    my $pmn = sub { my $i = shift; my $d = substr($sig, 104 * $ns + 8 * $i, 8); $d =~ s/\s+$//; $d };
    close $fh;

    is $dim->($eeg), 'uV', 'EDF signal header: EEG dimension is uV';
    is $dim->($dc),  'mV', 'EDF signal header: DC dimension is mV';
    ok length($pmn->($dc)) <= 8,
        'EDF signal header: the DC physical_min fits the 8-char field '
        . "(got '" . $pmn->($dc) . "'; in uV it would be -12002913, which is 9)";

    my $back = read_edf($edf);
    is scalar @{ $back->{labels} }, scalar @lab, 'read_edf: channel count survives';
    is $back->{units}[$dc], 'mV', 'read_edf reports the dimension the signal was stored in';

    # --- THE assertion: the values came back in uV, both kinds ---------------
    for my $c ($eeg, $dc) {
        my $a = $rec->{data}->slice("($c),:");
        my $b = $back->{data}->slice("($c),:");
        my $tol = 2 * ($rec->{gains}->at($c));      # one quantisation step, doubled
        my $err = ($a - $b)->abs->max;
        ok $err < $tol,
            sprintf('round trip: %s survives in uV (max err %.3g, tol %.3g)',
                    $lab[$c], $err, $tol);
    }
}

done_testing();
