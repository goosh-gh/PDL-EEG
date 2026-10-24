#!/usr/bin/env perl
# t/34_brainvision.t -- PDL::EEG::IO::BrainVision::read_brainvision on small
# files written here from the BrainVision Core Data Format description
# (header INI sections, multiplexed / vectorized binary, ChN=name,ref,resolution,unit,
# MkN=type,description,position,size,channel with 1-based positions).

use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use PDL;
use PDL::EEG::IO::BrainVision qw(read_brainvision);

my $dir = tempdir(CLEANUP => 1);

# write one recording; $vals is a Perl array [channel][sample] in file units
sub write_bv {
    my (%o) = @_;
    my $base = "$dir/$o{name}";
    my @ch   = @{ $o{channels} };                       # [name, resolution, unit]
    my $n    = @{ $o{vals}[0] };
    open my $h, '>:raw', "$base.vhdr" or die;
    print $h "Brain Vision Data Exchange Header File Version 1.0\r\n; a comment\r\n\r\n";
    print $h "[Common Infos]\r\nCodepage=UTF-8\r\nDataFile=$o{name}.eeg\r\nMarkerFile=$o{name}.vmrk\r\n";
    print $h "DataFormat=BINARY\r\nDataOrientation=$o{orient}\r\nNumberOfChannels=" . @ch . "\r\n";
    print $h "SamplingInterval=$o{si}\r\n\r\n[Binary Infos]\r\nBinaryFormat=$o{fmt}\r\n\r\n[Channel Infos]\r\n";
    for my $k (0 .. $#ch) {
        my ($nm, $res, $u) = @{ $ch[$k] };
        print $h "Ch" . ($k + 1) . "=$nm,,$res,$u\r\n";
    }
    close $h;
    my $pk = { INT_16 => 's<', IEEE_FLOAT_32 => 'f<', INT_32 => 'l<' }->{ $o{fmt} };
    open my $d, '>:raw', "$base.eeg" or die;
    if ($o{orient} eq 'MULTIPLEXED') {
        for my $i (0 .. $n - 1) { print $d pack("$pk*", map { $o{vals}[$_][$i] } 0 .. $#ch) }
    } else {
        for my $c (0 .. $#ch) { print $d pack("$pk*", @{ $o{vals}[$c] }) }
    }
    close $d;
    open my $m, '>:raw', "$base.vmrk" or die;
    print $m "Brain Vision Data Exchange Marker File, Version 1.0\r\n\r\n[Common Infos]\r\nDataFile=$o{name}.eeg\r\n\r\n";
    print $m "[Marker Infos]\r\n";
    my $k = 1;
    print $m "Mk" . $k++ . "=$_->[0],$_->[1],$_->[2],1,0\r\n" for @{ $o{markers} // [] };
    close $m;
    return "$base.vhdr";
}

# 1. INT_16, multiplexed, resolutions and units
{
    my @v = ([ 1, -2, 3, 1000, -32768 ], [ 10, 20, 30, 40, 50 ], [ 7, 7, 7, 7, 7 ]);
    my $f = write_bv(name => 'a', orient => 'MULTIPLEXED', fmt => 'INT_16', si => 2000,
                     channels => [ [ 'Fp1', 0.1, "\xB5V" ], [ 'ECG', 0.5, 'mV' ], [ 'X\1Y', 2, 'V' ] ],
                     vals => \@v,
                     markers => [ [ 'Stimulus', 'S  1', 1 ], [ 'Comment', 'eyes\1 closed', 4 ] ]);
    my $r = read_brainvision($f);
    is($r->{fs}, 500, 'fs from SamplingInterval (us)');
    is_deeply($r->{labels}, [ 'Fp1', 'ECG', 'X,Y' ], 'labels, \\1 -> comma');
    is_deeply([ $r->{data}->dims ], [ 3, 5 ], 'dims (n_ch, n_samp)');
    my $want = pdl([ map { 0.1 * $_ } @{ $v[0] } ], [ map { 0.5e3 * $_ } @{ $v[1] } ], [ map { 2e6 * $_ } @{ $v[2] } ])->transpose;
    ok(($r->{data} - $want)->abs->max->sclr < 1e-6, 'every sample scaled to uV (uV, mV, V)');
    is(scalar @{ $r->{events} }, 2, 'two markers');
    is($r->{events}[0]{onset}, 0, 'marker at position 1 -> 0 s');
    is($r->{events}[1]{onset}, 3 / 500, 'marker at position 4 -> 3 samples');
    is($r->{events}[1]{label}, 'Comment/eyes, closed', 'marker label type/description');
}

# 2. IEEE_FLOAT_32, vectorized, against the same numbers multiplexed
{
    my @v = map { my $c = $_; [ map { sin($_ * 0.37 + $c) * 100 } 0 .. 99 ] } 0 .. 3;
    my @ch = map { [ "C$_", 1, 'uV' ] } 1 .. 4;
    my $fv = write_bv(name => 'v', orient => 'VECTORIZED',  fmt => 'IEEE_FLOAT_32', si => 1000, channels => \@ch, vals => \@v);
    my $fm = write_bv(name => 'm', orient => 'MULTIPLEXED', fmt => 'IEEE_FLOAT_32', si => 1000, channels => \@ch, vals => \@v);
    my ($rv, $rm) = map { read_brainvision($_) } $fv, $fm;
    is_deeply([ $rv->{data}->dims ], [ 4, 100 ], 'vectorized: dims (n_ch, n_samp)');
    ok(($rv->{data} - $rm->{data})->abs->max->sclr == 0, 'vectorized and multiplexed give the same data');
    my $want = pdl(float, \@v)->transpose->double;           # float32 round trip
    ok(($rm->{data} - $want)->abs->max->sclr == 0, 'every float32 sample read exactly');
    is($rm->{fs}, 1000, 'fs 1000 Hz');
    is(scalar @{ $rm->{events} }, 0, 'no markers');

    # channel subset, in the order asked, for both orientations
    for my $pair ([ 'multiplexed', $fm ], [ 'vectorized', $fv ]) {
        my $rs = read_brainvision($pair->[1], channels => [ 'C3', 'C1' ]);
        is_deeply($rs->{labels}, [ 'C3', 'C1' ], "$pair->[0] subset: labels in the order asked");
        ok(($rs->{data} - $want->dice_axis(0, pdl(long, [ 2, 0 ])))->abs->max->sclr == 0,
           "$pair->[0] subset: every sample of the channels kept");
    }
    eval { read_brainvision($fm, channels => [ 'C9' ]) };
    like($@, qr/no channel 'C9'/, 'subset: unknown channel croaks');
}

# 2b. multiplexed file longer than one read block (65536 samples), subset
{
    my $n  = 70000;
    my @v  = map { my $c = $_; [ map { (($_ * 7 + $c * 13) % 2001) - 1000 } 0 .. $n - 1 ] } 0 .. 2;
    my $f  = write_bv(name => 'long', orient => 'MULTIPLEXED', fmt => 'INT_16', si => 1000,
                      channels => [ map { [ "L$_", 1, 'uV' ] } 1 .. 3 ], vals => \@v);
    my $r  = read_brainvision($f, channels => [ 'L3', 'L2' ]);
    my $wn = pdl(\@v)->transpose->dice_axis(0, pdl(long, [ 2, 1 ]));
    is_deeply([ $r->{data}->dims ], [ 2, $n ], 'long file subset: dims');
    ok(($r->{data} - $wn)->abs->max->sclr == 0, 'long file subset: every sample, across the block boundary');
}

# 3. negative control: an unsupported format is refused
{
    my $f = write_bv(name => 'bad', orient => 'MULTIPLEXED', fmt => 'INT_16', si => 1000,
                     channels => [ [ 'A', 1, 'uV' ] ], vals => [ [ 1, 2 ] ]);
    open my $h, '<:raw', $f or die; my $txt = do { local $/; <$h> }; close $h;
    $txt =~ s/BinaryFormat=INT_16/BinaryFormat=IEEE_FLOAT_64/;
    open $h, '>:raw', $f or die; print $h $txt; close $h;
    eval { read_brainvision($f) };
    like($@, qr/not supported/, 'unsupported BinaryFormat croaks');
}

done_testing;
