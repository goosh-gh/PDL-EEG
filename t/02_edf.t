use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::IO::EDF qw(write_edf read_edf);
use File::Temp qw(tempfile);
# ---------------------------------------------------------------------------
# Build a small in-memory record (no reader dependency), write EDF+C, then
# parse the header + annotations back with a tiny pure-Perl reader.
# ---------------------------------------------------------------------------
my $fs   = 100;
my $n_s  = 250;                                  # 2.5 s -> 3 records (last padded)
my $t    = sequence($n_s) / $fs;
my $data = pdl(
    (100 * sin($t * 6))->list,                   # Fp1  ~ +/-100 uV
    (($t * 40) - 20)->list,                      # Cz   ramp
    (5 * ones($n_s))->list,                      # DC01 flat 5 uV
)->reshape($n_s, 3)->xchg(0, 1)->sever;          # [3, n_s]
my $rec = {
    data    => $data->float,
    fs      => $fs,
    labels  => [qw(Fp1 Cz DC01)],
    events  => [ { time => 0.5, label => 'S1' }, [ 2.1, 'Resp' ] ],
    t_start => '2025-12-21 16:43:30',            # exercises date/time conversion
};
my (undef, $path) = tempfile(SUFFIX => '.edf', UNLINK => 1);
write_edf($rec, $path, phys => 'gain');
# --- tiny EDF reader --------------------------------------------------------
open my $fh, '<:raw', $path or die $!;
local $/; my $buf = <$fh>; close $fh;
my $p = 0;
my $g = sub { my $s = substr($buf, $p, $_[0]); $p += $_[0]; $s =~ s/\s+$//; $s };
my %h;
$h{ver}   = $g->(8);
$g->(80) for 1 .. 2;                             # subject or patient, recording
$h{sdate} = $g->(8);
$h{stime} = $g->(8);
$h{hbytes}   = $g->(8) + 0;
$h{reserved} = $g->(44);
$h{nrec}     = $g->(8) + 0;
$h{rdur}     = $g->(8) + 0;
$h{ns}       = $g->(4) + 0;
my $ns = $h{ns};
my @lab  = map { $g->(16) } 1 .. $ns;
$g->(80) for 1 .. $ns;
$g->(8)  for 1 .. $ns;                           # dim
my @pmn  = map { $g->(8) + 0 } 1 .. $ns;
my @pmx  = map { $g->(8) + 0 } 1 .. $ns;
$g->(8)  for 1 .. $ns; $g->(8) for 1 .. $ns;      # dmin, dmax
$g->(80) for 1 .. $ns;
my @nspr = map { $g->(8) + 0 } 1 .. $ns;
$g->(32) for 1 .. $ns;
is($p, $h{hbytes}, 'header length matches declared byte count');
is($h{hbytes}, 256 * (1 + $ns), 'header = 256*(1+ns)');
is($ns, 4, '3 signals + 1 annotation channel');
is($h{nrec}, 3, '250 samples / 100 spr -> 3 padded records');
is($h{reserved}, 'EDF+C', 'reserved field = EDF+C');
is($lab[3], 'EDF Annotations', 'last signal is annotation channel');
ok(abs($pmx[0] - 3199.902) < 0.01, 'gain-mode physical max ~3199.902');
# start date/time from t_start
is($h{sdate}, '21.12.25', 'startdate dd.mm.yy from t_start');
is($h{stime}, '16.43.30', 'starttime hh.mm.ss from t_start');
# --- data records + annotations --------------------------------------------
my (@sig, @annot);
for my $r (0 .. $h{nrec} - 1) {
    for my $c (0 .. $ns - 1) {
        my $raw = substr($buf, $p, $nspr[$c] * 2); $p += $nspr[$c] * 2;
        if ($c == $ns - 1) { push @annot, $raw }
        else { push @{ $sig[$c] }, unpack('s<*', $raw) }
    }
}
is($p, length($buf), 'consumed whole file');
my $gain = 0.09765625;
ok(abs($sig[0][30] * $gain - 100 * sin(0.30 * 6)) < 2 * $gain,
   'Fp1 sample 30 round-trips within ~1 LSB');
ok(abs($sig[2][10] * $gain - 5) < 2 * $gain, 'flat DC01 preserved at 5 uV');
ok(abs($sig[0][299] * $gain) < 2 * $gain, 'padded tail is ~0 uV');
my %ev;
for my $a (@annot) {
    for my $tal (split /\x00/, $a) {
        next unless length $tal;
        my ($onhdr, @txts) = split /\x14/, $tal, -1;
        (my $on = (split /\x15/, $onhdr)[0]) =~ s/^\+//;
        $ev{$_} = 0 + $on for grep { length } @txts;
    }
}
ok(abs(($ev{S1}   // -9) - 0.5) < 1e-6, 'event S1 @ 0.5 s');
ok(abs(($ev{Resp} // -9) - 2.1) < 1e-6, 'event Resp @ 2.1 s');

# ---------------------------------------------------------------------------
# read_edf round-trip: parse the file we just wrote with the REAL reader and
# confirm it honours the read_nk contract that edf_to_mul.pl depends on:
#   data => PDL[n_ch,n_samp] float uV, fs, labels (EDF Annotations excluded),
#   t_start "YYYY-MM-DD HH:MM:SS", events => [{onset,label},...].
# This is the piece the tiny in-line reader above does NOT cover.
# ---------------------------------------------------------------------------
{
    my $back = read_edf($path);
    is($back->{fs}, 100, 'read_edf: fs = spr / record_dur');
    is_deeply($back->{labels}, [qw(Fp1 Cz DC01)],
              'read_edf: data labels only (EDF Annotations excluded)');
    is($back->{data}->dim(0), 3, 'read_edf: 3 data channels (annotation dropped)');
    is($back->{data}->dim(1), 300,
       'read_edf: 3 records x 100 spr = 300 samples (includes padded tail)');
    is($back->{t_start}, '2025-12-21 16:43:30',
       'read_edf: t_start round-trips (full year from Startdate field)');
    is($back->{edf_type}, 'EDF+C', 'read_edf: edf_type = EDF+C');

    my $lsb = 0.09765625;
    ok(abs($back->{data}->at(0, 30) - 100 * sin(0.30 * 6)) < 2 * $lsb,
       'read_edf: Fp1 sample 30 within ~1 LSB');
    ok(abs($back->{data}->at(2, 10) - 5) < 2 * $lsb,
       'read_edf: flat DC01 preserved at 5 uV');
    ok(abs($back->{data}->at(0, 299)) < 2 * $lsb,
       'read_edf: padded tail ~0 uV');

    my %on = map { $_->{label} => $_->{onset} } @{ $back->{events} };
    ok(abs(($on{S1}   // -9) - 0.5) < 1e-6, 'read_edf: event S1 @ 0.5 s');
    ok(abs(($on{Resp} // -9) - 2.1) < 1e-6, 'read_edf: event Resp @ 2.1 s');
}
done_testing();
