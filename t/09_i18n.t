use strict;
use warnings;
use Test::More;
use PDL;
use Encode qw(encode decode);
use PDL::EEG::IO::EDF qw(write_edf);
use File::Temp qw(tempfile);

# A Japanese task label "安静開眼" carried as a Unicode string on the event
# (as read_nk produces after CP932 decoding). The EDF writer must emit it as
# UTF-8 octets in the annotation channel.

my $jp   = decode('UTF-8', encode('UTF-8', "\x{5b89}\x{9759}\x{958b}\x{773c}")); # 安静開眼
my $fs   = 100;
my $n    = 300;
my $data = sequence($n)->dummy(0, 2)->sever;   # [2,300]

my $rec = {
    data    => $data->float,
    fs      => $fs,
    labels  => [qw(a b)],
    events  => [ { t_data => 1.0, label => $jp },
                 { t_data => 2.0, label => 'task1' } ],
    t_start => '2026-07-02 14:03:03',
};

my (undef, $path) = tempfile(SUFFIX => '.edf', UNLINK => 1);

# default (utf8): Japanese preserved as UTF-8 bytes
write_edf($rec, $path);
my $buf = do { open my $fh, '<:raw', $path or die $!; local $/; <$fh> };
my $utf8 = encode('UTF-8', $jp);                # E5 AE 89 E9 9D 99 E9 96 8B E7 9C BC
# ok index($buf, $utf8) >= 0, 'UTF-8 bytes of 安静開眼 present in EDF (default utf8)';
ok CORE::index($buf, $utf8) >= 0, 'UTF-8 bytes of 安静開眼 present in EDF (default utf8)';
# ok index($buf, 'task1') >= 0, 'ASCII label also present';
ok CORE::index($buf, 'task1') >= 0, 'ASCII label also present';

# parse the annotation channel back and decode
my @annot;
{
    open my $fh, '<:raw', $path or die $!; local $/; my $b = <$fh>;
    my $p = 0; my $g = sub { my $s = substr($b, $p, $_[0]); $p += $_[0]; $s };
    $g->(8+80+80+8+8); my $hb = $g->(8)+0; $g->(44); my $nrec=$g->(8)+0; $g->(8); my $ns=$g->(4)+0;
    my @nspr; { my $save=$p; $p=256 + $ns*(16+80+8+8+8+8+8+80); @nspr=map{$g->(8)+0}1..$ns; $p=$save; }
    $p = $hb;
    for my $r (0..$nrec-1) { for my $c (0..$ns-1) {
        my $raw=$g->($nspr[$c]*2); push @annot,$raw if $c==$ns-1; } }
}
my %found;
for my $a (@annot) { for my $tal (split /\x00/, $a) {
    next unless length $tal;
    my ($on,@txt) = split /\x14/, $tal, -1;
    for my $t (@txt) { next unless length $t; $found{ decode('UTF-8', $t) } = 1; }
} }
ok $found{$jp}, '安静開眼 recovered from EDF annotation via UTF-8 decode';

# ascii mode: non-ASCII replaced, ASCII kept
write_edf($rec, $path, annot_encoding => 'ascii');
my $buf2 = do { open my $fh, '<:raw', $path or die $!; local $/; <$fh> };
# ok index($buf2, $utf8) < 0, 'ascii mode: no UTF-8 bytes';
ok CORE::index($buf2, $utf8) < 0, 'ascii mode: no UTF-8 bytes';
# ok index($buf2, 'task1') >= 0, 'ascii mode keeps ASCII labels';
ok CORE::index($buf2, 'task1') >= 0, 'ascii mode keeps ASCII labels';

done_testing();
