package PDL::EEG::IO::BrainVision;

# Read a BrainVision recording (.vhdr header, binary .eeg data, .vmrk markers)
# into the record layout of read_edf / read_nk.

use strict;
use warnings;
use Carp qw(croak carp);
use File::Basename qw(dirname);
use File::Spec;
use PDL;
use Exporter 'import';
our @EXPORT_OK = qw(read_brainvision);

# ---------------------------------------------------------------------------
# $rec = read_brainvision($vhdr_path, %opt)
#
#   channels => [names]   keep only these channels, in this order (unknown names
#                         croak); the data are read in blocks of samples, so a long
#                         recording needs memory only for the channels kept
#   skip_missing => 1     with channels: leave out names that are not in the file
#
# Returns { data => PDL (n_ch, n_samp) double in uV, fs, labels => [...],
#           units => [...] (as written in the header), events => [ {onset => s,
#           label => 'Type/Description'}, ... ], format => 'brainvision' }
#
# Supported: DataFormat=BINARY, DataOrientation MULTIPLEXED or VECTORIZED,
# BinaryFormat INT_16, UINT_16, INT_32 or IEEE_FLOAT_32, little-endian.  Each
# channel's samples are multiplied by its resolution (Channel Infos field 3,
# default 1) and converted to uV from the unit of field 4 (uV/µV, mV, V, nV;
# any other unit, e.g. that of a PPG channel, is left as it is and carped
# about once).  Marker positions are 1-based samples in the file.
# ---------------------------------------------------------------------------
sub read_brainvision {
    my ($vhdr, %o) = @_;
    my $ini  = _read_ini($vhdr);
    my $ci   = $ini->{'Common Infos'} or croak "$vhdr: no [Common Infos]";
    my $bi   = $ini->{'Binary Infos'} // {};
    my $ch   = $ini->{'Channel Infos'} or croak "$vhdr: no [Channel Infos]";
    ($ci->{DataFormat} // 'BINARY') eq 'BINARY'
        or croak "$vhdr: DataFormat $ci->{DataFormat} is not supported (BINARY only)";
    my $orient = $ci->{DataOrientation} // 'MULTIPLEXED';
    my $nch    = $ci->{NumberOfChannels} or croak "$vhdr: NumberOfChannels missing";
    my $si     = $ci->{SamplingInterval} or croak "$vhdr: SamplingInterval missing";
    my $fs     = 1e6 / $si;
    my %type   = (INT_16 => [short, 2], UINT_16 => [ushort, 2], INT_32 => [long, 4],
                  IEEE_FLOAT_32 => [float, 4]);
    my $bf     = $bi->{BinaryFormat} // 'INT_16';
    my ($ptype, $bytes) = @{ $type{$bf} // croak "$vhdr: BinaryFormat $bf is not supported" };

    # channels: ChN=<Name>,<Reference>,<Resolution>,<Unit>  ("\1" stands for a comma)
    my (@lab, @res, @unit);
    for my $k (1 .. $nch) {
        my $v = $ch->{"Ch$k"} // croak "$vhdr: Ch$k missing in [Channel Infos]";
        my @f = map { s/\\1/,/gr } split /,/, $v, -1;
        push @lab,  $f[0];
        push @res,  (defined $f[2] && length $f[2]) ? $f[2] : 1;
        push @unit, (defined $f[3] && length $f[3]) ? $f[3] : 'µV';
    }

    my $dir   = dirname($vhdr);
    my $dfile = File::Spec->catfile($dir, $ci->{DataFile} // croak "$vhdr: DataFile missing");
    my $size  = -s $dfile // croak "$dfile: $!";
    $size % ($bytes * $nch) == 0
        or carp "$dfile: size $size is not a multiple of $nch channels x $bytes bytes; the tail is ignored";
    my $ns    = int($size / ($bytes * $nch));

    # channels kept
    my @keep = 0 .. $nch - 1;
    if (my $want_ch = $o{channels}) {
        my %pos; @pos{@lab} = 0 .. $#lab;
        @keep = $o{skip_missing} ? map { exists $pos{$_} ? $pos{$_} : () } @$want_ch
                                 : map { $pos{$_} // croak "$vhdr: no channel '$_'" } @$want_ch;
        croak "$vhdr: none of the channels asked for is in the file" unless @keep;
    }
    my $kp = pdl(long, \@keep);
    my $raw;
    open my $fh, '<:raw', $dfile or croak "$dfile: $!";
    if ($orient eq 'VECTORIZED') {                               # one channel after another
        $raw = zeroes($ptype, scalar @keep, $ns);
        for my $j (0 .. $#keep) {
            my $row = zeroes($ptype, $ns);
            sysseek($fh, $keep[$j] * $ns * $bytes, 0) or croak "$dfile: seek: $!";
            _read_into($fh, $row, $ns * $bytes, $dfile);
            $raw->slice("($j),:") .= $row;
        }
    } else {                                                     # multiplexed: blocks of samples
        $raw = zeroes($ptype, scalar @keep, $ns);
        my $blk = 65536;
        for (my $i = 0; $i < $ns; $i += $blk) {
            my $m = $i + $blk <= $ns ? $blk : $ns - $i;
            my $b = zeroes($ptype, $nch, $m);
            _read_into($fh, $b, $m * $nch * $bytes, $dfile);
            $raw->slice(":," . $i . ":" . ($i + $m - 1)) .= $b->dice_axis(0, $kp);
        }
    }
    close $fh;
    @lab  = @lab[@keep];
    @res  = @res[@keep];
    @unit = @unit[@keep];
    $nch  = @keep;

    # scale to uV
    my %to_uv = ('µV' => 1, 'uV' => 1, "\xB5V" => 1, 'microV' => 1, 'mV' => 1e3, 'V' => 1e6, 'nV' => 1e-3);
    my $warned = 0;
    my @scale = map {
        my $u = $unit[$_];
        utf8::decode(my $ud = $u);
        my $f = $to_uv{$u} // $to_uv{$ud};
        if (!defined $f) { carp "$vhdr: unit '$u' of $lab[$_] left as it is" unless $warned++; $f = 1 }
        $res[$_] * $f
    } 0 .. $nch - 1;
    my $data = $raw->double * pdl(\@scale)->dummy(1);

    # markers
    my @ev;
    if (my $mk = $ci->{MarkerFile}) {
        my $mpath = File::Spec->catfile($dir, $mk);
        if (-e $mpath) {
            my $m = _read_ini($mpath)->{'Marker Infos'} // {};
            for my $k (sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] } grep { /^Mk\d+$/ } keys %$m) {
                my ($typ, $desc, $pos) = map { s/\\1/,/gr } split /,/, $m->{$k}, -1;
                next unless defined $pos && $pos =~ /^\d+$/;
                push @ev, { onset => ($pos - 1) / $fs,
                            label => join('/', grep { defined && length } $typ, $desc) };
            }
        } else {
            carp "$vhdr: marker file $mk not found";
        }
    }
    return { data => $data, fs => $fs, labels => \@lab, units => \@unit,
             events => \@ev, format => 'brainvision' };
}

# fill the piddle $p with $want bytes from $fh
sub _read_into {
    my ($fh, $p, $want, $name) = @_;
    my $ref = $p->get_dataref;
    my $got = 0;
    while ($got < $want) {
        my $n = sysread($fh, $$ref, $want - $got, $got);
        croak "$name: read error: $!" unless defined $n;
        last if $n == 0;
        $got += $n;
    }
    $p->upd_data;
    croak "$name: short read ($got of $want bytes)" if $got < $want;
}

# minimal INI reader: { section => { key => value } }; ';' comments; the
# first line ("Brain Vision Data Exchange Header File Version x") is skipped
sub _read_ini {
    my $path = shift;
    open my $fh, '<:raw', $path or croak "$path: $!";
    my (%ini, $sec);
    while (my $l = <$fh>) {
        $l =~ s/\r?\n\z//;
        $l =~ s/^\x{EF}\x{BB}\x{BF}//;                       # UTF-8 BOM (bytes)
        next if $l =~ /^\s*(;|$)/;
        if ($l =~ /^\s*\[(.+?)\]\s*$/) { $sec = $1; next }
        next unless defined $sec;
        my ($k, $v) = split /=/, $l, 2;
        next unless defined $v;
        $ini{$sec}{$k} = $v;
    }
    close $fh;
    return \%ini;
}

1;
