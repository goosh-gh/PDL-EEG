package PDL::EEG::IO::NihonKohden::PTN;

use strict;
use warnings;
use Carp qw(croak);
use Exporter 'import';

our @EXPORT_OK = qw(parse_ptn list_montages find_montage_file);
our $VERSION   = '0.02';

=head1 NAME

PDL::EEG::IO::NihonKohden::PTN - Read Nihon Kohden Neurofax pattern (montage) files

=head1 SYNOPSIS

  use PDL::EEG::IO::NihonKohden::PTN qw(parse_ptn find_montage_file list_montages);

  my $m = parse_ptn('Pattern_032.PTN');
  print "$m->{name}\n";                         # "21A"
  for my $ch (@{ $m->{channels} }) {
      printf "%-10s g1=%d g2=%d sens=0x%02x%s\n",
          ($ch->{inline} // '-'), $ch->{g1}, $ch->{g2}, $ch->{sens},
          $ch->{trigger} ? " (trigger)" : "";
  }

  my $path = find_montage_file('/path/to/subject.PTN', 'IIA');

=head1 DESCRIPTION

Both EEG-1100C and EEG-1200A store display montages as "EEG-1000/9000 Pattern
Info File" (.PTN). This module decodes the display-channel list independently of
recording format.

Each display channel is an 80-byte record beginning at offset 0x410:

  byte0   G1 electrode index (0-based; 0 for special DC/trigger)
  byte1   G2 reference electrode index (0 for special)
  byte2   SENS code (0x0f EEG, 0x05 DC/trigger, others special)
  b12-13  display x-position (LE); 0xFFFF marks the end of the used list
  b14..   inline display name (custom montages only; empty on stock ones)

Custom montages (e.g. "21A") give every channel an inline name including the
trigger channels ("TrigBit0/2/4/8"). Stock montages (e.g. "IA") leave EEG names
empty, so the caller supplies electrode names from the .21e / defaults.

B<Important:> trigger channels carry G1=0 (no electrode index), so a .PTN gives
the trigger B<count and names> but not which recorded channel each one is. Bind
those to real channels by signal (see L<PDL::EEG::IO::NihonKohden::Montage>).

=cut

use constant {
    REC_START => 1040,      # 0x410
    REC_SIZE  => 80,
    NAME_OFF  => 0x80,
    END_POS   => 0xFFFF,
};

sub parse_ptn {
    my ($file) = @_;
    open my $fh, '<:raw', $file or croak "parse_ptn: $file: $!";
    local $/; my $buf = <$fh>; close $fh;

    my $sig = unpack 'Z*', substr($buf, 0, 31);
    croak "parse_ptn: not a Neurofax pattern file: $file"
        unless $sig =~ /Pattern Info File/;

    my %m = (
        file     => $file,
        header   => $sig,
        name     => (unpack 'Z*', substr($buf, NAME_OFF, 16)),
        channels => [],
        triggers => [],
    );

    my $off = REC_START;
    my $n   = 0;
    while ($off + REC_SIZE <= length $buf) {
        my $r = substr($buf, $off, REC_SIZE); $off += REC_SIZE;
        my ($g1, $g2, $sens) = unpack 'CCC', $r;
        my $pos = unpack 'v', substr($r, 12, 2);
        last if $pos == END_POS;                              # end-of-list sentinel
        next if $sens == 0 && substr($r, 3, 7) eq "\0" x 7;   # stray blank slot

        my $inline = unpack 'Z*', substr($r, 14, 32);
        $inline =~ s/[^\x20-\x7e]//g;
        my $special = ($g1 == 0 && $g2 == 0);
        my $trigger = ($special && $sens == 0x05);

        push @{ $m{channels} }, {
            index   => $n,
            g1      => $g1,
            g2      => $g2,
            sens    => $sens,
            xpos    => $pos,
            inline  => (length $inline ? $inline : undef),
            special => $special ? 1 : 0,
            trigger => $trigger ? 1 : 0,
            ch_idx  => $special ? undef : $g1 + 1,   # EEG: recorded ch_idx = G1+1
        };
        push @{ $m{triggers} }, $n if $trigger;
        $n++;
    }
    $m{n} = $n;
    return \%m;
}

# List all Pattern_*.PTN in a directory -> { NAME => path }
sub list_montages {
    my ($dir) = @_;
    opendir my $dh, $dir or croak "list_montages: $dir: $!";
    my @files = sort grep { /\.PTN$/i } readdir $dh;
    closedir $dh;
    my %by_name;
    for my $f (@files) {
        my $path = "$dir/$f";
        my $m = eval { parse_ptn($path) } or next;
        $by_name{ $m->{name} } = $path if defined $m->{name} && length $m->{name};
    }
    return \%by_name;
}

# Find the .PTN whose NAME matches $name (case-insensitive, whitespace-insensitive)
sub find_montage_file {
    my ($dir, $name) = @_;
    return undef unless defined $name;
    (my $want = uc $name) =~ s/\s+//g;
    my $by = list_montages($dir);
    for my $nm (keys %$by) {
        (my $have = uc $nm) =~ s/\s+//g;
        return $by->{$nm} if $have eq $want;
    }
    return undef;
}

1;
