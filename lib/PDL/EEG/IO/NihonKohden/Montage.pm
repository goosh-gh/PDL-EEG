package PDL::EEG::IO::NihonKohden::Montage;

use strict;
use warnings;
use Carp qw(croak carp);
use Exporter 'import';
use PDL::EEG::IO::NihonKohden::PTN qw(parse_ptn find_montage_file);
use PDL::EEG::Signal qw(detect_square_pulses);

our @EXPORT_OK = qw(montage_from_log resolve_labels);
our $VERSION   = '0.01';

=head1 NAME

PDL::EEG::IO::NihonKohden::Montage - Resolve channel labels from .LOG montage + .PTN + signal

=head1 SYNOPSIS

  use PDL::EEG::IO::NihonKohden qw(read_nk);
  use PDL::EEG::IO::NihonKohden::Montage qw(resolve_labels);

  my $rec = read_nk('YJ0394VB.EEG', all_blocks => 1);

  my $r = resolve_labels($rec,
      ptn_dir => 'YJ0394VB.PTN',   # dir of Pattern_0NN.PTN (optional)
      apply   => 1,                # rewrite $rec->{labels} in place
  );
  # $r->{montage}   e.g. "IIA"  (recording montage, from .LOG)
  # $r->{triggers}  [ {ch_idx=>45,name=>'TrigBit0',range=>...}, ... ]
  # $r->{label_map} { 45=>'TrigBit0', 46=>'TrigBit2', ... }  (feed back to read_nk)

  # then, for a byte-clean result, re-read with the map:
  my $rec2 = read_nk('YJ0394VB.EEG', all_blocks=>1, label_map => $r->{label_map});

=head1 DESCRIPTION

Nihon Kohden trigger/DC channel names are B<not> reliably derivable from the
recording format: the same trigger appears as DC03-06 on one headbox and DC01-04
on another, and the .21e is a generic template. The authoritative naming lives in
the display/montage layer (.PTN) and the physical wiring, while which recorded
channel actually carries a trigger is only visible in the B<signal>.

This module combines all three sources, none of which alone is sufficient:

=over

=item * B<.LOG> gives the recording montage name (C<REC START IIA EEG> -> "IIA").

=item * B<.PTN> for that montage gives the trigger B<count and display names>
(e.g. TrigBit0/2/4/8) but not their recorded channel index.

=item * B<signal> (L<PDL::EEG::Signal>) identifies which recorded channels
actually carry TTL/square pulses.

=back

Detected trigger channels (sorted by ascending ch_idx) are zipped, in order,
onto the montage's trigger names, producing a C<label_map> keyed by 1-based
ch_idx suitable for C<< read_nk(..., label_map => ...) >>. A manual C<label_map>
always wins, so exceptions can be pinned by hand.

=head2 montage_from_log(\@events)

Return the recording montage name from read_nk's C<events> (first
C<REC START E<lt>NAMEE<gt> EEG|CAL> marker), or undef.

=head2 resolve_labels($rec, %opt)

  ptn_dir  => $dir     # directory containing Pattern_0NN.PTN (to look up names)
  ptn      => $file    # or a specific .PTN file (overrides montage lookup)
  montage  => $name    # override the montage name (else taken from .LOG events)
  n        => $n       # expected trigger count (if no .PTN available)
  names    => \@names  # override trigger names (else from .PTN, else Trig1..N)
  rel      => 4        # detector sensitivity (see PDL::EEG::Signal)
  skip_sec => 5        # skip leading calibration seconds
  apply    => 0        # if true, rewrite $rec->{labels} in place

Returns a hashref: { montage, ptn, n_expected, triggers=>[...], label_map, notes }.

=cut

sub montage_from_log {
    my ($events) = @_;
    return undef unless $events && @$events;
    for my $e (@$events) {
        my $lab = ref $e eq 'HASH' ? ($e->{label} // '') : "$e";
        return $1 if $lab =~ /REC\s+START\s+(\S+)\s+(?:EEG|CAL)\b/i;
    }
    return undef;
}

sub resolve_labels {
    my ($rec, %opt) = @_;
    croak "resolve_labels: need read_nk record hashref" unless ref $rec eq 'HASH';
    croak "resolve_labels: \$rec->{data} is not a PDL"
        unless eval { $rec->{data}->isa('PDL') };

    my @notes;
    my $montage = defined $opt{montage} ? $opt{montage}
                : montage_from_log($rec->{events});
    push @notes, "montage name not found in .LOG events" unless defined $montage;

    # locate + parse the .PTN for this montage
    my ($ptn_path, $ptn);
    if ($opt{ptn}) {
        $ptn_path = $opt{ptn};
    } elsif ($opt{ptn_dir} && defined $montage) {
        $ptn_path = find_montage_file($opt{ptn_dir}, $montage);
        push @notes, "no .PTN named '$montage' in $opt{ptn_dir}" unless $ptn_path;
    }
    if ($ptn_path) {
        $ptn = eval { parse_ptn($ptn_path) };
        push @notes, "failed to parse $ptn_path: $@" if $@;
    }

    # trigger names (montage slot order) and expected count
    my @names;
    if ($opt{names}) {
        @names = @{ $opt{names} };
    } elsif ($ptn && @{ $ptn->{triggers} }) {
        @names = map { $ptn->{channels}[$_]{inline} // "Trig" } @{ $ptn->{triggers} };
    }
    my $n_expected = @names ? scalar(@names)
                   : defined $opt{n} ? $opt{n}
                   : undef;

    # signal-based detection (device independent)
    my $cands = detect_square_pulses(
        $rec->{data},
        fs       => $rec->{fs},
        skip_sec => (defined $opt{skip_sec} ? $opt{skip_sec} : 5),
        rel      => (defined $opt{rel} ? $opt{rel} : 4),
        (defined $n_expected ? (n => $n_expected) : ()),
    );

    # map detector positions -> 1-based ch_idx, sort ascending
    my $chidx = $rec->{ch_indices};
    my @trig;
    for my $c (@$cands) {
        my $pos    = $c->{pos};
        my $ch_idx = ($chidx && defined $chidx->[$pos]) ? $chidx->[$pos] : ($pos + 1);
        push @trig, { pos => $pos, ch_idx => $ch_idx,
                      range => $c->{range}, med => $c->{med}, score => $c->{score} };
    }
    @trig = sort { $a->{ch_idx} <=> $b->{ch_idx} } @trig;

    # zip names (montage order) onto detected triggers (ch_idx order)
    if (@names) {
        if (@names != @trig) {
            push @notes, sprintf(
                "trigger count mismatch: montage lists %d (%s) but signal found %d",
                scalar @names, join(',', @names), scalar @trig);
        }
        for my $i (0 .. $#trig) {
            $trig[$i]{name} = $i < @names ? $names[$i] : sprintf("Trig%d", $i + 1);
        }
    } else {
        $trig[$_]{name} = sprintf("Trig%d", $_ + 1) for 0 .. $#trig;
        push @notes, "no montage trigger names; used generic Trig1..N";
    }

    my %label_map = map { $_->{ch_idx} => $_->{name} } @trig;

    if ($opt{apply}) {
        my $labels = $rec->{labels};
        if ($chidx && $labels) {
            for my $i (0 .. $#$labels) {
                my $ci = $chidx->[$i];
                $labels->[$i] = $label_map{$ci} if defined $ci && exists $label_map{$ci};
            }
        } else {
            push @notes, "apply: cannot rewrite labels (missing ch_indices/labels)";
        }
    }

    return {
        montage    => $montage,
        ptn        => $ptn_path,
        n_expected => $n_expected,
        triggers   => \@trig,
        label_map  => \%label_map,
        notes      => \@notes,
    };
}

1;
