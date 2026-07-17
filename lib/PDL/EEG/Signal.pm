package PDL::EEG::Signal;

use strict;
use warnings;
use Carp qw(croak);
use PDL;
use Exporter 'import';

our @EXPORT_OK = qw(detect_square_pulses);
our $VERSION   = '0.02';

=head1 NAME

PDL::EEG::Signal - Device-independent signal heuristics for EEG channels

=head1 SYNOPSIS

  use PDL::EEG::Signal qw(detect_square_pulses);

  # $data : PDL [n_ch, n_samples], any units
  my $cands = detect_square_pulses($data, fs => 1000, skip_sec => 5, n => 4);
  for my $c (@$cands) {
      printf "ch %d  range=%.0f  median=%.1f  score=%.1f\n",
          $c->{pos}, $c->{range}, $c->{med}, $c->{score};
  }

=head1 DESCRIPTION

C<detect_square_pulses> finds channels that behave like TTL / square-pulse
trigger lines: they sit at a baseline most of the time and jump to a level far
outside the range of the ordinary (EEG) channels. It is B<device independent> —
it takes only a PDL of samples and knows nothing about Nihon Kohden, montages,
or file formats. Vendor-specific naming/wiring is handled elsewhere (e.g.
L<PDL::EEG::IO::NihonKohden::Montage>).

The heuristic, per channel over an analysis window (start skipped to avoid
calibration signals):

  range   = max - min
  med     = median (the resting baseline)

A channel is flagged pulse-like when its C<range> is much larger than the
typical channel's range (C<< range > rel * median_range >>) AND its baseline
sits away from the rails (C<< |med| < baseline_frac * range >>), which rejects
channels stuck near saturation (constant markers) as well as ordinary EEG
(whose range is small). Results are ranked by range; if C<n> is given the top
C<n> candidates are returned.

=head2 detect_square_pulses($data, %opt)

  fs            => $hz     # sampling rate (enables skip_sec)
  skip_sec      => 5       # seconds to skip at the start (calibration guard)
  skip          => $n      # or skip this many samples directly
  rel           => 4       # range must exceed rel * median-of-ranges
  baseline_frac => 0.4     # max |median|/range to still count as pulse-like
  n             => undef   # if set, return only the top n candidates

Returns an arrayref of hashrefs (sorted by range, descending):

  { pos => $ch, range => , med => , maxabs => , score => range/median_range }

=cut

sub detect_square_pulses {
    my ($data, %opt) = @_;
    croak "detect_square_pulses: need a PDL" unless eval { $data->isa('PDL') };
    my ($n_ch, $n_samp) = $data->dims;
    croak "detect_square_pulses: expected 2-D [n_ch,n_samp]" unless defined $n_samp;

    my $rel   = defined $opt{rel}           ? $opt{rel}           : 4;
    my $bfrac = defined $opt{baseline_frac} ? $opt{baseline_frac} : 0.4;

    my $skip = defined $opt{skip} ? $opt{skip}
             : ($opt{skip_sec} && $opt{fs}) ? int($opt{skip_sec} * $opt{fs})
             : 0;
    $skip = 0 if $skip < 0 || $skip >= $n_samp;
    my $hi = $n_samp - 1;

    my $w  = $data->slice(":,$skip:$hi");   # [n_ch, m]
    my $wt = $w->transpose;                 # [m, n_ch] -> reductions collapse dim0
    my $mx = $wt->maximum;                  # [n_ch]
    my $mn = $wt->minimum;                  # [n_ch]
    my $md = $wt->medover;                  # [n_ch] per-channel median (baseline)

    my @info;
    for my $c (0 .. $n_ch - 1) {
        my $max = $mx->at($c);
        my $min = $mn->at($c);
        my $med = $md->at($c);
        my $range  = $max - $min;
        my $maxabs = (abs($max) > abs($min)) ? abs($max) : abs($min);
        push @info, { pos => $c, range => $range, med => $med, maxabs => $maxabs };
    }

    # robust scale = median of per-channel ranges (dominated by ordinary channels)
    my @sorted = sort { $a <=> $b } map { $_->{range} } @info;
    my $med_range = @sorted ? $sorted[int($#sorted / 2)] : 0;
    $med_range = 1e-9 if $med_range <= 0;

    my @cand;
    for my $h (@info) {
        $h->{score} = $h->{range} / $med_range;
        next unless $h->{range} > $rel * $med_range;         # much larger swing
        next unless abs($h->{med}) < $bfrac * $h->{range};   # baseline off the rails
        push @cand, $h;
    }
    @cand = sort { $b->{range} <=> $a->{range} } @cand;

    if (defined $opt{n} && $opt{n} >= 0 && @cand > $opt{n}) {
        @cand = @cand[0 .. $opt{n} - 1];
    }
    return \@cand;
}

1;
