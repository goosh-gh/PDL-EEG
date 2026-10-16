package PDL::EEG::SynthNoise;

# Noise for synthetic EEG test data: white or AR(1), Gaussian or uniform,
# rescaled to an exact SD, from a seed that gives the same numbers on any
# machine (see seed_rng and gaussian below).

use strict;
use warnings;
use Carp qw(croak);
use Exporter 'import';

use PDL;

our @EXPORT_OK = qw(seed_rng gaussian uniform synth_noise);

# ---------------------------------------------------------------------
# seed_rng($seed): seed PDL's random/grandom and Perl's rand together.
#   PDL renamed its seed function from srand to srandom in 2.090; the
#   name that exists is used. CORE::srand is seeded too, so scripts that
#   mix PDL draws and CORE::rand need only this one call.
# ---------------------------------------------------------------------
sub seed_rng {
    my ($seed) = @_;
    croak 'seed_rng: a seed is required' unless defined $seed;
    croak "seed_rng: seed must be a non-negative integer (got '$seed')"
        unless $seed =~ /\A\d+\z/;
    if    (defined &PDL::srandom) { PDL::srandom($seed) }
    elsif (defined &PDL::srand)   { PDL::srand($seed) }
    else  { croak 'seed_rng: this PDL has neither PDL::srandom nor PDL::srand' }
    CORE::srand($seed);
    return;
}

# ---------------------------------------------------------------------
# PDL keeps one random state per CPU. When an operation is split across
# threads (piddles above the autothreading size, 1M elements by default),
# each part draws from a different state, so the numbers depend on the
# number of threads. Drawing with autothreading off uses state 0 only,
# which depends on the seed alone.
# ---------------------------------------------------------------------
sub _single_thread {
    my ($code) = @_;
    my $old = get_autopthread_targ();
    set_autopthread_targ(0);
    my $x = eval { $code->() };
    my $err = $@;
    set_autopthread_targ($old);
    die $err if $err;
    return $x;
}

sub gaussian { my @d = @_; return _single_thread(sub { grandom(@d) }) }
sub uniform  { my @d = @_; return _single_thread(sub { random(@d) }) }

# ---------------------------------------------------------------------
# AR(1) along dim 0 (broadcast over the other dims):
#   x[0] = e[0] / sqrt(1 - a^2)        (stationary from the first sample)
#   x[i] = a * x[i-1] + e[i]
# Within a block starting at s, with c = x[s-1]:
#   x[s+j] = a^j * ( a*c + sum_{k<=j} e[s+k] * a^(-k) )
# a^(-k) grows, so blocks are kept short enough that |a|^(-L) <= 1e150.
# ---------------------------------------------------------------------
sub _ar1 {
    my ($e, $phi) = @_;
    my $nt = $e->dim(0);
    my $x  = $e->copy;
    return $x if $phi == 0;
    my $L = int(150 * log(10) / -log(abs($phi)));
    $L = $nt if $L > $nt;
    $L = 1   if $L < 1;
    my $pw  = $phi ** (-sequence(double, $L));     # a^(-k), k = 0..L-1
    my $pwi = $phi **   sequence(double, $L);      # a^j
    my $first = $x->slice('(0)');
    $first /= sqrt(1 - $phi * $phi);
    my $carry;                                   # x at the end of the previous block
    for (my $s = 0; $s < $nt; $s += $L) {
        my $e1 = ($s + $L > $nt) ? $nt - 1 : $s + $L - 1;
        my $n  = $e1 - $s + 1;
        my $blk = $x->slice("$s:$e1");
        my $w   = $pw->slice('0:' . ($n - 1));
        my $z   = ($blk * $w)->cumusumover;
        $z += $phi * $carry->dummy(0) if defined $carry;
        $blk .= $z * $pwi->slice('0:' . ($n - 1));
        $carry = $blk->slice('(-1)')->copy;
    }
    return $x;
}

my %OPT = map { $_ => 1 } qw(nch kind ar dist sd seed);

# ---------------------------------------------------------------------
# synth_noise($nt, %opt) -> (nt) or (nt, nch)
#   nch  => number of channels (default: none -> 1-D result)
#   kind => 'white' (default) | 'ar1'
#   ar   => AR(1) coefficient, |ar| < 1 (default 0.95)
#   dist => 'gauss' (default) | 'uniform'   distribution of the innovations
#   sd   => target SD per channel (default 1)
#   seed => seed_rng(seed) before drawing (default: continue the current state)
# Each channel is made zero-mean and rescaled to SD sd exactly (divisor nt).
# ---------------------------------------------------------------------
sub synth_noise {
    my ($nt, %o) = @_;
    croak 'synth_noise: nt must be an integer >= 2'
        unless defined $nt && $nt =~ /\A\d+\z/ && $nt >= 2;
    my @bad = grep { !$OPT{$_} } sort keys %o;
    croak "synth_noise: unknown option(s): @bad" if @bad;

    my $kind = $o{kind} // 'white';
    my $dist = $o{dist} // 'gauss';
    my $sd   = $o{sd}   // 1;
    my $phi    = $o{ar}   // 0.95;
    croak "synth_noise: kind must be 'white' or 'ar1' (got '$kind')"
        unless $kind eq 'white' || $kind eq 'ar1';
    croak "synth_noise: dist must be 'gauss' or 'uniform' (got '$dist')"
        unless $dist eq 'gauss' || $dist eq 'uniform';
    croak 'synth_noise: sd must be > 0' unless $sd > 0;
    croak 'synth_noise: ar must satisfy |ar| < 1' if $kind eq 'ar1' && abs($phi) >= 1;
    my $nch = $o{nch};
    croak 'synth_noise: nch must be an integer >= 1'
        if defined $nch && !($nch =~ /\A\d+\z/ && $nch >= 1);

    seed_rng($o{seed}) if defined $o{seed};
    my @dims = defined $nch ? ($nt, $nch) : ($nt);
    my $e = $dist eq 'gauss' ? gaussian(double, @dims)
                             : (uniform(double, @dims) - 0.5) * sqrt(12);
    my $x = $kind eq 'ar1' ? _ar1($e, $phi) : $e;

    $x -= $x->average->dummy(0);
    my $s = sqrt(($x * $x)->average);
    croak 'synth_noise: a channel came out constant' if any($s == 0);
    $x *= $sd / $s->dummy(0);
    return $x;
}

1;

__END__

=encoding utf8

=head1 NAME

PDL::EEG::SynthNoise - reproducible noise for synthetic EEG test data

=head1 SYNOPSIS

  use PDL;
  use PDL::EEG::SynthNoise qw(seed_rng gaussian uniform synth_noise);

  # 60 s at 1 kHz, 19 channels of AR(1) noise, SD 3 per channel
  my $bg = synth_noise(60_000, nch => 19, kind => 'ar1', ar => 0.95,
                       sd => 3, seed => 11);

  # bounded white noise: never crosses a 2.5 x MAD threshold by itself
  my $v = synth_noise(20_000, dist => 'uniform', sd => 1.7, seed => 1);

  seed_rng(7);                  # PDL random/grandom and CORE::rand
  my $m = gaussian(19, 6);      # N(0,1), same numbers on any machine

=head1 DESCRIPTION

=head2 seed_rng($seed)

Seeds PDL's random number generator and Perl's C<rand> with the same
non-negative integer. PDL 2.090 renamed its seed function from C<srand> to
C<srandom>; C<seed_rng> calls whichever exists. A bare C<srand> in a script
seeds PDL in older PDL and only Perl's C<rand> in newer PDL, so use
C<seed_rng> instead.

=head2 gaussian(@dims), uniform(@dims)

C<grandom(@dims)> and C<random(@dims)> drawn with PDL's automatic threading
turned off (and restored afterwards). PDL keeps one random state per CPU,
and an array above the autothreading size (1M elements by default) is
filled by several threads from different states, so C<grandom> gives
different numbers for the same seed on machines with different numbers of
cores. With threading off, the numbers depend only on the seed, and a
longer draw starts with the same numbers as a shorter one.

=head2 synth_noise($nt, %opt)

Returns C<$nt> samples of noise, or an (nt, nch) array when C<nch> is
given. Options:

=over

=item nch

Number of channels. Channels are independent.

=item kind

C<'white'> (default) or C<'ar1'>: C<x[i] = ar * x[i-1] + e[i]>, started
from the stationary distribution so there is no transient at the start.

=item ar

AR(1) coefficient, |ar| E<lt> 1. Default 0.95.

=item dist

Distribution of the innovations C<e>: C<'gauss'> (default) or
C<'uniform'>.

=item sd

Target SD (default 1). Each channel is made zero-mean and rescaled so that
its SD (divisor nt) is exactly C<sd>.

=item seed

If given, C<seed_rng(seed)> is called first. Otherwise the draw continues
from the current state. Because C<seed_rng> also seeds Perl's C<rand>, a
C<seed> here restarts the numbers that C<CORE::rand> returns afterwards.

=back

=head2 Noise and robust thresholds

Detectors such as C<PDL::EEG::GED::detect_blinks> mark samples above
C<k * 1.4826 * MAD>. Gaussian noise crosses that threshold for k = 2.5 in
about 1.2% of the samples whatever its SD, and AR(1) noise with Gaussian
innovations is Gaussian too: it crosses in fewer, longer runs. White
uniform noise on [-c, c] has MAD c/2, so its largest value is about
1.35 x (1.4826 x MAD) and it does not cross the threshold for any k above
about 1.35. Use C<< kind => 'white', dist => 'uniform' >> when a test needs
noise that stays below such a threshold by itself. AR(1) noise with uniform
innovations is close to Gaussian and does not have this property.

=head1 SEE ALSO

L<PDL::EEG::GED>

=cut
