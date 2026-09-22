package PDL::EEG::TFA;

use strict;
use warnings;
use PDL;
use PDL::FFT;
use Carp qw(croak);

use Exporter 'import';
our @EXPORT_OK = qw(
    morlet_wavelet
    cwt_morlet
    tfr_morlet
    tfr_superlet
    tfr_stat
    apply_baseline
);

# ------------------------------------------------------------------
# _next_pow2($n) -> smallest power of two >= $n
# ------------------------------------------------------------------
sub _next_pow2 {
    my ($n) = @_;
    my $p = 1;
    $p <<= 1 while $p < $n;
    return $p;
}

# ------------------------------------------------------------------
# morlet_wavelet($freq, $sfreq, %opt)
#
# Complex Morlet wavelet in the "number of cycles" parameterisation
# (the convention used by MNE-Python / Tallon-Baudry):
#
#     sigma_t = n_cycles / (2*pi*f)                 (Gaussian SD, seconds)
#     w(t)    = exp(2i*pi*f*t) * exp(-t^2/(2 sigma_t^2))
#
# The kernel is truncated at +/- sigma_win * sigma_t (default 5) and
# normalised to unit L2 energy, so the absolute magnitude scale is a
# fixed convention.  Baseline z-scoring (the usual final step) is
# invariant to this constant, so it does not affect the maps.
#
# Options:
#   n_cycles   => 7.0     wavelet width in cycles
#   sigma_win  => 5.0     truncation half-width in units of sigma_t
#
# Returns ($wr, $wi): real and imaginary parts, 1-D piddles of odd
# length centred on t = 0.
# ------------------------------------------------------------------
sub morlet_wavelet {
    my ($freq, $sfreq, %opt) = @_;
    croak "morlet_wavelet: freq must be > 0"  unless $freq  > 0;
    croak "morlet_wavelet: sfreq must be > 0" unless $sfreq > 0;

    my $n_cycles  = $opt{n_cycles}  // 7.0;
    my $sigma_win = $opt{sigma_win} // 5.0;

    my $sigma_t = $n_cycles / (2 * 3.14159265358979 * $freq);
    my $half    = int($sigma_win * $sigma_t * $sfreq);       # samples each side
    $half = 1 if $half < 1;
    my $n       = 2 * $half + 1;                              # odd length

    my $t   = (sequence($n) - $half) / $sfreq;               # centred time axis
    my $env = exp(-($t**2) / (2 * $sigma_t**2));
    my $ang = 2 * 3.14159265358979 * $freq * $t;
    my $wr  = $env * cos($ang);
    my $wi  = $env * sin($ang);

    # unit L2 energy
    my $norm = sqrt( ($wr**2 + $wi**2)->sum );
    $norm = 1 if $norm == 0;
    return ($wr / $norm, $wi / $norm);
}

# ------------------------------------------------------------------
# _fftconvolve_same($sr, $si, $wr, $wi) -> ($cr, $ci)
#
# Linear convolution of a (possibly complex) signal with a complex
# kernel via zero-padded FFT, returning the central 'same'-length
# window (length == signal length).
# ------------------------------------------------------------------
sub _fftconvolve_same {
    my ($sr, $si, $wr, $wi) = @_;
    my $ns   = $sr->dim(0);
    my $nw   = $wr->dim(0);
    my $L    = $ns + $nw - 1;
    my $nfft = _next_pow2($L);

    my $Sr = zeroes(double, $nfft); $Sr->slice("0:@{[$ns-1]}") .= $sr;
    my $Si = zeroes(double, $nfft); $Si->slice("0:@{[$ns-1]}") .= $si;
    my $Wr = zeroes(double, $nfft); $Wr->slice("0:@{[$nw-1]}") .= $wr;
    my $Wi = zeroes(double, $nfft); $Wi->slice("0:@{[$nw-1]}") .= $wi;

    fft($Sr, $Si);
    fft($Wr, $Wi);

    my $Pr = $Sr * $Wr - $Si * $Wi;
    my $Pi = $Sr * $Wi + $Si * $Wr;
    ifft($Pr, $Pi);

    my $start = int(($nw - 1) / 2);
    my $cr = $Pr->slice("$start:@{[$start + $ns - 1]}")->sever;
    my $ci = $Pi->slice("$start:@{[$start + $ns - 1]}")->sever;
    return ($cr, $ci);
}

# ------------------------------------------------------------------
# cwt_morlet($signal, $sfreq, $freqs, %opt) -> ($cr, $ci)
#
# Continuous Morlet wavelet transform of one real 1-D signal.
# $freqs is a 1-D piddle of centre frequencies (Hz).
#
# Returns two (nfreq, ntime) piddles: real and imaginary parts of the
# complex coefficients.  n_cycles may be a scalar or a per-frequency
# piddle.
# ------------------------------------------------------------------
sub cwt_morlet {
    my ($signal, $sfreq, $freqs, %opt) = @_;
    $signal = $signal->flat;
    my $nt = $signal->dim(0);
    my @f  = list($freqs);
    my $nf = scalar @f;

    my $ncyc = $opt{n_cycles} // 7.0;
    my @ncyc = ref($ncyc) && eval { $ncyc->isa('PDL') } ? list($ncyc)
             : (($ncyc) x $nf);
    croak "n_cycles length mismatch" unless @ncyc == $nf;

    my $si0 = zeroes(double, $nt);
    my $cr  = zeroes(double, $nf, $nt);
    my $ci  = zeroes(double, $nf, $nt);

    for my $k (0 .. $nf - 1) {
        my ($wr, $wi) = morlet_wavelet($f[$k], $sfreq,
                            n_cycles  => $ncyc[$k],
                            sigma_win => ($opt{sigma_win} // 5.0));
        my ($rr, $ii) = _fftconvolve_same($signal, $si0, $wr, $wi);
        $cr->slice("($k),:") .= $rr;
        $ci->slice("($k),:") .= $ii;
    }
    return ($cr, $ci);
}

# ------------------------------------------------------------------
# tfr_morlet($data, $sfreq, $freqs, %opt) -> $power   (or %result)
#
# Trial-averaged Morlet power.  $data is either a 1-D piddle (ntime)
# for a single trial, or a 2-D piddle (nepoch, ntime).
#
# Options:
#   n_cycles  => 7.0        scalar or (nfreq) piddle
#   sigma_win => 5.0
#   output    => 'power'    'power' | 'itc' | 'evoked' | 'induced' | 'complex' | 'all'
#
# Returns, by output:
#   'power'   : (nfreq, ntime) mean |coef|^2 across trials (total power)
#   'itc'     : (nfreq, ntime) inter-trial coherence in [0,1]
#   'evoked'  : (nfreq, ntime) phase-locked power = |mean coef|^2
#   'induced' : (nfreq, ntime) non-phase-locked power = total - evoked
#   'complex' : ($cr, $ci) trial-averaged complex (evoked) coefficients
#   'all'     : hash { power, itc, evoked, induced, evoked_r, evoked_i }
# ------------------------------------------------------------------
sub tfr_morlet {
    my ($data, $sfreq, $freqs, %opt) = @_;
    my $out = $opt{output} // 'power';

    $data = $data->dummy(0) if $data->ndims == 1;   # -> (1, ntime)
    my $ne = $data->dim(0);
    my $nt = $data->dim(1);
    my $nf = $freqs->nelem;

    my $power    = zeroes(double, $nf, $nt);        # sum |c|^2
    my $itc_r    = zeroes(double, $nf, $nt);        # sum c/|c| (real)
    my $itc_i    = zeroes(double, $nf, $nt);
    my $evoked_r = zeroes(double, $nf, $nt);        # sum c (real)
    my $evoked_i = zeroes(double, $nf, $nt);

    for my $e (0 .. $ne - 1) {
        my $sig = $data->slice("($e),:");
        my ($cr, $ci) = cwt_morlet($sig, $sfreq, $freqs, %opt);
        my $mag2 = $cr**2 + $ci**2;
        $power    += $mag2;
        $evoked_r += $cr;
        $evoked_i += $ci;
        if ($out eq 'itc' or $out eq 'all') {
            my $mag = sqrt($mag2);
            $mag->where($mag == 0) .= 1;
            $itc_r += $cr / $mag;
            $itc_i += $ci / $mag;
        }
    }

    $power    /= $ne;
    $evoked_r /= $ne;
    $evoked_i /= $ne;

    return $power if $out eq 'power';

    if ($out eq 'itc') {
        return sqrt($itc_r**2 + $itc_i**2) / $ne;
    }
    if ($out eq 'complex') {
        return ($evoked_r, $evoked_i);
    }
    if ($out eq 'evoked') {                     # phase-locked power = |mean c|^2
        return $evoked_r**2 + $evoked_i**2;
    }
    if ($out eq 'induced') {                    # non-phase-locked = total - evoked
        my $ind = $power - ($evoked_r**2 + $evoked_i**2);
        $ind->where($ind < 0) .= 0;
        return $ind;
    }
    # 'all'
    my $evoked_power = $evoked_r**2 + $evoked_i**2;
    my $induced      = $power - $evoked_power;
    $induced->where($induced < 0) .= 0;
    return (
        power    => $power,
        itc      => sqrt($itc_r**2 + $itc_i**2) / $ne,
        evoked   => $evoked_power,
        induced  => $induced,
        evoked_r => $evoked_r,
        evoked_i => $evoked_i,
    );
}

# ------------------------------------------------------------------
# tfr_superlet($data, $sfreq, $freqs, %opt) -> $power
#
# Trial-averaged (multiplicative, adaptive) superlet power, after
# Moca et al. 2021.  A superlet at frequency f is a set of Morlet
# wavelets sharing the centre frequency but with increasing cycle
# counts base_cycles*o (o = 1..order); the superlet response is the
# geometric mean of the individual wavelet powers.  Higher orders add
# frequency resolution; the geometric mean keeps the sharp time
# resolution of the low-order wavelet.  "Adaptive" ramps the order
# linearly with frequency, giving fine time resolution for fast HFO
# bursts and fine frequency resolution where oscillations are slow.
#
# Options:
#   base_cycles => 3        cycle count of the order-1 wavelet (c1)
#   order_min   => 1        order at the lowest frequency
#   order_max   => 5        order at the highest frequency (adaptive)
#   order       => undef    fixed order for all freqs (overrides min/max)
#   sigma_win   => 5.0
#
# Returns (nfreq, ntime) power.  Slower than tfr_morlet by ~mean(order).
# ------------------------------------------------------------------
sub tfr_superlet {
    my ($data, $sfreq, $freqs, %opt) = @_;
    $data = $data->dummy(0) if $data->ndims == 1;
    my $ne = $data->dim(0);
    my $nt = $data->dim(1);
    my @f  = list($freqs);
    my $nf = scalar @f;

    my $base = $opt{base_cycles} // 3;
    my $swin = $opt{sigma_win}   // 5.0;

    # per-frequency integer order
    my @order;
    if (defined $opt{order}) {
        @order = (int($opt{order})) x $nf;
    } else {
        my $omin = $opt{order_min} // 1;
        my $omax = $opt{order_max} // 5;
        my $span = ($f[-1] > $f[0]) ? ($f[-1] - $f[0]) : 1;
        for my $k (0 .. $nf - 1) {
            my $o = $omin + ($omax - $omin) * ($f[$k] - $f[0]) / $span;
            push @order, int($o + 0.5);
        }
    }
    $_ < 1 and $_ = 1 for @order;

    my $si0   = zeroes(double, $nt);
    my $power = zeroes(double, $nf, $nt);

    for my $e (0 .. $ne - 1) {
        my $sig = $data->slice("($e),:");
        for my $k (0 .. $nf - 1) {
            my $logsum = zeroes(double, $nt);       # sum of log-powers
            for my $o (1 .. $order[$k]) {
                my ($wr, $wi) = morlet_wavelet($f[$k], $sfreq,
                        n_cycles => $base * $o, sigma_win => $swin);
                my ($rr, $ii) = _fftconvolve_same($sig, $si0, $wr, $wi);
                my $p = $rr**2 + $ii**2;
                $p->where($p <= 0) .= 1e-300;        # geometric-mean guard
                $logsum += log($p);
            }
            $power->slice("($k),:") += exp($logsum / $order[$k]);   # geo mean
        }
    }
    return $power / $ne;
}

# ------------------------------------------------------------------
# tfr_stat($data, $sfreq, $freqs, $times, $baseline, %opt) -> %result
#
# Across-trial reliability of the evoked time-frequency power, computed
# streaming (no per-trial storage).  Each trial's power is baseline-
# corrected per frequency (mean over the baseline time window
# subtracted), then the mean and its standard error across trials are
# accumulated.  The returned z = mean / SEM is a per-(freq,time)
# one-sample statistic against the baseline: "is the averaged power at
# this point reliably different from baseline, given trial-to-trial
# variance?".
#
# This is distinct from apply_baseline('zscore'), which normalises by
# the TEMPORAL SD within the baseline window of the already-averaged
# map.  tfr_stat normalises by the ACROSS-TRIAL SD, so with many trials
# the Central Limit Theorem makes z approximately Gaussian even though
# single-trial power is skewed.  Suggested thresholds: |z|>=3 (liberal),
# |z|>=5 (~Bonferroni for ~1e5 time-frequency points).
#
# $data is (nepoch, ntime); $times is (ntime) seconds; $baseline is
# [t0,t1].  %opt is passed through to cwt_morlet (n_cycles, sigma_win).
#
# Returns hash: z (nfreq,ntime), mean (baseline-corrected mean power),
# sem, n.
# ------------------------------------------------------------------
sub tfr_stat {
    my ($data, $sfreq, $freqs, $times, $baseline, %opt) = @_;
    $data = $data->dummy(0) if $data->ndims == 1;
    my $ne = $data->dim(0);
    my $nt = $data->dim(1);
    my $nf = $freqs->nelem;

    my ($b0, $b1) = @$baseline;
    $b0 = $times->min->sclr unless defined $b0;
    $b1 = $times->max->sclr unless defined $b1;
    my $bidx = which(($times >= $b0) & ($times <= $b1));
    croak "tfr_stat: empty baseline window" unless $bidx->nelem;

    my $S1 = zeroes(double, $nf, $nt);      # sum of baseline-corrected power
    my $S2 = zeroes(double, $nf, $nt);      # sum of squares

    for my $e (0 .. $ne - 1) {
        my ($cr, $ci) = cwt_morlet($data->slice("($e),:"), $sfreq, $freqs, %opt);
        my $p  = $cr**2 + $ci**2;                               # (nf,nt)
        my $b  = $p->dice_axis(1, $bidx)->xchg(0,1)->avgover;   # (nf) baseline mean
        my $pc = $p - $b->dummy(1);
        $S1 += $pc;
        $S2 += $pc * $pc;
    }

    my $mean = $S1 / $ne;
    my $var  = ($S2 - $S1 * $S1 / $ne) / ($ne - 1);
    $var->where($var <= 0) .= 1e-300;
    my $sem  = sqrt($var / $ne);
    my $z    = $mean / $sem;
    return (z => $z, mean => $mean, sem => $sem, n => $ne);
}

# ------------------------------------------------------------------
# apply_baseline($power, $times, $baseline, %opt) -> $normalised
#
# Per-frequency baseline normalisation of a (nfreq, ntime) power map.
# $times is a 1-D piddle (seconds) aligned with the time axis.
# $baseline is [tmin, tmax] (seconds); either bound may be undef to
# extend to the edge.
#
# Options:
#   mode => 'zscore'   'zscore' | 'logratio' | 'percent' | 'ratio' | 'mean'
# ------------------------------------------------------------------
sub apply_baseline {
    my ($power, $times, $baseline, %opt) = @_;
    my $mode = $opt{mode} // 'zscore';
    my ($b0, $b1) = @$baseline;
    $b0 = $times->min->sclr unless defined $b0;
    $b1 = $times->max->sclr unless defined $b1;

    my $mask = ($times >= $b0) & ($times <= $b1);
    my $idx  = which($mask);
    croak "apply_baseline: empty baseline window" if $idx->nelem == 0;

    my $bl   = $power->dice_axis(1, $idx);        # (nfreq, n_bl)
    my $mean = $bl->xchg(0,1)->avgover;           # (nfreq)
    my $out  = $power->copy;

    if ($mode eq 'zscore') {
        my $sd = (($bl - $mean->dummy(1))**2)->xchg(0,1)->avgover->sqrt; # (nfreq)
        $sd->where($sd == 0) .= 1;
        $out = ($power - $mean->dummy(1)) / $sd->dummy(1);
    }
    elsif ($mode eq 'ratio') {
        $out = $power / $mean->dummy(1);
    }
    elsif ($mode eq 'logratio') {
        $out = log10($power / $mean->dummy(1));
    }
    elsif ($mode eq 'percent') {
        $out = ($power - $mean->dummy(1)) / $mean->dummy(1);
    }
    elsif ($mode eq 'mean') {
        $out = $power - $mean->dummy(1);
    }
    else {
        croak "apply_baseline: unknown mode '$mode'";
    }
    return $out;
}

1;

__END__

=head1 NAME

PDL::EEG::TFA - continuous complex-Morlet time-frequency analysis in PDL

=head1 SYNOPSIS

    use PDL;
    use PDL::EEG::TFA qw(tfr_morlet apply_baseline);

    # $epochs : (nepoch, ntime) real piddle, $sfreq in Hz
    my $freqs = 70 + 5*sequence(167);            # 70..900 Hz
    my $power = tfr_morlet($epochs, $sfreq, $freqs, n_cycles => 7);
    my $z     = apply_baseline($power, $times, [-0.05, -0.004], mode => 'zscore');

=head1 DESCRIPTION

Renderer- and IO-independent time-frequency engine: complex Morlet
continuous wavelet transform, trial-averaged power and inter-trial
coherence, and per-frequency baseline normalisation.  Convolution is
carried out in the frequency domain via C<PDL::FFT>.

=head1 FUNCTIONS

All functions are exported on request only (C<@EXPORT_OK>).

=head2 morlet_wavelet($freq, $sfreq, %opt)

    my ($wr, $wi) = morlet_wavelet(100, 1000, n_cycles => 7);

Complex Morlet wavelet in the "number of cycles" parameterisation:

    sigma_t = n_cycles / (2*pi*freq)            # Gaussian SD, seconds
    w(t)    = exp(2i*pi*freq*t) * exp(-t^2 / (2*sigma_t^2))

C<$freq> is the centre frequency and C<$sfreq> the sampling rate, both in
Hz; both must be E<gt> 0 or the function croaks.

Options:

=over 4

=item C<n_cycles> => number (default 7.0)

Wavelet width in cycles.

=item C<sigma_win> => number (default 5.0)

Truncation half-width in units of C<sigma_t>.  The kernel spans
C<int(sigma_win * sigma_t * sfreq)> samples on each side of t = 0 (at
least 1).

=back

Returns C<($wr, $wi)>: real and imaginary parts as 1-D piddles of odd
length, centred on t = 0, normalised to unit L2 energy
(C<sum($wr**2 + $wi**2) == 1>).  Kernel length shrinks as C<$freq>
increases.

=head2 cwt_morlet($signal, $sfreq, $freqs, %opt)

    my ($cr, $ci) = cwt_morlet($signal, $sfreq, $freqs, n_cycles => 7);

Continuous Morlet wavelet transform of one real signal.  C<$signal> is
flattened to 1-D (C<ntime>); C<$freqs> is a 1-D piddle of centre
frequencies in Hz.

Options:

=over 4

=item C<n_cycles> => number | piddle (nfreq) (default 7.0)

A scalar applies to every frequency; a piddle gives one cycle count per
frequency and must have C<nfreq> elements (croaks otherwise).

=item C<sigma_win> => number (default 5.0)

Passed to C<morlet_wavelet>.

=back

Each frequency is convolved with its wavelet by zero-padded FFT (length
rounded up to a power of two) and the central C<ntime> samples of the
linear convolution are kept ("same" mode).  Because the padding is with
zeros, coefficients within half a kernel length of either edge are
affected by the signal boundary.

Returns C<($cr, $ci)>: real and imaginary parts of the complex
coefficients, each a C<(nfreq, ntime)> double piddle.  A pure tone
produces its maximum power in the bin at its own frequency
(F<t/16_tfa.t>).

=head2 tfr_morlet($data, $sfreq, $freqs, %opt)

    my $power = tfr_morlet($epochs, $sfreq, $freqs, n_cycles => 7);
    my $itc   = tfr_morlet($epochs, $sfreq, $freqs, output => 'itc');
    my %r     = tfr_morlet($epochs, $sfreq, $freqs, output => 'all');

Trial-averaged Morlet power and related quantities.  C<$data> is a 1-D
piddle C<(ntime)> for a single trial or a 2-D piddle C<(nepoch, ntime)>.
Each trial is transformed with C<cwt_morlet>; C<n_cycles> and
C<sigma_win> are passed through to it.

Options:

=over 4

=item C<n_cycles> => number | piddle (nfreq) (default 7.0)

=item C<sigma_win> => number (default 5.0)

=item C<output> => string (default C<'power'>)

One of the following.  With C<c> the complex coefficient of one trial
and C<mean> the mean across trials:

=over 4

=item C<'power'>

C<(nfreq, ntime)>: C<mean(|c|^2)>, total power.

=item C<'itc'>

C<(nfreq, ntime)>: C<|mean(c / |c|)|>, inter-trial coherence.  A
coefficient with zero magnitude contributes 0 to the sum.

=item C<'evoked'>

C<(nfreq, ntime)>: C<|mean(c)|^2>, phase-locked power.

=item C<'induced'>

C<(nfreq, ntime)>: C<power - evoked>, non-phase-locked power; negative
values (rounding) are clamped to 0.

=item C<'complex'>

C<($er, $ei)>: real and imaginary parts of C<mean(c)>, each
C<(nfreq, ntime)>.

=item C<'all'>

A flat hash (list) with keys C<power>, C<itc>, C<evoked>, C<induced>,
C<evoked_r>, C<evoked_i>, the last two being the parts of C<mean(c)>.

=back

=back

Verified in F<t/16_tfa.t>: a Gaussian burst peaks at its own frequency
and time; ITC exceeds 0.9 for phase-locked trials and stays below 0.5
for random-phase trials; C<power == evoked + induced> to rounding
error.

=head2 tfr_superlet($data, $sfreq, $freqs, %opt)

    my $power = tfr_superlet($epochs, $sfreq, $freqs,
                             base_cycles => 3, order_min => 1, order_max => 6);

Trial-averaged multiplicative superlet power (Moca et al. 2021).
C<$data> is a 1-D piddle C<(ntime)> for a single trial or a 2-D piddle
C<(nepoch, ntime)>.

For each frequency C<f> and order C<o = 1 .. order(f)>, a Morlet wavelet
with C<n_cycles = base_cycles * o> is applied (via C<morlet_wavelet> and
the same FFT convolution as C<cwt_morlet>); the superlet power at C<f>
is the geometric mean of the C<order(f)> individual powers.  Individual
powers are floored at 1e-300 before taking the logarithm.  The result is
averaged across trials.

Options:

=over 4

=item C<base_cycles> => number (default 3)

Cycle count of the order-1 wavelet.

=item C<order> => integer (default undef)

Fixed order for every frequency.  When given, C<order_min>/C<order_max>
are ignored.

=item C<order_min>, C<order_max> => integer (defaults 1, 5)

Adaptive order: interpolated linearly from C<order_min> at
C<$freqs-E<gt>at(0)> to C<order_max> at the last element of C<$freqs>,
rounded to the nearest integer.  C<$freqs> is therefore expected to be
ascending; if the last element is not greater than the first, the span
is taken as 1 Hz.

=item C<sigma_win> => number (default 5.0)

Passed to C<morlet_wavelet>.

=back

Every order is clamped to at least 1.  Returns a C<(nfreq, ntime)>
power piddle.  Only power is produced (no ITC / evoked / induced
variants).  Cost is roughly C<mean(order)> times that of
C<tfr_morlet> with C<output =E<gt> 'power'>.

=head2 tfr_stat($data, $sfreq, $freqs, $times, $baseline, %opt)

    my %s = tfr_stat($epochs, $sfreq, $freqs, $times, [-0.05, -0.004],
                     n_cycles => 7);
    # $s{z}, $s{mean}, $s{sem} : (nfreq, ntime);  $s{n} : nepoch

Across-trial reliability of Morlet power, accumulated per trial without
storing single-trial maps.  C<$data> is C<(nepoch, ntime)> (a 1-D piddle
is treated as one trial); C<$times> is a C<(ntime)> piddle in seconds
aligned with the time axis; C<$baseline> is C<[t0, t1]> in seconds,
either bound may be C<undef> to extend to the edge of C<$times>.
Croaks if the baseline window selects no samples.

For each trial the power C<|c|^2> from C<cwt_morlet> is computed, the
mean over the baseline samples is subtracted per frequency, and the sum
and sum of squares of the corrected power are accumulated.  Then

    mean = S1 / n
    var  = (S2 - S1^2 / n) / (n - 1)      # floored at 1e-300
    sem  = sqrt(var / n)
    z    = mean / sem

Because the variance divides by C<n - 1>, at least two trials are
required for a meaningful result.

C<%opt> is passed through to C<cwt_morlet> (C<n_cycles>, C<sigma_win>).

Returns a hash with keys C<z>, C<mean>, C<sem> (each C<(nfreq, ntime)>)
and C<n> (the trial count).

This differs from C<apply_baseline(..., mode =E<gt> 'zscore')>: that
normalises the already trial-averaged map by the standard deviation
B<over time> within the baseline window, whereas C<tfr_stat> divides the
trial-averaged, baseline-corrected power by its standard error B<across
trials>.

=head2 apply_baseline($power, $times, $baseline, %opt)

    my $z = apply_baseline($power, $times, [-0.05, -0.004], mode => 'zscore');

Per-frequency baseline normalisation of a C<(nfreq, ntime)> power map.
C<$times> is a C<(ntime)> piddle in seconds aligned with the time axis;
C<$baseline> is C<[t0, t1]> in seconds, either bound may be C<undef> to
extend to the edge of C<$times>.  Croaks if the window selects no
samples.

Let C<m> be the mean of C<$power> over the baseline samples, computed
separately for each frequency.

Options:

=over 4

=item C<mode> => string (default C<'zscore'>)

=over 4

=item C<'zscore'>

C<(power - m) / sd>, where C<sd> is the standard deviation over the
baseline samples (divisor C<n_baseline>, not C<n_baseline - 1>); a zero
C<sd> is replaced by 1.

=item C<'ratio'>

C<power / m>

=item C<'logratio'>

C<log10(power / m)>

=item C<'percent'>

C<(power - m) / m>

=item C<'mean'>

C<power - m>

=back

Any other value croaks.

=back

Returns a new C<(nfreq, ntime)> piddle; C<$power> is not modified.
Verified in F<t/16_tfa.t>: the C<zscore> output averages to ~0 over
the baseline window, and every mode runs on a burst map without
producing BAD values.

=head1 SEE ALSO

F<t/16_tfa.t>, F<examples/sep_hfo_tfa.pl>.

=cut
