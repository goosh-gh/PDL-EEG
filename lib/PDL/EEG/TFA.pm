package PDL::EEG::TFA;

use strict;
use warnings;
use PDL;
use PDL::EEG::FFT::Dispatcher qw(fft_len fft_inplace ifft_inplace);
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

my $PI_TFA = 3.14159265358979;

# ------------------------------------------------------------------
# _half_len($freq, $sfreq, $n_cycles, $sigma_win) -> samples each side
# of t = 0 in the truncated kernel (kernel length = 2*half + 1).
# ------------------------------------------------------------------
sub _half_len {
    my ($freq, $sfreq, $n_cycles, $sigma_win) = @_;
    my $sigma_t = $n_cycles / (2 * $PI_TFA * $freq);
    my $half    = int($sigma_win * $sigma_t * $sfreq);
    return $half < 1 ? 1 : $half;
}

# ------------------------------------------------------------------
# _ncyc_list($n_cycles, $nf) -> list of nf cycle counts
# ------------------------------------------------------------------
sub _ncyc_list {
    my ($ncyc, $nf) = @_;
    $ncyc //= 7.0;
    my @ncyc = ref($ncyc) && eval { $ncyc->isa('PDL') } ? list($ncyc)
             : (($ncyc) x $nf);
    croak "n_cycles length mismatch" unless @ncyc == $nf;
    return @ncyc;
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

    my $sigma_t = $n_cycles / (2 * $PI_TFA * $freq);
    my $half    = _half_len($freq, $sfreq, $n_cycles, $sigma_win);  # samples each side
    my $n       = 2 * $half + 1;                              # odd length

    my $t   = (sequence($n) - $half) / $sfreq;               # centred time axis
    my $env = exp(-($t**2) / (2 * $sigma_t**2));
    my $ang = 2 * $PI_TFA * $freq * $t;
    my $wr  = $env * cos($ang);
    my $wi  = $env * sin($ang);

    # unit L2 energy
    my $norm = sqrt( ($wr**2 + $wi**2)->sum );
    $norm = 1 if $norm == 0;
    return ($wr / $norm, $wi / $norm);
}

# ------------------------------------------------------------------
# _kernel_fft($freq, $sfreq, $n_cycles, $sigma_win, $nfft)
#     -> ($Wr, $Wi, $half)
#
# FFT of the truncated, unit-energy Morlet kernel, zero-padded to
# $nfft with the kernel occupying samples 0 .. 2*half.
# ------------------------------------------------------------------
sub _kernel_fft {
    my ($freq, $sfreq, $n_cycles, $sigma_win, $nfft) = @_;
    my ($wr, $wi) = morlet_wavelet($freq, $sfreq,
                        n_cycles => $n_cycles, sigma_win => $sigma_win);
    my $nw = $wr->dim(0);
    my $Wr = zeroes(double, $nfft); $Wr->slice("0:@{[$nw-1]}") .= $wr;
    my $Wi = zeroes(double, $nfft); $Wi->slice("0:@{[$nw-1]}") .= $wi;
    fft_inplace($Wr, $Wi);
    return ($Wr, $Wi, ($nw - 1) / 2);
}

# Upper bound on elements of one (nfft, nsignal) work array; signals
# beyond this are processed in successive blocks.
our $BLOCK_ELEMS = 2**18;

# ------------------------------------------------------------------
# _cwt_run($sig, $sfreq, $specs, $sigma_win, $cb)
#
# Shared convolution engine.  $sig is real, (ntime, nsignal), time on
# dim 0.  $specs is a list of [freq, n_cycles].  For each block of
# signals and each spec k the callback is called as
#
#     $cb->($k, $cr, $ci, $first)
#
# with $cr, $ci the (ntime, nblock) real and imaginary coefficients
# ('same'-length linear convolution, zero-padded edges) and $first
# the index of the block's first signal.  Specs are visited in order
# within a block.  The piddles are only valid during the call.
#
# The signal block is transformed once and reused for every spec.
# One FFT length serves all specs: aliasing of the circular
# convolution only reaches the first (L - nfft) output samples of the
# full linear convolution, and those lie before the 'same' window as
# long as nfft >= ntime + half.
# ------------------------------------------------------------------
sub _cwt_run {
    my ($sig, $sfreq, $specs, $sigma_win, $cb) = @_;
    my ($nt, $ns) = $sig->dims;
    my $nspec = scalar @$specs;

    my $maxhalf = 0;
    for my $s (@$specs) {
        croak "cwt: freq must be > 0" unless $s->[0] > 0;
        my $h = _half_len($s->[0], $sfreq, $s->[1], $sigma_win);
        $maxhalf = $h if $h > $maxhalf;
    }
    my $need = $nt + $maxhalf;
    $need = 2 * $maxhalf + 1 if $need < 2 * $maxhalf + 1;
    my $nfft = fft_len($need);

    my $block = int($BLOCK_ELEMS / $nfft);
    $block = 1   if $block < 1;
    $block = $ns if $block > $ns;

    # keep kernel FFTs across blocks only when there is more than one
    # block and the table stays small
    my $keep = ($block < $ns && $nspec * $nfft <= 4 * $BLOCK_ELEMS);
    my @kern;

    for (my $s0 = 0; $s0 < $ns; $s0 += $block) {
        my $s1 = $s0 + $block - 1;
        $s1 = $ns - 1 if $s1 > $ns - 1;
        my $nb = $s1 - $s0 + 1;

        my $Sr = zeroes(double, $nfft, $nb);
        $Sr->slice("0:@{[$nt-1]}") .= $sig->slice(":,$s0:$s1");
        my $Si = zeroes(double, $nfft, $nb);
        fft_inplace($Sr, $Si);

        for my $k (0 .. $nspec - 1) {
            my ($Wr, $Wi, $half);
            if ($keep && $kern[$k]) {
                ($Wr, $Wi, $half) = @{ $kern[$k] };
            } else {
                ($Wr, $Wi, $half) = _kernel_fft($specs->[$k][0], $sfreq,
                                        $specs->[$k][1], $sigma_win, $nfft);
                $kern[$k] = [$Wr, $Wi, $half] if $keep;
            }
            my $Pr = $Sr * $Wr - $Si * $Wi;
            my $Pi = $Sr * $Wi + $Si * $Wr;
            ifft_inplace($Pr, $Pi);
            my $win = "$half:@{[$half + $nt - 1]}";
            $cb->($k, $Pr->slice($win), $Pi->slice($win), $s0);
        }
    }
    return;
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
    my @ncyc = _ncyc_list($opt{n_cycles}, $nf);

    my $cr  = zeroes(double, $nf, $nt);
    my $ci  = zeroes(double, $nf, $nt);

    _cwt_run($signal->dummy(1), $sfreq,
             [ map { [ $f[$_], $ncyc[$_] ] } 0 .. $nf - 1 ],
             ($opt{sigma_win} // 5.0),
             sub {
                 my ($k, $r, $i) = @_;
                 $cr->slice("($k),:") .= $r->slice(":,(0)");
                 $ci->slice("($k),:") .= $i->slice(":,(0)");
             });
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
    my @f  = list($freqs);
    my $nf = scalar @f;
    my @ncyc = _ncyc_list($opt{n_cycles}, $nf);

    my $want_pow = ($out eq 'power' or $out eq 'induced' or $out eq 'all');
    my $want_itc = ($out eq 'itc'   or $out eq 'all');
    my $want_evo = ($out ne 'power' and $out ne 'itc');

    my $power    = zeroes(double, $nf, $nt);        # sum |c|^2
    my $itc_r    = zeroes(double, $nf, $nt);        # sum c/|c| (real)
    my $itc_i    = zeroes(double, $nf, $nt);
    my $evoked_r = zeroes(double, $nf, $nt);        # sum c (real)
    my $evoked_i = zeroes(double, $nf, $nt);

    _cwt_run($data->xchg(0,1), $sfreq,
             [ map { [ $f[$_], $ncyc[$_] ] } 0 .. $nf - 1 ],
             ($opt{sigma_win} // 5.0),
             sub {
                 my ($k, $r, $i) = @_;                # (ntime, nblock)
                 my $row = "($k),:";
                 if ($want_evo) {
                     $evoked_r->slice($row) += $r->xchg(0,1)->sumover;
                     $evoked_i->slice($row) += $i->xchg(0,1)->sumover;
                 }
                 return unless $want_pow or $want_itc;
                 my $mag2 = $r * $r + $i * $i;         # (**2 is ~7x slower)
                 $power->slice($row) += $mag2->xchg(0,1)->sumover if $want_pow;
                 if ($want_itc) {
                     my $mag = sqrt($mag2);
                     $mag->where($mag == 0) .= 1;
                     $itc_r->slice($row) += ($r / $mag)->xchg(0,1)->sumover;
                     $itc_i->slice($row) += ($i / $mag)->xchg(0,1)->sumover;
                 }
             });

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

    my $power = zeroes(double, $nf, $nt);

    # one spec per (frequency, order); the engine visits them in order,
    # so the log-powers of one frequency arrive consecutively
    my (@specs, @which);
    for my $k (0 .. $nf - 1) {
        for my $o (1 .. $order[$k]) {
            push @specs, [ $f[$k], $base * $o ];
            push @which, [ $k, $o ];
        }
    }

    my $logsum;                                      # (ntime, nblock)
    _cwt_run($data->xchg(0,1), $sfreq, \@specs, $swin,
             sub {
                 my ($s, $r, $i) = @_;
                 my ($k, $o) = @{ $which[$s] };
                 my $p = $r * $r + $i * $i;
                 $p->where($p <= 0) .= 1e-300;       # geometric-mean guard
                 if ($o == 1) { $logsum = log($p) } else { $logsum += log($p) }
                 $power->slice("($k),:")
                     += exp($logsum / $order[$k])->xchg(0,1)->sumover   # geo mean
                     if $o == $order[$k];
             });
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
# [t0,t1].  %opt: n_cycles, sigma_win (as cwt_morlet).
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

    my @f    = list($freqs);
    my @ncyc = _ncyc_list($opt{n_cycles}, $nf);

    _cwt_run($data->xchg(0,1), $sfreq,
             [ map { [ $f[$_], $ncyc[$_] ] } 0 .. $nf - 1 ],
             ($opt{sigma_win} // 5.0),
             sub {
                 my ($k, $r, $i) = @_;                           # (ntime, nblock)
                 my $p    = $r * $r + $i * $i;
                 my $base = $p->dice_axis(0, $bidx)->avgover;    # (nblock) baseline mean
                 my $pc   = ($p - $base->dummy(0))->xchg(0,1);   # (nblock, ntime)
                 $S1->slice("($k),:") += $pc->sumover;
                 $S2->slice("($k),:") += ($pc * $pc)->sumover;
             });

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
carried out in the frequency domain via L<PDL::EEG::FFT::Dispatcher>.

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

The signal is transformed once by zero-padded FFT and the result is
reused for every frequency; each frequency multiplies it by the FFT of
its wavelet and transforms back, and the central C<ntime> samples of the
linear convolution are kept ("same" mode).  One FFT length is used for
all frequencies: C<fft_len> of L<PDL::EEG::FFT::Dispatcher> of C<ntime> plus the
half-length of the longest kernel, the smallest product of powers of 2, 3
and 5 at least that long that PDL::FFT transforms correctly.  Because the
padding is with zeros, coefficients within half a kernel length of either
edge are affected by the signal boundary.

The coefficients equal those of a direct time-domain convolution with
the kernel of C<morlet_wavelet> (F<t/16_tfa.t>).

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
All trials are transformed together (the FFT runs along time for a block
of trials at once; see L</BLOCK SIZE>), with the same kernels and edge
handling as C<cwt_morlet>; the result equals transforming each trial
with C<cwt_morlet> and combining (F<t/16_tfa.t>).  Only the sums needed
for the requested C<output> are accumulated.

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

Across-trial reliability of Morlet power, accumulated without storing
single-trial maps (one frequency of one block of trials is held at a
time; see L</BLOCK SIZE>).  C<$data> is C<(nepoch, ntime)> (a 1-D piddle
is treated as one trial); C<$times> is a C<(ntime)> piddle in seconds
aligned with the time axis; C<$baseline> is C<[t0, t1]> in seconds,
either bound may be C<undef> to extend to the edge of C<$times>.
Croaks if the baseline window selects no samples.

For each trial the power C<|c|^2> (same coefficients as C<cwt_morlet>) is computed, the
mean over the baseline samples is subtracted per frequency, and the sum
and sum of squares of the corrected power are accumulated.  Then

    mean = S1 / n
    var  = (S2 - S1^2 / n) / (n - 1)      # floored at 1e-300
    sem  = sqrt(var / n)
    z    = mean / sem

Because the variance divides by C<n - 1>, at least two trials are
required for a meaningful result.

C<%opt> takes C<n_cycles> and C<sigma_win>, as in C<cwt_morlet>.

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

=head1 BLOCK SIZE

C<tfr_morlet>, C<tfr_superlet> and C<tfr_stat> transform the trials in
blocks.  A block holds as many trials as fit in
C<$PDL::EEG::TFA::BLOCK_ELEMS> elements (default 2**18) of one
C<(nfft, ntrial)> work array, and at least one.  The results do not
depend on the block size beyond rounding (F<t/16_tfa.t>); it only trades
memory against speed.

=head1 SEE ALSO

F<t/16_tfa.t>, F<examples/sep_hfo_tfa.pl>.

=cut
