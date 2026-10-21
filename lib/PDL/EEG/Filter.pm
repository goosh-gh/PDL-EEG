package PDL::EEG::Filter;

use strict;
use warnings;
use Carp qw(croak);
use PDL;
use PDL::EEG::FFT::Dispatcher qw(fft_len fft_inplace ifft_inplace);
use Exporter 'import';

our @EXPORT_OK = qw(zp_butter analytic_signal envelope inst_phase);
our %EXPORT_TAGS = (all => \@EXPORT_OK);

my $PI = 3.14159265358979323846;

# ---------------------------------------------------------------------------
# zp_butter($x, fs => $hz, low => $f1, high => $f2, order => $n [, pad => $sec])
#
# Zero-phase Butterworth filter, applied in the frequency domain.
#
#   $x : (nt) or (nch, nt)
#   low only          -> high-pass at `low`
#   high only         -> low-pass  at `high`
#   low and high      -> band-pass  low..high
#
# The gain applied to each FFT bin is the squared magnitude of the digital
# (bilinear-transform) Butterworth filter of the given order,
#
#     G(f) = 1 / (1 + q(f)^(2*order)),   W = tan(pi f / fs)
#     q = Wc/W (high-pass),  W/Wc (low-pass),  (W^2 - W1 W2)/((W2 - W1) W) (band-pass)
#
# which is what a forward-backward pass (scipy.signal.filtfilt of
# scipy.signal.butter(order, ...)) has in steady state. The record is extended
# at both ends by odd reflection (`pad` seconds, default 3/low or 3/high,
# at most nt-1 samples) and the two extension ends are joined by a straight line
# up to an FFT-friendly length, so the circular FFT sees no step.
# ---------------------------------------------------------------------------
sub zp_butter {
    my ($x, %o) = @_;
    my $fs = $o{fs} or croak "zp_butter: fs required";
    my ($lo, $hi) = @o{qw(low high)};
    croak "zp_butter: give low and/or high" unless defined $lo || defined $hi;
    my $N  = $o{order} // 2;
    my $nd = $x->ndims;
    croak "zp_butter: x must be (nt) or (nch,nt)" unless $nd == 1 || $nd == 2;

    my $X  = ($nd == 1 ? $x->dummy(0) : $x)->mv(1, 0)->double;   # (nt, nch)
    my $nt = $X->dim(0);
    my $fc = defined $lo ? $lo : $hi;
    my $pad = int((defined $o{pad} ? $o{pad} : 3 / $fc) * $fs);
    $pad = $nt - 1 if $pad > $nt - 1;
    my ($re, $L) = _extend($X, $pad);
    my $im = zeroes(double, $re->dims);
    fft_inplace($re, $im);

    my $k = sequence($L);
    my $f = $k * ($fs / $L);
    $f = $f->where($k <= $L / 2)->append(($fs - $f)->where($k > $L / 2));
    my $W = tan($PI * $f / $fs);
    $W = $W->clip(1e-300, 1e300);
    my $q;
    if (defined $lo && defined $hi) {
        my ($W1, $W2) = map { tan($PI * $_ / $fs) } $lo, $hi;
        $q = ($W * $W - $W1 * $W2) / (($W2 - $W1) * $W);
    } elsif (defined $lo) {
        $q = tan($PI * $lo / $fs) / $W;
    } else {
        $q = $W / tan($PI * $hi / $fs);
    }
    my $g = 1 / (1 + ($q * $q) ** $N);                            # (L)
    $re *= $g; $im *= $g;
    ifft_inplace($re, $im);

    my $y = $re->slice("$pad:" . ($pad + $nt - 1))->mv(0, 1)->sever;   # (nch, nt)
    return $nd == 1 ? $y->slice("(0)")->sever : $y;
}

# ---------------------------------------------------------------------------
# ($re, $im) = analytic_signal($x, %opt)
#
#   $x        : (nt) or (nch, nt), real
#   $re, $im  : same shape; $re is $x, $im its Hilbert transform
#
# The spectrum is taken, the negative frequencies are set to zero and the
# positive ones doubled, and the result is transformed back.
#
# The channel mean is removed first (its Hilbert transform is 0) and the
# record is followed by zeros up to fft_len(nt) of PDL::EEG::FFT::Dispatcher, the next
# length whose only prime factors are 2, 3 and 5 and that PDL::FFT transforms
# correctly; nothing is added when nt is such a length already.
#
#   exact => 1    no zeros: one FFT of length nt, the record taken as one
#                 period. This is scipy.signal.hilbert. The transform is
#                 PDL::FFTW3 when it is installed; without it PDL::FFT, which
#                 needs time proportional to nt^2 when nt has a large prime
#                 factor and croaks for a length it cannot transform
#                 (PDL::EEG::FFT::Dispatcher).
# ---------------------------------------------------------------------------
sub analytic_signal {
    my ($x, %o) = @_;
    my $nd = $x->ndims;
    croak "analytic_signal: x must be (nt) or (nch,nt)" unless $nd == 1 || $nd == 2;
    my $X  = ($nd == 1 ? $x->dummy(0) : $x)->mv(1, 0)->double;   # (nt, nch)
    my $nt = $X->dim(0);
    croak "analytic_signal: at least 2 samples needed" if $nt < 2;

    my $L    = $o{exact} ? $nt : fft_len($nt);
    my $mean = $X->average;                                      # (nch)
    my $re   = zeroes(double, $L, $X->dim(1));
    $re->slice("0:" . ($nt - 1)) .= $X - $mean->dummy(0);
    my $im = zeroes(double, $re->dims);
    fft_inplace($re, $im);

    my $h = zeroes(double, $L);                                  # 1, 2, ..., 2, (1), 0, ..., 0
    $h->set(0, 1);
    my $half = int(($L + 1) / 2);                                # first index past the positive frequencies
    $h->slice("1:" . ($half - 1)) .= 2 if $half > 1;
    $h->set($L / 2, 1) if $L % 2 == 0;
    $re *= $h; $im *= $h;
    ifft_inplace($re, $im);

    my @out = map { $_->slice("0:" . ($nt - 1))->mv(0, 1)->sever } $re, $im;
    $out[0] += $mean->dummy(1);
    @out = map { $_->slice("(0)")->sever } @out if $nd == 1;
    return @out;
}

# envelope($x, %opt): magnitude of the analytic signal, same shape as $x
sub envelope {
    my ($x, %o) = @_;
    my ($re, $im) = analytic_signal($x, %o);
    return sqrt($re * $re + $im * $im);
}

# inst_phase($x, %opt): angle of the analytic signal in radians, -pi..pi,
# 0 at a positive maximum of a narrow-band $x; same shape as $x
sub inst_phase {
    my ($x, %o) = @_;
    my ($re, $im) = analytic_signal($x, %o);
    return atan2($im, $re);
}

# ---------------------------------------------------------------------------
# ($ext, $L) = _extend($X, $pad)
#
#   $X   : (nt, nch) double
#   $ext : (L, nch): $pad samples of odd reflection, $X, $pad samples of odd
#          reflection, then a straight line from the last extended sample back
#          to the first. $X sits at $pad .. $pad+nt-1.
# ---------------------------------------------------------------------------
sub _extend {
    my ($X, $pad) = @_;
    my $nt = $X->dim(0);
    my $L  = fft_len($nt + 2 * $pad + 16);
    my $ng = $L - $nt - 2 * $pad;                                # join length

    my $re = zeroes(double, $L, $X->dim(1));
    $re->slice("$pad:" . ($pad + $nt - 1)) .= $X;
    if ($pad > 0) {
        my $x0 = $X->slice("0:0");
        my $xn = $X->slice("-1:-1");
        $re->slice("0:" . ($pad - 1)) .= 2 * $x0 - $X->slice("$pad:1");
        $re->slice(($pad + $nt) . ":" . (2 * $pad + $nt - 1))
            .= 2 * $xn - $X->slice("-2:" . ($nt - 1 - $pad));
    }
    my $last  = $re->slice("(" . (2 * $pad + $nt - 1) . ")");    # (nch)
    my $first = $re->slice("(0)");
    my $w = (sequence($ng) + 1) / ($ng + 1);                     # (ng)
    $re->slice((2 * $pad + $nt) . ":-1")
        .= $last->dummy(0) + ($first - $last)->dummy(0) * $w->dummy(1);
    return ($re, $L);
}

1;

__END__

=head1 NAME

PDL::EEG::Filter - temporal filters for multichannel EEG: zero-phase Butterworth, analytic signal

=head1 SYNOPSIS

    use PDL::EEG::Filter qw(zp_butter analytic_signal envelope inst_phase);

    my $y   = zp_butter($X, fs => 1000, low => 0.5, high => 45, order => 4);   # (nch, nt)
    my $alpha = zp_butter($X, fs => 1000, low => 8, high => 13, order => 3);
    my $env = envelope($alpha);                # (nch, nt)
    my $ph  = inst_phase($alpha);              # radians, -pi..pi
    my ($re, $im) = analytic_signal($alpha);

=head1 DESCRIPTION

Filters that act along the time axis. Input is C<(nt)> or C<(nch, nt)>, the
layout of C<read_edf>; the output has the shape of the input. Filters that
combine channels are in L<PDL::EEG::Derivation>, L<PDL::EEG::GED>,
L<PDL::EEG::ICA> and L<PDL::EEG::Regress>.

All transforms go through L<PDL::EEG::FFT::Dispatcher>. The running time of L<PDL::FFT>
grows with the square of the largest prime factor of the length, so the record
is extended to C<fft_len> of L<PDL::EEG::FFT::Dispatcher>: a length whose only prime
factors are 2, 3 and 5 and that PDL::FFT transforms correctly (it fails on
some such lengths, e.g. 354294).

=over 4

=item zp_butter($x, fs=>, low=>, high=>, order=>, pad=>)

Zero-phase Butterworth high-pass (C<low> only), low-pass (C<high> only) or
band-pass (both), applied in the frequency domain with the squared-magnitude
response of the digital Butterworth filter of the given order (default 2).
This is the steady-state response of a forward-backward pass
(C<scipy.signal.filtfilt> of C<scipy.signal.butter(order, ...)>). The record
is extended at both ends by odd reflection: C<pad> seconds, default C<3/low>
or C<3/high>, at most C<nt-1> samples.

=item analytic_signal($x, exact=>)

Returns C<($re, $im)>: C<$re> is C<$x> and C<$im> its Hilbert transform.

The channel mean is removed before the transform and the record is followed by
zeros up to C<fft_len(nt)> (see above); when C<nt> is such a length nothing
is added. C<< exact => 1 >> transforms the record at its own length, as one
period of a periodic signal; this is C<scipy.signal.hilbert>. The transform
is PDL::FFTW3 when it is installed; without it PDL::FFT, which takes long when
C<nt> has a large prime factor and croaks for a length it cannot transform.

The Hilbert transform is not local. The values near both ends of the record
depend on what follows the record (zeros, or the record again), in either
mode.

=item envelope($x, %opt)

C<sqrt(re**2 + im**2)> of the analytic signal. Options as for
C<analytic_signal>.

=item inst_phase($x, %opt)

C<atan2(im, re)> of the analytic signal in radians, C<-pi..pi>, not unwrapped.
For a narrow-band signal the phase is 0 at a positive maximum and C<+-pi> at a
negative one. Options as for C<analytic_signal>.

=back

=cut
