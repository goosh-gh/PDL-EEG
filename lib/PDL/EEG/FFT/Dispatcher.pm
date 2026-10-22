package PDL::EEG::FFT::Dispatcher;

use strict;
use warnings;
use Carp qw(croak);
use PDL;
use PDL::FFT ();
use Exporter 'import';

our @EXPORT_OK = qw(fft_len fft_ok fft_backend fft_inplace ifft_inplace);
our %EXPORT_TAGS = (all => \@EXPORT_OK);

# croak reports the line of the caller of these modules, not a line inside them
our @CARP_NOT = qw(PDL::EEG::Filter PDL::EEG::TFA PDL::EEG::Heartbeat::Respiration);

# 'auto'      : PDL::FFTDispatcher (PDL::PocketFFT) when it is installed; otherwise
#               PDL::FFTW3 when it is installed; otherwise PDL::FFT for the lengths
#               it transforms correctly, and a croak for the others
# 'pocketfft' : PDL::PocketFFT through PDL::FFTDispatcher (croaks without it)
# 'fftw3'     : PDL::FFTW3 for every length
# 'pdl'       : PDL::FFT only; a length it cannot transform croaks
our $BACKEND = 'auto';

my $NFACTOR = 11;     # size of factor[] in PDL's lib/PDL/FFT/fftn.c
my %NF;               # length -> number of factors
my ($HAVE_FFTW3, $HAVE_FD);

sub _fftw3 {
    $HAVE_FFTW3 //= eval { require PDL::FFTW3; 1 } ? 1 : 0;
    return $HAVE_FFTW3;
}

sub _fd {
    $HAVE_FD //= eval { require PDL::FFTDispatcher; 1 } ? 1 : 0;
    return $HAVE_FD;
}

# ---------------------------------------------------------------------------
# _nfactor($n): the number of entries factorize() in PDL's lib/PDL/FFT/fftn.c
# writes into factor[] for a transform of length $n. Same counts as the C
# code, without its loop over every integer up to the last prime factor:
#   - each factor 16 gives a 4 and each square p^2 of an odd prime gives a p;
#     these are written twice (the second half mirrors the first);
#   - if what is left is at most 4, it is one more entry (unless it is 1);
#   - otherwise a remaining factor 4 gives a 2 (also mirrored), and each prime
#     left (each appears once) is one entry.
# ---------------------------------------------------------------------------
sub _distinct_primes {
    my $n = shift;
    my $c = 0;
    for (my $p = 2; $p * $p <= $n; $p += ($p == 2 ? 1 : 2)) {
        next if $n % $p;
        $c++;
        $n /= $p while $n % $p == 0;
    }
    $c++ if $n > 1;
    return $c;
}

sub _nfactor {
    my $n  = shift;
    my $sq = 0;
    while ($n % 16 == 0) { $sq++; $n /= 16 }
    my $j = 3;
    do {
        while ($n % ($j * $j) == 0) { $sq++; $n /= $j * $j }
        $j += 2;
    } while ($j * $j <= $n);
    return 2 * $sq + ($n != 1 ? 1 : 0) if $n <= 4;
    my $two = 0;
    if ($n % 4 == 0) { $two = 1; $n /= 4 }
    return 2 * ($sq + $two) + _distinct_primes($n);
}

sub _is235 {
    my $r = shift;
    for my $p (2, 3, 5) { $r = int($r / $p) while $r % $p == 0 }
    return $r == 1;
}

# fft_ok($n): true when PDL::FFT transforms length $n correctly
sub fft_ok {
    my $n = shift;
    croak "fft_ok: length must be a positive integer" unless $n >= 1 && $n == int $n;
    return ($NF{$n} //= _nfactor($n)) <= $NFACTOR;
}

# fft_len($n): smallest length >= $n whose only prime factors are 2, 3 and 5
# and that PDL::FFT transforms correctly
sub fft_len {
    my $n = shift;
    croak "fft_len: length must be a positive integer" unless $n >= 1 && $n == int $n;
    for (my $m = $n; ; $m++) {
        return $m if _is235($m) && fft_ok($m);
    }
}

# fft_backend($n): 'PDL::PocketFFT', 'PDL::FFTW3' or 'PDL::FFT', the library
# that transforms length $n under $BACKEND; croaks when none can
sub fft_backend {
    my $n = shift;
    croak "PDL::EEG::FFT::Dispatcher: \$BACKEND must be 'auto', 'pocketfft', 'fftw3' or 'pdl', not '$BACKEND'"
        unless $BACKEND =~ /\A(?:auto|pocketfft|fftw3|pdl)\z/;
    if ($BACKEND eq 'pocketfft') {
        croak "PDL::EEG::FFT::Dispatcher: \$BACKEND is 'pocketfft' but PDL::FFTDispatcher is not installed"
            unless _fd();
        return 'PDL::PocketFFT';
    }
    if ($BACKEND eq 'fftw3') {
        croak "PDL::EEG::FFT::Dispatcher: \$BACKEND is 'fftw3' but PDL::FFTW3 is not installed" unless _fftw3();
        return 'PDL::FFTW3';
    }
    if ($BACKEND eq 'auto') {
        return 'PDL::PocketFFT' if _fd();
        return 'PDL::FFTW3'     if _fftw3();
    }
    return 'PDL::FFT' if fft_ok($n);
    croak sprintf "PDL::EEG::FFT::Dispatcher: PDL::FFT cannot transform length %d (it splits it into %d factors, "
        . "its limit is %d) and %s", $n, $NF{$n}, $NFACTOR,
        $BACKEND eq 'pdl' ? "it is the only library allowed (\$PDL::EEG::FFT::Dispatcher::BACKEND = 'pdl')"
                          : 'neither PDL::FFTDispatcher nor PDL::FFTW3 is installed';
}

# ---------------------------------------------------------------------------
# fft_inplace($re, $im), ifft_inplace($re, $im)
#
# Complex FFT along dim 0, in place, broadcasting over the other dims; the
# conventions of PDL::FFT::fft and ifft (forward unscaled, inverse scaled by
# 1/n). The library is chosen (and a length PDL::FFT cannot transform is
# refused) here; the transform itself is done by PDL::FFTDispatcher when it is
# installed.
# ---------------------------------------------------------------------------
my %FD_BACKEND = ('PDL::PocketFFT' => 'pocketfft', 'PDL::FFTW3' => 'fftw3', 'PDL::FFT' => 'pdl');

sub fft_inplace  { _transform(@_, 0) }
sub ifft_inplace { _transform(@_, 1) }

sub _transform {
    my ($re, $im, $inv) = @_;
    croak "fft_inplace: re and im must have the same dims"
        unless join(',', $re->dims) eq join(',', $im->dims);
    my $be = fft_backend($re->dim(0));
    if (_fd()) {
        local $PDL::FFTDispatcher::BACKEND = $FD_BACKEND{$be};
        $inv ? PDL::FFTDispatcher::ifft_inplace($re, $im) : PDL::FFTDispatcher::fft_inplace($re, $im);
        return;
    }
    if ($be eq 'PDL::FFT') {
        $inv ? PDL::FFT::ifft($re, $im) : PDL::FFT::fft($re, $im);
        return;
    }
    my $z = czip($re, $im);
    $z = $inv ? PDL::FFTW3::ifft1($z) : PDL::FFTW3::fft1($z);
    $re .= $z->re;
    $im .= $z->im;
    return;
}

1;

__END__

=head1 NAME

PDL::EEG::FFT::Dispatcher - picks the FFT library: PDL::FFTDispatcher (PDL::PocketFFT) or PDL::FFTW3 when installed, otherwise PDL::FFT for the lengths it transforms correctly

=head1 SYNOPSIS

    use PDL::EEG::FFT::Dispatcher qw(fft_len fft_inplace ifft_inplace);

    my $L = fft_len($nt);                 # padded length
    my ($re, $im) = (zeroes(double, $L, $nch), zeroes(double, $L, $nch));
    $re->slice("0:" . ($nt - 1)) .= $X;   # $X : (nt, nch)
    fft_inplace($re, $im);
    ...
    ifft_inplace($re, $im);

=head1 DESCRIPTION

The transforms of L<PDL::EEG::Filter>, L<PDL::EEG::TFA> and
L<PDL::EEG::Heartbeat::Respiration> go through this module.

L<PDL::FFT> splits the length into at most 11 factors (C<NFACTOR> in
F<lib/PDL/FFT/fftn.c>). For a length with more factors it prints
C<Error: FFTRADIX() - exceeded number of factors> and returns the input
unchanged (Linux) or aborts the process (macOS). Lengths of the form
2^a 3^b 5^c are affected too: 51 of the 493 such lengths between 1000 and
2,000,000, the smallest 157464 = 2^3 3^9. Checked with PDL 2.085, 2.104 and
2.106.

L<PDL::FFTDispatcher> and L<PDL::PocketFFT> (github.com/goosh-gh) and
L<PDL::FFTW3> are optional; with neither installed, PDL::FFT is used for the
lengths it transforms correctly and other lengths croak.

=over 4

=item fft_len($n)

The smallest length E<gt>= C<$n> whose only prime factors are 2, 3 and 5 and
that PDL::FFT transforms correctly.

=item fft_ok($n)

True when PDL::FFT transforms length C<$n> correctly, i.e. its factorization
in F<fftn.c> has at most 11 factors. The count gives the same result as the
C code.

=item fft_inplace($re, $im), ifft_inplace($re, $im)

Complex FFT along dim 0, in place, broadcasting over the other dims, with the
conventions of C<PDL::FFT::fft> and C<ifft>: forward unscaled, inverse scaled
by 1/n. C<$re> and C<$im> must have the same dims.

The library is the one C<fft_backend> names. When L<PDL::FFTDispatcher> is
installed, it carries out the transform with that library. When PDL::FFT is
the library chosen, a length it cannot transform croaks with a message that
names the length; it does not return a wrong result.

=item fft_backend($n)

C<'PDL::PocketFFT'>, C<'PDL::FFTW3'> or C<'PDL::FFT'>: the library that
transforms length C<$n>; croaks as C<fft_inplace> does.

=back

=head2 $PDL::EEG::FFT::Dispatcher::BACKEND

C<'auto'> (default): PDL::PocketFFT through L<PDL::FFTDispatcher> when it is
installed, otherwise PDL::FFTW3 when it is installed, otherwise PDL::FFT for
the lengths C<fft_ok> accepts and a croak for the others. C<'pocketfft'> uses
PDL::PocketFFT and croaks when PDL::FFTDispatcher is not installed.
C<'fftw3'> uses PDL::FFTW3 and croaks when it is not installed. C<'pdl'> uses
PDL::FFT only (a croak for the lengths C<fft_ok> rejects).

=cut
