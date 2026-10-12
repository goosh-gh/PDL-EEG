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

# 'auto'  : PDL::FFT when it can transform the length, PDL::FFTW3 otherwise
# 'pdl'   : PDL::FFT only; a length it cannot transform croaks
# 'fftw3' : PDL::FFTW3 for every length
our $BACKEND = 'auto';

my $NFACTOR = 11;     # size of factor[] in PDL's lib/PDL/FFT/fftn.c
my %NF;               # length -> number of factors
my $HAVE_FFTW3;

# ---------------------------------------------------------------------------
# _nfactor($n): the number of factors that factorize() in PDL's
# lib/PDL/FFT/fftn.c produces for a transform of length $n, step by step as in
# the C code. Square factors (16, and squares of odd primes) are counted twice:
# once in the first half and once in the mirrored second half.
# ---------------------------------------------------------------------------
sub _nfactor {
    my $n = shift;
    my @f;
    my $kt;
    while ($n % 16 == 0) { push @f, 4; $n = int($n / 16) }
    my ($j, $jj) = (3, 9);
    do {
        while ($n % $jj == 0) { push @f, $j; $n = int($n / $jj) }
        $j += 2;
        $jj = $j * $j;
    } while ($jj <= $n);
    if ($n <= 4) {
        $kt = @f;
        push @f, $n if $n != 1;
    } else {
        if ($n % 4 == 0) { push @f, 2; $n = int($n / 4) }
        $kt = @f;
        $j = 2;
        do {
            if ($n % $j == 0) { push @f, $j; $n = int($n / $j) }
            $j = (int(($j + 1) / 2) << 1) + 1;
        } while ($j <= $n);
    }
    return @f + $kt;
}

sub _is235 {
    my $r = shift;
    for my $p (2, 3, 5) { $r = int($r / $p) while $r % $p == 0 }
    return $r == 1;
}

sub _fftw3 {
    $HAVE_FFTW3 //= eval { require PDL::FFTW3; 1 } ? 1 : 0;
    return $HAVE_FFTW3;
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

# fft_backend($n): 'PDL::FFT' or 'PDL::FFTW3', the library that transforms
# length $n under $BACKEND; croaks when neither can
sub fft_backend {
    my $n = shift;
    croak "PDL::EEG::FFT::Dispatcher: \$BACKEND must be 'auto', 'pdl' or 'fftw3', not '$BACKEND'"
        unless $BACKEND =~ /\A(?:auto|pdl|fftw3)\z/;
    if ($BACKEND eq 'fftw3') {
        croak "PDL::EEG::FFT::Dispatcher: \$BACKEND is 'fftw3' but PDL::FFTW3 is not installed" unless _fftw3();
        return 'PDL::FFTW3';
    }
    return 'PDL::FFT' if fft_ok($n);
    return 'PDL::FFTW3' if $BACKEND eq 'auto' && _fftw3();
    croak sprintf "PDL::EEG::FFT::Dispatcher: PDL::FFT cannot transform length %d (it splits it into %d factors, "
        . "its limit is %d) and PDL::FFTW3 is %s", $n, $NF{$n}, $NFACTOR,
        $BACKEND eq 'pdl' ? "switched off (\$PDL::EEG::FFT::Dispatcher::BACKEND = 'pdl')" : "not installed";
}

# ---------------------------------------------------------------------------
# fft_inplace($re, $im), ifft_inplace($re, $im)
#
# Complex FFT along dim 0, in place, broadcasting over the other dims; the
# conventions of PDL::FFT::fft and ifft (forward unscaled, inverse scaled by
# 1/n).
# ---------------------------------------------------------------------------
sub fft_inplace  { _transform(@_, 0) }
sub ifft_inplace { _transform(@_, 1) }

sub _transform {
    my ($re, $im, $inv) = @_;
    croak "fft_inplace: re and im must have the same dims"
        unless join(',', $re->dims) eq join(',', $im->dims);
    if (fft_backend($re->dim(0)) eq 'PDL::FFT') {
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

PDL::EEG::FFT::Dispatcher - picks the FFT library for each length: PDL::FFT where it works, PDL::FFTW3 where it does not

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
C<Error: FFTRADIX() - exceeded number of factors> and returns a wrong result
(Linux) or aborts the process (macOS). Lengths of the form 2^a 3^b 5^c are
affected too: 51 of the 493 such lengths between 1000 and 2,000,000, the
smallest 157464 = 2^3 3^9. Checked with PDL 2.085, 2.104 and 2.106.

=over 4

=item fft_len($n)

The smallest length E<gt>= C<$n> whose only prime factors are 2, 3 and 5 and
that PDL::FFT transforms correctly.

=item fft_ok($n)

True when PDL::FFT transforms length C<$n> correctly, i.e. its factorization
in F<fftn.c> has at most 11 factors. The count follows the C code step by
step.

=item fft_inplace($re, $im), ifft_inplace($re, $im)

Complex FFT along dim 0, in place, broadcasting over the other dims, with the
conventions of C<PDL::FFT::fft> and C<ifft>: forward unscaled, inverse scaled
by 1/n. C<$re> and C<$im> must have the same dims.

The transform uses PDL::FFT when C<fft_ok> holds for the length, and
L<PDL::FFTW3> otherwise. When neither can transform the length, the call
croaks with a message that names the length; it does not return a wrong
result.

=item fft_backend($n)

C<'PDL::FFT'> or C<'PDL::FFTW3'>: the library that transforms length C<$n>;
croaks as C<fft_inplace> does.

=back

=head2 $PDL::EEG::FFT::Dispatcher::BACKEND

C<'auto'> (default) as above; C<'pdl'> never uses PDL::FFTW3; C<'fftw3'> uses
PDL::FFTW3 for every length and croaks when it is not installed.

=cut
