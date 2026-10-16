package PDL::EEG::Linalg;

use strict;
use warnings;
use Exporter 'import';

use PDL;
use PDL::MatrixOps qw(eigens_sym);

our @EXPORT_OK = qw(eigh_sym linalg_backend);

# PDL::LinearAlgebra (LAPACK) when it is installed in version 0.436 or later;
# otherwise PDL::MatrixOps, with one warning per process.
our $MIN_LA_VERSION = '0.436';
our $HAVE_LAPACK;           # set below; tests may localise it
our $USE_LAPACK = 1;        # 0: always PDL::MatrixOps, without a warning
our $QUIET      = 0;        # 1: no warning when falling back
my  $LA_STATUS;             # why LAPACK is not used (for the warning)
my  $WARNED     = 0;

if (eval { require PDL::LinearAlgebra; 1 }) {
    if (eval { PDL::LinearAlgebra->VERSION($MIN_LA_VERSION); 1 }) {
        $HAVE_LAPACK = 1;
    }
    else {
        $HAVE_LAPACK = 0;
        $LA_STATUS = "PDL::LinearAlgebra $PDL::LinearAlgebra::VERSION is older than $MIN_LA_VERSION";
    }
}
else {
    $HAVE_LAPACK = 0;
    $LA_STATUS = 'PDL::LinearAlgebra is not installed';
}

sub _reset_warning { $WARNED = 0; return }

sub _use_lapack {
    return 1 if $USE_LAPACK && $HAVE_LAPACK;
    if ($USE_LAPACK && !$QUIET && !$WARNED) {
        $WARNED = 1;
        warn "PDL::EEG: " . ($LA_STATUS // 'PDL::LinearAlgebra is not available')
           . "; symmetric eigendecompositions use PDL::MatrixOps::eigens_sym"
           . " (slower, residual about 1e-12 instead of 1e-16)."
           . " Install PDL::LinearAlgebra $MIN_LA_VERSION or later to use LAPACK.\n";
    }
    return 0;
}

# which solver eigh_sym uses now
sub linalg_backend {
    return ($USE_LAPACK && $HAVE_LAPACK)
        ? "PDL::LinearAlgebra $PDL::LinearAlgebra::VERSION (LAPACK syevd)"
        : 'PDL::MatrixOps eigens_sym';
}

# eigh_sym($A) -> ($values, $vectors)
#   $A        symmetric (n,n)
#   $values   (n), ascending
#   $vectors  (n,n), one eigenvector per column: $A x $vectors = $vectors x stretcher($values)
sub eigh_sym {
    my ($A) = @_;
    my ($n0, $n1) = $A->dims;
    die "eigh_sym: need a square matrix, got (" . join(',', $A->dims) . ")\n"
        unless $A->ndims == 2 && $n0 == $n1;
    $A = $A->double;

    if (_use_lapack()) {
        my ($val, $vec, $info) = PDL::LinearAlgebra::msymeigen($A, 0, 1, 'syevd');
        my $bad = (ref $info ? $info->abs->max->sclr : abs $info);
        die "eigh_sym: LAPACK syevd failed (info = $info)\n" if $bad;
        return ($val, $vec);
    }

    my ($vec, $val) = eigens_sym($A);
    my $ord = $val->qsorti;
    return ($val->index($ord), $vec->dice_axis(0, $ord));
}

1;

__END__

=head1 NAME

PDL::EEG::Linalg - symmetric eigendecomposition for PDL::EEG, with or without PDL::LinearAlgebra

=head1 SYNOPSIS

    use PDL::EEG::Linalg qw(eigh_sym linalg_backend);

    my ($values, $vectors) = eigh_sym($C);   # values ascending, vectors in columns
    print linalg_backend(), "\n";

=head1 DESCRIPTION

C<eigh_sym> uses LAPACK C<syevd> through L<PDL::LinearAlgebra> (C<msymeigen>)
when PDL::LinearAlgebra 0.436 or later is installed, and
L<PDL::MatrixOps> C<eigens_sym> otherwise. Both give the same form: eigenvalues
in ascending order, one eigenvector per column, so that

    $A x $vectors == $vectors x stretcher($values)

When PDL::LinearAlgebra is missing or older than 0.436, the first call prints
one warning per process, saying why and that C<eigens_sym> is used.
C<$PDL::EEG::Linalg::QUIET = 1> suppresses the warning;
C<$PDL::EEG::Linalg::USE_LAPACK = 0> selects C<eigens_sym> without a warning.

Loading PDL::LinearAlgebra also makes PDL's C<x> use BLAS for the whole
process.

=head2 linalg_backend()

The solver C<eigh_sym> uses: C<"PDL::LinearAlgebra 0.436 (LAPACK syevd)"> or
C<"PDL::MatrixOps eigens_sym">.

=head1 SEE ALSO

L<PDL::EEG::ICA>, L<PDL::EEG::Inverse::MinimumNorm>, t/27_linalg.t

=cut
