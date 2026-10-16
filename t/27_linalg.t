use strict;
use warnings;
use Test::More;
use PDL;
use PDL::MatrixOps qw(stretcher identity);
use PDL::EEG::Linalg qw(eigh_sym linalg_backend);

# Ground truth: A = Q diag(lambda) Q^T, with Q a product of random Givens
# rotations (orthogonal by construction) and lambda chosen here, ascending.
sub givens_q {
    my ($n, $nrot) = @_;
    my $Q = identity($n)->copy;
    for (1 .. $nrot) {
        my $p = int(CORE::rand($n));
        my $q = int(CORE::rand($n - 1)); $q++ if $q >= $p;
        my $t = 2 * 3.14159265358979 * CORE::rand();
        my ($c, $s) = (cos($t), sin($t));
        my $G = identity($n)->copy;            # rotation in the (p, q) plane
        $G->set($p, $p, $c);  $G->set($q, $q, $c);
        $G->set($q, $p, -$s); $G->set($p, $q, $s);
        $Q = $Q x $G;
    }
    return $Q;
}

sub max_abs { return $_[0]->abs->max->sclr }

CORE::srand(27);
my @cases;
for my $n (5, 19, 64) {
    my $lam = 10 ** (sequence($n) / ($n - 1) * 6 - 3);    # 1e-3 .. 1e3
    $lam->set(0, -0.5);                                     # one negative value
    $lam = $lam->index($lam->qsorti);
    my $Q = givens_q($n, 4 * $n * $n);
    my $A = $Q x stretcher($lam) x $Q->transpose;
    $A = 0.5 * ($A + $A->transpose);                       # exactly symmetric
    push @cases, [$n, $A, $Q, $lam];
}

# check one result against the ground truth
sub check {
    my ($tag, $tol, $n, $A, $Q, $lam, $val, $vec) = @_;
    is_deeply([$val->dims], [$n], "$tag: values (n)");
    is_deeply([$vec->dims], [$n, $n], "$tag: vectors (n,n)");
    my $scale = $lam->abs->max->sclr;
    ok(max_abs($val - $lam) / $scale < $tol, "$tag: all values equal the known ones, ascending");
    my $resid = max_abs($A x $vec - $vec x stretcher($val)) / $scale;
    ok($resid < $tol, "$tag: A V = V diag(values) for every column (resid $resid)");
    ok(max_abs($vec->transpose x $vec - identity($n)) < $tol, "$tag: vectors orthonormal");
    # column j is +- column j of Q (all known values distinct)
    my $dot = ($Q->transpose x $vec)->diagonal(0, 1)->abs;
    ok(max_abs($dot - 1) < $tol * 10, "$tag: every vector matches the known one up to sign");
    # negative control: rows instead of columns breaks the eigen equation
    my $wrong = max_abs($A x $vec->transpose - $vec->transpose x stretcher($val)) / $scale;
    ok($wrong > 1e-3, "$tag: negative control, vectors read as rows fail ($wrong)");
}

# ---- PDL::MatrixOps path (always available) -------------------------------------
{
    local $PDL::EEG::Linalg::USE_LAPACK = 0;
    my @w; local $SIG{__WARN__} = sub { push @w, @_ };
    is(linalg_backend(), 'PDL::MatrixOps eigens_sym', 'USE_LAPACK = 0: backend is eigens_sym');
    for my $c (@cases) {
        my ($n, $A, $Q, $lam) = @$c;
        my $A0 = $A->copy;
        check("eigens_sym n=$n", 1e-9, $n, $A, $Q, $lam, eigh_sym($A));
        ok(($A == $A0)->all, "eigens_sym n=$n: input unchanged");
    }
    is(scalar @w, 0, 'USE_LAPACK = 0: no warning');
}

# ---- LAPACK path ------------------------------------------------------------------
SKIP: {
    skip 'PDL::LinearAlgebra 0.436 or later not installed', 3 * 8 + 2
        unless $PDL::EEG::Linalg::HAVE_LAPACK;
    like(linalg_backend(), qr/^PDL::LinearAlgebra .* \(LAPACK syevd\)$/, 'backend is LAPACK');
    for my $c (@cases) {
        my ($n, $A, $Q, $lam) = @$c;
        my $A0 = $A->copy;
        check("syevd n=$n", 1e-12, $n, $A, $Q, $lam, eigh_sym($A));
        ok(($A == $A0)->all, "syevd n=$n: input unchanged");
    }
    # both paths give the same values
    my ($n, $A) = @{ $cases[1] };
    my ($va) = eigh_sym($A);
    my ($vb) = do { local $PDL::EEG::Linalg::USE_LAPACK = 0; eigh_sym($A) };
    ok(max_abs($va - $vb) / $va->abs->max->sclr < 1e-9, 'syevd and eigens_sym give the same values');
}

# ---- falling back: one warning per process ------------------------------------------
{
    local $PDL::EEG::Linalg::HAVE_LAPACK = 0;
    PDL::EEG::Linalg::_reset_warning();
    my @w; local $SIG{__WARN__} = sub { push @w, @_ };
    my ($n, $A, $Q, $lam) = @{ $cases[0] };
    check('fallback n=5', 1e-9, $n, $A, $Q, $lam, eigh_sym($A));
    eigh_sym($A);
    is(scalar @w, 1, 'without PDL::LinearAlgebra: warned once for two calls');
    like($w[0] // '', qr/eigens_sym/, 'the warning names the solver used');
    is(linalg_backend(), 'PDL::MatrixOps eigens_sym', 'without PDL::LinearAlgebra: backend is eigens_sym');

    PDL::EEG::Linalg::_reset_warning();
    @w = ();
    { local $PDL::EEG::Linalg::QUIET = 1; eigh_sym($A); }
    is(scalar @w, 0, 'QUIET = 1: no warning');
    PDL::EEG::Linalg::_reset_warning();
}

eval { eigh_sym(zeroes(3, 4)) };
like($@, qr/square/, 'non-square input is refused');

done_testing;
