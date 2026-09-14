use strict; use warnings;
use PDL;
use Test::More;

# GED eye-blink removal (PDL::EEG::GED) validation on a synthetic ground truth.
#
# We build a known artifact and check that GED isolates and removes it:
#   * a fixed frontal scalp topography  a  (unit vector over nch channels)
#     carrying a time course  b           -> the "blink" (rank-1: a x b)
#   * a background of nsrc sources mixed by M, with M made ORTHOGONAL to a,
#     so the blink occupies a channel direction the background never uses.
# The blink is injected into a "blink block" (-> S) but not the "clean block"
# (-> R). A correct GED must then (1) be biorthogonal, (2) place the blink in a
# single dominant component, (3) recover a as that component's pattern, and
# (4) remove a x b while leaving the background untouched.
# Fully deterministic (no RNG / no seeding).

use_ok('PDL::EEG::GED');
use PDL::EEG::GED qw(ged_cov ged_operator apply_ged);

my ($nch,$nsrc,$nt) = (32,6,6000);
my $pi = 4*atan2(1,1);

# deterministic pseudo-uniform (-0.5,0.5) matrix (name avoids PDL's det())
sub synth {
    my ($d0,$d1,$seed) = @_;
    my $x = sin((sequence($d0*$d1)+1)*0.12345 + $seed*7.13) * 43758.5453;
    return ($x - floor($x) - 0.5)->reshape($d0,$d1);
}

my $ax = sequence($nch);
my $a  = exp(-(($ax-2)**2)/50); $a /= sqrt(sum($a**2));          # frontal blink topography

my $M = synth($nch,$nsrc,1);                                     # background mixing
$M = $M - $a->dummy(1) * ($a->dummy(1)*$M)->sumover->dummy(0);   # orthogonalise to a

my $t = sequence($nt);                                           # blink time course
my $b = zeroes($nt);
$b += sin(2*$pi*$_*$t/$nt + $_) for (1.5,2.7,4.1,6.3);
$b -= $b->avg; $b *= 8/sqrt(sum($b**2)/$nt);

my $Xb = (synth($nsrc,$nt,2) x $M) + $a->dummy(1)*$b->dummy(0);  # blink block -> S
my $Xc =  synth($nsrc,$nt,3) x $M;                              # clean block -> R
my $op = ged_operator(ged_cov($Xb), ged_cov($Xc), reg=>0.05);
my $k  = $nch-1;

my $err = sqrt(sum((($op->{W}->transpose x $op->{A})-identity($nch))**2));
ok($err < 1e-9, "biorthogonality ||W'A - I|| < 1e-9 (got ".sprintf('%.1e',$err).")");

my $gap = $op->{evals}->at($k)/$op->{evals}->at($k-1);
ok($gap > 20, "blink is one dominant component: lambda_max/lambda_2 > 20 (got ".sprintf('%.0f',$gap).")");

my $p = $op->{A}->slice("($k),")->flat->sever; $p = $p/sqrt(sum($p**2));
my $cos = abs(sum($p*$a));
ok($cos > 0.99, "top pattern recovers blink topography: cos > 0.99 (got ".sprintf('%.4f',$cos).")");

my $Xp  = $a->dummy(1)*$b->dummy(0);
my $resid = sqrt(sum(apply_ged($op,$Xp,n=>1)**2))/sqrt(sum($Xp**2));
ok($resid < 0.05, "pure blink removed to < 5% residual (got ".sprintf('%.1f%%',100*$resid).")");

my $ct = synth($nsrc,$nt,4) x $M; my $Xt = $ct + $a->dummy(1)*$b->dummy(0);
my $rec = sqrt(sum((apply_ged($op,$Xt,n=>1)-$ct)**2))/sqrt(sum($ct**2));
ok($rec < 0.05, "background recovered from blink+background: error < 5% (got ".sprintf('%.1f%%',100*$rec).")");

done_testing();
