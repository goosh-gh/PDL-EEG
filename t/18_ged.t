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
# Random numbers come from PDL::EEG::SynthNoise with fixed seeds, so the
# data are the same on every run and every machine.

eval { require PDL::EEG::GED; PDL::EEG::GED->import(qw(ged_cov ged_operator apply_ged)); 1 }
    or plan skip_all => "PDL::EEG::GED not loadable (needs PDL::LinearAlgebra): $@";
use PDL::EEG::SynthNoise qw(seed_rng gaussian);

my ($nch,$nsrc,$nt) = (32,6,6000);
my $pi = 4*atan2(1,1);

# (d0,d1) of N(0, 1/12) from a fixed seed: the SD (0.289) of uniform(-0.5,0.5),
# which this test used before, so the blink-to-background ratio is unchanged
sub synth {
    my ($d0,$d1,$seed) = @_;
    seed_rng($seed);
    return gaussian($d0,$d1) * sqrt(1/12);
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
