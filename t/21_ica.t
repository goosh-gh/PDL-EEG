use strict; use warnings;
use Test::More;
use PDL;
use FindBin; use lib "$FindBin::Bin/../lib";
use PDL::EEG::ICA qw(ica_decompose apply_ica);

sub sd  { my $x=shift; my $x0=$x-$x->avg; sqrt(($x0*$x0)->avg) }
sub cor { my($a,$b)=@_; my $a0=$a-$a->avg; my $b0=$b-$b->avg;
          (($a0*$b0)->avg)/(sd($a)*sd($b)+1e-30) }

srand(7);
my $N = 4000; my $k = 4;
my $t = sequence($N)/$N;
# 4 independent sources: laplacian(super-G), uniform(sub-G), sinusoid, sparse spikes
my $u  = random($N); my $sgn=($u<0.5)*2-1;
my $s1 = $sgn*log(1-2*($u-0.5)->abs+1e-12);
my $s2 = random($N)-0.5;
my $s3 = sin(2*3.14159*7*$t);
my $s4 = (random($N)<0.05)*grandom($N)*5;
my $Sd = zeroes($N,$k);
$Sd->slice(":,(0)").=$s1; $Sd->slice(":,(1)").=$s2;
$Sd->slice(":,(2)").=$s3; $Sd->slice(":,(3)").=$s4;
for my $i (0..$k-1){ my $c=$Sd->slice(":,($i)"); $c-=$c->avg; $c/=sd($c); }
my $A = grandom($k,$k);
my $X = $A x $Sd;                                  # observed (nt,nch)=(N,k)

my $ica = ica_decompose($X, nl=>'tanh', seed=>1);
ok($ica->{conv} < 1e-6, "FastICA converged (conv=".sprintf('%.2e',$ica->{conv}).")");
is($ica->{k}, $k, "kept all $k components");

# each true source recovered by some estimated component
my $Se = $ica->{sources};
for my $i (0..$k-1){
    my $best = 0;
    for my $j (0..$k-1){
        my $c = abs(cor($Sd->slice(":,($i)"), $Se->slice(":,($j)")));
        $best = $c if $c > $best;
    }
    ok($best > 0.99, sprintf("true source %d recovered (|corr|=%.4f)", $i, $best));
}

# removing no components is the identity
my $X0 = apply_ica($X, $ica, []);
my $err0 = ($X0-$X)->abs->max->sclr;
ok($err0 < 1e-6, "remove-none reconstruction is identity (max err=".sprintf('%.2e',$err0).")");

# mix x sources + mean reproduces X
my $Xrec = ($ica->{mix} x $ica->{sources}) + $ica->{mean}->dummy(0);
my $err1 = ($Xrec-$X)->abs->max->sclr;
ok($err1 < 1e-6, "mix x sources + mean reconstructs X (max err=".sprintf('%.2e',$err1).")");

done_testing();
