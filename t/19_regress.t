use strict; use warnings;
use Test::More;
use PDL;
use FindBin; use lib "$FindBin::Bin/../lib";
use PDL::EEG::Regress qw(regress_out);

srand(7);
my $nt = 4000; my $nch = 6;
my $t  = sequence($nt) / 200;
# two artifact sources + independent brain per channel
my $a1 = sin($t * 2.1);                       # "vertical"
my $a2 = sin($t * 0.7 + 1)->rotate(13);       # "horizontal"
my $mix1 = pdl(1.0, 0.9, 0.5, 0.4, 0.1, 0.05);
my $mix2 = pdl(0.2, -0.3, 1.0, -0.9, 0.1, 0.0);
my $brain = grandom($nch, $nt) * 0.05;
my $X = $mix1->dummy(1) * $a1->dummy(0) + $mix2->dummy(1) * $a2->dummy(0) + $brain; # (nch,nt)

# --- remove both true sources: artifact should be gone, brain ~preserved ---
my ($clean, $removed) = regress_out($X, $a1, $a2);
sub corr { my($x,$y)=@_; my $a=$x-$x->avg; my $b=$y-$y->avg;
           (sum($a*$b))/(sqrt(sum($a*$a)*sum($b*$b))+1e-30) }
my $ok_art = 1;
for my $c (0..$nch-1) {
    my $r = $clean->slice("($c),")->flat;
    $ok_art &&= (abs(corr($r,$a1)) < 0.05 && abs(corr($r,$a2)) < 0.05);
}
ok($ok_art, "both artifact sources removed from every channel (|corr|<0.05)");

# brain residual should still correlate with the injected brain
my $b0 = $brain->slice("(4),")->flat;         # channel 4: mostly brain
ok(corr($clean->slice("(4),")->flat, $b0) > 0.9, "brain signal preserved where artifact weak");

# --- span-invariance: any basis of span(a1,a2) gives the same removal ---
my $s1 = $a1 + 0.4*$a2;
my $s2 = $a2 - 0.25*$a1;
my $clean2 = regress_out($X, $s1, $s2);
ok(max(abs($clean - $clean2)) < 1e-8, "removal depends only on the source span (label-free)");

# --- DC preserved ---
my $Xdc = $X + 3.0;
my $cdc = regress_out($Xdc, $a1, $a2);
ok(abs($cdc->slice("(0),")->avg - $Xdc->slice("(0),")->avg) < 1e-6, "channel DC level preserved");

# --- removed + cleaned == input ---
ok(max(abs(($clean + $removed) - $X)) < 1e-8, "cleaned + removed reconstructs input");

done_testing();
