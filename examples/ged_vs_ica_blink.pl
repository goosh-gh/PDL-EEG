use strict; use warnings;
use PDL; use PDL::MatrixOps;
use FindBin; use lib "$FindBin::Bin/../lib";
use PDL::EEG::ICA qw(ica_decompose identify_by_reference inv_sqrt_sym);

sub sd  { my $x=shift; my $x0=$x-$x->avg; sqrt(($x0*$x0)->avg) }
sub cor { my($a,$b)=@_; my $a0=$a-$a->avg; my $b0=$b-$b->avg;
          (($a0*$b0)->avg)/(sd($a)*sd($b)+1e-30) }
sub fro { my $x=shift; sqrt(($x*$x)->sum) }         # Frobenius norm

# ============ synthetic EEG with a KNOWN blink ==================
srand(11);
my $fs = 1000; my $N = 60*$fs; my $nch = 19;
my $t = sequence($N)/$fs;

# --- background: 6 independent AR(1) colored sources, smooth spatial maps ---
my $nsrc = 6;
my $BG = zeroes($N,$nch);
for my $s (0..$nsrc-1){
    my $e = grandom($N); my $a = 0.95;             # AR(1)
    my $x = $e->copy;
    my $xx = $x->list ? [$x->list] : [];
    for my $i (1..$N-1){ $xx->[$i] += $a*$xx->[$i-1]; }
    my $src = pdl($xx); $src -= $src->avg; $src /= sd($src);
    my $map = grandom($nch); $map -= $map->avg;    # random spatial pattern
    $BG += $src->dummy(1) * $map->dummy(0);
}
$BG -= $BG->average->dummy(0);                     # zero-mean per channel
$BG /= sd($BG->flat);                              # scale overall

# --- known blink: frontal topography x sparse variable-amplitude bumps ---
my $atopo = pdl(1.0,1.0,0.8,0.8, 0.5,0.5,0.35,0.35, 0.2,0.2,0.15,0.15,
                0.1,0.1,0.08,0.08,0.05,0.05,0.05);  # 19 ch, frontal-heavy
$atopo /= fro($atopo);
my $bump_w = 300;                                   # 300 ms
my $bump = 0.5*(1-cos(2*3.14159*sequence($bump_w)/($bump_w-1))); # raised cosine
my $btc = zeroes($N);                               # blink time course
my $nblink = 45;
my @onsets;
for (1..$nblink){
    my $on = int(rand($N-$bump_w));
    my $amp = exp(0.6*grandom(1)->sclr) * 6;        # lognormal amps -> big outliers
    $btc->slice($on.":".($on+$bump_w-1)) += $bump*$amp;
    push @onsets, $on;
}
$btc -= $btc->avg;
my $BLINK = $btc->dummy(1) * $atopo->dummy(0);      # (N,nch)

# scale blink so it's a strong but not total contributor
$BLINK *= 3.0 / sd($BLINK->flat);
my $noise = grandom($N,$nch)*0.15;
my $X = $BG + $BLINK + $noise;

print "== data ==  N=$N nch=$nch  blinks=$nblink\n";
printf "SNR-ish  ||BG||=%.1f  ||BLINK||=%.1f  ||noise||=%.1f\n\n",
       fro($BG),fro($BLINK),fro($noise);

# peri-blink vs blink-free masks (from known onsets)
my $peri = zeroes($N);
for my $on (@onsets){ my $lo=$on-50; $lo=0 if $lo<0; my $hi=$on+$bump_w+50; $hi=$N-1 if $hi>$N-1; $peri->slice($lo.":".$hi).=1; }
my $free = 1-$peri;
my $ip = which($peri>0); my $if = which($free>0);
printf "peri samples=%d  free samples=%d\n\n", $ip->nelem, $if->nelem;

sub chan_cov { my ($M)=@_; my $Mc=$M-$M->average->dummy(0); ($Mc x $Mc->transpose)/$M->dim(0) }

# ============ GED (R-whitening route, eigens_sym only) ==========
my $Xperi = $X->dice_axis(0,$ip);
my $Xfree = $X->dice_axis(0,$if);
my $S = chan_cov($Xperi);
my $R = chan_cov($Xfree);
# shrinkage on R for stability
my $lam = 0.05; my $trace = $R->diagonal(0,1)->avg;
$R = (1-$lam)*$R + $lam*$trace*identity($nch);
my $Rih = inv_sqrt_sym($R);
my $St  = $Rih x $S x $Rih;                         # symmetric
my ($V,$th) = eigens_sym($St);
my $ord = $th->qsorti->slice('-1:0');               # descending
$th = $th->index($ord); $V = $V->dice_axis(0,$ord);
my $Wged = $Rih x $V;                               # filters, columns (dim0=comp)
my $Aged = $R x $Wged;                              # patterns  A = R W
printf "GED generalized eigenvalues (blink/clean var ratio), top5: %s\n",
       join(", ", map{sprintf"%.1f",$_} $th->slice('0:4')->list);

# component time courses = W^T X  ; check which correlate with known blink
my $Xc = $X - $X->average->dummy(0);
my $ged_tc = $Wged->transpose x $Xc;                # (N, nch) comp courses
print "GED comp |corr| with true blink tc, top5: ",
      join(", ", map{sprintf"%.3f",abs(cor($btc,$ged_tc->slice(":,($_)")))} 0..4), "\n\n";

# linear removal operator for a set of GED comps, applied to arbitrary Y (centered)
sub ged_remove {
    my ($Y,$idx)=@_;
    my $Abad = $Aged->dice_axis(0,$idx);            # (nb,nch) patterns
    my $sc   = ($Wged->transpose x $Y)->dice_axis(1,$idx); # (N,nb) comp courses
    return $Abad x $sc;                             # (N,nch)
}

# ============ ICA (module) ======================================
my $ica = ica_decompose($X, nl=>"tanh", seed=>3, maxit=>500, tol=>1e-6);
printf "ICA converged iters=%d conv=%.1e\n", $ica->{iters},$ica->{conv};
my @rank = identify_by_reference($ica->{sources}, $btc);   # proxy for vEOG
print "ICA comp |corr| with true blink tc, top5: ",
      join(", ", map{sprintf"%.3f",abs($rank[$_][1])} 0..4), "\n\n";

sub ica_remove {
    my ($Y,$idx)=@_;                                # Y (N,nch) centered
    my $Mi = $ica->{mix}->dice_axis(0,pdl(long,@$idx));   # (nb,nch)
    my $Ui = $ica->{unmix}->dice_axis(1,pdl(long,@$idx)); # (nch,nb)
    my $sc = $Ui x $Y;                              # (N,nb)
    my $rm = $Mi x $sc;                             # (N,nch)
    return $rm;
}

# ============ metrics: apply FIXED operator to BLINK and BG =====
my $BLc = $BLINK - $BLINK->average->dummy(0);
my $BGc = $BG    - $BG->average->dummy(0);
sub report {
    my ($name,$rmref)=@_;
    my $rm_blink = &$rmref($BLc);
    my $rm_bg    = &$rmref($BGc);
    my $resid = fro($BLc - $rm_blink)/fro($BLc);    # blink left behind (want low)
    my $dist  = fro($rm_bg)/fro($BGc);              # clean signal removed (want low)
    printf "  %-14s blink_residual=%.3f   clean_distortion=%.3f\n",$name,$resid,$dist;
}

print "== removal comparison (ground-truth separated) ==\n";
report("GED top-1", sub { ged_remove($_[0], pdl(long,0)) });
report("GED top-2", sub { ged_remove($_[0], pdl(long,0,1)) });
report("GED top-3", sub { ged_remove($_[0], pdl(long,0,1,2)) });
my $b1=[$rank[0][0]];
my $b2=[$rank[0][0],$rank[1][0]];
report("ICA blink-1", sub { ica_remove($_[0],$b1) });
report("ICA blink-2", sub { ica_remove($_[0],$b2) });
