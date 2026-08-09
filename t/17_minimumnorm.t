use strict; use warnings;
use Test::More;
use PDL;
use PDL::EEG::Inverse::MinimumNorm qw(forward_project centering inverse_operator
                         apply_inverse source_estimate source_power);

# ---- deterministic, well-conditioned synthetic leadfield (Ne=8, Ns=20) ----
# sources & electrodes on interleaved 1-D lines; K_es = exp(-d^2/2s^2).
# Ne != Ns so any accidental transpose is caught.
my ($Ne,$Ns) = (8,20);
my $es = sequence($Ne)/($Ne-1);                 # electrode x in [0,1]
my $ss = sequence($Ns)/($Ns-1);                 # source   x in [0,1]
my $sig = 0.18;
# K math (Ne x Ns) -> PDL (Ns,Ne): entry (s,e)
my $K = zeroes($Ns,$Ne);
for my $e (0..$Ne-1){ for my $s (0..$Ns-1){
    my $d = $es->at($e)-$ss->at($s);
    $K->set($s,$e, exp(-$d*$d/(2*$sig*$sig)));
}}
ok( ($K->dims)[0]==$Ns && ($K->dims)[1]==$Ne, 'K stored (Ns,Ne)' );

# ---------- forward_project ----------
my $j0 = zeroes($Ns); $j0->set(7,1); $j0->set(13,-0.6);
my $b0 = forward_project($K,$j0);
is_deeply([$b0->dims],[$Ne], 'forward_project -> (Ne)');
# multi-time forward
my $J = zeroes($Ns,3); $J->slice(':,0').=$j0; $J->slice(':,1').=$j0*2;
my $B = forward_project($K,$J);
is_deeply([$B->dims],[$Ne,3], 'forward_project -> (Ne,Nt)');
ok( (($B->slice(':,1') - 2*$b0->dummy(1))->abs->max) < 1e-10, 'forward linear in j' );

# ---------- centering / ref ----------
my $H = centering($Ne);
ok( (($H x $H) - $H)->abs->max < 1e-12, 'centering idempotent (projection)' );

# ---------- transpose-safety ----------
eval { inverse_operator($K->transpose, method=>'mne') };
ok( $@ =~ /transposed/, 'transposed leadfield rejected' );

# ---------- MNE operator matches direct formula ----------
{
    my $al=1e-3;
    my $op = inverse_operator($K, method=>'mne', ref=>'none', alpha=>$al);
    my $Kt=$K->transpose; my $C=($K x $Kt)+$al*identity($Ne);
    my $Tdirect = $Kt x PDL::MatrixOps::inv($C);
    ok( ($op->{T}-$Tdirect)->abs->max < 1e-6, 'MNE operator == K^T (KK^T+aI)^-1' );
}

# ---------- sLORETA: exact single-source localization (each source) ----------
{
    my $op = inverse_operator($K, method=>'sloreta', ref=>'none', alpha=>1e-6);
    my $hit=0;
    for my $src (0..$Ns-1){
        my $j=zeroes($Ns); $j->set($src,1);
        my $b=forward_project($K,$j);
        my $pw=source_power($op,$b);            # standardized^2
        $hit++ if $pw->maximum_ind == $src;
    }
    is($hit,$Ns,"sLORETA exact single-source localization $hit/$Ns");
}

# ---------- eLORETA: exact single-source localization + convergence ----------
{
    my $op = inverse_operator($K, method=>'eloreta', ref=>'none',
                              alpha=>1e-6, max_iter=>200, tol=>1e-12);
    ok( $op->{rel} < 1e-10, "eLORETA converged (iters=$op->{iters}, rel=".sprintf('%.1e',$op->{rel}).")" );
    ok( ($op->{w} > 0)->all, 'eLORETA weights positive' );
    my $hit=0;
    for my $src (0..$Ns-1){
        my $j=zeroes($Ns); $j->set($src,1);
        my $b=forward_project($K,$j);
        my $pw=source_power($op,$b);
        $hit++ if $pw->maximum_ind == $src;
    }
    is($hit,$Ns,"eLORETA exact single-source localization $hit/$Ns");
}

# ---------- car reference: build once, apply to many time points ----------
{
    my $op = inverse_operator($K, method=>'sloreta', ref=>'car', alpha=>1e-4);
    my $J = zeroes($Ns,4);
    $J->set(3,0,1); $J->set(9,1,1); $J->set(15,2,1); $J->set(19,3,1);
    my $B = forward_project($K,$J);             # (Ne,4)
    my $S = apply_inverse($op,$B);              # (Ns,4)
    is_deeply([$S->dims],[$Ns,4], 'apply_inverse multi-time -> (Ns,Nt)');
    my $peaks = ($S*$S)->maximum_ind;             # per-time argmax over sources
    ok( ($peaks - pdl(3,9,15,19))->abs->max <= 1, 'sLORETA(car) tracks 4 single sources across time');
}

# ---------- source_estimate convenience ----------
{
    my $b=forward_project($K,$j0);
    my $s1 = source_estimate($K,$b,method=>'eloreta',ref=>'none',alpha=>1e-4);
    is_deeply([$s1->dims],[$Ns], 'source_estimate one-shot -> (Ns)');
}

done_testing();
