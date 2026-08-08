use strict;
use warnings;
use Test::More;
use PDL;

BEGIN { use_ok('PDL::EEG::TFA', qw(morlet_wavelet cwt_morlet tfr_morlet tfr_superlet tfr_stat apply_baseline)); }

my $PI = 4 * atan2(1, 1);

# ---- morlet_wavelet: unit energy, odd length, freq-scaling ----------------
{
    my ($wr, $wi) = morlet_wavelet(100, 1000, n_cycles => 7);
    ok($wr->nelem % 2 == 1, 'wavelet has odd length (centred)');
    my $en = sqrt(($wr**2 + $wi**2)->sum);
    ok(abs($en - 1) < 1e-9, 'wavelet normalised to unit L2 energy');
    my ($wr2) = morlet_wavelet(200, 1000, n_cycles => 7);
    ok($wr2->nelem < $wr->nelem, 'higher frequency -> shorter kernel');
}

# ---- cwt_morlet: a pure tone peaks at its own frequency -------------------
{
    my $sfreq = 1000; my $nt = 1000;
    my $t = sequence($nt) / $sfreq;
    my $freqs = pdl(10, 20, 40, 80, 120, 160);
    my $tone  = sin(2 * $PI * 80 * $t);
    my ($cr, $ci) = cwt_morlet($tone, $sfreq, $freqs, n_cycles => 7);
    is_deeply([$cr->dims], [$freqs->nelem, $nt], 'cwt output shape (nfreq,ntime)');
    my $pw   = $cr**2 + $ci**2;
    my $prof = $pw->slice(":,200:800")->xchg(0,1)->maximum;   # per-freq over central time
    is($freqs->at($prof->maximum_ind->sclr), 80, 'tone peaks at 80 Hz bin');
}

# ---- tfr_morlet: Gaussian burst localises in time and frequency -----------
{
    my $sfreq = 1000; my $nt = 1000;
    my $t = sequence($nt) / $sfreq;
    my $freqs = pdl(40, 80, 120, 160, 200, 300);
    my $burst = exp(-(($t-0.5)**2)/(2*0.01**2)) * sin(2*$PI*200*$t);
    my $pw = tfr_morlet($burst, $sfreq, $freqs, n_cycles => 7);
    is_deeply([$pw->dims], [$freqs->nelem, $nt], 'tfr power shape');
    my $per_freq = $pw->xchg(0,1)->maximum;
    my $fi = $per_freq->maximum_ind->sclr;
    is($freqs->at($fi), 200, 'burst peaks at 200 Hz');
    my $ti = $pw->slice("($fi),:")->maximum_ind->sclr;
    ok(abs($ti/$sfreq - 0.5) < 0.01, 'burst peaks near t = 0.5 s');
}

# ---- tfr_morlet ITC: phase-locked ~ 1, random ~ low -----------------------
{
    my $sfreq = 1000; my $nt = 1000; my $ne = 40;
    my $t = sequence($nt) / $sfreq;
    my $freqs = pdl(60, 120, 200);
    my $env = exp(-(($t-0.5)**2)/(2*0.02**2));
    srand(1);
    my $locked = zeroes($ne,$nt); my $rand = zeroes($ne,$nt);
    for my $e (0..$ne-1) {
        $locked->slice("($e),:") .= $env * sin(2*$PI*120*$t);
        $rand->slice("($e),:")   .= $env * sin(2*$PI*120*$t + rand()*2*$PI);
    }
    my $bin = which($freqs == 120)->at(0);
    my $il = tfr_morlet($locked,$sfreq,$freqs,n_cycles=>7,output=>'itc');
    my $ir = tfr_morlet($rand,  $sfreq,$freqs,n_cycles=>7,output=>'itc');
    ok($il->at($bin,500) > 0.9,  'ITC high for phase-locked trials');
    ok($ir->at($bin,500) < 0.5,  'ITC low for random-phase trials');
}

# ---- apply_baseline: zscore baseline mean ~ 0; all modes run --------------
{
    my $sfreq = 1000; my $nt = 1000;
    my $t = sequence($nt) / $sfreq;
    my $freqs = pdl(40, 80, 160);
    my $burst = exp(-(($t-0.5)**2)/(2*0.01**2)) * sin(2*$PI*80*$t);
    my $pw = tfr_morlet($burst, $sfreq, $freqs, n_cycles => 7);
    my $z = apply_baseline($pw, $t, [0.0, 0.2], mode => 'zscore');
    ok(abs($z->slice(":,0:200")->avg) < 1e-6, 'zscore baseline mean ~ 0');
    for my $m (qw(ratio logratio percent mean)) {
        my $b = eval { apply_baseline($pw, $t, [0.0, 0.2], mode => $m) };
        ok(defined $b && !$b->isbad->any, "baseline mode '$m' runs cleanly");
    }
}

# ---- tfr_superlet: sharp burst localises ----------------------------------
{
    my $sfreq = 2000; my $nt = 600;
    my $t = (sequence($nt) - 100) / $sfreq;
    my $burst = exp(-(($t-0.020)**2)/(2*0.003**2)) * sin(2*$PI*275*$t);
    my $freqs = 100 + 10 * sequence(40);
    my $ps = tfr_superlet($burst, $sfreq, $freqs,
                          base_cycles => 3, order_min => 1, order_max => 6);
    is_deeply([$ps->dims], [$freqs->nelem, $nt], 'superlet power shape');
    my $fi = $ps->xchg(0,1)->maximum->maximum_ind->sclr;
    ok(abs($freqs->at($fi) - 275) <= 10, 'superlet peaks near 275 Hz');
    my $ti = $ps->slice("($fi),:")->maximum_ind->sclr;
    ok(abs($t->at($ti) - 0.020) < 0.004, 'superlet peaks near t = 20 ms');
}

# ---- tfr_stat: across-trial reliability separates signal / consistent / noise
{
    my $sfreq = 2000; my $nt = 400; my $ne = 200;
    my $t = (sequence($nt) - 100) / $sfreq;             # -0.05 .. 0.15 s
    my $freqs = pdl(100, 150, 200, 275, 400);
    srand(3);
    my $ep = zeroes($ne, $nt);
    for my $e (0 .. $ne-1) {
        my $burst   = 3   * exp(-(($t-0.040)**2)/(2*0.004**2)) * sin(2*$PI*275*$t);
        my $prebump = 1.2 * exp(-(($t+0.012)**2)/(2*0.003**2)) * sin(2*$PI*275*$t);
        $ep->slice("($e),:") .= $burst + $prebump + grandom($nt)*1.5;
    }
    my %s = tfr_stat($ep, $sfreq, $freqs, $t, [-0.05,-0.004], n_cycles=>7);
    is_deeply([$s{z}->dims], [$freqs->nelem, $nt], 'tfr_stat z shape');
    is($s{n}, $ne, 'tfr_stat reports trial count');
    my $fi = which($freqs == 275)->at(0);
    my $ti = which(abs($t-0.040) < 1/$sfreq)->at(0);
    my $ni = which(abs($t-0.040) < 1/$sfreq)->at(0);
    ok($s{z}->at($fi,$ti) > 8, 'genuine burst is highly reliable (z>8)');
    ok($s{z}->at(which($freqs==150)->at(0),$ti) < 3, 'empty band is not reliable (z<3)');
}

# ---- evoked / induced decomposition: total = evoked + induced ------------
{
    my $sfreq = 2000; my $nt = 400; my $ne = 300;
    my $t = (sequence($nt) - 80) / $sfreq;
    my $freqs = pdl(100, 150, 200, 250, 300);
    srand(9);
    my $ep = zeroes($ne, $nt);
    for my $e (0 .. $ne-1) {
        my $env = exp(-(($t-0.040)**2)/(2*0.004**2));
        $ep->slice("($e),:") .= 1.0*$env*sin(2*$PI*150*$t)                 # phase-locked
                              + 1.0*$env*sin(2*$PI*250*$t + rand()*2*$PI)  # random phase
                              + grandom($nt)*1.0;
    }
    my $tot = tfr_morlet($ep,$sfreq,$freqs,n_cycles=>7,output=>'power');
    my $evo = tfr_morlet($ep,$sfreq,$freqs,n_cycles=>7,output=>'evoked');
    my $ind = tfr_morlet($ep,$sfreq,$freqs,n_cycles=>7,output=>'induced');
    my $ti = which(abs($t-0.040) < 1/$sfreq)->at(0);
    my $f150 = which($freqs==150)->at(0);
    my $f250 = which($freqs==250)->at(0);
    ok($evo->at($f150,$ti) > $ind->at($f150,$ti), '150 Hz phase-locked -> evoked > induced');
    ok($ind->at($f250,$ti) > $evo->at($f250,$ti), '250 Hz random-phase -> induced > evoked');
    my $diff = ($tot - ($evo + $ind))->abs->max;
    ok($diff < 1e-9, 'total power = evoked + induced');
}

done_testing();
