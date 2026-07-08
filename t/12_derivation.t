use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::Derivation qw(bne derive rereference);

# ---------------------------------------------------------------------------
# Recorded data = true potential minus the acquisition reference Avr(C3,C4),
# exactly as Nihon Kohden stores it (SystemReference=C3,C4). If BNE is correct,
# re-referencing the RECORDED data must reproduce the TRUE balanced-non-cephalic
# signal, with Avr(C3,C4) cancelling out (weights sum to 1).
# ---------------------------------------------------------------------------
my @lab = qw(Fp1 Cz C3 C4 BN1 BN2 DC01 Trigger);
my %s   = (Fp1=>12, Cz=>8, C3=>4, C4=>-2, BN1=>3, BN2=>-5, DC01=>100, Trigger=>1);
my $sysref = ($s{C3} + $s{C4}) / 2;                 # Avr(C3,C4) = 1

# data[n_ch, n_samp], 2 identical samples; at(ch,samp) = s_ch - sysref
my @rows = map { my $c = $_; [ $s{$c} - $sysref, $s{$c} - $sysref ] } @lab;
my $data = pdl(\@rows)->xchg(0, 1)->sever;          # -> [n_ch=8, n_samp=2]
is($data->dim(0), 8, 'data has 8 channels');
is($data->dim(1), 2, 'data has 2 samples');

my $rec = { data => $data, fs => 1000, labels => [@lab],
            t_start => '2026-07-02 14:03:03' };

# ----- bne() ---------------------------------------------------------------
my $bn = bne($rec, prop => 0.6, suffix => '-BN');   # BN1=V, BN2=S

is_deeply($bn->{labels}, [qw(Fp1-BN Cz-BN C3-BN C4-BN DC01 Trigger)],
          'labels: EEG suffixed, DC/Trigger bare, BN1/BN2 dropped');
is($bn->{data}->dim(0), 6, '6 output channels (8 minus dropped BN1/BN2)');
is($bn->{data}->dim(1), 2, '2 samples preserved');

my $r = 0.6 * $s{BN1} + 0.4 * $s{BN2};              # true BNE reference = -0.2
ok(abs($bn->{data}->at(0,0) - ($s{Fp1} - $r)) < 1e-6, 'Fp1-BN = s_Fp1 - r  (C3,C4 cancels)');
ok(abs($bn->{data}->at(2,0) - ($s{C3}  - $r)) < 1e-6, 'C3-BN  = s_C3  - r');
ok(abs($bn->{data}->at(3,0) - ($s{C4}  - $r)) < 1e-6, 'C4-BN  = s_C4  - r');
ok(abs($bn->{data}->at(4,0) - ($s{DC01} - $sysref)) < 1e-6, 'DC01 passed through unchanged');
ok(abs($bn->{data}->at(5,0) - ($s{Trigger} - $sysref)) < 1e-6, 'Trigger passed through');
is($bn->{reference}, 'BNE', 'meta reference = BNE');
is($bn->{bne_prop}, 0.6,   'meta bne_prop = 0.6');
is($bn->{fs}, 1000,        'fs carried over');

# independence from the acquisition reference value
my @rows2 = map { my $c = $_; [ $s{$c} - 99.9, $s{$c} - 99.9 ] } @lab;
my $rec2  = { %$rec, data => pdl(\@rows2)->xchg(0,1)->sever };
my $bn2   = bne($rec2, prop => 0.6);
ok(abs($bn2->{data}->at(0,0) - ($s{Fp1} - $r)) < 1e-6,
   'BNE result is independent of the acquisition reference');

# drop_ref => 0 keeps the reference electrodes (themselves re-referenced)
my $bnk = bne($rec, prop => 0.6, suffix => '-BN', drop_ref => 0);
is($bnk->{data}->dim(0), 8, 'drop_ref=0 keeps all 8 channels');
ok((grep { $_ eq 'BN1-BN' } @{$bnk->{labels}}), 'BN1 retained and re-referenced');

# ----- general derive(): bipolar Fp1-Cz ------------------------------------
my $mini = { data => pdl([[10,10],[4,4],[1,1]])->xchg(0,1)->sever,
             fs => 1, labels => [qw(Fp1 Cz X)] };
my $bip = derive($mini, [[1,-1,0]], ['Fp1-Cz']);
is($bip->{data}->dim(0), 1, 'derive: one output row');
ok(abs($bip->{data}->at(0,0) - 6) < 1e-6, 'derive: Fp1-Cz = 10 - 4 = 6');

# ----- rereference() to a single channel -----------------------------------
my $rr = rereference($rec, 'Cz');
my $xCz = $s{Cz} - $sysref;
ok(abs($rr->{data}->at(0,0) - (($s{Fp1} - $sysref) - $xCz)) < 1e-6,
   'rereference to Cz: y_i = x_i - x_Cz');

done_testing();
