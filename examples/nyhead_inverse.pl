#!/usr/bin/env perl
# nyhead_inverse.pl - source localization on the New York Head leadfield.
#
#   Seed a known cortical source -> forward-project through the bundled
#   leadfield K -> recover it with sLORETA/eLORETA -> report the peak
#   vertex, its distance to the seed, and its Harvard-Oxford area.
#   --montage lets you add/remove electrodes and watch localization move.
#
# Real run (needs sa_nyhead.mat + PDL::IO::NYHead):
#   perl -I../lib examples/nyhead_inverse.pl --mat sa_nyhead.mat \
#        --electrodes 19 --method eloreta --seed-area "Precentral Gyrus" \
#        --snr 6 --out powers.dat
#
# Anywhere (synthetic, no .mat):
#   perl -I lib examples/nyhead_inverse.pl --demo --method eloreta
#
use strict; use warnings;
use PDL;
use Getopt::Long;
use PDL::EEG::Inverse::MinimumNorm qw(forward_project inverse_operator apply_inverse source_power);

my %o = (method=>'eloreta', ref=>'car', electrodes=>'19', snr=>0,
         reg_frac=>0.05, max_iter=>100, tol=>1e-8);
GetOptions(\%o,
  'mat=s','demo','method=s','ref=s','electrodes=s','seed=i','seed_area|seed-area=s',
  'snr=f','alpha=f','reg_frac|reg-frac=f','max_iter|max-iter=i','tol=f',
  'out=s','montage=s@','quiet') or die;

# ---- data source: real NYHead or synthetic demo -------------------------
my ($Kfull,$vc,$Ne_all,$Ns,$area_of,$elec_idx_of);  # Kfull: (Ns, Ne_all)
if ($o{demo}) {
    # synthetic cortex patch: Ns sources on a jittered grid, Ne electrodes above
    $Ns = 300; $Ne_all = 32;
    my $g = 20; my $gx = (sequence($Ns)%$g)/($g-1); my $gy = (sequence($Ns)/$g)->floor/($Ns/$g-1);
    $vc = ($gx->cat($gy, zeroes($Ns)))*80 - pdl(40,40,0)->dummy(0);   # (Ns,3) mm-ish
    my $ex = (sequence($Ne_all)%8)/7*80-40; my $ey=(sequence($Ne_all)/8)->floor/3*80-40;
    my $ez = ones($Ne_all)*40;
    my $ep = $ex->cat($ey,$ez);                             # (Ne_all,3) electrodes above sheet
    $Kfull = zeroes($Ns,$Ne_all);
    for my $e (0..$Ne_all-1){
        my $d = sqrt((($vc - $ep->slice("($e),:")->dummy(0))**2)->xchg(0,1)->sumover);  # (Ns)
        $Kfull->slice(":,($e)") .= 1.0/(1+($d/22)**2);      # smooth falloff
    }
    my $emean = $Kfull->xchg(0,1)->avgover;  $Kfull = $Kfull - $emean;  # per-source CAR
    $area_of     = sub { sprintf "region%d", int($_[0]/30) };
    $elec_idx_of = sub { my $n=shift; $n eq 'all' ? sequence($Ne_all)
                          : $n eq '19' ? pdl(0..18)               # first 19 as a stand-in
                          : pdl(split /,/,$n) };
} else {
    die "need --mat sa_nyhead.mat (or use --demo)\n" unless $o{mat};
    require PDL::IO::NYHead;
    my $ny = PDL::IO::NYHead->new($o{mat});
    # leadfield() is (electrode, source) = (231,74382); module wants (Ns,Ne)
    $Kfull = $ny->leadfield->transpose->sever;              # (74382, 231)
    ($Ns,$Ne_all) = $Kfull->dims;
    my $surf = $ny->surface('cortex75K');                   # {vc=>(Nv,3), tri=>...}
    $vc = ref($surf) eq 'HASH' ? $surf->{vc} : $surf;
    $area_of     = sub { eval { $ny->area_of_vertex($_[0]) } // '?' };
    $elec_idx_of = sub {
        my $n=shift;
        return sequence($Ne_all) if $n eq 'all';
        return $ny->idx19 if $n eq '19';
        return pdl(split /,/,$n);
    };
}

# ---- pick a seed vertex -------------------------------------------------
my $seed;
if (defined $o{seed}) { $seed = $o{seed}; }
elsif ($o{seed_area} && !$o{demo}) {
    # centroid vertex of the named HO area (nearest vertex to area centroid)
    my @in; for my $v (0..$Ns-1){ push @in,$v if $area_of->($v) =~ /\Q$o{seed_area}\E/i; }
    die "no vertices in area '$o{seed_area}'\n" unless @in;
    my $idx = pdl(@in); my $ctr = $vc->dice_axis(0,$idx)->avgover;   # (3) mean over vertices
    my $d = sqrt((($vc - $ctr->dummy(0))**2)->xchg(0,1)->sumover);
    $seed = $d->minimum_ind;
} else { $seed = int($Ns/2); }

# ---- build a montage list ----------------------------------------------
my @montages = $o{montage} && @{$o{montage}} ? @{$o{montage}} : ($o{electrodes});

printf "seed vertex %d  area '%s'  pos [%.1f %.1f %.1f]\n",
       $seed, $area_of->($seed), $vc->slice("($seed),:")->list unless $o{quiet};
printf "%-10s %5s  %-8s  %5s  %-s\n",'montage','#el','method','dist','peak area' unless $o{quiet};

my $last_pw;
for my $m (@montages) {
    my $eidx = $elec_idx_of->($m);
    my $K    = $Kfull->dice_axis(1,$eidx)->sever;           # (Ns, Ne_sub)
    $K = PDL::EEG::Inverse::MinimumNorm::avg_reference($K) if $o{ref} eq q{car};  # re-reference over the chosen electrodes
    my $Ne   = ($K->dims)[1];

    # forward: delta source at seed -> scalp topo, optional noise by SNR (dB)
    my $j = zeroes($Ns); $j->set($seed,1);
    my $b = forward_project($K,$j);                         # (Ne)
    if ($o{snr} > 0) {
        my $sp = ($b*$b)->avg; my $np = $sp/(10**($o{snr}/10));
        $b = $b + grandom($Ne)*sqrt($np);
    }

    my %iopt = (method=>$o{method}, ref=>$o{ref},
                max_iter=>$o{max_iter}, tol=>$o{tol});
    defined $o{alpha} ? ($iopt{alpha}=$o{alpha}) : ($iopt{reg_frac}=$o{reg_frac});
    my $op = inverse_operator($K,%iopt);
    my $pw = source_power($op,$b);                          # (Ns) >=0
    my $peak = $pw->maximum_ind;
    my $dist = sqrt((($vc->slice("($peak),:") - $vc->slice("($seed),:"))**2)->sum);
    my $conv = $op->{method} eq 'eloreta' ? sprintf(" [elor %diter]",$op->{iters}) : "";
    printf "%-10s %5d  %-8s  %5.1f  %s%s\n",
           $m,$Ne,$o{method},$dist,$area_of->($peak),$conv unless $o{quiet};
    $last_pw = $pw;
}

if ($o{out}) {
    $last_pw->wcols($o{out});      # one value per source (74382 rows) for GS3D colors
    print "wrote source power -> $o{out}  ($Ns rows)\n" unless $o{quiet};
    # GS3D: read these into scene->{colors} via your colormap, same path as pick.
}
