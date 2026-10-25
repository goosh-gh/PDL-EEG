#!/usr/bin/env perl
#
# nyhead_alpha_sim.pl — a pure posterior alpha rhythm (one sine, no noise, no fluctuation) on the
# New York Head, projected through the leadfield to the NYHead electrodes and, by interpolation of
# the electrode values, to own points (neck / chest electrode sites): how the alpha field reaches
# the neck, and what each candidate reference (Cz, Cv12, Cv45, ...) does to the scalp channels.
#
#   perl examples/nyhead_alpha_sim.pl --points neck_points.tsv --out alpha_sim/alpha
#   perl examples/nyhead_alpha_sim.pl --points neck_points.tsv --orient 0,-1,0.3
#   perl examples/nyhead_alpha_sim.pl --points neck_points.tsv --source patch --w-calc 0.5
#   perl examples/nyhead_alpha_sim.pl --selftest             # synthetic, no .mat needed
#
# Sources (--source):
#   dipole (default)  equivalent dipoles at --sites, one per site and hemisphere, at the cortical vertex
#                     nearest to the centroid of the site (free-orientation leadfield
#                     /sa/cortex75K/V_fem, about 0.4 GB); the right dipoles are the mirror images of the
#                     left ones (x -> -x).  With one site its orientation is swept over --ndir
#                     directions; with several sites the moments (orientation uniform, size uniform in
#                     0..1, per site) are drawn in batches of --batch until --target-kept configurations
#                     are kept or --max-draws are drawn.  Kept are the configurations whose map on the 19
#                     10-20 channels, against --crit-ref, looks like a posterior alpha rhythm:
#                       1. the largest amplitude at O1 or O2;
#                       2. O1, O2, P3, P4 and Pz of one sign;
#                       3. no increase towards the front along O1 P3 C3 F3 Fp1, O2 P4 C4 F4 Fp2,
#                          Pz Cz Fz, P7 T7 F7 and P8 T8 F8 (absolute values, strict).
#                     The values at the neck are summarised over the kept configurations, also for two
#                     halves of them (when the halves agree, enough were kept).  The EDF files are
#                     written for the kept configuration whose map is closest to the mean of the kept
#                     maps, or with one site for --orient X,Y,Z (left hemisphere; the right mirrored).
#   patch             cortical patches with the normal current spread evenly (V_fem_normal): the
#                     vertices of the --pos areas within --gap mm of each other (parieto-occipital
#                     sulcus, weight --w-pos) and the --calc area (calcarine sulcus, weight --w-calc).
#                     The two walls of a sulcus face each other, so the net field depends on how the
#                     patch is cut; kept for comparison.
#
#   --mat FILE       sa_nyhead.mat                                     (~/src/NYHead/sa_nyhead.mat)
#   --points TSV     own points: columns name, x, y, z[, kind] (MNI mm); with a kind column only
#                    the rows of kind 'point' (the --points-out file of PDL-IO-NYHead
#                    examples/nyhead_neck_electrodes.pl)
#   --point N=X,Y,Z  an own point (repeatable)
#   --out PREFIX     output: PREFIX_ref-<REF>.edf for every --refs, in the directory of PREFIX
#                    (it must exist)                                                      (alpha)
#   --refs LIST      references: AVG (the average of the 231 NYHead electrodes, as the leadfield
#                    is stored) or any channel name                (AVG,Cz,Cv12,Cv45,Nk2,BN1,Th12)
#   --far NAME       the channel farthest from the sources, used instead of an unknown zero: the
#                    values of the candidate references are given against it              (Th12)
#   --channels SET   1020 (the 19 10-20 electrodes, Oz, Iz, Exxz, Nk1..Nk4) or all (231)  (1020)
#   --source S       dipole | patch                                                     (dipole)
#   --sites LIST     ';'-separated: POS (the parieto-occipital patch below) or a Harvard-Oxford area
#                    regex, each hemisphere separately
#                    (POS;Intracalcarine Cortex;Occipital Pole;Lateral Occipital Cortex, inferior division;
#                     Lateral Occipital Cortex, superior division)
#   --crit-ref NAME  reference of the alpha-like test: a channel or own point; the average of the 19
#                    when it is not there                                                    (BN1)
#   --ndir N         orientations swept with one site                                       (2000)
#   --orient X,Y,Z   with one site: this orientation of the left dipole for the EDF files
#   --batch N, --target-kept N, --max-draws N, --seed N   drawing with several sites
#                                                                          (2000, 500, 100000, 1)
#   --freq HZ        alpha frequency                                                         (10)
#   --dur S          length                                                                   (30)
#   --fs HZ          sampling rate                                                           (500)
#   --amp-uv UV      amplitude of --amp-ch against the average reference                      (50)
#   --amp-ch NAME                                                                              (O1)
#   --pos REGEX      Harvard-Oxford areas meeting at the parieto-occipital sulcus, as A|B
#                                                                 (Cuneal Cortex|Precuneous Cortex)
#   --calc REGEX     the calcarine patch (--source patch)                  (Intracalcarine Cortex)
#   --gap MM                                                                                   (5)
#   --w-pos W, --w-calc W   patch weights (--source patch)                              (1, 0.5)
#   --phase-rl DEG   phase of the right hemisphere against the left                            (0)
#   --method K:N     interpolation of the own points (PDL::IO::NYHead::Interp)             (r3:all)
#   --strings S      auto | h5dump | lowlevel
#
# The leadfield holds cortical sources only and is stored against the average of the 231
# electrodes, so the true zero (a reference at infinity) is not known in this model; an alpha
# field that is spread widely moves the average too.  The model's neck ends below Nk1, and points
# below the lowest electrode are extrapolated (PDL-IO-NYHead examples/nyhead_interp_rows.pl gives
# the expected errors of the interpolated values).  The signal is not an EEG.

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use PDL;
use PDL::EEG::IO::EDF qw(write_edf read_edf);
use PDL::EEG::SynthNoise qw(seed_rng gaussian uniform);

my $PI = 3.14159265358979;
my %o = (mat => "$ENV{HOME}/src/NYHead/sa_nyhead.mat", points => undef, point => [], out => 'alpha',
         refs => 'AVG,Cz,Cv12,Cv45,Nk2,BN1,Th12', far => 'Th12', channels => '1020', source => 'dipole',
         ndir => 2000, orient => undef, sites => 'POS;Intracalcarine Cortex;Occipital Pole;Lateral Occipital Cortex, inferior division;Lateral Occipital Cortex, superior division',
         crit_ref => 'BN1', batch => 2000, target_kept => 500, max_draws => 100000, seed => 1, freq => 10, dur => 30, fs => 500, amp_uv => 50, amp_ch => 'O1',
         pos => 'Cuneal Cortex|Precuneous Cortex', calc => 'Intracalcarine Cortex', gap => 5, w_pos => 1,
         w_calc => 0.5, phase_rl => 0, method => 'r3:all', strings => 'auto', selftest => 0);
GetOptions('mat=s' => \$o{mat}, 'points=s' => \$o{points}, 'point=s@' => $o{point}, 'out=s' => \$o{out},
           'refs=s' => \$o{refs}, 'far=s' => \$o{far}, 'channels=s' => \$o{channels}, 'source=s' => \$o{source},
           'ndir=i' => \$o{ndir}, 'orient=s' => \$o{orient}, 'sites=s' => \$o{sites}, 'crit-ref=s' => \$o{crit_ref},
           'batch=i' => \$o{batch}, 'target-kept=i' => \$o{target_kept}, 'max-draws=i' => \$o{max_draws}, 'seed=i' => \$o{seed}, 'freq=f' => \$o{freq}, 'dur=f' => \$o{dur},
           'fs=i' => \$o{fs}, 'amp-uv=f' => \$o{amp_uv}, 'amp-ch=s' => \$o{amp_ch}, 'pos=s' => \$o{pos},
           'calc=s' => \$o{calc}, 'gap=f' => \$o{gap}, 'w-pos=f' => \$o{w_pos}, 'w-calc=f' => \$o{w_calc},
           'phase-rl=f' => \$o{phase_rl}, 'method=s' => \$o{method}, 'strings=s' => \$o{strings},
           'selftest' => \$o{selftest})
    or die "bad options (see the header of $0)\n";
if ($o{selftest}) { exit(selftest() ? 0 : 1) }
-f $o{mat} or die "not found: $o{mat} (give it with --mat)\n";
$o{source} =~ /^(dipole|patch)$/ or die "--source dipole | patch\n";

require PDL::IO::NYHead;
require PDL::IO::NYHead::Interp;
my $ny  = PDL::IO::NYHead->new($o{mat}, strings => $o{strings});
my @lab = @{ $ny->electrode_labels };
my $E   = $ny->electrode_pos->double;                           # (231,3)
my $vc  = $ny->surface('cortex75K')->{vc}->double;               # (74382,3)
my $hoi = $ny->atlas_index;                                      # (74382) 1-based, 0 = none
my $ho  = $ny->ho_labels;
my @i19 = $ny->idx19->list;
printf "%s: %d electrodes, %d cortical vertices, %d Harvard-Oxford areas; strings read by %s\n",
    $o{mat}, scalar @lab, $vc->dim(0), scalar @$ho, $ny->strings_method;

# ---- own points: interpolation weights (the interpolation is linear in the electrode values) ----
my @pts = read_points(\%o);
my ($W, %flag);
if (@pts) {
    my ($k, $n) = $o{method} =~ /^(r3?):(all|\d+)$/ or die "--method KERNEL:N, not '$o{method}'\n";
    my $P = pdl(map { [ @$_[1 .. 3] ] } @pts)->transpose->copy;  # (np,3)
    $W = PDL::IO::NYHead::Interp::interp_weights($E, $P, kernel => $k, ($n eq 'all' ? () : (k => $n)));   # (ne,np)
    my $D  = PDL::IO::NYHead::Interp::distances($E, $P);
    my $zlow = $E->slice(':,(2)')->min;
    print "\n-- own points (interpolated from the electrode values, $o{method})\n";
    for my $i (0 .. $#pts) {
        my $dp = $D->slice(":,($i)");
        my $j  = $dp->minimum_ind->sclr;
        my $z  = $P->at($i, 2);
        $flag{ $pts[$i][0] } = $z < $zlow ? 'extrapolated' : $dp->at($j) > 45 ? 'extrapolated (far)' : '';
        printf "  %-6s (%.1f %.1f %.1f)  nearest %s %.1f mm  %s\n", $pts[$i][0], $P->slice("($i),")->list, $lab[$j], $dp->at($j),
            $z < $zlow ? sprintf('extrapolated: %.0f mm below the lowest electrode', $zlow - $z) : $flag{ $pts[$i][0] };
    }
}
my @pname = map { $_->[0] } @pts;

# ---- sources -> maps (ne,2): left, right ---------------------------------------------
my @patch = patches(\%o, $vc, $hoi, $ho);                        # [name, hemi, vertices, weight]
my $map;
if ($o{source} eq 'patch') {
    my $L = $ny->leadfield;                                      # (231,74382)
    $map = zeroes(double, scalar @lab, 2);
    for my $p (@patch) {
        my ($name, $h, $v, $w) = @$p;
        $map->slice(":,($h)") += ($L->dice_axis(1, $v)->double * ($w / $v->nelem))->mv(1, 0)->sumover;   # over the vertices (dim 1)
    }
} else {
    my @sites = split /;/, $o{sites};
    my (@GL, @GR, @sname);
    print "\n-- equivalent dipoles (at the vertex nearest to the centroid of each site)\n";
    my $Vf = $ny->leadfield_full;                                # (231,74382,3)
    for my $site (@sites) {
        my @hv;
        if ($site eq 'POS') { @hv = map { my $h = $_; (grep { $_->[1] == $h && $_->[0] =~ /^parieto/ } @patch)[0][2] } 0, 1 }
        else {
            my ($v) = area_vertices($site, $hoi, $ho);
            my $x = $vc->slice(':,(0)');
            @hv = ($v->where($x->index($v) < 0), $v->where($x->index($v) > 0));
        }
        for my $h (0, 1) {
            my $hv = $hv[$h];
            $hv->nelem or die "site '$site': no vertex in the ", qw(left right)[$h], " hemisphere\n";
            my $c  = $vc->dice_axis(0, $hv)->average;            # (3): over the vertices (dim 0)
            my $d2 = zeroes(double, $hv->nelem);
            $d2 += ($vc->dice_axis(0, $hv)->slice(":,($_)") - $c->at($_))**2 for 0 .. 2;
            my $v  = $hv->at($d2->minimum_ind->sclr);
            printf "  %-5s %-45s vertex %5d (%6.1f %6.1f %6.1f), %s\n", qw(left right)[$h], $site, $v, $vc->slice("($v),")->list,
                $ny->area_of_vertex($v) // '?';
            push @{ $h ? \@GR : \@GL }, $Vf->slice(":,($v),:")->double->copy;   # (231,3)
        }
        push @sname, $site;
    }
    undef $Vf;
    my $ex = explore(\%o, \@lab, \@i19, \@GL, \@GR, \@sname, drawer(\%o, scalar @sites), $W, \@pname, \%flag);
    $map = cat($ex->{mapL}, $ex->{mapR})->copy;                  # (231,2)
}

my %ch = (map { $lab[$_] => $map->slice("($_),") } 0 .. $#lab);
if (@pts) {
    my $pm = PDL::IO::NYHead::Interp::interp_rows($map, $W);     # (np,2)
    $ch{ $pname[$_] } = $pm->slice("($_),") for 0 .. $#pname;
}
my @chan = $o{channels} eq 'all' ? @lab
         : ((map { $lab[$_] } @i19), grep { exists $ch{$_} } qw(Oz Iz Exxz Nk1 Nk2 Nk3 Nk4));
push @chan, @pname;
run(\%o, \@chan, \%ch, { AVG => zeroes(double, 2) }, \%flag);
exit 0;

# =====================================================================================
# drawer: a closure giving batches of left-hemisphere moments, an arrayref of (n,3) per site.
#   One site: the --ndir orientations of a Fibonacci sphere, once.  Several sites: random batches.
sub drawer {
    my ($opt, $ns) = @_;
    if ($ns == 1) {
        my $done = 0;
        return sub {
            return undef if $done++;
            my $nd = $opt->{ndir};
            my $k  = sequence(double, $nd);
            my $z  = 1 - 2 * ($k + 0.5) / $nd;
            my $r  = sqrt(1 - $z**2);
            my $ph = $k * $PI * (3 - sqrt(5));
            return [ cat($r * cos($ph), $r * sin($ph), $z) ];   # (nd,3): unit orientations
        };
    }
    seed_rng($opt->{seed});
    return sub {
        my $n = $opt->{batch};
        return [ map {
            my $d = gaussian($n, 3);
            $d / sqrt(($d**2)->mv(1, 0)->sumover)->dummy(1) * uniform($n)->dummy(1);   # direction x size
        } 1 .. $ns ];
    };
}

# alpha_like($V (19,n) against the test reference, \%j19 name -> row) -> (n) 1 = kept
sub alpha_like {
    my ($V, $j19) = @_;
    my $row = sub { $V->slice("($j19->{$_[0]}),") };
    my $amax = $V->abs->maximum_ind;                             # (n): largest over the 19 (dim 0)
    my $s1 = $row->('O1');
    my $keep = (($amax == $j19->{O1}) | ($amax == $j19->{O2}))->long;
    $keep *= ($s1 * $row->($_) > 0) for qw(O2 P3 P4 Pz);
    for my $chain ([qw(O1 P3 C3 F3 Fp1)], [qw(O2 P4 C4 F4 Fp2)], [qw(Pz Cz Fz)], [qw(P7 T7 F7)], [qw(P8 T8 F8)]) {
        $keep *= ($row->($chain->[$_])->abs <= $row->($chain->[ $_ - 1 ])->abs) for 1 .. $#$chain;
    }
    return $keep;
}

# explore: source configurations -> the alpha-like ones kept -> the neck values summarised.
#   $GL, $GR: arrayrefs of (ne,3) per site; $draw: from drawer(); $W (ne,np) or undef.
sub explore {
    my ($opt, $lab, $i19, $GL, $GR, $sname, $draw, $W, $pname, $flag) = @_;
    my %j19; @j19{ map { $lab->[$_] } @$i19 } = 0 .. $#$i19;
    for (qw(O1 O2 P3 P4 Pz C3 C4 Cz F3 F4 Fz Fp1 Fp2 P7 P8 T7 T8 F7 F8)) { exists $j19{$_} or die "no $_ among the 19 channels\n" }
    my %at; @at{@$lab} = 0 .. $#$lab;
    my %pi; @pi{@$pname} = 0 .. $#$pname;
    my $cref = $opt->{crit_ref} // '';
    my $zero_of = exists $at{$cref} ? sub { $_[0]->slice("($at{$cref}),") }
                : (defined $W && exists $pi{$cref}) ? sub { $_[1]->slice("($pi{$cref}),") }
                : do { print "  --crit-ref: no channel or point '$cref', the average of the 19 is used\n"; $cref = 'avg19';
                       sub { $_[0]->dice_axis(0, pdl(long, $i19))->average } };
    my $ns = @$GL;
    my (@kM, @kP, @kMom);
    my ($ndraw, $nkept) = (0, 0);
    while (defined(my $mom = $draw->())) {
        my $n = $mom->[0]->dim(0);
        # M(e,d) = sum_s sum_c G_s(e,c) m_s(d,c), dims (ne,n): (A x B)(col j, row i) = sum_k A(k,i) B(j,k),
        # so A = m^T (dims (3,n)) and B = G (dims (ne,3)); checked against the dipole formula in --selftest
        my $M = zeroes(double, scalar @$lab, $n);
        for my $s (0 .. $ns - 1) {
            my $mr = $mom->[$s] * pdl(-1, 1, 1)->dummy(0);       # mirrored for the right hemisphere
            $M += ($mom->[$s]->transpose x $GL->[$s]) + ($mr->transpose x $GR->[$s]);
        }
        my $P = defined $W ? PDL::IO::NYHead::Interp::interp_rows($M, $W) : undef;   # (np,n)
        my $V = $M->dice_axis(0, pdl(long, $i19)) - $zero_of->($M, $P)->dummy(0);    # (19,n)
        my $ki = which(alpha_like($V, \%j19));
        $ndraw += $n;
        if ($ki->nelem) {
            push @kM, $M->dice_axis(1, $ki);
            push @kP, $P->dice_axis(1, $ki) if defined $P;
            push @kMom, [ map { $_->dice_axis(0, $ki) } @$mom ];
            $nkept += $ki->nelem;
        }
        last if $nkept >= $opt->{target_kept} || $ndraw >= $opt->{max_draws};
    }
    printf "\n-- source configurations (%s): %d drawn, %d look like a posterior alpha rhythm against %s (%.1f %%)\n",
        join(', ', @$sname), $ndraw, $nkept, $cref, 100 * $nkept / ($ndraw || 1);
    $nkept or die "no configuration gives an alpha-like map; more sites or other sites are needed\n";
    my $glue = sub { my $g = shift; $g = $g->glue(1, $_) for @_; $g };
    my $Mk = $glue->(@kM);                                        # (ne,nk)
    my $Pk = @kP ? $glue->(@kP) : undef;                          # (np,nk)
    my @momk = map { my $s = $_; my $m = $kMom[0][$s]; $m = $m->glue(0, $_->[$s]) for @kMom[1 .. $#kMom]; $m } 0 .. $ns - 1;   # (nk,3)
    my %val = map { $lab->[$_] => $Mk->slice("($_),") } 0 .. $#$lab;
    if (defined $Pk) { $val{ $pname->[$_] } = $Pk->slice("($_),") for 0 .. $#$pname }

    my $far = $opt->{far};
    my $zero = exists $val{$far} ? $val{$far} : do { print "  --far: no channel $far, the average reference is used\n"; $far = 'AVG'; zeroes(double, $nkept) };
    my $o1 = $val{O1} - $zero;                                   # (nk)
    my $half = sequence(long, $nkept) % 2;                       # two halves of the kept ones
    my $q = sub {
        my $x = shift;
        my $pc = sub { my ($y, $f) = @_; my $sv = $y->qsort; 100 * $sv->at(int($f * ($sv->nelem - 1) + 0.5)) };
        my $str = sprintf '%+6.1f %% (%+6.1f .. %+6.1f)', $pc->($x, 0.5), $pc->($x, 0.1), $pc->($x, 0.9);
        $str .= sprintf '   halves %+6.1f / %+6.1f', $pc->($x->where($half == 0), 0.5), $pc->($x->where($half == 1), 0.5) if $nkept >= 20;
        return $str;
    };
    printf "  over the kept configurations: median (10 .. 90 %%), as a share of O1 - %s\n", $far;
    print  "  the channel itself against $far (how much alpha a reference there carries):\n";
    for my $c (grep { exists $val{$_} && $_ ne $far } qw(Fp1 Fz Cz C3 T7 Pz Oz Iz Exxz Nk1 Nk2 Nk3 Nk4), @$pname) {
        printf "    %-6s %s%s\n", $c, $q->(($val{$c} - $zero) / $o1), $flag->{$c} ? "   ($flag->{$c})" : '';
    }
    print "  pairs (the difference of two candidate references, as a share of O1 - $far):\n";
    my @cand = grep { exists $val{$_} } (@$pname, 'Nk2', 'Cz');
    for my $ia (0 .. $#cand) {
        for my $ib ($ia + 1 .. $#cand) {
            my ($ca, $cb) = @cand[$ia, $ib];
            next if $ca eq $far || $cb eq $far;
            printf "    %-14s %s\n", "$ca - $cb", $q->(($val{$ca} - $val{$cb}) / $o1);
        }
    }
    if ($ns > 1) {                                               # moment sizes per site
        my $sz = cat(map { sqrt(($_**2)->mv(1, 0)->sumover) } @momk);   # (nk,ns)
        my $sh = $sz / $sz->mv(1, 0)->sumover->dummy(1);
        print "  share of the dipole moment per site (kept configurations; deeper sites need more moment for the same scalp field):\n";
        for my $si (0 .. $ns - 1) {
            my $sv = $sh->slice(":,($si)")->qsort;
            printf "    %-45s median %4.2f (%4.2f .. %4.2f)\n", $sname->[$si], map { $sv->at(int($_ * ($nkept - 1) + 0.5)) } 0.5, 0.1, 0.9;
        }
    }
    # representative: the kept map (on the 19 against the test reference, unit RMS) closest to the mean of the kept ones
    my $Vk = $Mk->dice_axis(0, pdl(long, $i19)) - $zero_of->($Mk, $Pk)->dummy(0);
    my $s1 = $Vk->slice("($j19{O1}),");
    $Vk = $Vk / sqrt(($Vk**2)->average)->dummy(0) * (2 * ($s1 > 0) - 1)->dummy(0);   # one sign: O1 positive
    my $mean = $Vk->mv(1, 0)->average;                           # (19): over the configurations (dim 1)
    my $best = ((($Vk - $mean->dummy(1))**2)->average)->minimum_ind->sclr;
    my @mb = map { $_->slice("($best),")->copy } @momk;          # (3) per site
    if ($ns == 1 && defined $opt->{orient}) { $mb[0] = pdl(split /,/, $opt->{orient}); $mb[0] /= sqrt(($mb[0]**2)->sum) }
    print "  EDF configuration (left dipoles; the right ones mirrored)", ($ns == 1 && defined $opt->{orient}) ? ' (--orient)' : ' (the kept map closest to their mean)', ":\n";
    printf "    %-45s moment (%+.2f %+.2f %+.2f)\n", $sname->[$_], $mb[$_]->list for 0 .. $ns - 1;
    my $mapL = zeroes(double, scalar @$lab);
    my $mapR = zeroes(double, scalar @$lab);
    for my $s (0 .. $ns - 1) {
        $mapL += ($mb[$s]->dummy(0)->transpose x $GL->[$s])->slice(':,(0)');                    # (3,1) x (ne,3) -> (ne)
        $mapR += (($mb[$s] * pdl(-1, 1, 1))->dummy(0)->transpose x $GR->[$s])->slice(':,(0)');
    }
    return { mapL => $mapL, mapR => $mapR, val => \%val, far => $far, nkept => $nkept, ndraw => $ndraw, mom => \@momk, crit => $cref };
}

# run: the alpha signal in every reference, the table and the EDF files
sub run {
    my ($opt, $chan, $ch, $refmap, $flag) = @_;
    exists $ch->{ $opt->{amp_ch} } or die "--amp-ch: no channel $opt->{amp_ch}\n";
    my $ph = $opt->{phase_rl} * $PI / 180;
    my $amp = sub {                                              # peak amplitude of l sin(wt) + r sin(wt + ph)
        my $m = shift;
        my ($l, $r) = ($m->at(0), $m->at(1));
        return sqrt(($l + $r * cos($ph))**2 + ($r * sin($ph))**2) * ($l + $r * cos($ph) < 0 ? -1 : 1);
    };
    my $scale = $opt->{amp_uv} / abs($amp->($ch->{ $opt->{amp_ch} }));   # uV per unit of the map
    my $nt = int($opt->{dur} * $opt->{fs} + 0.5);
    my $t  = sequence(double, $nt) / $opt->{fs};
    my $wt = 2 * $PI * $opt->{freq} * $t;
    my $S  = cat(sin($wt), sin($wt + $ph));                        # (nt,2): left, right time courses
    printf "\nalpha %.1f Hz, %.0f s at %d Hz; %s %.0f uV peak against the average reference; right hemisphere phase %+.0f deg\n",
        $opt->{freq}, $opt->{dur}, $opt->{fs}, $opt->{amp_ch}, $opt->{amp_uv}, $opt->{phase_rl};

    my @refs = grep { exists $refmap->{$_} || exists $ch->{$_} or do { print "  --refs: no channel $_, left out\n"; 0 } } split /,/, $opt->{refs};
    my %seen;
    my @show = grep { exists $ch->{$_} && !$seen{$_}++ } (qw(Fp1 Fp2 Fz F3 F4 Cz C3 C4 T7 T8 Pz P3 P4 O1 O2 Oz Iz Exxz Nk1 Nk2 Nk3 Nk4), @$chan);
    print "\n-- peak alpha amplitude (uV; the sign is the phase against the left-hemisphere source) in each reference\n";
    printf "  %-6s%s\n", '', join('', map { sprintf '%9s', $_ } @refs);
    for my $c (@show) {
        printf "  %-6s%s%s\n", $c, join('', map { my $r = $refmap->{$_} // $ch->{$_}; $_ eq $c ? sprintf('%9s', '-') : sprintf('%9.1f', $scale * $amp->($ch->{$c} - $r)) } @refs),
            $flag->{$c} ? "   ($flag->{$c})" : '';
    }
    my $far = exists $ch->{ $opt->{far} } ? $opt->{far} : undef;
    print "\n-- per reference: frontal / occipital in that reference",
          ($far ? "; the reference site against $far as a share of O1 - $far" : ''), "\n";
    for my $rn (@refs) {
        my $r = $refmap->{$rn} // $ch->{$rn};
        my $ratio = sub { my ($x, $y) = @_; my $dy = $amp->($ch->{$y} - $r); abs($dy) > 1e-12 ? $amp->($ch->{$x} - $r) / $dy : 0 };
        my $site = $far ? sprintf('   site-%s %+6.1f %%', $far, 100 * $amp->($r - $ch->{$far}) / $amp->($ch->{O1} - $ch->{$far})) : '';
        printf "  %-5s Fz/O1 %+6.2f   Fp1/O1 %+6.2f   Cz/O1 %+6.2f   Pz/O1 %+6.2f%s%s\n", $rn,
            (exists $ch->{Fz} && $rn ne 'O1') ? map { $ratio->($_, 'O1') } qw(Fz Fp1 Cz Pz) : (0, 0, 0, 0),
            $site, $flag->{$rn} ? "   ($flag->{$rn})" : '';
    }

    # ---- EDF, one per reference
    my @files;
    for my $rn (@refs) {
        my $r  = $refmap->{$rn} // $ch->{$rn};
        my @cl = grep { $_ ne $rn } @$chan;
        my $M  = cat(map { $ch->{$_} - $r } @cl) * $scale;         # dims (2,nc): the matrix nc x 2
        my $X  = ($M x $S)->float;                                 # (nc x 2)(2 x nt): dims (nt,nc)
        my $rec = { data => $X->transpose->copy, fs => $opt->{fs}, labels => \@cl };   # (nc,nt)
        (my $rs = $rn) =~ s{[^\w.-]}{_}g;
        my $f = "$opt->{out}_ref-$rs.edf";
        write_edf($rec, $f, plus => 1, unit => 'uV', start => [2026, 1, 1, 0, 0, 0], recording => "NYHead_alpha_sim_ref_$rs");
        push @files, $f;
    }
    print "\nwrote ", join(', ', @files), "\n",
          "The signal is one sine per hemisphere with no noise and no fluctuation: it is not an EEG, it shows how the\n",
          "alpha field of the chosen sources spreads to each channel and what each reference does to it.\n";
    return \@files;
}

sub patches {
    my ($opt, $vc, $hoi, $ho) = @_;
    my $area = sub { area_vertices($_[0], $hoi, $ho) };
    my $hname = sub {                                            # the names of one hemisphere
        my ($names, $h) = @_;
        my @n = grep { !/^(Left|Right)\b/ || ($h ? /^Right\b/ : /^Left\b/) } @$names;
        return join(' + ', @n ? @n : @$names);
    };
    my ($ra, $rb) = split /\|/, $opt->{pos};
    defined $rb or die "--pos needs two areas, A|B\n";
    my ($va, $na) = $area->($ra);
    my ($vb, $nb) = $area->($rb);
    my ($vcal, $ncal) = $area->($opt->{calc});
    my $x = $vc->slice(':,(0)');
    my @p;
    for my $h (0, 1) {                                           # 0 = left (x < 0), 1 = right
        my $side = sub { my $v = shift; $v->where($h ? $x->index($v) > 0 : $x->index($v) < 0) };
        my ($ha, $hb) = ($side->($va), $side->($vb));
        my $D = dist2($vc->dice_axis(0, $ha), $vc->dice_axis(0, $hb)); # (na,nb) squared mm
        my $g2 = $opt->{gap}**2;
        # A vertices near B: minimum over B (dim 1); B vertices near A: minimum over A (dim 0)
        my $pos = $ha->where($D->mv(1, 0)->minimum < $g2)->append($hb->where($D->minimum < $g2));
        my $cal = $side->($vcal);
        $pos->nelem or die "the areas /$ra/ and /$rb/ do not meet within $opt->{gap} mm in the ", qw(left right)[$h], " hemisphere\n";
        push @p, [ 'parieto-occipital (' . $hname->($na, $h) . ' | ' . $hname->($nb, $h) . ')', $h, $pos, $opt->{w_pos} ],
                 [ 'calcarine (' . $hname->($ncal, $h) . ')', $h, $cal, $opt->{w_calc} ];
    }
    print "\n-- patches", $opt->{source} eq 'patch' ? ' (normal current, spread evenly over the vertices)' : ' (to place the dipoles)', "\n";
    for (@p) {
        next if $opt->{source} ne 'patch' && $_->[0] !~ /^parieto/;
        printf "  %-5s %-55s %5d vertices, centroid (%.1f %.1f %.1f)%s\n", qw(left right)[ $_->[1] ], $_->[0], $_->[2]->nelem,
            $vc->dice_axis(0, $_->[2])->average->list, $opt->{source} eq 'patch' ? ", weight $_->[3]" : '';
    }
    return @p;
}

sub area_vertices {                                              # vertices of the areas matching a regex, and their names
    my ($re, $hoi, $ho) = @_;
    my @id = grep { $ho->[$_] =~ /$re/i } 0 .. $#$ho;
    @id or die "no Harvard-Oxford area matches /$re/\n";
    my $m = zeroes(byte, $hoi->nelem);
    $m |= ($hoi == $_ + 1) for @id;
    return (which($m), [ map { $ho->[$_] } @id ]);
}

sub dist2 {                                                      # (na,3), (nb,3) -> (na,nb), axis by axis
    my ($A, $B) = @_;
    my $d = zeroes(double, $A->dim(0), $B->dim(0));
    $d += ($A->slice(":,($_)")->dummy(1) - $B->slice(":,($_)")->dummy(0))**2 for 0 .. 2;
    return $d;
}

sub read_points {
    my ($opt) = @_;
    my @p;
    for my $s (@{ $opt->{point} }) {
        my ($n, $xyz) = $s =~ /^([^=]+)=(.+)$/ or die "--point NAME=X,Y,Z, not '$s'\n";
        my @c = split /,/, $xyz; @c == 3 or die "--point $n: three coordinates\n";
        push @p, [ $n, @c ];
    }
    if (defined $opt->{points}) {
        open my $fh, '<', $opt->{points} or die "$opt->{points}: $!\n";
        my @col;
        while (my $l = <$fh>) {
            chomp $l;
            next if $l =~ /^\s*(#|$)/;
            my @f = split /\t/, $l;
            if (!@col) { @col = map { lc } @f; next }
            my %r; @r{@col} = @f;
            next if defined $r{kind} && $r{kind} ne 'point';
            push @p, [ @r{qw(name x y z)} ];
        }
    }
    return @p;
}

# =====================================================================================
sub selftest {
    my $ok = 1;
    my $chk = sub { my ($c, $m) = @_; printf "  [%s] %s\n", $c ? 'PASS' : 'FAIL', $m; $ok &&= $c };
    my $quiet = sub { my $code = shift; my $out = ''; open my $fh, '>', \$out or die; my $old = select $fh; my @r = $code->(); select $old; close $fh; return ($out, @r) };

    # ---- run(): maps given directly, EDF written and read back
    my %ch = (Fp1 => pdl(-0.2, -0.2), Fz => pdl(-0.3, -0.3), Cz => pdl(-0.1, -0.1), O1 => pdl(1.0, 0.4),
              O2 => pdl(0.4, 1.0), Pz => pdl(0.5, 0.5), Cv12 => pdl(0.6, 0.6), Th12 => pdl(0.1, 0.1));
    require File::Temp;
    my $dir = File::Temp::tempdir(CLEANUP => 1);
    my %opt = (%o, out => "$dir/t", refs => 'AVG,Cz,Cv12', far => 'Th12', dur => 2, fs => 250, amp_uv => 50, phase_rl => 0);
    my ($out, $files) = $quiet->(sub { run(\%opt, [qw(Fp1 Fz Cz O1 O2 Pz Cv12 Th12)], \%ch, { AVG => zeroes(double, 2) }, {}) });
    print $out;
    $chk->((@$files == 3 && !grep { !-s $_ } @$files), 'three EDF files written');
    my $rec = read_edf("$dir/t_ref-Cv12.edf");
    my %ix; @ix{ @{ $rec->{labels} } } = 0 .. $#{ $rec->{labels} };
    $chk->(!exists $ix{Cv12} && exists $ix{O1}, 'the reference channel is left out of its own file');
    my $o1 = $rec->{data}->slice("($ix{O1}),");
    my $want = 50 * (1.4 - 1.2) / 1.4;                             # O1 - Cv12, both hemispheres in phase
    $chk->(abs($o1->abs->max - $want) < 0.02 * $want + 0.05, sprintf 'O1 - Cv12 peak %.2f uV (expected %.2f)', $o1->abs->max, $want);
    my $ra = read_edf("$dir/t_ref-AVG.edf");
    my %ia; @ia{ @{ $ra->{labels} } } = 0 .. $#{ $ra->{labels} };
    my $d = ($ra->{data}->slice("($ia{O1}),") - $ra->{data}->slice("($ia{Cv12}),") - $o1)->abs->max;
    $chk->($d < 0.05, sprintf 'AVG file: O1 - Cv12 equals the Cv12 file (max diff %.3f uV)', $d);
    $chk->(abs($ra->{data}->slice("($ia{O1}),")->abs->max - 50) < 0.1, 'O1 is 50 uV peak against the average');
    # Cv12 - Th12 = 1.0, O1 - Th12 = 1.2 -> 83.3 %
    $chk->($out =~ /Cv12\s+Fz\/O1.*site-Th12\s+\+83\.3 %/, 'reference site against --far: Cv12 - Th12 = 83.3 % of O1 - Th12');

    # ---- patches(): two areas meeting at y = 0 on a 2 mm grid, and a third area elsewhere
    my $gx = (xvals(10, 10) * 2 - 9)->flat; my $gy = (yvals(10, 10) * 2 - 9)->flat;
    my $vcs = cat($gx, $gy, zeroes(100))->append(pdl([[-5, 30, 0], [5, 30, 0]])->transpose)->copy;   # (100,3) + (2,3)
    my $hoi = (1 + 2 * ($gy > 0)->long)->append(pdl(long, 2, 2));   # 1 = A (y < 0), 3 = B (y > 0), 2 = C
    my (undef, @pp) = $quiet->(sub { patches({ %o, source => 'patch', pos => '^Aaa$|^Bbb$', calc => '^Ccc$', gap => 3 }, $vcs, $hoi, [ 'Aaa', 'Ccc', 'Bbb' ]) });
    my $sel = $pp[0][2];
    $chk->($sel->nelem == 10 && all(abs($gy->index($sel)) == 1) && all($gx->index($sel) < 0),
           sprintf 'left parieto-occipital patch: the 10 grid vertices next to the border (got %d)', $sel->nelem);
    $chk->($pp[1][2]->nelem == 1 && $pp[1][2]->at(0) == 100 && $pp[2][1] == 1, 'calcarine patch and hemispheres split by x');

    # ---- explore(): dipoles in a sphere (infinite medium, V = m.(p - q)/|p - q|^3), electrodes named
    #      at 10-20-like positions on it; maps checked against the formula, and the kept orientations
    my @nm = qw(Fp1 Fp2 F3 F4 C3 C4 P3 P4 O1 O2 F7 F8 T7 T8 P7 P8 Fz Cz Pz Nk2 Th12);
    my %xyz = (Fp1 => [-25, 80, 20], Fp2 => [25, 80, 20], F3 => [-40, 50, 60], F4 => [40, 50, 60], C3 => [-55, 0, 70],
               C4 => [55, 0, 70], P3 => [-40, -55, 60], P4 => [40, -55, 60], O1 => [-25, -85, 15], O2 => [25, -85, 15],
               F7 => [-65, 50, 10], F8 => [65, 50, 10], T7 => [-85, 0, 10], T8 => [85, 0, 10], P7 => [-65, -55, 10],
               P8 => [65, -55, 10], Fz => [0, 55, 70], Cz => [0, 0, 90], Pz => [0, -55, 70], Nk2 => [0, -60, -80], Th12 => [0, -55, -140]);
    my $Es = pdl(map { my $v = pdl($xyz{$_}); [ ($v * 90 / sqrt(($v**2)->sum))->list ] } @nm)->transpose->copy;   # (21,3) on the sphere
    $Es->slice('(20),') .= pdl(0, -55, -140);                    # Th12 left far below
    my $Gof = sub {                                              # (ne,3): the field of unit dipoles along x, y, z at q
        my $q = pdl(@_);
        my $dd = $Es - $q->dummy(0);
        my $r3 = (($dd**2)->mv(1, 0)->sumover)**1.5;              # (ne): over x, y, z (dim 1)
        return $dd / $r3->dummy(1);
    };
    require PDL::IO::NYHead::Interp;                            # an own point, to go through the interpolation path
    my $Wt = PDL::IO::NYHead::Interp::interp_weights($Es, pdl(0, -75, -50)->dummy(0));   # (21,1)
    my @chains = ([qw(O1 P3 C3 F3 Fp1)], [qw(O2 P4 C4 F4 Fp2)], [qw(Pz Cz Fz)], [qw(P7 T7 F7)], [qw(P8 T8 F8)]);
    my $indep = sub {                                            # the alpha-like test written out channel by channel
        my %v = @_;
        my ($mx) = sort { abs($v{$b}) <=> abs($v{$a}) } @nm[0 .. 18];
        return 0 unless $mx eq 'O1' || $mx eq 'O2';
        for (qw(O2 P3 P4 Pz)) { return 0 unless $v{O1} * $v{$_} > 0 }
        for my $c (@chains) { for (1 .. $#$c) { return 0 if abs($v{ $c->[$_] }) > abs($v{ $c->[ $_ - 1 ] }) } }
        return 1;
    };
    my $field = sub {                                            # channel values of moments m_s (3) per site, via the formula
        my ($sites, $mom) = @_;
        my %v;
        for my $i (0 .. $#nm) {
            my $p = $Es->slice("($i),");
            my $sum = 0;
            for my $s (0 .. $#$sites) {
                for my $hm ([ pdl(@{ $sites->[$s] }), $mom->[$s] ], [ pdl(@{ $sites->[$s] }) * pdl(-1, 1, 1), $mom->[$s] * pdl(-1, 1, 1) ]) {
                    my $d = $p - $hm->[0];
                    $sum += ($hm->[1] * $d)->sum / (($d**2)->sum)**1.5;
                }
            }
            $v{ $nm[$i] } = $sum;
        }
        return %v;
    };
    my $check_run = sub {
        my ($tag, $sites, %oo) = @_;
        my @GLs = map { $Gof->(@$_) } @$sites;
        my @GRs = map { $Gof->($_->[0] * -1, @$_[1, 2]) } @$sites;
        my %o3 = (%o, far => 'Th12', crit_ref => 'Th12', orient => undef, %oo);
        my ($eout, $ex) = $quiet->(sub { explore(\%o3, \@nm, [0 .. 18], \@GLs, \@GRs, [ map { "s$_" } 0 .. $#$sites ],
                                                drawer(\%o3, scalar @$sites), $Wt, ['Cv12'], {}) });
        print $eout;
        my $nk = $ex->{nkept};
        $chk->($nk > 0, sprintf '%s: %d of %d configurations kept', $tag, $nk, $ex->{ndraw});
        my ($fok, $aok) = (1, 1);
        for my $j (0 .. $nk - 1) {
            my %v = $field->($sites, [ map { $_->slice("($j),") } @{ $ex->{mom} } ]);
            $fok &&= abs($ex->{val}{$_}->at($j) - $v{$_}) < 1e-9 * abs($v{$_}) + 1e-15 for @nm;
            my %r = map { $_ => $v{$_} - $v{Th12} } @nm;
            $aok &&= $indep->(%r);
        }
        $chk->($fok, "$tag: the kept maps equal the dipole formula at every channel");
        $chk->($aok, "$tag: every kept map passes the test written out channel by channel");
        my $cw = ($Wt->slice(':,(0)') * pdl(map { $ex->{val}{$_}->at(0) } @nm))->sum;
        $chk->(abs($ex->{val}{Cv12}->at(0) - $cw) < 1e-12 + 1e-9 * abs($cw) && $eout =~ /Cv12 - Nk2/,
               "$tag: an own point goes through (interpolated values, pairs listed)");
        return ($ex, $eout);
    };
    my %j19; @j19{ @nm[0 .. 18] } = 0 .. 18;
    # one site: the orientation sweep
    my ($e1) = $check_run->('one site', [[-20, -72, 10]], ndir => 400);
    for my $neg (['an upward dipole at a deep parieto-occipital site', [-10, -65, 25], [0, 0, 1]],
                 ['a frontal dipole pointing forward', [-20, 50, 30], [0, 1, 0]]) {
        my %f = $field->([ $neg->[1] ], [ pdl(@{ $neg->[2] }) ]);
        my %r = map { $_ => $f{$_} - $f{Th12} } @nm;
        my $V1 = pdl(map { $r{$_} } @nm[0 .. 18])->dummy(1);    # (19,1)
        $chk->(!$indep->(%r) && alpha_like($V1, \%j19)->at(0) == 0, "$neg->[0] fails the alpha-like test");
    }
    # the test itself: a map with Fp1 larger than C3 is not kept, its correction is
    my %good = (O1 => 10, O2 => 10, P3 => 6, P4 => 6, Pz => 6, C3 => 3, C4 => 3, Cz => 3, F3 => 2, F4 => 2, Fz => 2,
                Fp1 => 1, Fp2 => 1, P7 => 5, P8 => 5, T7 => 2, T8 => 2, F7 => 1, F8 => 1);
    my %bad = (%good, Fp1 => -4, F3 => -3.5);
    my $V2 = pdl([ map { $good{$_} } @nm[0 .. 18] ], [ map { $bad{$_} } @nm[0 .. 18] ])->transpose->transpose;   # (19,2)
    my $k2 = alpha_like($V2, \%j19);
    $chk->($k2->at(0) == 1 && $k2->at(1) == 0, 'alpha_like keeps a decreasing map and drops one with |Fp1| > |F3|');
    # two sites: random draws, reproducible with the seed
    my ($e2) = $check_run->('two sites', [[-10, -65, 25], [-25, -80, 5]], batch => 300, target_kept => 60, max_draws => 30000, seed => 3);
    my ($e3) = $check_run->('two sites again', [[-10, -65, 25], [-25, -80, 5]], batch => 300, target_kept => 60, max_draws => 30000, seed => 3);
    $chk->($e2->{nkept} == $e3->{nkept} && ($e2->{mom}[1] - $e3->{mom}[1])->abs->max == 0, 'the same seed gives the same configurations');
    print $ok ? "selftest passed\n" : "selftest FAILED\n";
    return $ok;
}
