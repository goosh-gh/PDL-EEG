#!/usr/bin/env perl
#
# nyhead_dipole_gof.pl — how well does one cortical dipole of the New York Head explain
# a scalp topography?  Used to compare topographies of known non-brain origin (the
# cardiac field on the head) with candidate brain topographies (HEP) on the same
# electrodes.
#
#   perl examples/nyhead_dipole_gof.pl --mat ~/src/NYHead/sa_nyhead.mat PATTERNS.tsv
#   perl examples/nyhead_dipole_gof.pl --mat ... PATTERNS.tsv --columns QRS.PC1,T2.PC1 --free
#   perl examples/nyhead_dipole_gof.pl --selftest            # synthetic check, no .mat needed
#
#   PATTERNS.tsv     tab separated; first column = channel name, one topography per further
#                    column, first line = header (heartbeat_neck_field_rank.pl writes
#                    PREFIX_head_patterns.tsv in this form)
#   --mat FILE       sa_nyhead.mat                                  (~/src/NYHead/sa_nyhead.mat)
#   --columns LIST   topographies to fit                            (all columns)
#   --map LIST       OLD=NEW channel names for NYHead               (T3=T7,T4=T8,T5=P7,T6=P8)
#                    a channel that is still not an NYHead electrode is left out (listed)
#   --exclude LIST   channels to leave out                          (none)
#   --free           also fit free-orientation dipoles (/sa/cortex75K/V_fem, about 0.4 GB)
#   --nrand N        random (spatially white) topographies for the reference distribution (200)
#   --seed N         seed for them                                  (1)
#   --top N          best sources listed per topography             (5)
#   --strings S      how PDL::IO::NYHead reads the names: auto | h5dump | lowlevel (auto)
#
# The leadfield (V_fem_normal, normal to the cortex, 74382 sources) is restricted to the
# channels of the file and re-referenced to their average, as the topographies are.
# For a source with gain vector g (centred) and the centred topography y:
#   GOF = (g.y)^2 / ((g.g)(y.y))           (normal orientation)
#   GOF = b' M^-1 b / (y.y),  b = G'y, M = G'G, G = the 3 orientations   (--free)
# GOF is the share of the topography's power one dipole explains (1 = perfect).
# With few electrodes smooth topographies are explained well by deep dipoles, so a GOF
# means something only next to the GOF of topographies of known origin on the same
# channels, and next to the reference distribution of random topographies.

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use PDL;
use PDL::MatrixOps qw(inv);
use Scalar::Util qw(looks_like_number);

my %o = (mat => "$ENV{HOME}/src/NYHead/sa_nyhead.mat", columns => '', map => 'T3=T7,T4=T8,T5=P7,T6=P8',
         exclude => '', free => 0, nrand => 200, seed => 1, top => 5, strings => 'auto', selftest => 0);
GetOptions('mat=s' => \$o{mat}, 'columns=s' => \$o{columns}, 'map=s' => \$o{map}, 'exclude=s' => \$o{exclude},
           'free' => \$o{free}, 'nrand=i' => \$o{nrand}, 'seed=i' => \$o{seed}, 'top=i' => \$o{top},
           'strings=s' => \$o{strings}, 'selftest' => \$o{selftest})
    or die "bad options (see the header of $0)\n";

if ($o{selftest}) { selftest(); exit 0 }

my $tsv = shift @ARGV or die "usage: $0 --mat sa_nyhead.mat PATTERNS.tsv [options]\n";
-f $o{mat} or die "not found: $o{mat} (give it with --mat)\n";
require PDL::IO::NYHead;

# ---- topographies ---------------------------------------------------------------
open my $fh, '<', $tsv or die "$tsv: $!";
chomp(my $hdr = <$fh>);
my (undef, @cols) = split /\t/, $hdr;
my (@ch, @rows);
while (<$fh>) {
    chomp; next if /^\s*(#|$)/;
    my ($c, @v) = split /\t/;
    push @ch, $c; push @rows, [@v];
}
close $fh;
my @want = $o{columns} ne '' ? split(/,/, $o{columns}) : @cols;
my %colix; @colix{@cols} = 0 .. $#cols;
exists $colix{$_} or die "column '$_' not in $tsv (have: @cols)\n" for @want;
# every value used must be a number: a channel table or other text file stops here
for my $w (@want) {
    my $j = $colix{$w};
    for my $k (0 .. $#ch) {
        my $v = $rows[$k][$j];
        defined $v && looks_like_number($v)
            or die sprintf "%s: column '%s', channel '%s' is not a number (%s).\n"
                         . "  The file must be a topography table: channel name, then one number per column.\n",
                       $tsv, $w, $ch[$k], defined $v ? "'$v'" : 'empty';
    }
}

# ---- channels -> NYHead electrodes ---------------------------------------------
my $ny   = PDL::IO::NYHead->new($o{mat}, strings => $o{strings});
my @lab  = @{ $ny->electrode_labels };
my %nyix; @nyix{@lab} = 0 .. $#lab;
my %map  = map { split /=/, $_, 2 } grep { /=/ } split /,/, $o{map};
my %excl = map { $_ => 1 } split /,/, $o{exclude};
my (@use, @nyi, @miss);
for my $k (0 .. $#ch) {
    next if $excl{ $ch[$k] };
    my $n = $map{ $ch[$k] } // $ch[$k];
    if (exists $nyix{$n}) { push @use, $k; push @nyi, $nyix{$n} } else { push @miss, $ch[$k] }
}
@use >= 4 or die "fewer than 4 channels are NYHead electrodes\n";
printf "%s: %d channels used: %s\n", $tsv, scalar @use,
    join(' ', map { my $n = $map{ $ch[$use[$_]] } // $ch[$use[$_]]; $n eq $ch[$use[$_]] ? $n : "$ch[$use[$_]]=$n" } 0 .. $#use);
printf "left out (not an NYHead electrode name; give --map OLD=NEW): %s\n", join(' ', @miss) if @miss;

my $Y = pdl(map { my $j = $colix{$_}; [ map { $rows[$_][$j] } @use ] } @want);   # (ne, ntopo)
$Y = $Y->dummy(1) if $Y->ndims == 1;
my $ne = @use;

my $L  = $ny->leadfield->dice_axis(0, pdl(long, \@nyi))->double;           # (ne, ns)
my $ns = $L->dim(1);
printf "leadfield V_fem_normal: %d sources\n", $ns;
my $G3 = $o{free} ? $ny->leadfield_full->dice_axis(0, pdl(long, \@nyi))->double : undef;   # (ne, ns, 3)

my ($gof, $gof3) = fit_gof($L, $G3, $Y);                                   # (ns, ntopo)
my $vc = $ny->surface('cortex75K')->{vc};

# reference distribution: random spatially white topographies on the same channels
CORE::srand($o{seed});
my $R = pdl(map { [ map { gauss() } 1 .. $ne ] } 1 .. $o{nrand});           # (ne, nrand)
my (@rb, @rb3);                                                             # in chunks of 25 (memory)
for (my $c0 = 0; $c0 < $o{nrand}; $c0 += 25) {
    my $c1 = $c0 + 24 < $o{nrand} - 1 ? $c0 + 24 : $o{nrand} - 1;
    my ($rg, $rg3) = fit_gof($L, $G3, $R->slice(":,$c0:$c1"));
    push @rb, $rg->maximum->list;                                           # best GOF of each
    push @rb3, $rg3->maximum->list if defined $rg3;
}
my $rbest  = pdl(\@rb);
my $rbest3 = @rb3 ? pdl(\@rb3) : undef;
printf "random white topographies (%d): best-dipole GOF median %.3f, 95th pct %.3f%s\n",
    $o{nrand}, $rbest->median, $rbest->pct(0.95),
    defined $rbest3 ? sprintf(' (free orientation: median %.3f, 95th pct %.3f)', $rbest3->median, $rbest3->pct(0.95)) : '';

for my $t (0 .. $#want) {
    my $g = $gof->slice(":,($t)");
    my @ord = $g->qsorti->slice("-1:0")->list;
    printf "\n%s: best GOF %.3f (normal)%s\n", $want[$t], $g->at($ord[0]),
        defined $gof3 ? sprintf(', %.3f (free orientation)', $gof3->slice(":,($t)")->max) : '';
    for my $i (@ord[0 .. $o{top} - 1]) {
        printf "  source %6d  GOF %.3f  MNI (%6.1f %6.1f %6.1f)  %s\n", $i, $g->at($i),
            (map { $vc->at($i, $_) } 0 .. 2), $ny->area_of_vertex($i) // '(no HO area)';
    }
    my $frac = ($g >= 0.9 * $g->at($ord[0]))->sum / $ns;
    printf "  sources within 90%% of the best GOF: %.1f %% of the cortex\n", 100 * $frac;
}

# =====================================================================================
# GOF of every source for every topography.  $L (ne, ns), $G3 (ne, ns, 3) or undef,
# $Y (ne, ntopo).  Leadfields and topographies are re-referenced to the average of the
# ne channels here.  Returns (ns, ntopo) piddles; the second is undef without $G3.
sub fit_gof {
    my ($L, $G3, $Y) = @_;
    my $Lc = $L - $L->average->dummy(0);                          # (ne, ns)
    my $Yc = $Y - $Y->average->dummy(0);                          # (ne, nt)
    my $yy = ($Yc**2)->sumover;                                   # (nt)
    my $gg = ($Lc**2)->sumover;                                   # (ns)
    # PDL matrices are (columns, rows): Lc is ns x ne, Yc is ne x nt, so Lc x Yc is ns x nt
    my $gy = $Lc x $Yc->transpose;                                # (nt, ns): Yc' is ne x nt
    my $gof = ($gy**2 / $gg->dummy(0) / $yy->dummy(1))->transpose->sever;    # (ns, nt)
    return ($gof, undef) unless defined $G3;
    my $Gc = $G3 - $G3->average->dummy(0);                        # (ne, ns, 3)
    my $Gm = $Gc->mv(2, 1);                                       # (ne, 3, ns)
    my $M  = ($Gm->dummy(1) * $Gm->dummy(2))->sumover;            # (3, 3, ns): M[i,j] = sum_e G_ei G_ej
    my $Mi = inv($M);                                             # broadcast over sources
    # b_k = G_k' y for the 3 orientations, each a matrix product (ns x ne)(ne x nt)
    my $b  = cat(map { $Gc->slice(":,:,($_)") x $Yc->transpose } 0 .. 2);  # (nt, ns, 3)
    $b = $b->mv(2, 0)->mv(1, 2)->sever;                           # (3, ns, nt)
    my $Mb = zeroes(3, $b->dim(1), $b->dim(2));                   # (M^-1 b)_c = sum_a Mi[a,c] b_a (M symmetric)
    $Mb += $Mi->slice("($_),:,:")->dummy(2) * $b->slice("($_),:,:")->dummy(0) for 0 .. 2;
    my $g3 = ($b * $Mb)->sumover / $yy->dummy(0);                 # (ns, nt)
    return ($gof, $g3);
}

sub gauss {                                                       # standard normal, CORE::rand
    my ($u, $v) = (CORE::rand() || 1e-12, CORE::rand());
    return sqrt(-2 * log($u)) * cos(6.283185307179586 * $v);
}

# synthetic check: a random leadfield with planted dipoles
sub selftest {
    CORE::srand(3);
    my ($ne, $ns) = (17, 400);
    my $G3 = pdl([ map { gauss() } 1 .. $ne * $ns * 3 ])->reshape($ne, $ns, 3);   # (ne, ns, 3)
    my $L  = $G3->slice(':,:,(0)')->sever;                        # "normal" = first orientation
    my $y1 = $L->slice(':,(123)') * 2.5 + 7;                       # planted normal dipole + offset (reference)
    # planted free dipole: G (ne x 3 matrix = dims (3, ne)) times the moment (3 x 1 = dims (1, 3))
    my $y2 = ($G3->slice(':,(45),:')->transpose x pdl([[0.3], [-1.2], [0.8]]))->flat;
    my $Y  = cat($y1, $y2);                                       # (ne, 2)
    my ($g, $g3) = fit_gof($L, $G3, $Y);
    my $ok = 1;
    my $chk = sub { my ($c, $m) = @_; printf "  [%s] %s\n", $c ? 'PASS' : 'FAIL', $m; $ok &&= $c };
    $chk->($g->slice(':,(0)')->maximum_ind == 123 && abs($g->at(123, 0) - 1) < 1e-12, 'normal: planted source 123 found, GOF 1');
    $chk->($g3->slice(':,(1)')->maximum_ind == 45 && abs($g3->at(45, 1) - 1) < 1e-12, 'free: planted source 45 found, GOF 1');
    $chk->($g3->at(123, 0) >= $g->at(123, 0) - 1e-12, 'free GOF >= normal GOF at the same source');
    $chk->(all($g3 <= 1 + 1e-12) && all($g >= 0), 'GOF within [0, 1]');
    # brute force for one source, independent of the broadcast code
    my $G = $G3->slice(':,(123),:'); $G = $G - $G->average->dummy(0);        # (ne, 3), centred
    my $Gm = $G->transpose;                                                  # ne x 3 matrix
    my $yc = ($y2 - $y2->avg)->dummy(0);                                     # ne x 1
    my $fit = $Gm x inv($Gm->transpose x $Gm) x ($Gm->transpose x $yc);
    my $bf  = (($fit**2)->sum / ($yc**2)->sum)->sclr;
    $chk->(abs($bf - $g3->at(123, 1)) < 1e-10, sprintf('free GOF at another source equals the brute-force projection (%.6f)', $bf));
    # negative control: a topography that is pure reference offset plus noise is not GOF 1 anywhere
    my $ynoise = pdl([ map { gauss() } 1 .. $ne ]);
    my ($gn) = fit_gof($L, undef, $ynoise->dummy(1));
    $chk->($gn->max < 0.9, sprintf('random topography: best GOF %.3f < 0.9', $gn->max));
    print $ok ? "selftest passed\n" : "selftest FAILED\n";
    exit($ok ? 0 : 1);
}
