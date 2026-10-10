#!/usr/bin/env perl
# examples/heartbeat_hep_scalp.pl
#
#   perl examples/heartbeat_hep_scalp.pl REC.edf [options]
#
#     --bn1 L --bn2 L   the two non-cephalic electrodes (BN1, BN2); R peaks are
#                       detected on BN1-BN2
#     --skip S          ignore the first S seconds                    (10)
#     --band LO,HI      zero-phase Butterworth band, order 4 [Hz]     (0.5,45)
#     --reject UV       per-channel peak-to-peak limit, scalp and ears (150)
#     --nsurr N         surrogate (non-heartbeat) realisations         (200)
#     --seed N          seed for the surrogate offsets                 (1)
#     --windows LIST    latency windows [ms]      (100-200,200-300,300-450,450-600)
#     --out PREFIX      output prefix                 (REC without .edf + _hep)
#     --montage ELC     ASA .elc with the electrode positions
#     --no-fig          text and TSV only
#
# R-locked average of the scalp channels in the reference they were recorded in
# (no re-referencing), and what part of it the cardiac electric field can account
# for. For the R sample and for each latency window the script prints
#
#   z          amplitude of every scalp channel divided by that channel's SD over
#              the non-heartbeat averages (every R moved by up to +-0.5 s)
#   sum z^2    over scalp channels, with its rank among the non-heartbeat sets (p)
#   max |z|    with the 95th percentile of max |z| over channels at non-heartbeat
#              points, and the channels above it
#   halves     correlation over channels of the z maps from the first and the
#              second half of the beats
#   cos(R)     cosine between this map and the map at the R sample
#   BN(QRS)    share of the map (weighted by 1/SD^2) fitted by a_c*BN1 + b_c*BN2
#              with a_c, b_c estimated per channel in the QRS window (-40..60 ms)
#              of the averages; sum z^2 and p of what is left
#   BN(win)    the same with a_c, b_c estimated in --bnwin (150..450 ms); next to
#              it the share this fit takes from the non-heartbeat averages when
#              they are fitted to the real BN1, BN2 waves
#
# Outputs: PREFIX_maps.tsv (maps per window: uV, z, and uV after each removal),
#          PREFIX_time.tsv (per-sample RMS over scalp channels, non-heartbeat band, averages),
#          PREFIX_wave.png, PREFIX_maps.png

use strict;
use warnings;
use PDL;
use PDL::MatrixOps qw(inv);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use PDL::EEG::IO::EDF   qw(read_edf);
use PDL::EEG::Heartbeat qw(:all);

my ($bn1, $bn2, $skip, $band, $reject, $nsurr, $seed) = ('BN1', 'BN2', 10, '0.5,45', 150, 200, 1);
my ($wins, $bnwin, $out, $montage, $nofig) = ('100-200,200-300,300-450,450-600', '150,450');
GetOptions('bn1=s' => \$bn1, 'bn2=s' => \$bn2, 'skip=f' => \$skip, 'band=s' => \$band,
           'reject=f' => \$reject, 'nsurr=i' => \$nsurr, 'seed=i' => \$seed,
           'windows=s' => \$wins, 'bnwin=s' => \$bnwin, 'out=s' => \$out,
           'montage=s' => \$montage, 'no-fig' => \$nofig) or die;
my $file = shift @ARGV or die "usage: $0 REC.edf [options]\n";
($out = $file) =~ s/\.edf$//i, $out .= '_hep' unless defined $out;
$montage //= "$FindBin::Bin/../t/data/standard_1020_subset.elc";
my ($f_lo, $f_hi) = split /,/, $band;
my @win  = map { [ split /-/, $_, 2 ] } split /,/, $wins;
my @bnw  = split /,/, $bnwin;
$nsurr >= 20 or die "--nsurr must be at least 20\n";

# ---- read, select, filter --------------------------------------------------
my $e  = read_edf($file);
my $fs = $e->{fs};
my %ix; @ix{ @{ $e->{labels} } } = 0 .. $#{ $e->{labels} };
my %alt   = (T3 => 'T7', T4 => 'T8', T5 => 'P7', T6 => 'P8');
my @scalp = grep { exists $ix{$_} } map { exists $ix{$_} ? $_ : ($alt{$_} // $_) }
            qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 C3 Cz C4 T4 T5 P3 Pz P4 T6 O1 O2);
my @ear   = grep { exists $ix{$_} } qw(A1 A2);
exists $ix{$_} or die "channel '$_' not in $file\n" for $bn1, $bn2;
my @head  = (@scalp, @ear);
my @ch    = (@head, $bn1, $bn2);
my ($ns, $nh) = (scalar @scalp, scalar @head);
my $D  = zp_butter($e->{data}->dice_axis(0, pdl(long, [ @ix{@ch} ])),
                   fs => $fs, low => $f_lo, high => $f_hi, order => 4);
my $nt = $D->dim(1);

# ---- R peaks, real and surrogate averages ----------------------------------
my ($tmin, $tmax, @base) = (-0.3, 0.8, -0.30, -0.22);
my ($r) = detect_r($D->slice("($nh)") - $D->slice("(" . ($nh + 1) . ")"), $fs);
my ($n0, $n1) = (int($tmin * $fs - 0.5), int($tmax * $fs + 0.5));
my ($lo, $hi) = (int($skip * $fs) - $n0, $nt - $n1);
$r = $r->where(($r >= $lo) & ($r < $hi));
my $lim = ones($nh + 2) * 1e30;
$lim->slice("0:" . ($nh - 1)) .= $reject;
my %opt = (tmin => $tmin, tmax => $tmax, baseline => \@base, reject => $lim);

sub halves {                    # -> (full, first half, second half) averages, each (nch, nep)
    my $rr = shift;
    my $nr = $rr->nelem;
    my @h  = map { beat_average($D, $rr->slice($_), $fs, %opt) } '0:' . (int($nr / 2) - 1), int($nr / 2) . ':-1';
    my ($na, $nb) = map { $_->{n}->dummy(1) } @h;
    return (($h[0]{avg} * $na + $h[1]{avg} * $nb) / ($na + $nb), $h[0]{avg}, $h[1]{avg}, $h[0]);
}
my ($A, $Ao, $Ae, $h0) = halves($r);          # full, first half, second half
my $t   = $h0->{t_ms};
my $nep = $t->nelem;
my $SA  = zeroes(double, $nh + 2, $nep, $nsurr);       # surrogate averages
my $SO  = zeroes(double, $nh + 2, $nep, $nsurr);
my $SE  = zeroes(double, $nh + 2, $nep, $nsurr);
for my $k (0 .. $nsurr - 1) {
    my $sr = surrogate_r($r, $nt, $fs, min => 0, max => 0.5, lo => $lo, hi => $hi,
                         seed => $seed * 100000 + $k + 1);
    my ($a, $o, $ev) = halves($sr);
    $SA->slice(":,:,($k)") .= $a; $SO->slice(":,:,($k)") .= $o; $SE->slice(":,:,($k)") .= $ev;
}
printf "%s\n  %s, %d Hz, %.1f s, band %s-%s Hz, reference as recorded\n",
    $file, $e->{recording}, $fs, $nt / $fs, $f_lo, $f_hi;
printf "  %d beats (R on %s-%s), %d non-heartbeat sets, beats kept per scalp channel %d-%d\n",
    $r->nelem, $bn1, $bn2, $nsurr, beat_average($D, $r, $fs, %opt)->{n}->slice("0:" . ($ns - 1))->minmax;

# ---- BN regression coefficients (per head channel) -------------------------
sub widx { my ($a, $b) = @_; which(($t >= $a) & ($t <= $b)) }
sub bn_coef {                    # least squares  head(t) ~ a*bn1(t) + b*bn2(t)  in a window
    my ($a, $b, $Y) = @_;        # $Y: averages to fit (default the real ones), (nch, nep)
    my $k = widx($a, $b);
    my $B = $A->slice("$nh:" . ($nh + 1))->dice_axis(1, $k);                # 2 cols, nk rows
    $Y = ($Y // $A)->slice("0:" . ($nh - 1))->dice_axis(1, $k);             # nh cols, nk rows
    return inv($B->transpose x $B) x $B->transpose x $Y;                    # nh cols, 2 rows
}
sub bn_fit {                     # coefficients (nh,2), averages (nch, n...) -> fitted (nh, n...)
    my ($C, $V) = @_;
    return $C->slice(":,(0)") * $V->slice("($nh)")->dummy(0)
         + $C->slice(":,(1)") * $V->slice("(" . ($nh + 1) . ")")->dummy(0);
}
my $Cq = bn_coef(-40, 60);
my $Cw = bn_coef(@bnw);

# ---- maps: R sample and windows --------------------------------------------
sub wmean { my ($V, $a, $b) = @_; $V->dice_axis(1, widx($a, $b))->mv(1, 0)->average }   # drops the time dim
my $iR = widx(0, 0)->at(0);
my @map = (['R 0 ms', 0, 0], map { ["$_->[0]-$_->[1] ms", @$_] } @win);
my $S0   = "0:" . ($ns - 1);                             # scalp channels
my $H0   = "0:" . ($nh - 1);                             # scalp + ears
my $pct  = sub { my ($v, $p) = @_; $v->qsort->at(int($p * ($v->nelem - 1) + 0.5)) };
my $pval = sub { my ($obs, $null) = @_; (1 + sum($null >= $obs)) / (1 + $null->nelem) };
my $corr = sub { my ($a, $b) = map { $_ - $_->avg } @_; sum($a * $b) / sqrt(sum($a * $a) * sum($b * $b)) };
my $cos  = sub { my ($p, $q) = @_; sum($p * $q) / sqrt(sum($p * $p) * sum($q * $q)) };
my $mR   = $A->slice("$H0,($iR)");
my %keep;
print "\nMaps in the recorded reference; z = amplitude / SD over the non-heartbeat averages; p = rank among them.\n";
for my $m (@map) {
    my ($name, $w0, $w1) = @$m;
    my $full = wmean($A, $w0, $w1);                        # (nch)
    my $sful = wmean($SA, $w0, $w1);                       # (nch, nsurr)
    my ($v, $sv) = ($full->slice($H0), $sful->slice($H0));
    my $sd   = sqrt((($sv - $sv->mv(1, 0)->average->dummy(1)) ** 2)->mv(1, 0)->sumover / ($nsurr - 1));   # (nh)
    my ($z, $sz) = ($v / $sd, $sv / $sd->dummy(1));
    my $x2   = sub { ($_[0]->slice($S0) ** 2)->sumover };            # sum z^2 over scalp, broadcasts
    my $zs   = $z->slice($S0);
    my $thr  = $pct->($sz->slice($S0)->abs->maximum, 0.95);
    my @over = sort { abs($zs->at($b)) <=> abs($zs->at($a)) } grep { abs($zs->at($_)) > $thr } 0 .. $ns - 1;
    my ($vo, $ve) = map { (wmean($_, $w0, $w1)->slice($H0) / $sd)->slice($S0) } $Ao, $Ae;
    my $srh  = pdl(map { my $k = $_;
        $corr->(map { (wmean($_->slice(":,:,($k)"), $w0, $w1)->slice($H0) / $sd)->slice($S0) } $SO, $SE) } 0 .. $nsurr - 1);
    my $wsh  = sub { my ($u, $fit) = @_;                # weighted share over scalp channels
        1 - sum((($u - $fit) / $sd)->slice($S0) ** 2) / sum(($u / $sd)->slice($S0) ** 2) };
    # BN fits: QRS coefficients carried over; window coefficients (in-sample)
    my $qf   = bn_fit($Cq, $full->dummy(1))->slice(":,(0)");
    my $wf   = bn_fit($Cw, $full->dummy(1))->slice(":,(0)");
    my $sqf  = bn_fit($Cq, $sful);                                           # same a_c, b_c on non-hb BN
    my $swf  = cat(map { bn_fit(bn_coef(@bnw, $SA->slice(":,:,($_)")), $full->dummy(1))->slice(":,(0)") }
                   0 .. $nsurr - 1);                                         # non-hb scalp refitted to real BN
    my $swsh = pdl(map { $wsh->($sv->slice(":,($_)"), $swf->slice(":,($_)")) } 0 .. $nsurr - 1);
    printf "\n%s\n", $name;
    printf "  uV %s\n", join(' ', map { sprintf '%s %+.1f', $scalp[$_], $v->at($_) } 0 .. $ns - 1);
    printf "  z  %s\n", join(' ', map { sprintf '%s %+.1f', $scalp[$_], $z->at($_) } 0 .. $ns - 1);
    printf "  sum z^2 = %.0f (non-hb mean %.0f), p = %.3f\n", $x2->($z), $x2->($sz)->avg, $pval->($x2->($z), $x2->($sz));
    printf "  max |z| = %.1f at %s; non-hb 95th of max |z| = %.1f; above it: %s\n",
        $zs->abs->max, $scalp[ $zs->abs->maximum_ind->sclr ], $thr, (@over ? join(' ', @scalp[@over]) : 'none');
    printf "  %s %+.2f, %s %+.2f uV (mean in the window)\n", $bn1, $full->at($nh), $bn2, $full->at($nh + 1);
    printf "  halves r = %+.2f (non-hb mean %+.2f, 95th %+.2f)\n", $corr->($vo, $ve), $srh->avg, $pct->($srh, 0.95);
    printf "  cos(R) = %+.2f (non-hb 95th of |cos| %.2f)\n", $cos->($v, $mR),
        $pct->(pdl(map { abs($cos->($sv->slice(":,($_)"), $mR)) } 0 .. $nsurr - 1), 0.95) unless $w0 == 0 && $w1 == 0;
    printf "  BN(QRS) share %.2f; left sum z^2 = %.0f, p = %.3f\n", $wsh->($v, $qf),
        $x2->(($v - $qf) / $sd), $pval->($x2->(($v - $qf) / $sd), $x2->(($sv - $sqf) / $sd->dummy(1)));
    printf "  BN(%s..%s) share %.2f (non-hb fitted the same way: mean %.2f, 95th %.2f); left sum z^2 = %.0f, p = %.3f\n",
        @bnw, $wsh->($v, $wf), $swsh->avg, $pct->($swsh, 0.95),
        $x2->(($v - $wf) / $sd), $pval->($x2->(($v - $wf) / $sd), $x2->(($sv - $swf) / $sd->dummy(1)));
    $keep{$name} = { raw => $v, z => $z, bnq => $v - $qf, bnw => $v - $wf };
}

# ---- TSV -------------------------------------------------------------------
{
    open my $fh, '>', "${out}_maps.tsv" or die "${out}_maps.tsv: $!";
    print $fh join("\t", 'window', 'kind', @head), "\n";
    for my $m (@map) {
        for my $k (qw(raw z bnq bnw)) {
            print $fh join("\t", $m->[0], $k, map { sprintf '%.4f', $_ } $keep{ $m->[0] }{$k}->list), "\n";
        }
    }
}
my $Ah    = $A->slice("0:" . ($nh - 1));
my $rmsT  = sqrt(($Ah->slice("0:" . ($ns - 1)) ** 2)->average);                       # (nep)
my $srmsT = sqrt(($SA->slice("0:" . ($ns - 1)) ** 2)->average);                       # (nep, nsurr)
my $s95   = $srmsT->mv(1, 0)->qsort->slice("(" . int(0.95 * ($nsurr - 1) + 0.5) . ")");   # per sample
my $smaxT = $SA->slice("0:" . ($ns - 1))->abs->maximum->maximum;                      # (nsurr): max over ch, time
my $thrT  = $pct->($smaxT, 0.95);
{
    open my $fh, '>', "${out}_time.tsv" or die $!;
    print $fh join("\t", qw(time_ms rms_scalp nonhb_rms_mean nonhb_rms_95), @head, $bn1, $bn2), "\n";
    my $sm = $srmsT->mv(1, 0)->average;
    for my $i (0 .. $nep - 1) {
        print $fh join("\t", $t->at($i), map { sprintf '%.4f', $_ }
            $rmsT->at($i), $sm->at($i), $s95->at($i), $A->slice(":,($i)")->list), "\n";
    }
}
printf "\nwhole epoch: 95th percentile of the largest |amplitude| over all scalp channels and samples at non-heartbeat points = %.2f uV\n", $thrT;
print "wrote ${out}_maps.tsv, ${out}_time.tsv\n";
exit 0 if $nofig;

# ---- figures ---------------------------------------------------------------
eval { require PDL::Graphics::Cairo; require PDL::EEG::MAP2D; 1 }
    or do { print "PDL::Graphics::Cairo not available: no figures\n"; exit 0 };
PDL::Graphics::Cairo->import(qw(subplots figure));
{
    my ($fig, $p, $q, $g) = subplots(3, 1, width => 900, height => 940);
    $p->line($t, $Ah->slice("($_)"), color => '#1565C0', lw => 0.7) for 0 .. $ns - 1;
    $p->axhline($_ * $thrT, color => '#C62828', lw => 0.8) for -1, 1;
    $p->set_title(sprintf('R-locked average, reference as recorded: %d scalp channels, %d beats', $ns, $r->nelem));
    $p->invert_yaxis; $p->set_ylabel('uV (negative up)');
    $q->line($t, $rmsT, color => '#1565C0', lw => 1.4, label => 'RMS over scalp channels');
    $q->line($t, $s95,  color => '#C62828', lw => 0.9, label => 'non-heartbeat 95th percentile');
    $q->set_ylabel('uV'); $q->legend;
    my %j; @j{@head} = 0 .. $#head;
    my @show = grep { exists $j{$_} } qw(Fz Cz Pz);
    my @cc   = ('#2E7D32', '#1565C0', '#6A1B9A');
    $g->line($t, $Ah->slice("($j{$show[$_]})"), color => $cc[$_], lw => 1.4, label => $show[$_]) for 0 .. $#show;
    $g->line($t, ($A->slice("($nh)") - $A->slice("(" . ($nh + 1) . ")")) / 20, color => '#C62828', lw => 0.9,
             label => "($bn1-$bn2)/20");
    $g->axhline(0, color => '#999999', lw => 0.5);
    $g->set_ylim(-6, 6); $g->invert_yaxis; $g->set_ylabel('uV (negative up)'); $g->legend;
    $g->set_xlabel(sprintf('time from R peak [ms]   (top panel lines: non-heartbeat 95%% of max |amplitude|, %.2f uV)', $thrT));
    $_->axvline(0, color => '#999999', lw => 0.5) for $p, $q, $g;
    $p->axhline(0, color => '#999999', lw => 0.5);
    $fig->save("${out}_wave.png");
}
{
    my @kind = (['raw', 'uV as recorded'], ['z', 'z'],
                ['bnq', 'uV, BN1,BN2 fit (QRS window) removed'], ['bnw', "uV, BN1,BN2 fit ($bnw[0]..$bnw[1] ms) removed"]);
    my $nc  = scalar @map;
    my $fig = figure(width => 250 * $nc + 40, height => 270 * @kind + 30);
    my $gs  = $fig->add_gridspec(scalar(@kind), $nc, left => 0.02, right => 0.98, top => 0.95,
                                 bottom => 0.02, wspace => 0.10, hspace => 0.22);
    for my $c (0 .. $nc - 1) {
        my $nm = $map[$c][0];
        for my $rw (0 .. $#kind) {
            my $cl = $keep{$nm}{ $kind[$rw][0] eq 'z' ? 'z' : 'raw' }->slice("0:" . ($ns - 1))->abs->max->sclr;
            my $ax = $fig->add_subplot($gs->at($rw, $c));
            PDL::EEG::MAP2D::plot_topomap(
                values => $keep{$nm}{ $kind[$rw][0] }->slice("0:" . ($ns - 1)), labels => \@scalp,
                montage => $montage, clim => $cl, contours => 0, names => 0, title_size => 9,
                title => sprintf('%s, +-%.1f', $nm, $cl), _ax => $ax, _fig => $fig);
        }
    }
    $fig->suptitle('rows: ' . join(' / ', map { $_->[1] } @kind) . '   (red positive)', fontsize => 9);
    $fig->save("${out}_maps.png");
}
print "wrote ${out}_wave.png, ${out}_maps.png\n";
