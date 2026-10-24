#!/usr/bin/env perl
# examples/heartbeat_slowwave_tests.pl
#
#   perl examples/heartbeat_slowwave_tests.pl REC1.edf [REC2.edf ...] [options]
#
#     --bn1 L --bn2 L   non-cephalic electrodes; R peaks are detected on BN1-BN2
#     --skip S          ignore the first S seconds                        (10)
#     --band LO,HI      zero-phase Butterworth band, order 4 [Hz]         (0.5,40)
#     --exclude LIST    head channels left out of everything              (none)
#     --vessel LIST     channels left out for the second average reference
#                       (electrodes that may sit over an artery)  (Fp1,Fp2,F8,A2)
#     --reject UV       a beat is dropped when any head channel's peak-to-peak
#                       over -400..+1400 ms exceeds this                   (200)
#     --ngroup N        groups of beats by the RR interval to the next R    (3)
#     --rrwin S         group by the mean RR over +-S/2 seconds around each
#                       beat instead of the beat's own RR (0: the beat's own) (60)
#     --chunk N         beats per chunk when averaging (memory)            (500)
#     --out PREFIX      output prefix (first REC without .edf + _slowwave)
#     --no-fig          text only
#
# Is the slow, heartbeat-locked wave (rising ~150-200 ms after R, largest near
# 300 ms, Cz/Pz against Fp/F8/A2 in the average reference) a brain response
# (HEP) or a pulse (vascular) artefact?  Three tests, per record:
#
#  A. reference: the slow-wave map in the average reference of all head
#     channels, and in the average reference without the --vessel channels;
#     the slow-wave weight at Cz and Pz as recorded (system reference), in both
#     average references.  If the Cz/Pz positivity were only the mirror image
#     of a few electrodes over arteries, it would vanish without them.
#  B. RR dependence: beats in --ngroup groups by the interval to the next R
#     (--rrwin > 0: by the mean RR over the surrounding --rrwin seconds, so that
#     the groups follow the slow changes of heart rate and not respiratory sinus
#     arrhythmia, which would sort the beats by breathing phase as well).
#     For each group, the slow-wave time course aligned to its own R and to the
#     next R.  The fall of the wave (first time after 450 ms below half of its
#     250..450 ms mean) is printed with the group's median RR: a fall at a fixed
#     latency after R points to a response locked to the beat; a fall that
#     moves with RR (fixed time before the next R) points to the cardiac cycle.
#  D. spatial smoothness: every channel is predicted from all the others by a
#     spherical spline (Perrin et al. 1989, m = 4, 50 Legendre terms, positions
#     from the standard 10-20 .elc); the share of the map's power the
#     predictions miss.  A field volume-conducted from a brain source (or from
#     the heart) is smooth and is predicted well; a wave that lives at a few
#     electrodes (a vessel under an electrode) is not.  Done for the slow-wave
#     map (mean 250..450 ms), the map at the R peak (smooth, volume conducted),
#     and 20 non-heartbeat averages (noise) at the same lags; channels with the
#     largest misses are listed.
#  C. across records (two or more files): the cosine between the slow-wave maps
#     of every pair of records over the channels they share.  A pulse artefact
#     depends on where each electrode happens to sit, so it changes when the
#     electrodes are put on again; a brain response should not.
#
# The slow-wave map: the R-locked average (every accepted beat, baseline
# -75..-50 ms) in the 450..650 ms window, with the first two principal maps of
# its QRS window (-45..55 ms) projected out, first principal component; sign
# so that its weight at 300 ms is positive.
#
# Output: PREFIX_slowwave.tsv (the maps of every record), and per record
#         <REC>_slowwave_rr.png, <REC>_slowwave_maps.png

use strict;
use warnings;
use PDL;
use PDL::MatrixOps ();
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use PDL::EEG::IO::EDF   qw(read_edf);
use PDL::EEG::Filter    qw(zp_butter);
use PDL::EEG::Heartbeat qw(detect_r beat_epochs);
use PDL::EEG::IO::ASA   qw(read_elc);
use PDL::EEG::Linalg    qw(eigh_sym);

my ($bn1, $bn2, $skip, $band, $exclude, $vessel, $reject, $ngroup, $chunk, $out, $nofig)
    = ('BN1', 'BN2', 10, '0.5,40', '', 'Fp1,Fp2,F8,A2', 200, 3, 500);
my $rrwin = 60;
GetOptions('bn1=s' => \$bn1, 'bn2=s' => \$bn2, 'skip=f' => \$skip, 'band=s' => \$band,
           'exclude=s' => \$exclude, 'vessel=s' => \$vessel, 'reject=f' => \$reject,
           'ngroup=i' => \$ngroup, 'rrwin=f' => \$rrwin, 'chunk=i' => \$chunk, 'out=s' => \$out,
           'no-fig' => \$nofig) or die;
@ARGV or die "usage: $0 REC1.edf [REC2.edf ...] [options]\n";
($out = $ARGV[0]) =~ s/\.edf$//i, $out .= '_slowwave' unless defined $out;
my ($f_lo, $f_hi) = split /,/, $band;
my %excl = map { $_ => 1 } split /,/, $exclude;
my %vess = map { $_ => 1 } split /,/, $vessel;
my @ALL  = qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 C3 Cz C4 T4 T5 P3 Pz P4 T6 O1 O2 A1 A2);
my %alias = (T3 => 'T7', T4 => 'T8', T5 => 'P7', T6 => 'P8');   # names in the .elc
my $montage = "$FindBin::Bin/../t/data/standard_1020_subset.elc";
my $have_fig = !$nofig && eval { require PDL::Graphics::Cairo; require PDL::EEG::MAP2D; 1 };
PDL::Graphics::Cairo->import(qw(figure subplots)) if $have_fig;

my %maps;                                   # file => { label => value }
for my $file (@ARGV) {
    my $res = run_record($file);
    $maps{$file} = $res if $res;
}

# ---- C. across records -----------------------------------------------------
my @files = grep { $maps{$_} } @ARGV;
open my $fh, '>', "${out}.tsv" or die "${out}.tsv: $!";
print $fh join("\t", 'ch', map { s{.*/}{}r } @files), "\n";
for my $ch (@ALL) {
    next unless grep { exists $maps{$_}{$ch} } @files;
    print $fh join("\t", $ch, map { exists $maps{$_}{$ch} ? sprintf('%.4f', $maps{$_}{$ch}) : '' } @files), "\n";
}
close $fh;
print "\nwrote ${out}.tsv\n";
if (@files >= 2) {
    print "\n== C. cosine between the slow-wave maps of the records (average reference over the shared channels)\n";
    printf "%-24s %s\n", '', join(' ', map { sprintf '%8.8s', s{.*/}{}r } @files);
    for my $f1 (@files) {
        printf "%-24.24s %s\n", $f1 =~ s{.*/}{}r, join(' ', map { sprintf '%8.3f', _cos_maps($maps{$f1}, $maps{$_}) } @files);
    }
}

# ---------------------------------------------------------------------------
sub run_record {
    my $file = shift;
    my $e  = read_edf($file);
    my $fs = $e->{fs};
    my %ix; @ix{ @{ $e->{labels} } } = 0 .. $#{ $e->{labels} };
    for ($bn1, $bn2) { exists $ix{$_} or do { print "$file: no channel $_, skipped\n"; return } }
    my @head = grep { exists $ix{$_} && !$excl{$_} } @ALL;
    my $nh   = @head;
    my $X = zp_butter($e->{data}->dice_axis(0, pdl(long, [ @ix{@head, $bn1, $bn2} ])),
                      fs => $fs, low => $f_lo, high => $f_hi, order => 4);
    my $Hrec = $X->slice("0:" . ($nh - 1) . ",:");
    my $nt   = $Hrec->dim(1);
    my $ecg  = $X->slice("($nh),:") - $X->slice("(" . ($nh + 1) . "),:");
    my $s    = sub { my $v = shift() * 1e-3 * $fs; $v >= 0 ? int($v + 0.5) : -int(-$v + 0.5) };
    my ($n0, $n1) = ($s->(-400), $s->(1400));
    my %wi = map { my ($nm, $a0, $a1) = @$_; $nm => [ $s->($a0) - $n0, $s->($a1) - $n0 - 1 ] }
             (['base', -75, -50], ['QRS', -45, 55], ['T', 250, 450], ['Late', 450, 650]);
    my $tms = (sequence($n1 - $n0) + $n0) * (1000 / $fs);

    # beats with a following R, whose epochs around both fit in the record
    my ($rall) = detect_r($ecg, $fs);
    my $rcur = $rall->slice("0:-2");
    my $rnx  = $rall->slice("1:-1");
    my $ok   = ($rcur >= $s->(1000 * $skip) - $n0) & ($rnx + $n1 <= $nt) & ($rcur + $n0 >= 0);
    $rcur = $rcur->where($ok);
    $rnx  = $rnx->where($ok);

    # rejection over the head channels as recorded (one pass, in chunks)
    my $Havg = $Hrec - $Hrec->average->dummy(0);
    my $keep = _accept($Hrec, $rcur, $n0, $n1);
    $rcur = $rcur->where($keep);
    $rnx  = $rnx->where($keep);
    my $rr = ($rnx - $rcur)->double * (1000 / $fs);           # ms
    printf "\n==== %s: fs %g Hz, %d head channels, %d beats accepted (reject %g uV), RR median %.0f ms (5..95%%: %.0f..%.0f)\n",
        $file, $fs, $nh, $rcur->nelem, $reject, $rr->median, $rr->pct(0.05), $rr->pct(0.95);
    return if $rcur->nelem < 30;

    # ---- A. slow-wave map in the average reference, all head channels --------
    my $A   = _mean($Havg, $rcur, $n0, $n1, $wi{base});
    my $sw  = _slow_map($A, \%wi, $s->(300) - $n0);              # nh, unit
    my $tc  = ($A x $sw->dummy(0))->flat;
    my @vi  = grep { !$vess{ $head[$_] } } 0 .. $nh - 1;
    my $Hv  = $Hrec->dice_axis(0, pdl(long, \@vi));
    $Hv     = $Hv - $Hv->average->dummy(0);
    my $Av  = _mean($Hv, $rcur, $n0, $n1, $wi{base});
    my $swv = _slow_map($Av, \%wi, $s->(300) - $n0);
    my $Arec = _mean($Hrec, $rcur, $n0, $n1, $wi{base});

    my ($t0, $t1) = @{ $wi{T} };
    my %amp;                                                   # 250..450 ms mean per channel
    $amp{avg}  = $A->slice(":,$t0:$t1")->mv(1, 0)->average;
    $amp{avgv} = $Av->slice(":,$t0:$t1")->mv(1, 0)->average;
    $amp{rec}  = $Arec->slice(":,$t0:$t1")->mv(1, 0)->average;
    print "-- A. slow-wave map (unit vector, average reference of all head channels):\n   ",
        join(' ', map { sprintf '%s %.2f', $head[$_], $sw->at($_) } 0 .. $nh - 1), "\n";
    print "   without ", join(',', grep { $vess{$_} } @head), ":\n   ",
        join(' ', map { sprintf '%s %.2f', $head[$vi[$_]], $swv->at($_) } 0 .. $#vi), "\n";
    my %hpos; @hpos{@head} = 0 .. $nh - 1;
    my %vpos; @vpos{ map { $head[$_] } @vi } = 0 .. $#vi;
    print "   mean over 250..450 ms [uV] (baseline -75..-50 ms):\n";
    printf "   %-5s %10s %10s %14s\n", 'ch', 'recorded', 'avg ref', 'avg w/o vessel';
    for my $ch (grep { exists $hpos{$_} } qw(Fz Cz Pz C3 C4 P3 P4 O1 O2 Fp1 Fp2 F8 A1 A2)) {
        printf "   %-5s %10.2f %10.2f %14s\n", $ch, $amp{rec}->at($hpos{$ch}), $amp{avg}->at($hpos{$ch}),
            exists $vpos{$ch} ? sprintf('%.2f', $amp{avgv}->at($vpos{$ch})) : '-';
    }
    printf "   cosine between the two maps over the channels they share: %.3f\n",
        _cos_vec($sw->index(pdl(long, \@vi)), $swv);

    # ---- D. spatial smoothness (leave-one-out spherical spline) ----------------
    {
        my @ml = map { $alias{$_} // $_ } @head;
        my $P  = _unit_positions(\@ml);                               # (3, nh)
        my $T  = _loo_matrix($P);                                     # nh x nh: v -> predictions
        my $i0r = $s->(0) - $n0;
        my %mp = (slow => $amp{avg}, 'R peak' => $A->slice(":,($i0r)")->flat->copy);
        my @nz;
        for my $k (1 .. 20) {
            # time points uniform over the usable record (shifted R times would
            # stay at a similar cardiac phase and share the heartbeat-locked signal)
            CORE::srand($k);
            my ($lo, $hi) = ($s->(1000 * $skip) - $n0, $nt - $n1);
            my $sr = pdl(long, [ map { $lo + int(CORE::rand() * ($hi - $lo)) } 1 .. $rcur->nelem ])->qsort;
            my $S  = _mean($Havg, $sr, $n0, $n1, $wi{base});
            push @nz, $S->slice(":,$t0:$t1")->mv(1, 0)->average;
        }
        print "-- D. leave-one-out spherical spline: share of the map power missed\n";
        for my $nm ('R peak', 'slow') {
            my $v = $mp{$nm} - $mp{$nm}->avg;
            my $d = $v - ($T x $v->dummy(0))->flat;
            my @top = (sort { abs($d->at($b)) <=> abs($d->at($a)) } 0 .. $nh - 1)[0 .. 4];
            printf "   %-7s %.3f   largest misses [uV]: %s\n", $nm, (($d**2)->sum / ($v**2)->sum)->sclr,
                join(' ', map { sprintf '%s %+.2f (map %+.2f)', $head[$_], $d->at($_), $v->at($_) } @top);
        }
        my @nr = map { my $v = $_ - $_->avg; my $d = $v - ($T x $v->dummy(0))->flat;
                       (($d**2)->sum / ($v**2)->sum)->sclr } @nz;
        my $nrp = pdl(\@nr);
        printf "   noise   %.3f   (median of 20 non-heartbeat averages; 5..95%%: %.3f..%.3f)\n",
            $nrp->median, $nrp->pct(0.05), $nrp->pct(0.95);
        # what a perfectly smooth map would miss: only its noise, i.e. about
        # (noise miss share) x (noise power / map power)
        my $npw = pdl(map { my $v = $_ - $_->avg; ($v**2)->sum } @nz)->median;
        for my $nm ('R peak', 'slow') {
            my $v = $mp{$nm} - $mp{$nm}->avg;
            my $ratio = ($npw / ($v**2)->sum)->sclr;
            printf "   %-7s map power / noise map power %.1f: a smooth map plus this noise would miss about %.3f\n",
                $nm, 1 / $ratio, $ratio * $nrp->median;
        }
    }

    # ---- B. RR groups ---------------------------------------------------------
    print "-- B. groups by the RR to the next R (slow-wave weight, all-channel average reference)\n";
    my $gv  = $rrwin > 0 ? _runmean_rr($rcur, $rr, $rrwin * $fs) : $rr;
    printf "   grouping variable: %s\n", $rrwin > 0
        ? sprintf('mean RR over +-%g s (5..95%%: %.0f..%.0f ms)', $rrwin / 2, $gv->pct(0.05), $gv->pct(0.95))
        : "the beat's own RR";
    my $ord = $gv->qsorti;
    my $nb  = $rr->nelem;
    my @grp;
    for my $g (0 .. $ngroup - 1) {
        my $i0 = int($g * $nb / $ngroup);
        my $i1 = int(($g + 1) * $nb / $ngroup) - 1;
        my $sel = $ord->slice("$i0:$i1");
        my $Ag  = _mean($Havg, $rcur->index($sel), $n0, $n1, $wi{base});
        my $An  = _mean($Havg, $rnx->index($sel),  $n0, $n1, $wi{base});
        my $tcg = ($Ag x $sw->dummy(0))->flat;
        my $tcn = ($An x $sw->dummy(0))->flat;
        my $ecg_g = _mean($ecg->dummy(0), $rcur->index($sel), $n0, $n1, $wi{base})->flat;
        my $pk  = $tcg->slice("$t0:$t1")->avg;
        my $rrm = $rr->index($sel)->median;
        my $gvm = $gv->index($sel)->median;
        my $i450 = $s->(450) - $n0;
        my $iend = _min($s->($rrm) - $n0, $tcg->nelem - 1);
        my $below = which($tcg->slice("$i450:$iend") < 0.5 * $pk);
        my $fall = $below->nelem ? $tms->at($below->at(0) + $i450) : undef;
        printf "   group %d: %4d beats, grouping value median %4.0f ms, RR median %4.0f ms (%4.0f..%4.0f); mean 250..450 ms %.2f uV; "
             . "fall below half at %s ms after R%s\n",
            $g + 1, $sel->nelem, $gvm, $rrm, $rr->index($sel)->minmax, $pk,
            defined $fall ? sprintf('%.0f', $fall) : '(none before the next R)',
            defined $fall ? sprintf(' = %.0f ms before the next R', $rrm - $fall) : '';
        push @grp, { tc => $tcg, tcn => $tcn, ecg => $ecg_g, rr => $rrm, n => $sel->nelem };
    }

    if ($have_fig) {
        my $base = $file =~ s/\.edf$//ir;
        my ($fig, $p, $q, $c) = subplots(3, 1, width => 900, height => 1050);
        my @cc = ('#1565C0', '#2E7D32', '#C62828', '#6A1B9A', '#EF6C00');
        for my $g (0 .. $#grp) {
            my $col = $cc[$g % @cc];
            my $lab = sprintf('RR %.0f ms (n=%d)', $grp[$g]{rr}, $grp[$g]{n});
            $p->line($tms, $grp[$g]{tc},  color => $col, lw => 1.3, label => $lab);
            $p->axvline($grp[$g]{rr}, color => $col, lw => 0.8);
            $q->line($tms, $grp[$g]{tcn}, color => $col, lw => 1.3, label => $lab);
            $q->axvline(-$grp[$g]{rr}, color => $col, lw => 0.8);
            $c->line($tms, $grp[$g]{ecg}, color => $col, lw => 1.0, label => $lab);
        }
        $p->set_title('slow-wave weight, aligned to R (vertical lines: median RR of each group)');
        $q->set_title('slow-wave weight, aligned to the next R (vertical lines: minus median RR)');
        $c->set_title("$bn1-$bn2 (ECG), aligned to R");
        for my $ax ($p, $q, $c) {
            $ax->axvline(0, color => '#999999', lw => 0.5);
            $ax->axhline(0, color => '#999999', lw => 0.5);
            $ax->set_xlim($tms->at(0), $tms->at(-1));
            $ax->legend;
        }
        $p->set_ylabel('uV'); $q->set_ylabel('uV'); $c->set_ylabel('uV');
        $c->set_xlabel('time [ms]');
        $fig->save("${base}_slowwave_rr.png");

        my $fig2 = figure(width => 760, height => 330);
        my $gs   = $fig2->add_gridspec(1, 2, left => 0.02, right => 0.98, top => 0.86, bottom => 0.02, wspace => 0.10);
        my $amp  = $tc->slice("$t0:$t1")->avg;
        for my $pan ([0, $sw * $amp, [ @head ], 'all head channels'],
                     [1, $swv * (($Av x $swv->dummy(0))->flat->slice("$t0:$t1")->avg), [ @head[@vi] ],
                      'without ' . join(',', grep { $vess{$_} } @head)]) {
            my ($col, $vals, $lab, $ttl) = @$pan;
            my $ax = $fig2->add_subplot($gs->at(0, $col));
            PDL::EEG::MAP2D::plot_topomap(values => $vals, labels => [ map { $alias{$_} // $_ } @$lab ],
                montage => $montage, clim => $vals->abs->max->sclr, contours => 0, names => 1,
                name_size => 6, title_size => 9, title => sprintf('%s, +-%.1f uV', $ttl, $vals->abs->max->sclr),
                _ax => $ax, _fig => $fig2);
        }
        $fig2->suptitle(sprintf('%s: slow-wave map, weight x mean 250..450 ms (red positive)', $file =~ s{.*/}{}r),
                        fontsize => 10);
        $fig2->save("${base}_slowwave_maps.png");
        print "wrote ${base}_slowwave_rr.png, ${base}_slowwave_maps.png\n";
    }
    return { map { $head[$_] => $sw->at($_) } 0 .. $nh - 1 };
}

# beats whose head-channel peak-to-peak over the epoch stays within --reject
sub _accept {
    my ($Y, $rp, $n0, $n1) = @_;
    my @k;
    for (my $i = 0; $i < $rp->nelem; $i += $chunk) {
        my $j  = _min($i + $chunk, $rp->nelem) - 1;
        my $Em = beat_epochs($Y, $rp->slice("$i:$j"), $n0, $n1)->mv(1, 0);   # (nep, nh, nb)
        push @k, (($Em->maximum - $Em->minimum)->maximum <= $reject)->list;  # per beat
    }
    return pdl(byte, \@k);
}

# average of baseline-corrected epochs, accumulated in chunks of beats
sub _mean {
    my ($Y, $rp, $n0, $n1, $bw) = @_;
    my $sum = zeroes(double, $Y->dim(0), $n1 - $n0);
    for (my $i = 0; $i < $rp->nelem; $i += $chunk) {
        my $j = _min($i + $chunk, $rp->nelem) - 1;
        my $E = beat_epochs($Y, $rp->slice("$i:$j"), $n0, $n1)->double;
        $E -= $E->slice(":,$bw->[0]:$bw->[1]")->mv(1, 0)->average->dummy(1);
        $sum += $E->mv(2, 0)->sumover;
    }
    return $sum / $rp->nelem;
}

# slow-wave map: Late window with the first two QRS maps projected out, first PC
sub _slow_map {
    my ($A, $wi, $i300) = @_;
    my $Q  = (_pcs($A->slice(":,$wi->{QRS}[0]:$wi->{QRS}[1]")))[1]->slice("0:1,:");
    my $L  = $A->slice(":,$wi->{Late}[0]:$wi->{Late}[1]");
    $L     = $L - ($L x $Q) x $Q->transpose;
    my $v  = (_pcs($L))[1]->slice("(0),:")->copy;
    $v    *= -1 if ($A->slice(":,($i300)") * $v)->sum->sclr < 0;
    return $v;
}
sub _pcs {
    my $Z = shift;
    my ($v, $V) = eigh_sym($Z->transpose x $Z);
    my $o = sequence(long, $v->nelem)->slice("-1:0");
    return ($v->index($o), $V->dice_axis(0, $o));
}
sub _cos_vec { my ($x, $y) = @_; (($x * $y)->sum / sqrt(($x**2)->sum * ($y**2)->sum))->sclr }
# cosine between two maps given as {label => value}, over the shared labels,
# each re-referenced to the average of those labels
sub _cos_maps {
    my ($m1, $m2) = @_;
    my @sh = grep { exists $m2->{$_} } grep { exists $m1->{$_} } @ALL;
    my $x = pdl(@{$m1}{@sh}); $x -= $x->avg;
    my $y = pdl(@{$m2}{@sh}); $y -= $y->avg;
    return _cos_vec($x, $y);
}
sub _min { $_[0] < $_[1] ? $_[0] : $_[1] }

# unit vectors (3, n) of the labels from the .elc, after fitting a sphere
sub _unit_positions {
    my $lab = shift;
    my $mon = read_elc($montage);
    my @miss = grep { !exists $mon->{pos}{$_} } @$lab;
    die "no position for @miss in $montage\n" if @miss;
    my $X = pdl(double, [ map { [ @{ $mon->{pos}{$_} }{qw(x y z)} ] } @$lab ]);   # (3, n)
    # sphere fit |x|^2 = 2 c.x + (r^2 - |c|^2): least squares in (c, k)
    my $M = $X->glue(0, ones(1, $X->dim(1))) ;                                    # (4, n)
    $M->slice("0:2,:") *= 2;
    my $y = ($X**2)->sumover;                                                     # (n)
    my $sol = PDL::MatrixOps::inv($M->transpose x $M) x ($M->transpose x $y->dummy(0));
    my $c = $sol->flat->slice("0:2");
    my $U = $X - $c->dummy(1);
    return $U / sqrt(($U**2)->sumover)->dummy(0);
}

# leave-one-out spherical spline (Perrin 1989, m = 4, 50 terms): matrix T with
# predictions = T x v, prediction of channel i from all the other channels
sub _loo_matrix {
    my $U  = shift;                                          # (3, n)
    my $n  = $U->dim(1);
    my $cs = ($U->dummy(1) * $U->dummy(2))->sumover->clip(-1, 1);   # (n, n) cosines
    my $G  = zeroes(double, $n, $n);
    my ($p0, $p1) = (ones(double, $n, $n), $cs->copy);
    for my $k (1 .. 50) {
        $G += (2 * $k + 1) / ($k * ($k + 1))**4 * $p1;
        my $p2 = ((2 * $k + 1) * $cs * $p1 - $k * $p0) / ($k + 1);
        ($p0, $p1) = ($p1, $p2);
    }
    $G /= 4 * 3.141592653589793;
    my $T = zeroes(double, $n, $n);
    for my $i (0 .. $n - 1) {
        my $o  = pdl(long, [ grep { $_ != $i } 0 .. $n - 1 ]);
        my $m  = $o->nelem;
        my $Go = $G->dice_axis(0, $o)->dice_axis(1, $o);
        my $L  = zeroes(double, $m + 1, $m + 1);
        $L->slice("0:" . ($m - 1) . ",0:" . ($m - 1)) .= $Go;
        $L->slice("($m),0:" . ($m - 1)) .= 1;
        $L->slice("0:" . ($m - 1) . ",($m)") .= 1;
        my $Li = PDL::MatrixOps::inv($L);                   # (c, c0) = Li x (v_o, 0)
        my $g  = $G->slice("$i,:")->flat->index($o);        # G(i, o)
        my $row = $g->append(pdl(1)) x $Li;                  # 1 x (m+1): prediction weights
        $T->slice(":,($i)")->index($o) .= $row->flat->slice("0:" . ($m - 1));
    }
    return $T;
}

# mean RR of the beats whose R lies within +-$w/2 samples of each beat's R
# (R times ascending; two pointers over Perl arrays, cumulative sums)
sub _runmean_rr {
    my ($rp, $rr, $w) = @_;
    my @t = $rp->list;
    my @c = (0); push @c, $c[-1] + $_ for $rr->list;
    my ($lo, $hi) = (0, 0);
    my @m;
    for my $i (0 .. $#t) {
        $lo++ while $t[$lo] < $t[$i] - $w / 2;
        $hi++ while $hi + 1 <= $#t && $t[$hi + 1] <= $t[$i] + $w / 2;
        push @m, ($c[$hi + 1] - $c[$lo]) / ($hi - $lo + 1);
    }
    return pdl(double, \@m);
}
