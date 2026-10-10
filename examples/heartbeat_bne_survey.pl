#!/usr/bin/env perl
# examples/heartbeat_bne_survey.pl
#
#   perl examples/heartbeat_bne_survey.pl REC.edf [options]
#
#     --bn1 L --bn2 L   labels of the two balanced non-cephalic electrodes,
#                       each recorded against the system reference (BN1, BN2)
#     --skip S          ignore the first S seconds                  (10)
#     --band LO,HI      zero-phase Butterworth band, order 4 [Hz]   (0.5,45)
#     --reject UV       per-channel peak-to-peak limit for scalp channels (150)
#     --nsurr N         surrogate (non-heartbeat) realisations       (100)
#     --seed N          seed for the surrogate offsets               (1)
#     --out PREFIX      output prefix               (REC without .edf + _hb)
#     --montage ELC     ASA .elc for the topomaps   (t/data/standard_1020_subset.elc)
#     --no-fig          text and TSV only
#
# R peaks are detected on BN1-BN2 (a body-surface lead that contains no head
# electrode). Every channel is then averaged around R, and the script reports
#
#   1. P / R / T amplitude of every channel at the latencies where BN1-BN2 peaks
#   2. the weight r for which  r*BN1 + (1-r)*BN2  is smallest, separately for
#      the P, QRS and T windows, with the residual each r leaves in each window
#   3. for six references of the scalp channels (system, average, A1+A2,
#      BN r=0.5, BN r=r_QRS, BN r=r_T): the R-locked scalp amplitude in the QRS,
#      T and three later windows, next to the same measure at non-heartbeat
#      time points
#
# Outputs: PREFIX_avg_<ref>.tsv, PREFIX_bn.png, PREFIX_scalp.png,
#          PREFIX_topo_<sys|avg>_<P|R|T>.png

use strict;
use warnings;
use PDL;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use PDL::EEG::IO::EDF   qw(read_edf);
use PDL::EEG::Heartbeat qw(:all);

my ($bn1, $bn2, $skip, $band, $reject, $nsurr, $seed) = ('BN1', 'BN2', 10, '0.5,45', 150, 100, 1);
my ($out, $montage, $nofig);
GetOptions('bn1=s' => \$bn1, 'bn2=s' => \$bn2, 'skip=f' => \$skip, 'band=s' => \$band,
           'reject=f' => \$reject, 'nsurr=i' => \$nsurr, 'seed=i' => \$seed,
           'out=s' => \$out, 'montage=s' => \$montage, 'no-fig' => \$nofig) or die;
my $file = shift @ARGV or die "usage: $0 REC.edf [options]\n";
($out = $file) =~ s/\.edf$//i, $out .= '_hb' unless defined $out;
$montage //= "$FindBin::Bin/../t/data/standard_1020_subset.elc";
my ($f_lo, $f_hi) = split /,/, $band;

# ---- read, select, filter --------------------------------------------------
my $e  = read_edf($file);
my $fs = $e->{fs};
my %ix; @ix{ @{ $e->{labels} } } = 0 .. $#{ $e->{labels} };
my @want  = qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 C3 Cz C4 T4 T5 P3 Pz P4 T6 O1 O2);
my %alt   = (T3 => 'T7', T4 => 'T8', T5 => 'P7', T6 => 'P8');
my @scalp = grep { exists $ix{$_} } map { exists $ix{$_} ? $_ : ($alt{$_} // $_) } @want;
my @ear   = grep { exists $ix{$_} } qw(A1 A2);
exists $ix{$_} or die "channel '$_' not in $file\n" for $bn1, $bn2;
my @ch = (@scalp, @ear, $bn1, $bn2);
my %j;  @j{@ch} = 0 .. $#ch;
my $D  = zp_butter($e->{data}->dice_axis(0, pdl(long, [ @ix{@ch} ])),
                   fs => $fs, low => $f_lo, high => $f_hi, order => 4);   # (nch, nt)
my $nt = $D->dim(1);
my $ns = scalar @scalp;
my $S  = $D->slice("0:" . ($ns - 1));
my ($B1, $B2) = map { $D->slice("($j{$_})") } $bn1, $bn2;
my $ecg = $B1 - $B2;

# ---- R peaks ---------------------------------------------------------------
my ($tmin, $tmax, @base) = (-0.3, 0.8, -0.30, -0.22);
my ($r, $sign) = detect_r($ecg, $fs);
my ($n0, $n1) = (int($tmin * $fs - 0.5), int($tmax * $fs + 0.5));
my ($lo, $hi) = (int($skip * $fs) - $n0, $nt - $n1);
$r = $r->where(($r >= $lo) & ($r < $hi));
my $rr = ($r->slice("1:-1") - $r->slice("0:-2"))->double / $fs;
printf "%s\n  %s, %d Hz, %.1f s, band %s-%s Hz\n", $file, $e->{recording}, $fs, $nt / $fs, $f_lo, $f_hi;
printf "  R peaks on %s-%s: %d beats (polarity %s), RR median %.3f s, range %.3f-%.3f s\n",
    $bn1, $bn2, $r->nelem, ($sign > 0 ? '+' : '-'), $rr->median, $rr->min, $rr->max;

my %avgopt = (tmin => $tmin, tmax => $tmax, baseline => \@base);
my $t;                                                    # latency axis [ms]
sub win { my ($a, $b) = @_; which(($t >= $a) & ($t <= $b)) }

# ---- 1. P / R / T in every channel, system reference -------------------------
my $lim = ones(scalar(@ch) + 1) * 1e30;                   # no limit on BN / ear
$lim->slice("0:" . ($ns - 1)) .= $reject;
my $A = beat_average($D->glue(0, $ecg->dummy(0)), $r, $fs, %avgopt, reject => $lim);
$t = $A->{t_ms};
my $avg = $A->{avg};
my $se  = $A->{sd} / sqrt($A->{n}->dummy(1));
my @lab = (@ch, "$bn1-$bn2");
my $lead = $avg->slice("(-1)");
my %lat;
for (['P', -250, -70], ['R', -40, 60], ['T', 150, 450]) {
    my ($w, $a, $b) = @$_;
    my $k = win($a, $b);
    $lat{$w} = $k->at($lead->index($k)->abs->maximum_ind);
}
printf "\n1. Amplitude [uV] at the latencies where %s-%s peaks: P %d ms, R %d ms, T %d ms\n",
    $bn1, $bn2, map { $t->at($lat{$_}) } qw(P R T);
print "   system reference, baseline $base[0]..$base[1] s; +-SE over beats\n";
printf "   %-9s %15s %15s %15s %6s\n", 'channel', 'P', 'R', 'T', 'beats';
for my $c (0 .. $#lab) {
    printf "   %-9s %s %6d\n", $lab[$c],
        join(' ', map { sprintf '%8.2f +-%4.2f', $avg->at($c, $lat{$_}), $se->at($c, $lat{$_}) } qw(P R T)),
        $A->{n}->at($c);
}

# ---- 2. BN balance per wave ------------------------------------------------
my %W = (P   => win($t->at($lat{P}) - 40, $t->at($lat{P}) + 40),
         QRS => win(-40, 60),
         T   => win($t->at($lat{T}) - 80, $t->at($lat{T}) + 80));
my ($a1, $a2) = map { $avg->slice("($j{$_})") } $bn1, $bn2;
my $rms = sub { my $v = shift; sqrt(sum($v * $v) / $v->nelem) };
my %rw;
print "\n2. Weight r minimising r*$bn1 + (1-r)*$bn2 in each window of the R-locked average\n";
printf "   %-4s %14s %7s %10s %10s %10s   %s\n", 'wave', 'window [ms]', 'r',
    "rms $bn1", "rms $bn2", 'rms mix', 'r from odd / even beats';
my @half = map { my $k = $_;
    beat_average($D, $r->index(which(sequence($r->nelem) % 2 == $k)), $fs, %avgopt)->{avg} } 0, 1;
for my $w (qw(P QRS T)) {
    my $k = $W{$w};
    my ($rv, $m, $m1, $m2) = bn_balance($a1->index($k), $a2->index($k));
    $rw{$w} = $rv;
    my @h = map { (bn_balance($_->slice("($j{$bn1})")->index($k), $_->slice("($j{$bn2})")->index($k)))[0] } @half;
    printf "   %-4s %6d..%-6d %7.3f %10.2f %10.2f %10.2f   %.3f / %.3f\n",
        $w, $t->at($k->at(0)), $t->at($k->at(-1)), $rv, $m1, $m2, $m, @h;
}
print "\n   RMS [uV] left in each window by each choice of r\n";
printf "   %-12s %8s %8s %8s\n", 'r', 'P', 'QRS', 'T';
for (['0.5', 0.5], map { ["r_$_", $rw{$_}] } qw(P QRS T)) {
    my ($nm, $rv) = @$_;
    printf "   %-5s %6.3f %s\n", $nm, $rv,
        join(' ', map { sprintf '%8.2f', $rms->(($rv * $a1 + (1 - $rv) * $a2)->index($W{$_})) } qw(P QRS T));
}

# ---- 3. references of the scalp channels -----------------------------------
my $bn = sub { my $w = shift; $w * $B1 + (1 - $w) * $B2 };
my @ref = (['system', 'sys', undef],
           ['average', 'avg', $S->average]);
push @ref, ['A1+A2', 'ear', ($D->slice("($j{A1})") + $D->slice("($j{A2})")) / 2] if @ear == 2;
push @ref, ['BN r=0.5', 'bn50', $bn->(0.5)],
           [sprintf('BN r=%.3f (QRS)', $rw{QRS}), 'bnqrs', $bn->($rw{QRS})],
           [sprintf('BN r=%.3f (T)',   $rw{T}),   'bnt',   $bn->($rw{T})];
my @meas = (['QRS p-p', -40, 60, 'pp'], ['T', $t->at($lat{T}) - 40, $t->at($lat{T}) + 40, 'mean'],
            ['200-300', 200, 300, 'mean'], ['300-450', 300, 450, 'mean'], ['450-600', 450, 600, 'mean']);
sub measures {                       # (nscalp, nep) -> (nmeas) RMS over channels, (nmeas, nscalp)
    my $a = shift;
    my @v;
    for my $m (@meas) {
        my $seg = $a->dice_axis(1, win($m->[1], $m->[2]))->mv(1, 0);    # (nwin, nscalp)
        push @v, $m->[3] eq 'pp' ? $seg->maximum - $seg->minimum : $seg->average;
    }
    my $v = cat(@v);                                               # (nscalp, nmeas)
    return (sqrt(($v * $v)->average), $v);
}
my (%refavg, %refval);
print "\n3. R-locked scalp amplitude under each reference: RMS over the $ns scalp channels [uV]\n";
printf "   QRS p-p = max-min in -40..60 ms; other columns = mean amplitude in the window [ms]\n";
printf "   second line: the same at %d sets of non-heartbeat time points (R +- up to 0.5 s), mean and 95th percentile\n", $nsurr;
printf "   %-18s %s\n", 'reference', join(' ', map { sprintf '%13s', $_->[0] } @meas);
for my $rf (@ref) {
    my ($name, $tag, $sig) = @$rf;
    my $Xr = defined $sig ? $S - $sig->dummy(0) : $S;
    my $Ar = beat_average($Xr, $r, $fs, %avgopt, reject => $reject);
    my ($m, $v) = measures($Ar->{avg});
    $refavg{$tag} = $Ar; $refval{$tag} = $v;
    printf "   %-18s %s\n", $name, join(' ', map { sprintf '%13.2f', $_ } $m->list);
    if ($nsurr > 0) {
        my @sm;
        for my $k (1 .. $nsurr) {
            my $sr = surrogate_r($r, $nt, $fs, min => 0, max => 0.5, lo => $lo, hi => $hi,
                                 seed => $seed * 100000 + $k);
            push @sm, (measures(beat_average($Xr, $sr, $fs, %avgopt, reject => $reject)->{avg}))[0];
        }
        my $sm = cat(@sm);                                         # (nmeas, nsurr)
        printf "   %-18s %s\n", '  non-heartbeat', join(' ', map {
            sprintf '%5.2f (%5.2f)', $sm->slice("($_)")->avg, $sm->slice("($_)")->pct(0.95)->sclr } 0 .. $#meas);
    }
    open my $fh, '>', "${out}_avg_$tag.tsv" or die "${out}_avg_$tag.tsv: $!";
    print $fh join("\t", 'time_ms', @scalp), "\n";
    my $M = $Ar->{avg};
    for my $i (0 .. $t->nelem - 1) {
        print $fh join("\t", $t->at($i), map { sprintf '%.4f', $M->at($_, $i) } 0 .. $ns - 1), "\n";
    }
    close $fh;
}
print "\n   per channel, QRS p-p [uV]\n";
printf "   %-18s %s\n", 'reference', join(' ', map { sprintf '%5s', $_ } @scalp);
printf "   %-18s %s\n", $_->[0], join(' ', map { sprintf '%5.1f', $_ } $refval{ $_->[1] }->slice(":,(0)")->list) for @ref;
print "\n   per channel, mean 300-450 ms [uV]\n";
printf "   %-18s %s\n", 'reference', join(' ', map { sprintf '%5s', $_ } @scalp);
printf "   %-18s %s\n", $_->[0], join(' ', map { sprintf '%5.2f', $_ } $refval{ $_->[1] }->slice(":,(3)")->list) for @ref;
{
    open my $fh, '>', "${out}_avg_bn.tsv" or die $!;
    print $fh join("\t", 'time_ms', @lab), "\n";
    for my $i (0 .. $t->nelem - 1) {
        print $fh join("\t", $t->at($i), map { sprintf '%.4f', $avg->at($_, $i) } 0 .. $#lab), "\n";
    }
}
print "\nwrote ${out}_avg_{bn," . join(',', map { $_->[1] } @ref) . "}.tsv\n";
exit 0 if $nofig;

# ---- figures ---------------------------------------------------------------
eval { require PDL::Graphics::Cairo; require PDL::EEG::MAP2D; 1 }
    or do { print "PDL::Graphics::Cairo not available: no figures\n"; exit 0 };
PDL::Graphics::Cairo->import(qw(subplots));
my @col = ('#1565C0', '#C62828', '#2E7D32', '#6A1B9A', '#00838F', '#37474F');

{   # BN1, BN2, their difference, and the mixes (negative up)
    my ($fig, $p, $q) = subplots(2, 1, width => 900, height => 760);
    $p->line($t, $a1,   color => $col[0], lw => 1.4, label => $bn1);
    $p->line($t, $a2,   color => $col[1], lw => 1.4, label => $bn2);
    $p->line($t, $lead, color => $col[5], lw => 1.0, label => "$bn1-$bn2");
    $p->set_title("R-locked average, system reference (" . $r->nelem . " beats)");
    my @mix = (['r=0.5', 0.5], map { [sprintf('r=%.3f (%s)', $rw{$_}, $_), $rw{$_}] } qw(P QRS T));
    $q->line($t, $mix[$_][1] * $a1 + (1 - $mix[$_][1]) * $a2, color => $col[$_ + 2 > 5 ? 0 : $_ + 2],
             lw => 1.3, label => $mix[$_][0]) for 0 .. $#mix;
    $q->set_title("r*$bn1 + (1-r)*$bn2");
    for my $ax ($p, $q) {
        $ax->axhline(0, color => '#999999', lw => 0.5);
        $ax->axvline(0, color => '#999999', lw => 0.5);
        $ax->invert_yaxis; $ax->set_ylabel('uV (negative up)'); $ax->legend;
    }
    $q->set_xlabel('time from R peak [ms]');
    $fig->save("${out}_bn.png");
}
{   # scalp channels under each reference, common scale
    my $nr = scalar @ref;
    my ($fig, @axs) = subplots($nr, 1, width => 900, height => 230 * $nr);
    my $ymax = 0;
    for (@ref) { my $m = $refavg{ $_->[1] }{avg}->abs->max; $ymax = $m if $m > $ymax }
    for my $i (0 .. $nr - 1) {
        my $ax = $axs[$i];
        my $M  = $refavg{ $ref[$i][1] }{avg};
        $ax->line($t, $M->slice("($_)"), color => '#1565C0', lw => 0.6) for 0 .. $ns - 1;
        $ax->axhline(0, color => '#999999', lw => 0.5);
        $ax->axvline(0, color => '#999999', lw => 0.5);
        $ax->set_ylim(-$ymax * 1.05, $ymax * 1.05); $ax->invert_yaxis;
        $ax->set_ylabel('uV'); $ax->set_title("reference: $ref[$i][0]");
    }
    $axs[-1]->set_xlabel('time from R peak [ms]  (negative up, ' . $ns . ' scalp channels)');
    $fig->save("${out}_scalp.png");
}
for my $rf (grep { $_->[1] =~ /^(?:sys|avg)$/ } @ref) {   # scalp maps at the P, R, T latencies
    my $M = $refavg{ $rf->[1] }{avg};
    for my $w (qw(P R T)) {
        PDL::EEG::MAP2D::plot_topomap(
            values => $M->slice(":,($lat{$w})"), labels => \@scalp, montage => $montage,
            contours => 6, names => 1,
            title => sprintf('%s %d ms, ref %s', $w, $t->at($lat{$w}), $rf->[0]),
            outfile => "${out}_topo_$rf->[1]_$w.png");
    }
}
print "wrote ${out}_bn.png, ${out}_scalp.png, ${out}_topo_<sys|avg>_<P|R|T>.png\n";
