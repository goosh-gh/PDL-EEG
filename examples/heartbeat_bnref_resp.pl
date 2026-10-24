#!/usr/bin/env perl
# examples/heartbeat_bnref_resp.pl
#
#   perl examples/heartbeat_bnref_resp.pl REC.edf [options]
#
#     --bn1 L --bn2 L   non-cephalic electrodes                                (BN1, BN2)
#     --ratio R         BN = R*BN1 + (1-R)*BN2, the reference for the head      (0.65)
#     --skip S          ignore the first S seconds                              (10)
#     --band LO,HI      zero-phase Butterworth band of the head channels [Hz]   (0.5,40)
#     --exclude LIST    head channels left out                                  (none)
#     --reject UV       a beat is dropped when a head channel's peak-to-peak
#                       over -400..+700 ms exceeds this                          (200)
#     --resp KIND       per-beat respiratory variable: 'template' (the breath
#                       shape of phase_template at each beat's phase) or 'rsize'
#                       (the R size itself)                                     (template)
#     --chunk N         beats per chunk when accumulating (memory)              (500)
#     --out PREFIX      output prefix (REC without .edf + _bnresp)
#     --no-fig          text only
#
# R-locked averages of the head channels (19 scalp + A1, A2) re-referenced to the
# balanced non-cephalic electrode BN = R*BN1 + (1-R)*BN2, with the respiratory
# variation of every beat taken as a continuous variable.
#
# Respiration (from the ECG only): BN1-BN2 band-passed 0.5-45 Hz, R peaks
# (detect_r), per-beat features (PDL::EEG::Heartbeat::BeatFeatures); the R size
# (R amplitude against its running mean, %) gives the phase (resp_phase 0.15-0.5 Hz,
# then phase_template, window 5 s, K = 2, 3 rounds, as heartbeat_resp_phase.pl).
# rho = the template (the shape of one breath) at each beat's phase, or the R
# size itself; scaled to mean 0, SD 1 over the beats used.
#
# For every head channel and lag, the beats are fitted as
#     E_k(t) = a(t) + b(t) rho_k + e_k(t)
# a(t) is the average at rho = 0, b(t) the change per SD of rho [uV/SD], with its
# standard error from the residuals.  The same fit is made for the weights of two
# fixed maps of the BN-referenced average: the first principal map of the QRS
# window (-45..55 ms) and the slow-wave map (first principal map of 450..650 ms
# with the first two QRS maps projected out).
#
# Text: the respiratory variable; per channel the average at the R peak and over
# 250..450 ms in the BN reference and as recorded; the QRS principal maps; the
# slow-wave map; the fitted b with its t value for the QRS-map weight at R, the
# slow-wave weight over 250..450 ms, and Fz, Cz, Pz over 250..450 ms.
# Figure PREFIX.png: averages of some channels with the ECG; b(t) +- 2 SE of
# the map weights; b(t) +- 2 SE of Fz, Cz, Pz.

use strict;
use warnings;
use PDL;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use PDL::EEG::IO::EDF                 qw(read_edf);
use PDL::EEG::Filter                  qw(zp_butter);
use PDL::EEG::Heartbeat               qw(detect_r beat_epochs);
use PDL::EEG::Heartbeat::BeatFeatures qw(beat_features);
use PDL::EEG::Heartbeat::Respiration  qw(resp_phase phase_at phase_template phase_agreement fourier_eval);
use PDL::EEG::Linalg                  qw(eigh_sym);

my ($bn1, $bn2, $ratio, $skip, $band, $exclude, $reject, $resp, $chunk, $out, $nofig)
    = ('BN1', 'BN2', 0.65, 10, '0.5,40', '', 200, 'template', 500);
GetOptions('bn1=s' => \$bn1, 'bn2=s' => \$bn2, 'ratio=f' => \$ratio, 'skip=f' => \$skip,
           'band=s' => \$band, 'exclude=s' => \$exclude, 'reject=f' => \$reject, 'resp=s' => \$resp,
           'chunk=i' => \$chunk, 'out=s' => \$out, 'no-fig' => \$nofig) or die;
my $file = shift @ARGV or die "usage: $0 REC.edf [options]\n";
($out = $file) =~ s/\.edf$//i, $out .= '_bnresp' unless defined $out;
$resp =~ /^(template|rsize)$/ or die "--resp must be template or rsize\n";
my ($f_lo, $f_hi) = split /,/, $band;
my %excl = map { $_ => 1 } split /,/, $exclude;

my $e  = read_edf($file);
my $fs = $e->{fs};
my %ix; @ix{ @{ $e->{labels} } } = 0 .. $#{ $e->{labels} };
exists $ix{$_} or die "channel '$_' not in $file\n" for $bn1, $bn2;
my @head = grep { exists $ix{$_} && !$excl{$_} }
           qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 C3 Cz C4 T4 T5 P3 Pz P4 T6 O1 O2 A1 A2);
my $nh = @head;
my %hp; @hp{@head} = 0 .. $nh - 1;

# ---- respiration from BN1-BN2 ------------------------------------------------
my $raw = $e->{data};
my $ecg = zp_butter($raw->slice("($ix{$bn1})")->double - $raw->slice("($ix{$bn2})")->double,
                    fs => $fs, low => 0.5, high => 45, order => 4);
my ($rall, $sign) = detect_r($ecg, $fs);
$ecg = $ecg * $sign;
$rall = $rall->where(($rall >= $fs) & ($rall + $fs < $ecg->nelem));
my $B  = beat_features($ecg, $fs, r => $rall);
my $u  = $B->{use};
my $tk = $B->{tb}->index($u);
my $ysz = $B->{F}{R}->index($u);
my $pH = phase_at(resp_phase($tk, $ysz, band => [0.15, 0.5]), $tk);
my $PT = phase_template($tk, $ysz, $pH, win => 5, K => 2, iter => 3);
my $pL = $PT->{phi};
my $q  = phase_agreement($tk, $pH, $pL, win => 120, min => 40) * (180 / 3.141592653589793);
my $rho = $resp eq 'template' ? fourier_eval($PT->{c}, $pL)->flat : $ysz->copy;
my $rb  = $B->{r}->index($u);                               # R samples of these beats

# ---- head channels in the BN reference and as recorded -------------------------
my $Hrec = zp_butter($raw->dice_axis(0, pdl(long, [ @ix{@head} ]))->double,
                     fs => $fs, low => $f_lo, high => $f_hi, order => 4);
my $bnr  = zp_butter($ratio * $raw->slice("($ix{$bn1})")->double + (1 - $ratio) * $raw->slice("($ix{$bn2})")->double,
                     fs => $fs, low => $f_lo, high => $f_hi, order => 4);
my $Hbn  = $Hrec - $bnr->dummy(0);
my $nt   = $Hrec->dim(1);
my $s    = sub { my $v = shift() * 1e-3 * $fs; $v >= 0 ? int($v + 0.5) : -int(-$v + 0.5) };
my ($n0, $n1) = ($s->(-400), $s->(700));
my $nep = $n1 - $n0;
my %wi = map { my ($nm, $a0, $a1) = @$_; $nm => [ $s->($a0) - $n0, $s->($a1) - $n0 - 1 ] }
         (['base', -75, -50], ['QRS', -45, 55], ['T', 250, 450], ['Late', 450, 650]);
my $iR  = $s->(0) - $n0;
my $tms = (sequence($nep) + $n0) * (1000 / $fs);

my $fit_ok = ($rb >= $s->(1000 * $skip) - $n0) & ($rb + $n1 <= $nt) & ($rb + $n0 >= 0);
$rb  = $rb->where($fit_ok);
$rho = $rho->where($fit_ok);
$q   = $q->where($fit_ok);
my $keep = _accept($Hbn, $rb);
$rb  = $rb->where($keep);
$rho = $rho->where($keep);
$q   = $q->where($keep);
$rho = $rho - $rho->avg;
$rho = $rho / sqrt(($rho**2)->avg);
my $qf = $q->where($q->isfinite);
printf "%s: fs %g Hz, %d head channels, BN = %.2f*%s + %.2f*%s; %d beats used (reject %g uV)\n",
    $file, $fs, $nh, $ratio, $bn1, 1 - $ratio, $bn2, $rb->nelem, $reject;
printf "respiratory variable: %s; q_agree median %s deg; corr(rho, R size) %.2f\n",
    $resp eq 'template' ? 'breath template at phase_L' : 'R size',
    $qf->nelem ? sprintf('%.1f', $qf->median) : 'nan',
    _corr($rho, $ysz->where($fit_ok)->where($keep));

# ---- averages ---------------------------------------------------------------
my $Abn  = _mean($Hbn, $rb);
my $Arec = _mean($Hrec, $rb);
my ($t0, $t1) = @{ $wi{T} };
print "\nper channel [uV]:     at R peak (BN ref, as recorded)   mean 250..450 ms (BN ref, as recorded)\n";
for my $c (0 .. $nh - 1) {
    printf "  %-4s %12.2f %12.2f %18.2f %12.2f\n", $head[$c], $Abn->at($c, $iR), $Arec->at($c, $iR),
        $Abn->slice("($c),$t0:$t1")->avg, $Arec->slice("($c),$t0:$t1")->avg;
}

# ---- maps of the BN-referenced average ------------------------------------------
my ($vq, $VQ) = _pcs($Abn->slice(":,$wi{QRS}[0]:$wi{QRS}[1]"));
my $cum = $vq->cumusumover / $vq->sum;
printf "\nQRS window (BN ref): share of power in the first 1..4 principal maps: %s\n",
    join(' ', map { sprintf '%.3f', $cum->at($_) } 0 .. 3);
my $mq = $VQ->slice("(0),:")->copy;
$mq *= -1 if ($Abn->slice(":,($iR)") * $mq)->sum->sclr < 0;
my $Q2 = $VQ->slice("0:1,:");
my $L  = $Abn->slice(":,$wi{Late}[0]:$wi{Late}[1]");
$L     = $L - ($L x $Q2) x $Q2->transpose;
my $ms = (_pcs($L))[1]->slice("(0),:")->copy;
$ms *= -1 if ($Abn->slice(":,(" . ($s->(300) - $n0) . ")") * $ms)->sum->sclr < 0;
print "QRS map 1:  ", join(' ', map { sprintf '%s %.2f', $head[$_], $mq->at($_) } 0 .. $nh - 1), "\n";
print "slow map:   ", join(' ', map { sprintf '%s %.2f', $head[$_], $ms->at($_) } 0 .. $nh - 1), "\n";

# ---- regression on rho -----------------------------------------------------------
my @sel  = grep { exists $hp{$_} } qw(Fz Cz Pz);
my $Wm   = cat($mq, $ms)->transpose;                         # nh x 2 (columns = maps)
my $reg  = _regress($Hbn, $rb, $rho, $Wm, [ @hp{@sel} ]);     # rows: QRS, slow, @sel
my @rows = ('QRS map 1', 'slow map', @sel);
print "\nfit E = a + b rho (rho in SD):\n";
printf "  %-10s at R peak: a %8.2f  b %7.2f uV/SD (%+.1f%% of a), t %6.1f\n", $rows[0],
    $reg->{a}->at(0, $iR), $reg->{b}->at(0, $iR), 100 * $reg->{b}->at(0, $iR) / $reg->{a}->at(0, $iR),
    $reg->{b}->at(0, $iR) / $reg->{se}->at(0, $iR);
for my $k (1 .. $#rows) {
    my $am = $reg->{a}->slice("($k),$t0:$t1")->avg;
    my $bm = $reg->{b}->slice("($k),$t0:$t1")->avg;
    my $sm = _se_window($reg, $k, $t0, $t1);
    printf "  %-10s mean 250..450 ms: a %6.2f  b %6.2f uV/SD, t %5.1f (SE from the per-beat window means)\n",
        $rows[$k], $am, $bm, $bm / $sm;
}

if (!$nofig && eval { require PDL::Graphics::Cairo; 1 }) {
    PDL::Graphics::Cairo->import(qw(subplots));
    my ($fig, $p, $g, $h) = subplots(3, 1, width => 900, height => 1050);
    my @cc = ('#1565C0', '#C62828', '#2E7D32', '#6A1B9A', '#EF6C00', '#00838F');
    my @show = grep { exists $hp{$_} } qw(Fz Cz Pz Fp1 O1 T3);
    $p->line($tms, $Abn->slice("($hp{$show[$_]}),:"), color => $cc[$_ % @cc], lw => 1.2, label => $show[$_])
        for 0 .. $#show;
    $p->line($tms, _mean($ecg->dummy(0), $rb)->flat / 20, color => '#000000', lw => 0.8, label => "($bn1-$bn2)/20");
    $p->set_title(sprintf('R-locked average, reference BN = %.2f %s + %.2f %s (%d beats)', $ratio, $bn1, 1 - $ratio, $bn2, $rb->nelem));
    $p->set_ylim(-15, 15);
    for my $k (0, 1) {
        my $col = $cc[$k];
        $g->line($tms, $reg->{b}->slice("($k),:"), color => $col, lw => 1.3, label => "b, $rows[$k] weight");
        $g->line($tms, $reg->{b}->slice("($k),:") + 2 * $reg->{se}->slice("($k),:"), color => $col, lw => 0.5);
        $g->line($tms, $reg->{b}->slice("($k),:") - 2 * $reg->{se}->slice("($k),:"), color => $col, lw => 0.5);
    }
    $g->set_title('change per SD of the respiratory variable (thin lines: +-2 SE)');
    for my $k (2 .. $#rows) {
        my $col = $cc[$k];
        $h->line($tms, $reg->{b}->slice("($k),:"), color => $col, lw => 1.3, label => "b, $rows[$k]");
        $h->line($tms, $reg->{b}->slice("($k),:") + 2 * $reg->{se}->slice("($k),:"), color => $col, lw => 0.5);
        $h->line($tms, $reg->{b}->slice("($k),:") - 2 * $reg->{se}->slice("($k),:"), color => $col, lw => 0.5);
    }
    $h->set_title('change per SD of the respiratory variable, single channels (BN ref)');
    $h->set_xlabel('time from R peak [ms]');
    for my $ax ($p, $g, $h) {
        $ax->axvline(0, color => '#999999', lw => 0.5);
        $ax->axhline(0, color => '#999999', lw => 0.5);
        $ax->set_xlim($tms->at(0), $tms->at(-1));
        $ax->set_ylabel('uV');
        $ax->legend;
    }
    $fig->save("${out}.png");
    print "wrote ${out}.png\n";
}

# ---------------------------------------------------------------------------
sub _accept {
    my ($Y, $rp) = @_;
    my @k;
    for (my $i = 0; $i < $rp->nelem; $i += $chunk) {
        my $j  = _min($i + $chunk, $rp->nelem) - 1;
        my $Em = beat_epochs($Y, $rp->slice("$i:$j"), $n0, $n1)->mv(1, 0);
        push @k, (($Em->maximum - $Em->minimum)->maximum <= $reject)->list;
    }
    return pdl(byte, \@k);
}
sub _epochs {
    my ($Y, $rp) = @_;
    my $E = beat_epochs($Y, $rp, $n0, $n1)->double;
    $E -= $E->slice(":,$wi{base}[0]:$wi{base}[1]")->mv(1, 0)->average->dummy(1);
    return $E;
}
sub _mean {
    my ($Y, $rp) = @_;
    my $sum = zeroes(double, $Y->dim(0), $nep);
    for (my $i = 0; $i < $rp->nelem; $i += $chunk) {
        my $j = _min($i + $chunk, $rp->nelem) - 1;
        $sum += _epochs($Y, $rp->slice("$i:$j"))->mv(2, 0)->sumover;
    }
    return $sum / $rp->nelem;
}
# least squares of y = a + b x per row and lag; rows are the map weights
# (columns of $Wm) and the head channels @$chs.  Also keeps, per row, the
# per-beat means over the T window for a window-level standard error.
sub _regress {
    my ($Y, $rp, $x, $Wm, $chs) = @_;
    my $nr = $Wm->dim(0) + @$chs;
    my ($Sy, $Sxy, $Syy) = map { zeroes(double, $nr, $nep) } 1 .. 3;
    my ($Sx, $Sxx, $n) = (0, 0, 0);
    my @wm;                                                    # per-beat T-window means (nr, nb)
    my $ci = pdl(long, $chs);
    for (my $i = 0; $i < $rp->nelem; $i += $chunk) {
        my $j  = _min($i + $chunk, $rp->nelem) - 1;
        my $E  = _epochs($Y, $rp->slice("$i:$j"));             # (nh, nep, nb)
        my $xs = $x->slice("$i:$j");
        my $Wt = $E x $Wm;                                     # (k, nep, nb): map weights
        $Wt    = $Wt->glue(0, $E->dice_axis(0, $ci)) if $ci->nelem;
        $Sy  += $Wt->mv(2, 0)->sumover;
        $Sxy += ($Wt * $xs->dummy(0)->dummy(0))->mv(2, 0)->sumover;
        $Syy += ($Wt**2)->mv(2, 0)->sumover;
        $Sx  += $xs->sum; $Sxx += ($xs**2)->sum; $n += $xs->nelem;
        push @wm, $Wt->slice(":,$t0:$t1,:")->mv(1, 0)->average;  # (nr, nb)
    }
    my $den = $n * $Sxx - $Sx**2;
    my $b   = ($n * $Sxy - $Sx * $Sy) / $den;
    my $a   = ($Sy - $b * $Sx) / $n;
    my $rss = $Syy - $a * $Sy - $b * $Sxy;
    my $se  = sqrt(($rss / ($n - 2))->clip(0, undef) / ($Sxx - $Sx**2 / $n));
    return { a => $a, b => $b, se => $se, wm => (@wm == 1 ? $wm[0] : $wm[0]->glue(1, @wm[1 .. $#wm])), x => $x };
}
# standard error of the slope of the per-beat window means of row $k
sub _se_window {
    my ($reg, $k) = @_;
    my $y = $reg->{wm}->slice("($k),:");
    my $x = $reg->{x};
    my $n = $x->nelem;
    my $xc = $x - $x->avg;
    my $bb = ($xc * $y)->sum / ($xc**2)->sum;
    my $res = $y - $y->avg - $bb * $xc;
    return sqrt((($res**2)->sum / ($n - 2)) / ($xc**2)->sum)->sclr;
}
sub _pcs {
    my $Z = shift;
    my ($v, $V) = eigh_sym($Z->transpose x $Z);
    my $o = sequence(long, $v->nelem)->slice("-1:0");
    return ($v->index($o), $V->dice_axis(0, $o));
}
sub _corr { my ($x, $y) = @_; $x = $x - $x->avg; $y = $y - $y->avg; (($x * $y)->sum / sqrt(($x**2)->sum * ($y**2)->sum))->sclr }
sub _min { $_[0] < $_[1] ? $_[0] : $_[1] }
