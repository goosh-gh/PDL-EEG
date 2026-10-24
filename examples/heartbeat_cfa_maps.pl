#!/usr/bin/env perl
# examples/heartbeat_cfa_maps.pl
#
#   perl examples/heartbeat_cfa_maps.pl REC.edf [options]
#
#     --bn1 L --bn2 L   non-cephalic electrodes; R peaks are detected on BN1-BN2
#     --skip S          ignore the first S seconds                         (10)
#     --band LO,HI      zero-phase Butterworth band, order 4 [Hz]          (0.5,40)
#     --exclude LIST    head channels to leave out, comma separated        (none)
#     --k K             number of QRS maps used to describe the cardiac field (2)
#     --montage ELC     ASA .elc with the electrode positions
#     --out PREFIX      output prefix                   (REC without .edf + _cfamaps)
#
# The cardiac field on the head as a few fixed maps with time-varying weights.
# Head channels (19 scalp + A1, A2), average reference, every beat, baseline
# -75..-50 ms.  First the slow wave (first principal component of the Late
# window, 450..650 ms, of the average) is projected out; its map is drawn too.
# Then the first K principal components m_1..m_K of the QRS window (-45..55 ms)
# give fixed maps, and at every lag t the average map V(t) is fitted as
#     V(t) ~ a_1(t) m_1 + ... + a_K(t) m_K        (a_k(t) = m_k' V(t))
#
# Figures:
#   PREFIX_maps.png   the K QRS maps and the slow-wave map (unit vectors, x the
#                     largest weight they reach), then the measured map and the
#                     K-map fit at the R peak, at the P-wave and T-wave maxima of
#                     a_1..a_K, and at 300 ms
#   PREFIX_time.png   the BN1-BN2 average (ECG) on the same time axis;
#                     a_1(t)..a_K(t) and the slow-wave weight; the share of the
#                     map power the K maps explain at every lag, with the same
#                     share for random K-dimensional subspaces (5th..95th
#                     percentile) for comparison; the loop a_1 against a_2 for
#                     the P, QRS and T stretches (when K >= 2)
# Text: the share explained in the P, QRS and T windows.

use strict;
use warnings;
use PDL;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use PDL::EEG::IO::EDF   qw(read_edf);
use PDL::EEG::Filter    qw(zp_butter);
use PDL::EEG::Heartbeat qw(detect_r beat_epochs);
use PDL::EEG::Linalg    qw(eigh_sym);
use PDL::EEG::SynthNoise qw(seed_rng);

my ($bn1, $bn2, $skip, $band, $exclude, $K, $montage, $out) = ('BN1', 'BN2', 10, '0.5,40', '', 2);
GetOptions('bn1=s' => \$bn1, 'bn2=s' => \$bn2, 'skip=f' => \$skip, 'band=s' => \$band,
           'exclude=s' => \$exclude, 'k=i' => \$K, 'montage=s' => \$montage, 'out=s' => \$out) or die;
my $file = shift @ARGV or die "usage: $0 REC.edf [options]\n";
($out = $file) =~ s/\.edf$//i, $out .= '_cfamaps' unless defined $out;
$montage //= "$FindBin::Bin/../t/data/standard_1020_subset.elc";
my ($f_lo, $f_hi) = split /,/, $band;

my $e  = read_edf($file);
my $fs = $e->{fs};
my %ix; @ix{ @{ $e->{labels} } } = 0 .. $#{ $e->{labels} };
my %excl = map { $_ => 1 } split /,/, $exclude;
my @head = grep { exists $ix{$_} && !$excl{$_} }
           qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 C3 Cz C4 T4 T5 P3 Pz P4 T6 O1 O2 A1 A2);
my %alias = (T3 => 'T7', T4 => 'T8', T5 => 'P7', T6 => 'P8');   # names in the .elc
my @mlab  = map { $alias{$_} // $_ } @head;
my $nh = @head;
my $X  = zp_butter($e->{data}->dice_axis(0, pdl(long, [ @ix{@head, $bn1, $bn2} ])),
                   fs => $fs, low => $f_lo, high => $f_hi, order => 4);
my $H  = $X->slice("0:" . ($nh - 1) . ",:");
$H     = $H - $H->average->dummy(0);
my $nt = $H->dim(1);
my $ecg = $X->slice("($nh),:") - $X->slice("(" . ($nh + 1) . "),:");

my ($n0, $n1) = (_s(-300), _s(650));
my ($r) = detect_r($ecg, $fs);
$r = $r->where(($r >= _s(1000 * $skip) - $n0) & ($r + $n1 <= $nt));
my %w  = (base => [-75, -50], P => [-250, -90], QRS => [-45, 55], T => [150, 450], Late => [450, 650]);
my %wi = map { $_ => [ _s($w{$_}[0]) - $n0, _s($w{$_}[1]) - $n0 - 1 ] } keys %w;
my $tms = (sequence($n1 - $n0) + $n0) * (1000 / $fs);

my $A = _mean($H, $r);                                       # (nh, nep)
my (undef, $VL) = _pcs(_win($A, 'Late'));
my $ql = $VL->slice("0:0,:");                                # nh x 1
my $as = ($A x $ql)->flat;
if ($as->index($as->abs->maximum_ind)->sclr < 0) { $ql = -$ql; $as = -$as }
my $A1 = $A - ($A x $ql) x $ql->transpose;                  # slow wave removed
my (undef, $VQ) = _pcs(_win($A1, 'QRS'));
my $Q  = $VQ->slice("0:" . ($K - 1) . ",:")->copy;           # nh x K
my $wt  = $A1 x $Q;                                           # (K, nep): a_k(t)
for my $k (0 .. $K - 1) {                                    # sign: positive at the R peak
    my $i0 = _s(0) - $n0;
    if ($wt->at($k, $i0) < 0) { $Q->slice("$k,:") *= -1; $wt->slice("$k,:") *= -1 }
}
my $fit   = $wt x $Q->transpose;                              # (nh, nep)
my $pw    = ($A1**2)->sumover;                               # (nep)
my $share = ($wt**2)->sumover / $pw->clip(1e-30, undef);

# random K-dimensional subspaces of the data's row space (orthogonal to the
# constant vector and to the slow-wave map): what K maps take by chance
seed_rng(1);
my $B = _rowspace_basis($ql);                                # nh x (nh-2)
my @rs;
for (1 .. 200) {
    my $G = $B x grandom($K, $B->dim(0));                    # nh x K
    my ($v, $U) = eigh_sym($G->transpose x $G);
    my $Qr = $G x ($U / sqrt($v)->dummy(1)) x $U->transpose;  # orthonormal columns
    push @rs, (($A1 x $Qr)**2)->sumover / $pw->clip(1e-30, undef);
}
my $RS  = cat(@rs);                                          # (nep, 200)
my $r05 = $RS->mv(1, 0)->pctover(0.05);
my $r95 = $RS->mv(1, 0)->pctover(0.95);

printf "%s: %d beats, %d head channels (%s left out), K = %d QRS maps\n",
    $file, $r->nelem, $nh, ($exclude || 'none'), $K;
for my $wn (qw(P QRS T)) {
    my ($i0, $i1) = @{ $wi{$wn} };
    printf "  %-4s power share of the %d maps: %.3f (random %d-dim subspaces: median %.3f)\n", $wn, $K,
        (($wt->slice(":,$i0:$i1")**2)->sum / $pw->slice("$i0:$i1")->sum)->sclr, $K,
        (($RS->slice("$i0:$i1,:") * $pw->slice("$i0:$i1")->dummy(1))->sumover
            / $pw->slice("$i0:$i1")->sum)->median;
}

# ---- figures ---------------------------------------------------------------
eval { require PDL::Graphics::Cairo; require PDL::EEG::MAP2D; 1 }
    or do { print "PDL::Graphics::Cairo not available: no figures\n"; exit 0 };
PDL::Graphics::Cairo->import(qw(figure subplots));

# lags for the measured / fitted maps: R, the P and T maxima of |a|, 300 ms
my $an = sqrt(($wt**2)->sumover);                             # (nep)
my %lag = (R => _s(0) - $n0);
for my $wn (qw(P T)) {
    my ($i0, $i1) = @{ $wi{$wn} };
    $lag{$wn} = $i0 + $an->slice("$i0:$i1")->maximum_ind->sclr;
}
$lag{'300'} = _s(300) - $n0;
my @lags = sort { $lag{$a} <=> $lag{$b} } keys %lag;

{
    my $nc  = _max($K + 1, scalar @lags);
    my $fig = figure(width => 230 * $nc + 40, height => 3 * 250 + 40);
    my $gs  = $fig->add_gridspec(3, $nc, left => 0.02, right => 0.98, top => 0.93, bottom => 0.02,
                                 wspace => 0.10, hspace => 0.25);
    for my $k (0 .. $K) {
        my ($v, $amp, $nm) = $k < $K
            ? ($Q->slice("($k),:"), $wt->slice("($k),:")->abs->max->sclr, sprintf('QRS map %d', $k + 1))
            : ($ql->flat, $as->abs->max->sclr, 'slow wave');
        my $vals = $v * $amp;
        _topo($fig, $gs->at(0, $k), $vals, $vals->abs->max->sclr, sprintf('%s (x %.1f)', $nm, $amp));
    }
    for my $c (0 .. $#lags) {
        my $i  = $lag{ $lags[$c] };
        my $vm = $A1->slice(":,($i)");
        my $cl = $vm->abs->max->sclr;
        _topo($fig, $gs->at(1, $c), $vm, $cl,
              sprintf('measured %d ms, +-%.1f', $tms->at($i), $cl));
        _topo($fig, $gs->at(2, $c), $fit->slice(":,($i)"), $cl,
              sprintf('%d maps, share %.2f', $K, $share->at($i)));
    }
    $fig->suptitle(sprintf('%s: cardiac field as %d fixed maps (slow wave removed; red positive, uV)',
                           $file =~ s{.*/}{}r, $K), fontsize => 10);
    $fig->save("${out}_maps.png");
}
{
    my ($fig, $c, $p, $q, $g) = subplots(4, 1, width => 900, height => 1250);
    my @cc = ('#1565C0', '#C62828', '#2E7D32', '#6A1B9A');
    my $ecgA = _mean($ecg->dummy(0), $r)->flat;              # same beats and baseline
    $c->line($tms, $ecgA, color => '#000000', lw => 1.2, label => "$bn1-$bn2");
    $c->axhline(0, color => '#999999', lw => 0.5);
    $c->set_ylabel('uV'); $c->legend;
    $c->set_title("R-locked average of $bn1-$bn2 (ECG)");
    $p->line($tms, $wt->slice("($_),:"), color => $cc[$_ % 4], lw => 1.3, label => 'a_' . ($_ + 1))
        for 0 .. $K - 1;
    $p->line($tms, $as, color => '#777777', lw => 1.0, label => 'slow wave');
    $p->axhline(0, color => '#999999', lw => 0.5);
    $p->set_ylabel('weight [uV]'); $p->legend;
    $p->set_title('weights of the fixed maps in the R-locked average');
    $q->line($tms, $share, color => '#1565C0', lw => 1.3, label => "$K QRS maps");
    $q->line($tms, $r05, color => '#999999', lw => 0.8, label => "random $K-dim, 5th / 95th pct");
    $q->line($tms, $r95, color => '#999999', lw => 0.8);
    $q->set_ylim(0, 1.02); $q->set_ylabel('share of map power');
    $q->set_xlabel('time from R peak [ms]'); $q->legend;
    for my $ax ($c, $p, $q) {
        $ax->axvline(0, color => '#999999', lw => 0.5);
        $ax->set_xlim($tms->at(0), $tms->at(-1));
    }
    if ($K >= 2) {
        for my $wn (qw(P QRS T)) {
            my ($i0, $i1) = @{ $wi{$wn} };
            $g->line($wt->slice("(0),$i0:$i1"), $wt->slice("(1),$i0:$i1"),
                     color => { P => '#2E7D32', QRS => '#1565C0', T => '#C62828' }->{$wn}, lw => 1.2,
                     label => "$wn $w{$wn}[0]..$w{$wn}[1] ms");
        }
        $g->axhline(0, color => '#999999', lw => 0.5); $g->axvline(0, color => '#999999', lw => 0.5);
        $g->set_xlabel('a_1 [uV]'); $g->set_ylabel('a_2 [uV]'); $g->legend;
        $g->set_title('loop of the weights (a_1, a_2)');
    }
    $fig->save("${out}_time.png");
}
print "wrote ${out}_maps.png, ${out}_time.png\n";

# ---------------------------------------------------------------------------
sub _topo {
    my ($fig, $cell, $vals, $cl, $title) = @_;
    $cl = 1e-6 if $cl <= 0;
    my $ax = $fig->add_subplot($cell);
    PDL::EEG::MAP2D::plot_topomap(values => $vals, labels => \@mlab, montage => $montage,
        clim => $cl, contours => 0, names => 1, name_size => 6, title_size => 9,
        title => $title, _ax => $ax, _fig => $fig);
}
sub _s    { my $ms = shift; my $v = $ms * 1e-3 * $fs; $v >= 0 ? int($v + 0.5) : -int(-$v + 0.5) }
sub _max  { $_[0] > $_[1] ? $_[0] : $_[1] }
sub _win  { my ($Z, $n) = @_; $Z->slice(":,$wi{$n}[0]:$wi{$n}[1]") }
sub _mean {
    my ($Y, $rp) = @_;
    my $E = beat_epochs($Y, $rp, $n0, $n1)->double;
    $E -= $E->slice(":,$wi{base}[0]:$wi{base}[1]")->mv(1, 0)->average->dummy(1);
    return $E->mv(2, 0)->average;
}
sub _pcs {
    my $Z = shift;
    my ($v, $V) = eigh_sym($Z->transpose x $Z);
    my $ord = sequence(long, $v->nelem)->slice("-1:0");
    return ($v->index($ord), $V->dice_axis(0, $ord));
}
# orthonormal basis (columns) of the vectors orthogonal to the constant vector and to $q
sub _rowspace_basis {
    my $q = shift;
    my $C = identity($nh) - ones($nh, $nh) / $nh;            # removes the constant vector
    $C    = $C - $q x $q->transpose;
    my ($v, $V) = eigh_sym($C);
    my $keep = which($v > 0.5);
    return $V->dice_axis(0, $keep);
}
