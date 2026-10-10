package PDL::EEG::Heartbeat;

use strict;
use warnings;
use Carp qw(croak);
use PDL;
use PDL::FFT;
use Exporter 'import';

our @EXPORT_OK = qw(
    zp_butter detect_r beat_epochs beat_average surrogate_r bn_balance
);
our %EXPORT_TAGS = (all => \@EXPORT_OK);

# ---------------------------------------------------------------------------
# zp_butter($x, fs => $hz, low => $f1, high => $f2, order => $n [, pad => $sec])
#
# Zero-phase Butterworth filter, applied in the frequency domain.
#
#   $x : (nt) or (nch, nt)
#   low only          -> high-pass at `low`
#   high only         -> low-pass  at `high`
#   low and high      -> band-pass  low..high
#
# The gain applied to each FFT bin is the squared magnitude of the digital
# (bilinear-transform) Butterworth filter of the given order,
#
#     G(f) = 1 / (1 + q(f)^(2*order)),   W = tan(pi f / fs)
#     q = Wc/W (high-pass),  W/Wc (low-pass),  (W^2 - W1 W2)/((W2 - W1) W) (band-pass)
#
# which is what a forward-backward pass (scipy.signal.filtfilt of
# scipy.signal.butter(order, ...)) has in steady state. The record is extended
# at both ends by odd reflection (`pad` seconds, default 3/low or 3/high,
# at most nt-1 samples) and the two extension ends are joined by a straight line
# up to an FFT-friendly length, so the circular FFT sees no step.
# ---------------------------------------------------------------------------
sub zp_butter {
    my ($x, %o) = @_;
    my $fs = $o{fs} or croak "zp_butter: fs required";
    my ($lo, $hi) = @o{qw(low high)};
    croak "zp_butter: give low and/or high" unless defined $lo || defined $hi;
    my $N  = $o{order} // 2;
    my $nd = $x->ndims;
    croak "zp_butter: x must be (nt) or (nch,nt)" unless $nd == 1 || $nd == 2;

    my $X  = ($nd == 1 ? $x->dummy(0) : $x)->mv(1, 0)->double;   # (nt, nch)
    my $nt = $X->dim(0);
    my $fc = defined $lo ? $lo : $hi;
    my $pad = int((defined $o{pad} ? $o{pad} : 3 / $fc) * $fs);
    $pad = $nt - 1 if $pad > $nt - 1;
    my $L  = _fast_len($nt + 2 * $pad + 16);
    my $ng = $L - $nt - 2 * $pad;                                # join length

    my $re = zeroes(double, $L, $X->dim(1));
    $re->slice("$pad:" . ($pad + $nt - 1)) .= $X;
    if ($pad > 0) {
        my $x0 = $X->slice("0:0");
        my $xn = $X->slice("-1:-1");
        $re->slice("0:" . ($pad - 1)) .= 2 * $x0 - $X->slice("$pad:1");
        $re->slice(($pad + $nt) . ":" . (2 * $pad + $nt - 1))
            .= 2 * $xn - $X->slice("-2:" . ($nt - 1 - $pad));
    }
    {   # straight line from the last extended sample back to the first
        my $a = $re->slice("(" . (2 * $pad + $nt - 1) . ")");    # (nch)
        my $b = $re->slice("(0)");
        my $w = (sequence($ng) + 1) / ($ng + 1);                 # (ng)
        $re->slice((2 * $pad + $nt) . ":-1")
            .= $a->dummy(0) + ($b - $a)->dummy(0) * $w->dummy(1);
    }
    my $im = zeroes(double, $re->dims);
    fft($re, $im);

    my $k = sequence($L);
    my $f = $k * ($fs / $L);
    $f = $f->where($k <= $L / 2)->append(($fs - $f)->where($k > $L / 2));
    my $W = tan(3.14159265358979323846 * $f / $fs);
    $W = $W->clip(1e-300, 1e300);
    my $q;
    if (defined $lo && defined $hi) {
        my ($W1, $W2) = map { tan(3.14159265358979323846 * $_ / $fs) } $lo, $hi;
        $q = ($W * $W - $W1 * $W2) / (($W2 - $W1) * $W);
    } elsif (defined $lo) {
        $q = tan(3.14159265358979323846 * $lo / $fs) / $W;
    } else {
        $q = $W / tan(3.14159265358979323846 * $hi / $fs);
    }
    my $g = 1 / (1 + ($q * $q) ** $N);                            # (L)
    $re *= $g; $im *= $g;
    ifft($re, $im);

    my $y = $re->slice("$pad:" . ($pad + $nt - 1))->mv(0, 1)->sever;   # (nch, nt)
    return $nd == 1 ? $y->slice("(0)")->sever : $y;
}

sub _fast_len {                       # smallest 2^a 3^b 5^c >= n
    my $n = shift;
    for (my $m = $n; ; $m++) {
        my $r = $m;
        for my $p (2, 3, 5) { $r /= $p while $r % $p == 0 }
        return $m if $r == 1;
    }
}

# ---------------------------------------------------------------------------
# ($r, $sign) = detect_r($ecg, $fs, %opt)
#
#   $ecg : (nt) one cardiac channel
#   $r   : (nbeat) long, sample index of each R peak
#   $sign: +1 or -1, polarity of the R wave in $ecg
#
# 1. zero-phase band-pass (band => [5,20] Hz, order 3)
# 2. polarity: the larger of the 99.5th and 0.5th percentile magnitudes
# 3. local maxima above height_frac (0.35) x 99th percentile, at least
#    min_rr (0.30 s) apart; when two are closer the higher one is kept
# 4. each peak moved to the maximum of the sign-corrected, 0.5 Hz high-passed
#    ECG within +-refine_ms (40 ms)
#
# sign => +1 | -1 fixes the polarity instead of estimating it.
# ---------------------------------------------------------------------------
sub detect_r {
    my ($ecg, $fs, %o) = @_;
    croak "detect_r: ecg must be 1-D" unless $ecg->ndims == 1;
    my ($b1, $b2) = @{ $o{band} // [5, 20] };
    my $f = zp_butter($ecg, fs => $fs, low => $b1, high => $b2, order => 3);
    my $sign = $o{sign}
        // ($f->pct(0.995)->sclr >= -$f->pct(0.005)->sclr ? 1 : -1);
    $f = $f * $sign;
    my $thr  = ($o{height_frac} // 0.35) * $f->pct(0.99)->sclr;
    my $dist = int(($o{min_rr} // 0.30) * $fs);

    my $nt  = $f->nelem;
    my $mid = $f->slice("1:-2");
    my $pk  = which(($mid > $f->slice("0:-3")) & ($mid >= $f->slice("2:-1"))
                    & ($mid >= $thr)) + 1;
    return (null->long, $sign) if $pk->isempty;

    # keep the higher of any two peaks closer than $dist (highest first)
    my @p = $pk->list;
    my @h = $f->index($pk)->list;
    my @keep = (1) x @p;
    for my $j (sort { $h[$b] <=> $h[$a] } 0 .. $#p) {
        next unless $keep[$j];
        for (my $k = $j - 1; $k >= 0  && $p[$j] - $p[$k] < $dist; $k--) { $keep[$k] = 0 }
        for (my $k = $j + 1; $k <= $#p && $p[$k] - $p[$j] < $dist; $k++) { $keep[$k] = 0 }
    }
    $pk = pdl(long, [ @p[ grep { $keep[$_] } 0 .. $#p ] ]);

    my $x = $sign * zp_butter($ecg, fs => $fs, low => 0.5, order => 2);
    my $w = int(($o{refine_ms} // 40) * 1e-3 * $fs + 0.5);
    my $idx = ($pk->dummy(0) + sequence(long, 2 * $w + 1)->dummy(1) - $w)
                  ->clip(0, $nt - 1);                              # (2w+1, npk)
    my $am  = $x->index($idx)->maximum_ind;                        # (npk)
    my $r   = $idx->index($am);                                    # idx[am[j], j]
    return ($r->uniq->long, $sign);
}

# ---------------------------------------------------------------------------
# $E = beat_epochs($X, $r, $n0, $n1)
#
#   $X  : (nch, nt)
#   $r  : (nbeat) sample indices; every $r+$n0 .. $r+$n1-1 must lie in the record
#   $E  : (nch, nep, nbeat),  nep = $n1 - $n0,  sample k is at lag ($n0 + k)
# ---------------------------------------------------------------------------
sub beat_epochs {
    my ($X, $r, $n0, $n1) = @_;
    croak "beat_epochs: X must be (nch,nt)" unless $X->ndims == 2;
    my $nt = $X->dim(1);
    croak "beat_epochs: an epoch runs off the record"
        if $r->nelem && ($r->min + $n0 < 0 || $r->max + $n1 > $nt);
    my $idx = $r->long->dummy(0) + sequence(long, $n1 - $n0)->dummy(1) + $n0;
    return $X->mv(1, 0)->index($idx->dummy(0))->sever;             # (nch,nep,nbeat)
}

# ---------------------------------------------------------------------------
# $res = beat_average($X, $r, $fs, %opt)
#
#   tmin, tmax => seconds around R (default -0.3, 0.8)
#   baseline   => [t0, t1] seconds: the mean over this window is removed from
#                 every beat and channel before combining (default: none)
#   reject     => uV, a number or an (nch) piddle of per-channel limits: a beat
#                 is left out of a channel's average when that channel's
#                 peak-to-peak over the epoch exceeds the limit (default: none)
#   stat       => 'mean' (default) | 'median'
#
# Beats whose epoch does not fit in the record are dropped.
# Returns { avg (nch,nep), sd (nch,nep), n (nch), t_ms (nep), r (beats used),
#           n0, n1 }.  sd is the across-beat standard deviation (n-1).
# ---------------------------------------------------------------------------
sub beat_average {
    my ($X, $r, $fs, %o) = @_;
    my $n0 = _rnd(($o{tmin} // -0.3) * $fs);
    my $n1 = _rnd(($o{tmax} //  0.8) * $fs);
    my $nt = $X->dim(1);
    $r = $r->where(($r + $n0 >= 0) & ($r + $n1 <= $nt));
    croak "beat_average: no beat fits in the record" if $r->isempty;
    my $E = beat_epochs($X, $r, $n0, $n1)->double;                 # (nch,nep,nb)
    if (my $b = $o{baseline}) {
        my ($i0, $i1) = (_rnd($b->[0] * $fs) - $n0, _rnd($b->[1] * $fs) - $n0 - 1);
        croak "beat_average: baseline window outside the epoch"
            if $i0 < 0 || $i1 >= $n1 - $n0 || $i1 < $i0;
        $E -= $E->slice(":,$i0:$i1")->mv(1, 0)->average->dummy(1);
    }
    my $good = ones(byte, $E->dim(0), $E->dim(2));                 # (nch,nb)
    if (defined $o{reject}) {
        my $Em = $E->mv(1, 0);
        my $thr = ref $o{reject} ? $o{reject}->flat->dummy(1) : $o{reject};
        $good = ($Em->maximum - $Em->minimum) <= $thr;
    }
    my $n = $good->mv(1, 0)->sumover->double;                             # (nch)
    croak "beat_average: a channel has no accepted beat" if $n->min < 1;
    my ($avg, $sd);
    if (($o{stat} // 'mean') eq 'median') {
        $avg = $E->setbadif(!$good->dummy(1))->mv(2, 0)->medover->setbadtoval(0);
    } else {
        $avg = ($E * $good->dummy(1))->mv(2, 0)->sumover / $n->dummy(1);
    }
    my $dev = ($E - $avg->dummy(2)) * $good->dummy(1);
    $sd = sqrt(($dev * $dev)->mv(2, 0)->sumover / ($n->dummy(1) - 1)->clip(1, undef));
    return { avg => $avg, sd => $sd, n => $n,
             t_ms => (sequence($n1 - $n0) + $n0) * (1000 / $fs),
             r => $r, n0 => $n0, n1 => $n1 };
}

sub _rnd { my $v = shift; $v >= 0 ? int($v + 0.5) : -int(-$v + 0.5) }

# ---------------------------------------------------------------------------
# $s = surrogate_r($r, $nt, $fs, %opt)
#
# Non-heartbeat time points: every R index is moved by a random offset whose
# magnitude is uniform in [min, max] seconds (default 0.25 .. 0.60) and whose
# sign is random (min => 0, max => 0.5 gives a uniform +-0.5 s jitter). Points that would leave [lo, hi) samples are reflected to the
# other sign; lo/hi default to 0 and $nt.
#
#   seed => integer : the offsets come from CORE::rand after CORE::srand(seed)
# ---------------------------------------------------------------------------
sub surrogate_r {
    my ($r, $nt, $fs, %o) = @_;
    my ($mn, $mx) = (($o{min} // 0.25) * $fs, ($o{max} // 0.60) * $fs);
    my ($lo, $hi) = ($o{lo} // 0, $o{hi} // $nt);
    CORE::srand($o{seed}) if defined $o{seed};
    my @s;
    for my $p ($r->list) {
        my $d = int($mn + CORE::rand() * ($mx - $mn) + 0.5);
        $d = -$d if CORE::rand() < 0.5;
        $d = -$d if $p + $d < $lo || $p + $d >= $hi;
        push @s, $p + $d if $p + $d >= $lo && $p + $d < $hi;
    }
    return pdl(long, \@s)->qsort;
}

# ---------------------------------------------------------------------------
# ($w, $rms, $rms1, $rms2) = bn_balance($bn1, $bn2)
#
# Weight w that minimises the energy of  w*bn1 + (1-w)*bn2  over the samples
# given (closed form), the RMS of that mix, and the RMS of each input.
# $bn1, $bn2: same shape; all elements are used.
# ---------------------------------------------------------------------------
sub bn_balance {
    my ($a, $b) = map { $_->flat->double } @_;
    croak "bn_balance: shapes differ" unless $a->nelem == $b->nelem;
    my $d   = $a - $b;
    my $den = sum($d * $d);
    croak "bn_balance: the two inputs are identical" if $den <= 0;
    my $w   = -sum($b * $d) / $den;
    my $mix = $w * $a + (1 - $w) * $b;
    my $rms = sub { sqrt(sum($_[0] * $_[0]) / $_[0]->nelem) };
    return ($w, $rms->($mix), $rms->($a), $rms->($b));
}

1;

__END__

=head1 NAME

PDL::EEG::Heartbeat - R-peak detection and R-locked averaging of multichannel EEG

=head1 SYNOPSIS

    use PDL::EEG::IO::EDF   qw(read_edf);
    use PDL::EEG::Heartbeat qw(detect_r beat_average bn_balance);

    my $e  = read_edf('rec.edf');                    # data (nch, nt), uV
    my $X  = $e->{data};
    my $fs = $e->{fs};
    my ($r, $sign) = detect_r($X->slice("($i_ecg)"), $fs);

    my $res = beat_average($X, $r, $fs,
                           tmin => -0.3, tmax => 0.8, baseline => [-0.3, -0.2]);
    # $res->{avg} (nch, nep), $res->{t_ms} (nep)

=head1 DESCRIPTION

Data layout is C<(nch, nt)>, matching C<read_edf>, L<PDL::EEG::GED> and
L<PDL::EEG::Regress>.

=over 4

=item zp_butter($x, fs=>, low=>, high=>, order=>)

Zero-phase Butterworth high-, low- or band-pass, applied in the frequency
domain with the squared-magnitude response of the digital Butterworth filter
(the steady-state response of a forward-backward pass).

=item detect_r($ecg, $fs, %opt)

R-peak sample indices and R polarity from one cardiac channel.

=item beat_epochs($X, $r, $n0, $n1)

C<(nch, nep, nbeat)> block of samples around each R index.

=item beat_average($X, $r, $fs, %opt)

R-locked average and across-beat SD of every channel, with optional per-beat
baseline removal, per-channel peak-to-peak rejection, and median instead of
mean.

=item surrogate_r($r, $nt, $fs, %opt)

Time points displaced from every R peak, for the same average at
non-heartbeat times.

=item bn_balance($bn1, $bn2)

The weight C<w> for which C<w*bn1 + (1-w)*bn2> has least energy over the
samples given.

=back

=cut
