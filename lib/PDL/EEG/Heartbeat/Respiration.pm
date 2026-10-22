package PDL::EEG::Heartbeat::Respiration;

use strict;
use warnings;
use Carp qw(croak);
use PDL;
use PDL::EEG::FFT::Dispatcher qw(fft_inplace);
use PDL::EEG::Filter          qw(zp_butter analytic_signal);
use PDL::EEG::Heartbeat       qw(beat_epochs);
use Exporter 'import';

our @EXPORT_OK = qw(rr_series r_amplitude resp_phase phase_at phase_locking
                    phase_template phase_template_multi phase_agreement
                    fourier_fit fourier_eval wrap_phase unwrap_phase circ_sd);
our %EXPORT_TAGS = (all => \@EXPORT_OK);

my $PI = 3.14159265358979323846;

# ---------------------------------------------------------------------------
# ($t, $rr) = rr_series($r, $fs)
#
#   $r  : (nbeat) R-peak sample indices, ascending
#   $t  : (nbeat-1) time [s] of beats 1 .. nbeat-1
#   $rr : (nbeat-1) interval [s] that ends at each of those beats
# ---------------------------------------------------------------------------
sub rr_series {
    my ($r, $fs) = @_;
    croak "rr_series: fs required" unless $fs;
    croak "rr_series: at least 2 beats needed" if $r->nelem < 2;
    my $t = $r->double / $fs;
    return ($t->slice("1:-1")->sever, ($t->slice("1:-1") - $t->slice("0:-2"))->sever);
}

# ---------------------------------------------------------------------------
# $amp = r_amplitude($X, $r, $fs, %opt)
#
#   $X : (nt) or (nlead, nt)
#   $amp: (nbeat) or (nlead, nbeat): the sample at each R index minus the mean
#        of the window `pre` before it
#
#   pre => [t0, t1]   window in seconds relative to R, default [-0.060, -0.041]
# ---------------------------------------------------------------------------
sub r_amplitude {
    my ($X, $r, $fs, %o) = @_;
    croak "r_amplitude: fs required" unless $fs;
    my $nd = $X->ndims;
    croak "r_amplitude: X must be (nt) or (nlead,nt)" unless $nd == 1 || $nd == 2;
    my ($p0, $p1) = map { _rnd($_ * $fs) } @{ $o{pre} // [-0.060, -0.041] };
    croak "r_amplitude: pre must lie before R, t0 <= t1 < 0" unless $p0 <= $p1 && $p1 < 0;
    my $E   = beat_epochs($nd == 1 ? $X->dummy(0) : $X, $r, $p0, 1)->double;   # (nlead, -p0+1, nbeat)
    my $pre = $E->slice(":,0:" . ($p1 - $p0))->mv(1, 0)->average;             # (nlead, nbeat)
    my $amp = ($E->slice(":,(-1)") - $pre)->sever;
    return $nd == 1 ? $amp->slice("(0)")->sever : $amp;
}

# ---------------------------------------------------------------------------
# $ph = resp_phase($t, $y, %opt)
#
#   $t : (n) times [s] of a beat series, ascending
#   $y : (n) the series (RR interval, R amplitude, ...)
#
#   fs    => 4             rate [Hz] the series is resampled to (linear interpolation)
#   band  => [0.20, 0.45]  pass band [Hz]
#   order => 3             Butterworth order
#
# The series is put on a regular grid from $t(0) to $t(-1), its mean is
# removed, it is band-passed with zp_butter and its analytic signal is taken.
#
# Returns a hash reference:
#   t         (m)  grid times [s]
#   fs, band
#   series    (m)  band-passed series
#   re, im    (m)  analytic signal of `series`
#   peak_hz        frequency with the largest power of the resampled series
#                  inside the band (Hann window)
#   band_frac      share of the resampled series' power that lies in the band
# ---------------------------------------------------------------------------
sub resp_phase {
    my ($t, $y, %o) = @_;
    croak "resp_phase: t and y must be 1-D of equal length" unless $t->ndims == 1 && $y->ndims == 1 && $t->nelem == $y->nelem;
    my $fs = $o{fs} // 4;
    my ($lo, $hi) = @{ $o{band} // [0.20, 0.45] };
    croak "resp_phase: band must satisfy 0 < low < high < fs/2" unless $lo > 0 && $lo < $hi && $hi < $fs / 2;
    my $td = $t->double;
    croak "resp_phase: t must be ascending" unless all($td->slice("1:-1") > $td->slice("0:-2"));
    my $t0 = $td->at(0);
    my $m  = int(($td->at(-1) - $t0) * $fs) + 1;
    croak "resp_phase: the series spans fewer than 3 cycles of the lower band edge" if ($m - 1) / $fs < 3 / $lo;
    my $tg = sequence(double, $m) / $fs + $t0;
    my ($yg) = interpol($tg, $td, $y->double);
    $yg = $yg - $yg->avg;
    my $bp = zp_butter($yg, fs => $fs, low => $lo, high => $hi, order => $o{order} // 3);
    my ($re, $im) = analytic_signal($bp);

    my $win = 0.5 - 0.5 * cos(2 * $PI * sequence($m) / ($m - 1));
    my ($sr, $si) = ($yg * $win, zeroes(double, $m));
    fft_inplace($sr, $si);
    my $nh  = int($m / 2);
    my $pw  = ($sr * $sr + $si * $si)->slice("1:$nh");
    my $f   = (sequence($nh) + 1) * $fs / $m;
    my $in  = which(($f >= $lo) & ($f <= $hi));
    croak "resp_phase: no spectral line inside the band" if $in->isempty;
    my $pin = $pw->index($in);
    return {
        t => $tg, fs => $fs, band => [$lo, $hi], series => $bp, re => $re, im => $im,
        peak_hz   => $f->index($in)->at($pin->maximum_ind->sclr),
        band_frac => ($pin->sum / $pw->sum)->sclr,
    };
}

# ---------------------------------------------------------------------------
# $phi = phase_at($ph, $times)
#
#   $ph    : result of resp_phase
#   $times : times [s], any shape, within the grid of $ph (up to one grid step
#            beyond either end is taken as the end value)
#   $phi   : radians, -pi..pi; 0 where the band-passed series has a maximum,
#            +-pi where it has a minimum
#
# The real and the imaginary part of the analytic signal are interpolated
# linearly, then the angle is taken.
# ---------------------------------------------------------------------------
sub phase_at {
    my ($ph, $times) = @_;
    my $tg = $ph->{t};
    my $tt = pdl($times)->double;
    my ($g0, $g1, $step) = ($tg->at(0), $tg->at(-1), 1 / $ph->{fs});
    croak sprintf("phase_at: times must lie within %.3f .. %.3f s", $g0, $g1)
        if $tt->min < $g0 - $step || $tt->max > $g1 + $step;
    my $tc = $tt->clip($g0, $g1);
    my ($c) = interpol($tc, $tg, $ph->{re});
    my ($s) = interpol($tc, $tg, $ph->{im});
    return atan2($s, $c);
}

# ---------------------------------------------------------------------------
# ($plv, $diff) = phase_locking($p, $q)
#
#   $p, $q : phases [rad], equal length
#   $plv   : length of the mean of exp(i (p - q)), 0 .. 1
#   $diff  : angle of that mean [rad], -pi..pi: the mean of p - q
# ---------------------------------------------------------------------------
sub phase_locking {
    my ($p, $q) = @_;
    croak "phase_locking: p and q must have the same number of elements" unless $p->nelem == $q->nelem;
    my $d = $p->flat->double - $q->flat->double;
    my ($c, $s) = (cos($d)->avg, sin($d)->avg);
    return (sqrt($c * $c + $s * $s)->sclr, atan2($s, $c)->sclr);
}

# ===========================================================================
# Respiratory phase re-estimated against the shape of one breath
# ===========================================================================

# ---------------------------------------------------------------------------
# $w = wrap_phase($x)      angle(s) to -pi .. pi (pi itself goes to -pi)
# $u = unwrap_phase($p)    phases of successive beats made continuous
# $s = circ_sd($d)         circular SD [rad] of the angles $d:
#                          sqrt(-2 ln R), R = length of the mean of exp(i d)
# ---------------------------------------------------------------------------
sub wrap_phase {
    my $x = shift;
    return $x - 2 * $PI * (($x + $PI) / (2 * $PI))->floor;
}

sub unwrap_phase {
    my $p = shift;
    return $p->copy if $p->nelem < 2;
    my $d = wrap_phase($p->slice('1:-1') - $p->slice('0:-2'));
    return append($p->slice('0:0'), $p->at(0) + $d->cumusumover);
}

sub circ_sd {
    my $d = shift;
    my $R = sqrt(cos($d)->avg ** 2 + sin($d)->avg ** 2);
    $R = 1 if $R > 1;                                   # rounding may put R just above 1
    return sqrt(-2 * log($R));
}

# ---------------------------------------------------------------------------
# $c = fourier_fit($phi, $y, $K)
#
#   $phi : (n) phase of each beat [rad]
#   $y   : (n) or (n, nf) one or several per-beat series
#   $c   : (2K+1, nf) least-squares coefficients of
#          y = c0 + sum_k (c_{2k-1} cos k phi + c_{2k} sin k phi), k = 1 .. K
#
# $f = fourier_eval($c, $phi)   (n, nf) the series of $c at the phases $phi
# ---------------------------------------------------------------------------
sub _fourier_basis {
    my ($phi, $K) = @_;
    my @col = (ones($phi->nelem));
    push @col, cos($_ * $phi), sin($_ * $phi) for 1 .. $K;
    return cat(@col)->transpose;                        # (2K+1, n) as a matrix: n rows
}

sub fourier_fit {
    my ($phi, $y, $K) = @_;
    croak "fourier_fit: K must be >= 1" unless $K && $K >= 1;
    croak "fourier_fit: phi and y must have the same number of beats" unless $phi->nelem == $y->dim(0);
    my $X = _fourier_basis($phi, $K);
    my $Y = $y->ndims == 1 ? $y->dummy(1) : $y;
    my $c = inv($X->transpose x $X) x ($X->transpose x $Y->transpose);
    return $c->transpose;
}

sub fourier_eval {
    my ($c, $phi) = @_;
    my $K = ($c->dim(0) - 1) / 2;
    return (_fourier_basis($phi, $K) x $c->transpose)->transpose;
}

# For each beat, the index range of the beats within +-$half s (times ascending), found by two pointers.
# The same sets as which(abs($tb - $tb->at($k)) <= $half), without a whole-length temporary per beat: with
# PDL 2.085 such temporaries made a loop over 30000 beats grow to several GB.
sub _bounds {
    my ($tb, $half) = @_;
    my @t = $tb->list;
    my (@lo, @hi);
    my ($l, $h) = (0, 0);
    for my $k (0 .. $#t) {
        $l++ while $t[$l] < $t[$k] - $half;
        $h = $k if $h < $k;
        $h++ while $h + 1 <= $#t && $t[$h + 1] <= $t[$k] + $half;
        push @lo, $l;
        push @hi, $h;
    }
    return (\@lo, \@hi);
}

sub _win { my ($lo, $hi, $k) = @_; return sequence(long, $hi->[$k] - $lo->[$k] + 1) + $lo->[$k] }

# local rate [rad/s]: slope of the unwrapped starting phase over the window
sub _local_rate {
    my ($tt, $pp) = @_;
    return 0 if $tt->nelem < 2;
    my $tc = $tt - $tt->avg;
    return sum($tc * ($pp - $pp->avg)) / sum($tc ** 2);
}

# grid minimum refined by a parabola through it and its two neighbours (the grid is circular)
sub _grid_min {
    my ($ss, $grid) = @_;
    my $ng = $ss->nelem;
    my $i  = $ss->minimum_ind->sclr;
    my ($s0, $s1, $s2) = map { $ss->at(($i + $_) % $ng) } -1, 0, 1;
    my $den = $s0 - 2 * $s1 + $s2;
    my $dx  = $den > 0 ? 0.5 * ($s0 - $s2) / $den : 0;
    return $grid->at($i) + $dx * 2 * $PI / $ng;
}

# One pass of the local fit, one or several features.  $Y (n, nf), templates $C (2K+1, nf), weights $wt (nf).
sub _phase_local {
    my ($tb, $Y, $phi_init, $C, $wt, $win, $ng) = @_;
    my ($n, $nf) = $Y->dims;
    my $K  = ($C->dim(0) - 1) / 2;
    my $Cc = $C->copy;
    $Cc->slice('(0),:') .= 0;                           # the template's mean is not part of the shape
    my $pu   = unwrap_phase($phi_init);
    my $grid = sequence($ng) / $ng * 2 * $PI - $PI;
    my $out  = zeroes($n);
    my ($blo, $bhi) = _bounds($tb, $win / 2);
    for my $k (0 .. $n - 1) {
        my $wi = _win($blo, $bhi, $k);
        my $tt = $tb->index($wi) - $tb->at($k);
        my $m  = $wi->nelem;
        my $om = _local_rate($tt, $pu->index($wi));
        my $P  = ($grid->dummy(0) + $om * $tt->dummy(1))->flat;      # m*ng trial phases, beats fastest
        my $S  = (_fourier_basis($P, $K) x $Cc->transpose)->transpose; # (m*ng, nf) template values
        my $tot = zeroes($ng);
        for my $f (0 .. $nf - 1) {
            my $Sf = $S->slice(":,($f)")->copy->reshape($m, $ng);
            my $Sc = $Sf - $Sf->average->dummy(0);
            my $yc = $Y->slice(":,($f)")->index($wi);
            $yc = $yc - $yc->avg;
            my $gain = (($Sc * $yc->dummy(1))->sumover / (($Sc ** 2)->sumover + 1e-12))->clip(0, undef);
            $tot += $wt->at($f) * (($yc->dummy(1) - $gain->dummy(0) * $Sc) ** 2)->sumover;
        }
        $out->set($k, _grid_min($tot, $grid));
    }
    return wrap_phase($out);
}

# ---------------------------------------------------------------------------
# $res = phase_template($tb, $y, $phi0, %opt)
#
#   $tb   : (n) beat times [s], ascending
#   $y    : (n) per-beat series (e.g. R amplitude minus its slow mean)
#   $phi0 : (n) starting phase at the beats [rad] (e.g. phase_at of resp_phase)
#
#   win   => 5     window [s] around each beat
#   K     => 2     harmonics of the template
#   iter  => 3     rounds
#   ngrid => 120   phase grid of the search
#
# Each round: the template (the shape of one breath) is the Fourier fit of $y on
# the current phase; then, for every beat, the beats within +-win/2 s are fitted
# by y_j = a + g s(phi_c + w (t_j - t_k)) with the template s, the local rate w
# from the current phase, a free offset a and gain g >= 0; phi_c is the new phase
# of the beat.  Returns { phi => (n) phase [rad], c => (2K+1, 1) template }.
# ---------------------------------------------------------------------------
sub phase_template {
    my ($tb, $y, $phi0, %o) = @_;
    croak "phase_template: tb, y and phi0 must be 1-D of equal length"
        unless $tb->ndims == 1 && $y->ndims == 1 && $tb->nelem == $y->nelem && $phi0->nelem == $y->nelem;
    my $K = $o{K} // 2;
    my $phi = $phi0;
    my $c;
    for (1 .. ($o{iter} // 3)) {
        $c   = fourier_fit($phi, $y, $K);
        $phi = _phase_local($tb, $y->dummy(1), $phi, $c, ones(1), $o{win} // 5, $o{ngrid} // 120);
    }
    return { phi => $phi, c => fourier_fit($phi, $y, $K) };
}

# ---------------------------------------------------------------------------
# $res = phase_template_multi($tb, $Y, $phi0, %opt)
#
#   $Y : (n, nf) several per-beat series that share one phase
#
#   gate => 0.1   a series is used only if the fundamental on $phi0 explains at
#                 least this share of its variance
#   win, K, iter, ngrid as phase_template
#
# As phase_template, with one template per series, a free offset and gain >= 0
# per window and series, and the series weighted by 1 / their residual variance.
# Returns { phi, C (2K+1, nused), used (indices into the columns of $Y),
#           w (nused) weights, r2 (nf) share of variance on $phi0 }.
# ---------------------------------------------------------------------------
sub phase_template_multi {
    my ($tb, $Y, $phi0, %o) = @_;
    croak "phase_template_multi: Y must be (n, nf)" unless $Y->ndims == 2 && $Y->dim(0) == $tb->nelem && $phi0->nelem == $tb->nelem;
    my $K = $o{K} // 2;
    my ($n, $nf) = $Y->dims;
    my $r2 = pdl(map {
        my $y = $Y->slice(":,($_)");
        my $f = fourier_eval(fourier_fit($phi0, $y, 1), $phi0)->flat;
        (sum(($f - $y->avg) ** 2) / sum(($y - $y->avg) ** 2))->sclr
    } 0 .. $nf - 1);
    my $used = which($r2 >= ($o{gate} // 0.1));
    $used = pdl(long, [0]) if $used->isempty;
    my $Yu  = $Y->dice_axis(1, $used);
    my $phi = $phi0;
    my ($C, $wt);
    for (1 .. ($o{iter} // 3)) {
        $C = fourier_fit($phi, $Yu, $K);
        my $res = $Yu - fourier_eval($C, $phi);
        $wt  = 1 / (($res ** 2)->sumover / ($n - 2 * $K - 1));
        $phi = _phase_local($tb, $Yu, $phi, $C, $wt, $o{win} // 5, $o{ngrid} // 120);
    }
    return { phi => $phi, C => fourier_fit($phi, $Yu, $K), used => $used, w => $wt, r2 => $r2 };
}

# ---------------------------------------------------------------------------
# $q = phase_agreement($tb, $p, $q2, %opt)
#
#   $tb     : (n) beat times [s], ascending
#   $p, $q2 : (n) two estimates of the phase at the beats [rad]
#
#   win => 120   window [s] centred on each beat (shorter at the ends)
#   min => 40    fewest beats in a window; below it the value is NaN
#
# Circular SD [rad] of p - q2 over the beats within +-win/2 s of each beat.  A
# constant difference does not count.  With p the phase of resp_phase and q2 that
# of phase_template on the same series, this is the quality index q_agree.
# ---------------------------------------------------------------------------
sub phase_agreement {
    my ($tb, $p, $q2, %o) = @_;
    croak "phase_agreement: tb, p and q2 must have the same number of beats" unless $tb->nelem == $p->nelem && $p->nelem == $q2->nelem;
    my $half = ($o{win} // 120) / 2;
    my $min  = $o{min} // 40;
    my $d = wrap_phase($p->double - $q2->double);
    my @c = cos($d)->list;
    my @s = sin($d)->list;
    my @t = $tb->list;
    my $n = @t;
    my ($lo, $hi, $sc, $ss) = (0, -1, 0, 0);
    my @out;
    for my $k (0 .. $n - 1) {
        while ($hi + 1 < $n && $t[$hi + 1] < $t[$k] + $half) { $hi++; $sc += $c[$hi]; $ss += $s[$hi] }
        while ($t[$lo] < $t[$k] - $half) { $sc -= $c[$lo]; $ss -= $s[$lo]; $lo++ }
        my $m = $hi - $lo + 1;
        if ($m < $min) { push @out, 9**9**9 / 9**9**9; next }      # NaN
        my $R = sqrt($sc ** 2 + $ss ** 2) / $m;
        $R = 1 if $R > 1;
        push @out, sqrt(-2 * log($R));
    }
    return pdl(double, \@out);
}

sub _rnd { my $v = shift; $v >= 0 ? int($v + 0.5) : -int(-$v + 0.5) }

1;

__END__

=head1 NAME

PDL::EEG::Heartbeat::Respiration - respiratory phase from beat-to-beat series of the heartbeat (ECG-derived respiration)

=head1 SYNOPSIS

    use PDL::EEG::Heartbeat qw(detect_r);
    use PDL::EEG::Heartbeat::Respiration qw(:all);

    my ($r) = detect_r($ecg, $fs);

    my ($t, $rr) = rr_series($r, $fs);                 # RR interval series
    my $ph_rr    = resp_phase($t, $rr);

    my $amp = r_amplitude($X2, $r, $fs);               # (2, nbeat): two leads
    my $mag = sqrt(($amp ** 2)->sumover);              # (nbeat): size of the R vector
    my $ph_am = resp_phase($r->double / $fs, $mag);

    my $phi = phase_at($ph_rr, $r->double / $fs);      # phase at each R peak
    my ($plv, $diff) = phase_locking(phase_at($ph_am, $t), phase_at($ph_rr, $t));

=head1 DESCRIPTION

Breathing changes the interval between heartbeats and the size of the R wave
in a given lead. This module takes such a beat-to-beat series, isolates the
band of the breathing rate and returns the phase of that oscillation at any
time. This is known as ECG-derived respiration (EDR).

The series has one value per heartbeat. C<resp_phase> puts it on a regular
grid by linear interpolation; a breathing rate above half the heart rate is
not represented.

The phase is that of the series, not of the airflow: it is 0 where the
band-passed series has a maximum and C<+-pi> where it has a minimum. Which of
the two is inspiration is not determined here.

=over 4

=item rr_series($r, $fs)

C<($t, $rr)>: the time of every beat but the first, and the interval that ends
at it, both in seconds.

=item r_amplitude($X, $r, $fs, pre => [t0, t1])

The sample at each R index minus the mean of the window C<pre> before it
(seconds relative to R, default C<[-0.060, -0.041]>). C<$X> is C<(nt)> or
C<(nlead, nt)>; the result is C<(nbeat)> or C<(nlead, nbeat)>.

=item resp_phase($t, $y, fs => 4, band => [0.20, 0.45], order => 3)

Resamples the series C<$y> at times C<$t> to C<fs> Hz, removes its mean,
band-passes it (C<zp_butter> of L<PDL::EEG::Filter>) and takes the analytic
signal. Returns a hash reference with the grid times C<t>, the band-passed
C<series>, the analytic signal C<re>, C<im>, the frequency C<peak_hz> of the
largest spectral line of the resampled series inside the band, and
C<band_frac>, the share of the resampled series' power inside the band.

The spectrum for C<peak_hz> and C<band_frac> is taken at the length of the
grid, with the library L<PDL::EEG::FFT::Dispatcher> chooses: PDL::PocketFFT
(through PDL::FFTDispatcher) or PDL::FFTW3 when installed, otherwise PDL::FFT,
which croaks for a length it cannot transform.

=item phase_at($ph, $times)

Phase in radians (C<-pi..pi>) of the result C<$ph> of C<resp_phase> at the
given times.

=item phase_locking($p, $q)

C<($plv, $diff)> of two phase series: the length (0..1) and the angle of the
mean of C<exp(i(p - q))>. Two series derived from different quantities of the
same recording (RR interval, R amplitude) that follow the same breathing give
a C<$plv> near 1.

=item phase_template($tb, $y, $phi0, win => 5, K => 2, iter => 3, ngrid => 120)

Re-estimates the phase at the beats against the shape of one breath. Each
round, the template is the Fourier series (C<K> harmonics) of C<$y> fitted on
the current phase; then for every beat the beats within C<+-win/2> s are fitted
by the template at phases C<phi_c + w (t_j - t_k)>, with the local rate C<w>
from the current phase and a free offset and gain (gain E<gt>= 0) per window;
C<phi_c> (grid search, refined by a parabola) is the new phase of the beat.
The beat samples are used as they are, so a harmonic of the breath above half
the beat rate does not fold into the phase, as it does in C<resp_phase>.
C<$phi0> is the starting phase, usually C<phase_at> of C<resp_phase> on the
same series. Returns C<{ phi =E<gt> (n) phase [rad], c =E<gt> (2K+1, 1) template }>.
The phase origin is that of C<$phi0>.

=item phase_template_multi($tb, $Y, $phi0, gate => 0.1, win, K, iter, ngrid)

As C<phase_template> for several per-beat series C<$Y> C<(n, nf)> that share
one phase: one template, offset and gain per series, the series weighted by
1 / their residual variance. A series is used only if its fundamental on
C<$phi0> explains at least C<gate> of its variance. Returns
C<{ phi, C, used, w, r2 }>.

=item phase_agreement($tb, $p, $q, win => 120, min => 40)

At each beat, the circular SD [rad] of C<$p - $q> over the beats within
C<+-win/2> s (fewer at the ends of the series; NaN where fewer than C<min>).
A constant difference does not count. With C<$p> from C<resp_phase> and C<$q>
from C<phase_template> on the same series, this is an index of how reliable
the phase is in that stretch: the two estimates rest on different
assumptions (a sinusoid in the band; a repeating shape) and differ where the
series is not dominated by regular breathing.

=item fourier_fit($phi, $y, $K), fourier_eval($c, $phi)

Least-squares coefficients C<(2K+1, nf)> of
C<y = c0 + sum_k (c_(2k-1) cos k phi + c_(2k) sin k phi)> for each column of
C<$y> C<(n)> or C<(n, nf)>, and the series of given coefficients at C<$phi>.

=item wrap_phase($x), unwrap_phase($p), circ_sd($d)

Angle to C<-pi .. pi>; phases of successive beats made continuous; circular
SD C<sqrt(-2 ln R)> of angles, C<R> the length of the mean of C<exp(i d)>.

=back

=cut
