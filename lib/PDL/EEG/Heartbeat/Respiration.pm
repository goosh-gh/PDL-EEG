package PDL::EEG::Heartbeat::Respiration;

use strict;
use warnings;
use Carp qw(croak);
use PDL;
use PDL::EEG::FFT::Dispatcher qw(fft_inplace);
use PDL::EEG::Filter          qw(zp_butter analytic_signal);
use PDL::EEG::Heartbeat       qw(beat_epochs);
use Exporter 'import';

our @EXPORT_OK = qw(rr_series r_amplitude resp_phase phase_at phase_locking);
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
grid. When PDL::FFT cannot transform that length, PDL::FFTW3 is used; without
it the call croaks (L<PDL::EEG::FFT::Dispatcher>).

=item phase_at($ph, $times)

Phase in radians (C<-pi..pi>) of the result C<$ph> of C<resp_phase> at the
given times.

=item phase_locking($p, $q)

C<($plv, $diff)> of two phase series: the length (0..1) and the angle of the
mean of C<exp(i(p - q))>. Two series derived from different quantities of the
same recording (RR interval, R amplitude) that follow the same breathing give
a C<$plv> near 1.

=back

=cut
