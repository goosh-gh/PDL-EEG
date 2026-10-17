package PDL::EEG::Heartbeat::BeatFeatures;

# Per-beat shape features of one ECG lead (size, shift and stretch of QRS and
# T against the average beat, R size relative to its slow mean, RR intervals),
# as series for the respiratory phase (PDL::EEG::Heartbeat::Respiration).

use strict;
use warnings;
use Carp qw(croak);
use PDL;
use PDL::EEG::Heartbeat              qw(beat_epochs);
use PDL::EEG::Heartbeat::Respiration qw(r_amplitude);
use Exporter 'import';

our @EXPORT_OK = qw(beat_features);

my ($E0_MS, $E1_MS) = (-300, 500);                 # epoch around R [ms]
my @WINDOWS = (                                     # [from, to, centre of the stretch] ms, feature prefix
    [-60, 70,  0,   'Q'],
    [150, 480, 310, 'T'],
);

# mean of the beats within +-$half s of each beat; $E (ne, nb); running sums, no nb x nb matrix
sub _runmean_beats {
    my ($tb, $E, $half) = @_;
    my @t  = $tb->list;
    my $nb = @t;
    my (@lo, @hi);
    my ($l, $h) = (0, 0);
    for my $k (0 .. $nb - 1) {
        $l++ while $t[$l] < $t[$k] - $half;
        $h = $k if $h < $k;
        $h++ while $h + 1 < $nb && $t[$h + 1] <= $t[$k] + $half;
        push @lo, $l;
        push @hi, $h;
    }
    my $C = $E->double->xchg(0, 1)->cumusumover->xchg(0, 1);
    my $Z = zeroes($E->dim(0), 1)->glue(1, $C);     # (ne, nb+1): column k = sum of beats 0 .. k-1
    my ($lo, $hi) = (pdl(long, \@lo), pdl(long, \@hi));
    my $S = $Z->dice_axis(1, $hi + 1) - $Z->dice_axis(1, $lo);
    return $S / ($hi - $lo + 1)->dummy(0);
}

# ---------------------------------------------------------------------------
# $B = beat_features($x, $fs, r => $r, slow => 10)
#
#   $x    : (nt) one ECG lead, band-passed, R peak positive
#   $fs   : sampling rate [Hz]
#   r     : (nbeat) R-peak sample indices (detect_r); beats whose epoch
#           (-300 .. +500 ms) does not fit in $x are dropped
#   slow  : half width [s] of the running mean of the beats (default 10)
#
# Returns a hash reference:
#   tb      (nb) beat times [s]        r (nb) R samples kept
#   use     indices of the beats with all features (not the first / last)
#   F       { name => (nb) }, names [R Qsz Qwid Qlat Tsz Tlat Twid RRp RRn]
#   A       (nb) R amplitude (r_amplitude), units of $x
#   template  average beat;  lag_ms  its time axis [ms]
# ---------------------------------------------------------------------------
sub beat_features {
    my ($x, $fs, %o) = @_;
    croak "beat_features: x must be 1-D" unless $x->ndims == 1;
    croak "beat_features: fs required" unless $fs;
    croak "beat_features: r => R samples required" unless defined $o{r};
    my $ms  = 1000 / $fs;
    my $smp = sub { my $v = $_[0] / $ms; $v >= 0 ? int($v + 0.5) : -int(-$v + 0.5) };
    my $half = $o{slow} // 10;
    my ($l0, $l1) = ($smp->($E0_MS), $smp->($E1_MS));
    my $nt = $x->nelem;
    my $r  = $o{r}->where(($o{r} + $l0 >= 0) & ($o{r} + $l1 < $nt));
    croak "beat_features: fewer than 3 beats inside the signal" if $r->nelem < 3;
    my $tb = $r->double / $fs;
    my $nb = $tb->nelem;

    my $A = r_amplitude($x, $r, $fs);
    my $E = beat_epochs($x->dummy(0), $r, $l0, $l1)->slice('(0)')->double;      # (ne, nb)
    my ($b0, $b1) = ($smp->(-100) - $l0, $smp->(-60) - $l0);                     # baseline -100 .. -60 ms
    $E = $E - $E->slice("$b0:$b1")->average->dummy(0);
    my $T   = $E->mv(1, 0)->average;                                              # average beat
    my $D   = $E - _runmean_beats($tb, $E, $half);                                # beat minus its slow mean
    my $lag = sequence($l1 - $l0) + $l0;                                          # samples
    my $Td  = ($T->slice('2:-1') - $T->slice('0:-3')) / 2;
    $Td = append(append($Td->slice('0:0'), $Td), $Td->slice('-1:-1'));            # dT/dsample

    my %F;
    for my $wn (@WINDOWS) {
        my ($w0, $w1, $tc, $p) = @$wn;
        my $sl = ($smp->($w0) - $l0) . ':' . ($smp->($w1) - $l0);
        my $G  = cat($T->slice($sl), -$Td->slice($sl), -($lag->slice($sl) - $smp->($tc)) * $Td->slice($sl))->transpose;
        my $C  = $D->slice($sl . ',:') x $G x inv($G->transpose x $G)->transpose;   # (3, nb)
        $F{"${p}sz"}  = 100 * $C->slice('(0),:')->flat->copy;
        $F{"${p}lat"} = $ms * $C->slice('(1),:')->flat->copy;
        $F{"${p}wid"} = 100 * $C->slice('(2),:')->flat->copy;
    }
    my $rm = _runmean_beats($tb, $A->dummy(0), $half)->flat;
    $F{R} = ($A - $rm) / $A->avg * 100;
    my $rr = ($r->slice('1:-1') - $r->slice('0:-2'))->double * $ms;
    $F{RRp} = append(pdl($rr->at(0)), $rr);
    $F{RRn} = append($rr, pdl($rr->at(-1)));
    my $use = which((sequence($nb) > 0) & (sequence($nb) < $nb - 1));
    return {
        tb => $tb, r => $r, use => $use, F => \%F, names => [qw(R Qsz Qwid Qlat Tsz Tlat Twid RRp RRn)],
        A => $A, template => $T, lag_ms => $lag * $ms,
    };
}

1;

__END__

=head1 NAME

PDL::EEG::Heartbeat::BeatFeatures - per-beat shape features of one ECG lead

=head1 SYNOPSIS

    use PDL::EEG::Filter                 qw(zp_butter);
    use PDL::EEG::Heartbeat              qw(detect_r);
    use PDL::EEG::Heartbeat::BeatFeatures qw(beat_features);

    my $x   = zp_butter($bn1 - $bn2, fs => $fs, low => 0.5, high => 45, order => 4);
    my ($r) = detect_r($x, $fs);
    my $B   = beat_features($x, $fs, r => $r);
    my $u   = $B->{use};
    my ($tb, $rsize) = ($B->{tb}->index($u), $B->{F}{R}->index($u));

=head1 DESCRIPTION

Each beat (-300 .. +500 ms around R, baseline -100 .. -60 ms) minus the mean of
the beats within C<+-slow> s is projected, in the QRS window (-60 .. +70 ms)
and in the T window (+150 .. +480 ms), on the average beat C<T(t)>, its shift
C<-T'(t)> and its stretch C<-(t - tc) T'(t)> (C<tc> 0 and 310 ms). Features:

    R          R amplitude minus its running mean, % of the mean R amplitude
    Qsz, Tsz   size [% of the average beat]
    Qlat, Tlat shift [ms]
    Qwid, Twid stretch [%]
    RRp, RRn   RR interval ending / starting at the beat [ms]

The lead should have the R peak positive (multiply by the polarity that
C<detect_r> returns if needed). The running means use cumulative sums, so long
records (tens of thousands of beats) need no C<nb x nb> matrix.

=cut
