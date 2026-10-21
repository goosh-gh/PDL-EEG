use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::FFT::Dispatcher qw(fft_len fft_ok fft_backend fft_inplace ifft_inplace);
use PDL::EEG::Filter qw(zp_butter analytic_signal);
use PDL::EEG::TFA qw(cwt_morlet);
use PDL::EEG::Heartbeat qw(detect_r);
use PDL::EEG::Heartbeat::Respiration qw(resp_phase);

my $PI = 3.14159265358979323846;
my $HAVE_FFTW3 = eval { require PDL::FFTW3; 1 } ? 1 : 0;
my $BAD = 157464;                 # 2^3 3^9: 12 factors in PDL::FFT

sub is235 { my $r = shift; for my $p (2, 3, 5) { $r /= $p while $r % $p == 0 } $r == 1 }

sub dft {                         # direct DFT along dim 0: ($re, $im) of real (n, ...) input
    my $x = shift;
    my $n = $x->dim(0);
    my $k = sequence(double, $n);
    my $ang = -2 * $PI * $k->dummy(0) * $k->dummy(1) / $n;        # (j, k)
    my $re = (cos($ang)->dummy(2) * $x->dummy(1))->sumover;     # (k, ...)
    my $im = (sin($ang)->dummy(2) * $x->dummy(1))->sumover;
    return ($re, $im);
}

subtest 'factor count and lengths' => sub {
    ok(!fft_ok($_), "fft_ok($_) false")
        for 157464, 262440, 279936, 354294, 393660, 437400, 466560, 472392;
    ok(fft_ok($_), "fft_ok($_) true") for 1, 2, 12, 157500, 165888, 360000;

    my ($all, $bad) = (0, 0);
    for my $m (1000 .. 2_000_000) {
        next unless is235($m);
        $all++;
        $bad++ unless fft_ok($m);
    }
    is("$all $bad", "493 51", '1000..2e6: 51 of the 493 lengths 2^a 3^b 5^c exceed the limit');

    my $ok = 1;
    for my $n (1000 .. 500_000) {
        next unless is235($n) && !fft_ok($n);
        my $L = fft_len($n);
        $ok = 0, last unless $L > $n && is235($L) && fft_ok($L);
    }
    ok($ok, 'fft_len skips every length PDL::FFT cannot transform up to 500000');
    is(fft_len(351000 + 1216), 360000, 'fft_len(352216) = 360000 (not 354294)');
    is(fft_len(1000), 1000, 'fft_len keeps a length that is already fine');
};

subtest 'fft_inplace against a direct DFT' => sub {
    my @backends = ('auto', ($HAVE_FFTW3 ? 'fftw3' : ()));
    for my $be (@backends) {
        local $PDL::EEG::FFT::Dispatcher::BACKEND = $be;
        for my $n (12, 7) {
            my $x = random(double, $n, 3);
            my ($er, $ei) = dft($x);
            my ($re, $im) = ($x->copy, zeroes(double, $n, 3));
            fft_inplace($re, $im);
            ok(max(abs($re - $er)) < 1e-12 && max(abs($im - $ei)) < 1e-12, "$be: forward, n=$n");
            ifft_inplace($re, $im);
            ok(max(abs($re - $x)) < 1e-12 && max(abs($im)) < 1e-12, "$be: inverse returns the input, n=$n");
        }
    }
};

# 'auto' prefers PDL::FFTW3 for every length when it is installed; without it,
# PDL::FFT for the lengths it transforms, a croak for the others. The case
# without PDL::FFTW3 runs in a child perl where loading PDL::FFTW3 fails.
subtest "BACKEND 'auto': which library" => sub {
    my $child = sub {
        my ($code) = @_;
        my $hide = 'BEGIN { unshift @INC, sub { die "hidden\n" if $_[1] eq "PDL/FFTW3.pm"; return } }';
        my @inc = map { "-I$_" } grep { !ref } @INC;
        open my $fh, '-|', $^X, @inc, '-e', "$hide use PDL; use PDL::EEG::FFT::Dispatcher qw(fft_backend); $code"
            or die "child: $!";
        local $/;
        my $out = <$fh>;
        close $fh;
        return $out;
    };
    is($child->('print fft_backend(1000), " ", (eval { fft_backend(' . $BAD . ') } // "croak"), "\n"'),
       "PDL::FFT croak\n", "without PDL::FFTW3: PDL::FFT for 1000, croak for $BAD");
    SKIP: {
        skip 'PDL::FFTW3 not installed', 2 unless $HAVE_FFTW3;
        local $PDL::EEG::FFT::Dispatcher::BACKEND = 'auto';
        is(fft_backend(1000), 'PDL::FFTW3', 'with PDL::FFTW3: PDL::FFTW3 for 1000 too');
        local $PDL::EEG::FFT::Dispatcher::BACKEND = 'pdl';
        is(fft_backend(1000), 'PDL::FFT', "BACKEND 'pdl': PDL::FFT for 1000");
    }
};

subtest "length $BAD" => sub {
    {
        local $PDL::EEG::FFT::Dispatcher::BACKEND = 'pdl';
        my ($re, $im) = (sin(sequence(double, $BAD) * 0.37), zeroes(double, $BAD));
        eval { fft_inplace($re, $im) };
        like($@, qr/cannot transform length $BAD/, "BACKEND 'pdl': croaks instead of calling PDL::FFT");
    }
    SKIP: {
        skip 'PDL::FFTW3 not installed', 3 unless $HAVE_FFTW3;
        is(fft_backend($BAD), 'PDL::FFTW3', 'auto: PDL::FFTW3');
        my $x = sin(sequence(double, $BAD) * 0.37);
        my ($re, $im) = ($x->copy, zeroes(double, $BAD));
        fft_inplace($re, $im);
        my $q = (sum($re * $re + $im * $im) / $BAD / sum($x * $x))->sclr;
        ok(abs($q - 1) < 1e-9, 'auto: Parseval holds');
        ifft_inplace($re, $im);
        ok(max(abs($re - $x)) < 1e-9, 'auto: inverse returns the input');
    }
};

# The case found in HBrBP003: zp_butter on 351000 samples with 600 samples of
# padding used to choose length 354294 and return a broken result.
subtest 'zp_butter, nt = 351000, PDL::FFT only' => sub {
    local $PDL::EEG::FFT::Dispatcher::BACKEND = 'pdl';
    my $fs = 1000;
    my $x = sin(2 * $PI * 10 * sequence(double, 351000) / $fs);
    my $y = zp_butter($x, fs => $fs, low => 5, high => 20, order => 2, pad => 0.6);
    my $mid = "2000:-2001";
    ok(max(abs($y->slice($mid) - $x->slice($mid))) < 1e-3, '10 Hz passes a 5-20 Hz band-pass unchanged');
};

subtest 'detect_r, nt = 351000 (HBrBP003)' => sub {
    local $PDL::EEG::FFT::Dispatcher::BACKEND = 'pdl';
    # synthetic ECG at 1000 Hz: R waves about 0.9 s apart, a T wave 300 ms
    # after each, 0.25 Hz drift and 50 Hz. The 5-20 Hz band-pass inside
    # detect_r used to run at length 354294; 30 of the 389 beats were found.
    my ($fs, $nt) = (1000, 351000);
    my $t = sequence(double, $nt) / $fs;
    my $x = 0.3 * sin(2 * $PI * 0.25 * $t) + 0.05 * sin(2 * $PI * 50 * $t);
    my ($c, $n) = (0.5, 0);
    while ($c < $nt / $fs - 0.5) {
        my ($lo, $hi) = (int(($c - 0.1) * $fs), int(($c + 0.5) * $fs));
        my $s = $t->slice("$lo:$hi");
        $x->slice("$lo:$hi") += exp(-0.5 * (($s - $c) / 0.012) ** 2)
                              + 0.25 * exp(-0.5 * (($s - $c - 0.3) / 0.06) ** 2);
        $n++;
        $c += 0.9 + 0.05 * sin(0.7 * $n);
    }
    my ($r) = detect_r($x, $fs);
    is($r->nelem, $n, "all $n beats found");
};

subtest 'cwt_morlet, nfft that used to be 157464' => sub {
    local $PDL::EEG::FFT::Dispatcher::BACKEND = 'pdl';
    my ($sf, $f) = (1000, 10);
    my $half = int(5 * 7 / (2 * $PI * $f) * $sf);              # 557: longest kernel (lowest freq)
    my $nt = $BAD - $half;                                      # need = 157464
    my $x = sin(2 * $PI * $f * sequence(double, $nt) / $sf);
    my ($cr, $ci) = cwt_morlet($x, $sf, pdl(10, 20, 40), n_cycles => 7);
    my $p = ($cr * $cr + $ci * $ci)->slice(":,5000:-5001");    # (3, n)
    my $m = $p->slice("(0)")->avg;
    ok(abs($m / max($p->slice("(0)")) - 1) < 1e-3 && $m > 100 * max($p->slice("1:2")),
       'a 10 Hz tone: constant power at 10 Hz, little at 20 and 40 Hz');
};

subtest "exact lengths: analytic_signal(exact => 1), resp_phase" => sub {
    my $x = cos(2 * $PI * 1000 * sequence(double, $BAD) / $BAD);   # whole number of cycles
    my $nb = 39366;                                                 # beats about 1 s apart
    my $t  = sequence(double, $nb) * (39365.8 / ($nb - 1));         # span 39365.8 s
    my $y  = sin(2 * $PI * 0.3 * $t);                               # grid at 4 Hz: 157464 samples
    {
        local $PDL::EEG::FFT::Dispatcher::BACKEND = 'pdl';
        eval { analytic_signal($x, exact => 1) };
        like($@, qr/cannot transform length $BAD/, "analytic_signal(exact): croaks with BACKEND 'pdl'");
        eval { resp_phase($t, $y, fs => 4) };
        like($@, qr/cannot transform length $BAD/, "resp_phase: croaks with BACKEND 'pdl'");
    }
    SKIP: {
        skip 'PDL::FFTW3 not installed', 2 unless $HAVE_FFTW3;
        my ($re, $im) = analytic_signal($x, exact => 1);
        ok(max(abs(sqrt($re * $re + $im * $im) - 1)) < 1e-9, 'analytic_signal(exact): envelope 1');
        my $ph = resp_phase($t, $y, fs => 4);
        ok(abs($ph->{peak_hz} - 0.3) < 1e-3, 'resp_phase: peak at 0.3 Hz');
    }
};

done_testing;
