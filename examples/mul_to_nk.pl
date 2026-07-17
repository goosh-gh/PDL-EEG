#!/usr/bin/env perl
# examples/mul_to_nk.pl
#
#   perl -Ilib examples/mul_to_nk.pl vendor.m01                  # just inspect it
#   perl -Ilib examples/mul_to_nk.pl vendor.m01 --eeg subject.EEG [--bne[=0.5]]
#
# Closes the loop. The Nihon Kohden viewer exports .mul itself, so its own export
# is a ground truth we did not write: read the .mul back, read the .EEG with
# read_nk, line them up, and diff. If our reader has the block boundaries, the
# channel order or the gains wrong, this says so in microvolts.
#
# That matters because every bug found in this file so far -- the 442-byte block
# headers read as samples, the DC gain off by 1000x, the DC channel numbering --
# was invisible until something independent disagreed.
#
# The vendor's .mul is BNE re-referenced (labels carry "-BN"), so to compare like
# with like, pass --bne and read_nk's data goes through PDL::EEG::Derivation::bne
# first. Without --bne only the DC and Trigger columns are comparable.
#
# ALIGNMENT. The vendor's export is a RANGE THE OPERATOR SELECTED IN THE VIEWER,
# by hand, watching the trace until the amplifiers settle and the trigger line
# stops swinging. Both ends are arbitrary. It is NOT "the block, minus two
# seconds" -- across one recording the trims run 2 s, 3 s, 16 s, and one export
# is LONGER than the segment its Time= falls in, i.e. the operator dragged the
# selection straight across a recording break.
#
# So: read the whole recording (all_blocks => 1), take the .mul's Time= only as a
# starting hint, and then slide the .mul against the concatenation to find the
# offset that actually minimises the error. A .mul that spans a break is the most
# useful one of all -- it is the vendor telling us, sample by sample, what the
# concatenation across that break is supposed to look like.

use strict;
use warnings;
use PDL;
use PDL::EEG::IO::NihonKohden  qw(read_nk block_extents);
use PDL::EEG::IO::BESA::ASCII  qw(read_mul);
use Getopt::Long;
use POSIX ();

binmode(STDOUT, ':encoding(UTF-8)');

my ($eeg, $bne, $max_shift, $nshow, $solve);
$max_shift = 5;          # seconds to search either way when aligning
$nshow     = 12;         # channels to list
GetOptions('eeg=s' => \$eeg, 'bne:s' => \$bne, 'solve-bne' => \$solve,
           'max-shift=f' => \$max_shift, 'show=i' => \$nshow)
    or die "bad options\n";
my $mulf = shift @ARGV
    or die "usage: $0 file.mul [--eeg subject.EEG] [--bne[=0.5]] [--solve-bne]\n";

# An empty --eeg is a shell variable that did not expand, not "no --eeg". Say so,
# rather than sailing on and croaking 30 lines later with a blank filename.
if (defined $eeg && $eeg !~ /\S/) {
    die "--eeg was given but is EMPTY -- a shell variable that did not expand?\n"
      . "  e.g.  --eeg \"\$EEG\"   with EEG unset.\n";
}
die "no such .mul file: '$mulf'\n" unless -f $mulf;
die "--eeg: no such file: '$eeg'\n" if defined $eeg && !-f $eeg;

# ---------------------------------------------------------------------------
my $m = read_mul($mulf);
printf "mul    : %s\n", $mulf;
printf "  %d channels x %d samples @ %g Hz   (%.1f s)\n",
    $m->{n_ch}, $m->{data}->dim(1), $m->{fs}, $m->{data}->dim(1) / $m->{fs};
printf "  Channels= %d   columns present %d   -> %s\n",
    $m->{n_report}, $m->{n_ch},
    ($m->{n_report} == $m->{n_ch}
        ? 'the Trigger column IS counted'
        : sprintf('%d column(s) not counted (Trigger excluded?)',
                  $m->{n_ch} - $m->{n_report}));
printf "  Bins/uV=%g   BeginSweep=%s   Date/Time cols=%d\n",
    $m->{bins_per_uv}, ($m->{begin_ms} // '-'), $m->{date_time};
printf "  t_start: %s\n", ($m->{t_start} // '(none)');
printf "  labels : %s\n", join(',', @{ $m->{labels} });

# the vendor names all four DC columns "Fp1"; say so rather than let it confuse
{
    my %seen;
    my @dup = grep { $seen{$_}++ } @{ $m->{labels} };
    my %u; @dup = grep { !$u{$_}++ } @dup;
    printf "  !! duplicate label(s) in the .mul: %s\n"
         . "     (the vendor names every DC column after the first channel --\n"
         . "      match those columns by POSITION, never by name)\n",
        join(',', @dup) if @dup;
}
printf "  Trigger column: %s\n\n",
    defined $m->{trig_idx} ? "index $m->{trig_idx}" : 'none found';

exit 0 unless defined $eeg;

# ---------------------------------------------------------------------------
# Which segment of the .EEG does this .mul come from?
# ---------------------------------------------------------------------------
sub epoch {
    my ($Y,$M,$D,$h,$mi,$s) = ($_[0] // '') =~
        /^(\d{4})-(\d{2})-(\d{2})[ T](\d{2}):(\d{2}):(\d{2})/ or return undef;
    POSIX::mktime($s, $mi, $h, $D, $M - 1, $Y - 1900, 0, 0, -1);
}

my $ext = block_extents($eeg);
printf "eeg    : %s\n  %d segment%s\n", $eeg, scalar @$ext, (@$ext == 1 ? '' : 's');
printf "    seg %-2d  %8d samp  %s\n", @{$_}{qw(index n_samp t_start)} for @$ext;

my $me = epoch($m->{t_start});
my ($seg) = sort {
    abs(($me // 0) - (epoch($a->{t_start}) // 0))
        <=> abs(($me // 0) - (epoch($b->{t_start}) // 0))
} @$ext;
die "cannot match the .mul's Time= against any segment\n"
    unless $seg && defined $me;
my $lag_in_seg = $me - epoch($seg->{t_start});
my $hint = $seg->{start_samp} + $lag_in_seg * $m->{fs};   # in CONCATENATED samples

printf "\n  Time= lands %+.0f s into segment %d (%s)\n",
    $lag_in_seg, $seg->{index}, $seg->{t_start};
printf "  the .mul is %.1f s long; that segment is %.1f s%s\n",
    $m->{data}->dim(1) / $m->{fs}, $seg->{n_samp} / $m->{fs},
    (($lag_in_seg * $m->{fs} + $m->{data}->dim(1)) > $seg->{n_samp}
        ? "  -> the .mul RUNS PAST THE END of it: the operator dragged the\n"
        . "     selection across a recording break. This is the interesting case:\n"
        . "     it tells us what the concatenation across that break should be."
        : '');

# Read the WHOLE recording, not one segment: the .mul may straddle a break.
my $rec = read_nk($eeg, all_blocks => 1);
if (defined $bne) {
    require PDL::EEG::Derivation;
    my $prop = length($bne) ? $bne + 0 : 0.5;
    $rec = PDL::EEG::Derivation::bne($rec, prop => $prop)
        or die "--bne failed (need BN1/BN2)\n";
    printf "  read_nk: BNE re-referenced (prop=%.3g)\n", $prop;
}
my @nk_lab = @{ $rec->{labels} };
printf "  read_nk: %d ch x %d samp\n", $rec->{data}->dim(0), $rec->{data}->dim(1);

# ---------------------------------------------------------------------------
# Match columns. Names first (after stripping the montage suffix); the DC columns
# have to go by position, because the vendor calls them all "Fp1".
# ---------------------------------------------------------------------------
my %nk_ix;
for my $i (0 .. $#nk_lab) {
    (my $bare = $nk_lab[$i]) =~ s/-[A-Za-z0-9]+$//;
    $nk_ix{ lc $bare } //= $i;
    $nk_ix{ lc $nk_lab[$i] } //= $i;
}
my @pair;                       # [ mul_col, nk_col, label ]
for my $c (0 .. $m->{n_ch} - 1) {
    (my $bare = $m->{labels}[$c]) =~ s/-[A-Za-z0-9]+$//;
    my $j = $nk_ix{ lc $bare };
    push @pair, [ $c, $j, $m->{labels}[$c] ] if defined $j;
}
printf "  matched %d of %d .mul columns by name\n\n", scalar @pair, $m->{n_ch};

# ---------------------------------------------------------------------------
# Find the sample offset that actually minimises the error, instead of trusting
# the 2 s. Use one clean scalp channel.
# ---------------------------------------------------------------------------
# Align on a CHANNEL DIFFERENCE, not on a channel. Whatever reference the vendor
# used, it is common to every channel, so it cancels in (ch_a - ch_b) -- and the
# alignment stops depending on getting the reference right, which is the thing we
# are about to measure.
my @scalp = grep { $_->[2] =~ /-BN$/ } @pair;
die "no -BN scalp columns to align on\n" unless @scalp >= 2;
my ($pa, $pb) = @scalp[0, $#scalp];

my $nm  = $m->{data}->dim(1);
my $nn  = $rec->{data}->dim(1);
my $win = 20 * $m->{fs};
$win = $nm if $win > $nm;

my ($best_off, $best_err) = (undef, undef);
my $lo = int($hint - $max_shift * $m->{fs});
my $hi = int($hint + $max_shift * $m->{fs});
$lo = 0 if $lo < 0;
my $mdiff = $m->{data}->slice("($pa->[0]),0:" . ($win - 1))
          - $m->{data}->slice("($pb->[0]),0:" . ($win - 1));
for (my $off = $lo; $off <= $hi; $off += 1) {
    last if $off + $win > $nn;
    my $ndiff = $rec->{data}->slice("($pa->[1])," . $off . ":" . ($off + $win - 1))
              - $rec->{data}->slice("($pb->[1])," . $off . ":" . ($off + $win - 1));
    my $e = (($mdiff - $ndiff) ** 2)->sum;
    if (!defined $best_err || $e < $best_err) { ($best_err, $best_off) = ($e, $off) }
}
die "could not align (searched +/-${max_shift}s around " . ($hint / $m->{fs}) . " s;"
  . " try --max-shift)\n" unless defined $best_off;
printf "  alignment: the .mul starts at concatenated sample %d (%.3f s)%s\n",
    $best_off, $best_off / $m->{fs},
    (abs($best_off - $hint) > 0.5 * $m->{fs} ? '   <- not where Time= pointed' : '');

# does the aligned .mul cross any recording break?
my @crossed = grep { $_->{start_samp} > $best_off
                  && $_->{start_samp} < $best_off + $m->{data}->dim(1) } @$ext;
printf "  it spans %d recording break%s: %s\n",
    scalar @crossed, (@crossed == 1 ? '' : 's'),
    (@crossed ? join(', ', map { sprintf('into segment %d at %.3f s of the .mul',
                                         $_->{index},
                                         ($_->{start_samp} - $best_off) / $m->{fs})
                               } @crossed)
              : 'none (it sits inside one segment)')
    if @$ext > 1;
print "\n";

# ---------------------------------------------------------------------------
# Diff.
# ---------------------------------------------------------------------------
my $n = $nm;
$n = $nn - $best_off if $best_off + $n > $nn;
printf "diff over %d samples (%.1f s)\n", $n, $n / $m->{fs};
printf "  %-12s %-10s %12s %12s %10s  %s\n",
    'mul label', 'nk label', 'max |diff|', 'rms', 'nk gain', 'verdict';

my $bad = 0;
for my $p (@pair) {
    my ($c, $j, $lab) = @$p;
    my $a = $m->{data}->slice("($c),0:" . ($n - 1));
    my $b = $rec->{data}->slice("($j)," . $best_off . ":" . ($best_off + $n - 1));
    my $d = ($a - $b)->abs;
    my $mx  = $d->max;
    my $rms = sqrt((($a - $b) ** 2)->avg);
    my $g   = eval { $rec->{gains}->at($j) } // 1;
    # one ADC step is the floor: the .mul is printed to 2 decimals, we are not
    my $tol = 2 * $g + 0.01;
    my $ok  = $mx <= $tol;
    $bad++ unless $ok;
    printf "  %-12s %-10s %12.3f %12.3f %10.4g  %s\n",
        $lab, $nk_lab[$j], $mx, $rms, $g,
        $ok ? 'match' : '*** DIFFERS ***';
}

if (@crossed) {
    print "\n";
    printf "This .mul crosses %d break(s). If the channels above match, the vendor\n"
         . "and read_nk agree on what the concatenation across a recording break\n"
         . "looks like -- sample for sample. That is the strongest confirmation\n"
         . "available that the 442-byte block headers are being skipped correctly\n"
         . "and that nothing is padded or dropped at a boundary.\n", scalar @crossed;
}

# ---------------------------------------------------------------------------
# --solve-bne : measure the recorder's BN balance against the vendor's own export
#
# The residual above is the SAME SIGNAL on every scalp channel. That is not a
# waveform error -- it is a REFERENCE error: something common has been subtracted
# from all of them, and we subtracted a slightly different something.
#
#     vendor : mul_c = raw_c - (p*BN1 + (1-p)*BN2)      p unknown
#     ours   : nk_c  = raw_c - (0.5*BN1 + 0.5*BN2)
#     so     : raw_c - mul_c = p*BN1 + (1-p)*BN2   -- identical for every c
#
# So take D = raw_c - mul_c (from the UN-re-referenced read), and least-squares it
# onto BN1 and BN2. The coefficients ARE the balance. If they sum to 1, the BNE
# model is confirmed; if a != 0.5, the hardware pot is not centred -- which is
# exactly the number we could not get out of the file itself.
# ---------------------------------------------------------------------------
if ($solve) {
    print "\n--- solving for the BN balance against the vendor's export ---\n";
    my $raw = read_nk($eeg, all_blocks => 1);
    my @rl  = @{ $raw->{labels} };
    my %rx; $rx{ lc $rl[$_] } = $_ for 0 .. $#rl;
    my ($i1, $i2) = ($rx{bn1}, $rx{bn2});
    die "no BN1/BN2 channels in the .EEG\n" unless defined $i1 && defined $i2;

    my $N   = $n;
    my $sl  = sub { $raw->{data}->slice("($_[0])," . $best_off . ":"
                                        . ($best_off + $N - 1))->double };
    my $bn1 = $sl->($i1);
    my $bn2 = $sl->($i2);

    # D from every scalp channel; they must agree, so check that first
    my (@D, @names);
    for my $p (@scalp) {
        (my $bare = $p->[2]) =~ s/-BN$//;
        my $j = $rx{ lc $bare } // next;
        push @D, $sl->($j) - $m->{data}->slice("($p->[0]),0:" . ($N - 1))->double;
        push @names, $bare;
    }
    die "no scalp channels resolved for the fit\n" unless @D >= 2;

    my $spread = 0;
    for my $k (1 .. $#D) {
        my $d = ($D[$k] - $D[0])->abs->max->sclr;
        $spread = $d if $d > $spread;
    }
    printf "  D = raw - mul, over %d scalp channels: they agree to %.3f uV\n"
         . "  -> %s\n", scalar @D, $spread,
        ($spread < 1.0
            ? 'CONFIRMED: the difference is one common signal, i.e. a reference'
            : 'NOT a pure reference difference -- something else is going on too');

    my $Dm = $D[0];
    for my $k (1 .. $#D) { $Dm = $Dm + $D[$k] }
    $Dm = $Dm / scalar(@D);                     # average away the quantisation

    # ---- solve  D = a*BN1 + b*BN2 + c  --------------------------------------
    #
    # NOT in the BN1/BN2 basis. BN1 and BN2 are both non-cephalic electrodes and
    # they are strongly correlated, so the normal equations in that basis are
    # nearly singular -- solving them directly gives coefficients of order 1e9 and
    # a residual larger than the data, which is exactly what the first version of
    # this did.
    #
    # Re-span the same plane with the SUM and the DIFFERENCE:
    #
    #     u = (BN1 + BN2)/2      the common part -- large
    #     v = (BN1 - BN2)/2      the differential part -- small, but it is the ONLY
    #                            thing that carries information about the balance
    #
    #     a*BN1 + b*BN2  =  (a+b)*u + (a-b)*v  =  alpha*u + beta*v
    #     a = (alpha+beta)/2,  b = (alpha-beta)/2
    #
    # u and v are near-orthogonal, so the 2x2 closes in closed form and stays sane.
    my $u = ($bn1 + $bn2) / 2;
    my $v = ($bn1 - $bn2) / 2;

    my $cen = sub { my $x = shift; $x - $x->avg };
    my ($uc, $vc, $yc) = (map { $cen->($_) } $u, $v, $Dm);

    # ->sclr, every time. ->sum returns a PDL SCALAR PIDDLE, not a Perl number, and
    # "my $piv = $A[$k][$k]" then aliases the piddle rather than copying its value
    # -- so "$A[$k][$k] /= $piv" divides the pivot by itself IN PLACE, leaves 1
    # behind, and every later "/= $piv" divides by 1. That is what wrecked the
    # first version of this solve (coefficients of order 1e9). Get Perl numbers out
    # of PDL before doing scalar linear algebra.
    my $Suu = ($uc * $uc)->sum->sclr;
    my $Svv = ($vc * $vc)->sum->sclr;
    my $Suv = ($uc * $vc)->sum->sclr;
    my $Suy = ($uc * $yc)->sum->sclr;
    my $Svy = ($vc * $yc)->sum->sclr;
    my $det = $Suu * $Svv - $Suv * $Suv;

    my ($b1c, $b2c) = ($cen->($bn1), $cen->($bn2));
    my $r12 = ($b1c * $b2c)->sum->sclr
            / (sqrt((($b1c ** 2)->sum->sclr) * (($b2c ** 2)->sum->sclr)) + 1e-30);
    printf "\n  BN1 vs BN2 : r = %.4f   rms(BN1-BN2) = %.3f uV\n",
        $r12, sqrt((($bn1 - $bn2) ** 2)->avg->sclr);
    print "  (if BN1 and BN2 were identical the balance would be unmeasurable --\n"
        . "   only their DIFFERENCE carries it)\n";

    # BN1 == BN2 makes Svv zero, and then "det < 1e-9 * Suu*Svv" is "0 < 0", which
    # is false -- the guard never fires and the next line divides by zero. Test the
    # DIFFERENTIAL energy directly: it is the only thing that carries the balance.
    if ($Svv <= 0 || $Svv < 1e-10 * $Suu) {
        print "\n  *** BN1 and BN2 are the SAME SIGNAL in this recording ***\n"
            . "      rms(BN1-BN2) = 0, so p*BN1 + (1-p)*BN2 = BN1 for ANY p.\n"
            . "      The balance is not measurable here -- and it does not matter:\n"
            . "      every value of prop gives the identical reference. bne() with\n"
            . "      any prop is correct for this file.\n"
            . "      (Whether the two BN electrodes were physically the same point,\n"
            . "       or only one was connected, is a question for the recording.)\n";

        # Still worth checking the model: D should just BE BN1, plus an offset.
        my $b1c  = $cen->($bn1);
        my $S11  = ($b1c * $b1c)->sum->sclr;
        if ($S11 > 0) {
            my $a  = ($b1c * $yc)->sum->sclr / $S11;
            my $c0 = $Dm->avg->sclr - $a * $bn1->avg->sclr;
            my $res  = sqrt((($Dm - ($bn1 * $a + $c0)) ** 2)->avg->sclr);
            my $drms = sqrt((($cen->($Dm)) ** 2)->avg->sclr);
            printf "\n  D  =  %.4f * BN1  +  %.2f uV\n", $a, $c0;
            printf "  residual after the fit: %.3f uV rms   (D is %.3f uV rms)\n",
                $res, $drms;
            printf "  the fit explains %.3f%% of D  -> %s\n",
                100 * (1 - ($res / ($drms + 1e-30)) ** 2),
                (abs($a - 1) < 0.02 && $res < 0.3 * $drms
                    ? 'the reference IS the BN electrode. Everything checks out.'
                    : 'the reference is NOT simply the BN electrode -- look further.');
        }
    }
    elsif (abs($det) < 1e-12 * $Suu * $Svv) {
        print "\n  *** the normal equations are singular: cannot solve ***\n";
    }
    else {
        my $alpha = ($Svv * $Suy - $Suv * $Svy) / $det;
        my $beta  = ($Suu * $Svy - $Suv * $Suy) / $det;
        my $a     = ($alpha + $beta) / 2;
        my $b     = ($alpha - $beta) / 2;
        my $c0    = $Dm->avg->sclr - $a * $bn1->avg->sclr - $b * $bn2->avg->sclr;

        my $fit  = $bn1 * $a + $bn2 * $b + $c0;
        my $res  = sqrt((($Dm - $fit) ** 2)->avg->sclr);
        my $drms = sqrt((($cen->($Dm)) ** 2)->avg->sclr);

        printf "\n  D  =  %.4f * BN1  +  %.4f * BN2  +  %.2f uV\n", $a, $b, $c0;
        printf "  a + b = %.4f   %s\n", $a + $b,
            (abs($a + $b - 1) < 0.02
                ? '(== 1: the acquisition reference cancels, so the BNE model holds)'
                : '(NOT 1 -- the BNE model does not describe what the vendor did)');
        printf "  residual after the fit: %.3f uV rms   (D is %.3f uV rms)\n",
            $res, $drms;
        printf "  the fit explains %.3f%% of D\n",
            100 * (1 - ($res / ($drms + 1e-30)) ** 2);

        # a sanity cross-check that does not go through the regression at all
        printf "  cross-check: rms((p-0.5)*(BN1-BN2)) = %.3f uV -- this is what the\n"
             . "               residual should be if you re-run with --bne (prop=0.5)\n",
            abs($a - 0.5) * sqrt((($bn1 - $bn2) ** 2)->avg->sclr);

        if (abs($a + $b - 1) < 0.02 && $res < 0.3 * $drms) {
            printf "\n  *** BN BALANCE  prop = %.4f ***   (we have been assuming 0.500)\n", $a;
            printf "      y = x - (%.4f*BN1 + %.4f*BN2)\n", $a, 1 - $a;
            printf "      re-run with --bne=%.4f; the scalp channels should then match.\n", $a;
        } else {
            print "\n  The fit does not support a simple BNE reference. What the vendor\n"
                . "  subtracted is something else -- look at what D actually is.\n";
        }
    }
}

print "\n";
if ($bad) {
    printf "%d channel(s) disagree with the vendor's own export.\n", $bad;
    print  "A constant ratio between them is a GAIN error; a constant offset is a\n"
         . "REFERENCE error (try --bne); a channel that matches a DIFFERENT one of\n"
         . "ours is a CHANNEL ORDER error.\n";
} else {
    print "Every matched channel agrees with the vendor's export to within one ADC\n"
        . "step. The block boundaries, the channel order and the gains are right.\n";
}
