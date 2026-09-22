#!/usr/bin/env perl
# Apply pure-PDL FastICA (PDL::EEG::ICA) to a continuous EDF recording, identify
# the ocular components from the EOG channels, and write the recording back out
# as BESA .mul with those component(s) removed by back-projection.
#
# Ocular components are identified PHYSIOLOGICALLY, not by a single ranking:
#   #1 blink       = component most correlated with the vertical EOG   (vEOG)
#   #2 horizontal  = component most correlated with the horizontal EOG (hEOG)
#   #3 next ocular = next component by max(|corr vEOG|,|corr hEOG|)
# (Vertical blinks project onto vEOG; horizontal saccades onto hEOG, and the two
#  EOGs are ~orthogonal, so ranking by vEOG alone will NOT surface the
#  horizontal component -- it needs hEOG.)
#
# The ICA is fitted on the scalp EEG together with the EOG channels, so removing
# an ocular component subtracts it from the EOG channels too (the effect on the
# EOG is visible in the output). Every non-ICA channel (DC / trigger inputs,
# Mark*, Events, ...) is passed through UNCHANGED, so the output still contains
# everything you need to epoch and average.
#
# Three files are written (base = input name with .edf stripped):
#   <base>-ica-1.mul       #1 removed              (blink)
#   <base>-ica-1-2.mul     #1 and #2 removed       (blink + horizontal)
#   <base>-ica-1-2-3.mul   #1, #2 and #3 removed   (+ next ocular)
#
# If no hEOG is present the script falls back to ranking by vEOG alone and #2/#3
# are simply the next vEOG-correlated components (reported as such).
#
# Usage (from the PDL-EEG repo, with the ICA dist on @INC):
#   perl -Ilib -I/path/to/PDL-EEG-ICA/lib examples/task_ica_mul.pl task.edf
#
# Options:
#   --veog LABEL     vertical-EOG channel   (default: auto-detect vEOG/EOG)
#   --heog LABEL     horizontal-EOG channel (default: auto-detect hEOG)
#   --eeg  a,b,c     explicit scalp-EEG channels; replaces the whole default set
#                    (default: 10-20 19ch + any of E,A1,A2,X1,BN1,BN2 present)
#   --extra a,b,c    non-10-20 electrodes to add to the default EEG set when
#                    present (default: E,A1,A2,X1,BN1,BN2; --extra "" for none)
#   --eog  a,b,c     explicit EOG channels to include in the ICA mixture
#                    (default: auto-detect any label matching /eog/i)
#   --no-eog         do NOT include EOG channels in the ICA mixture
#   --nl tanh|gauss  FastICA contrast (default tanh; gauss is more robust to
#                    spiky super-Gaussian sources such as blinks)
#   --seed N         RNG seed (default 1)
#   --keep N         reduce to top-N PCs before ICA (default: full rank)
#   --maxit N        max FastICA iterations (default 500)
#   --tol X          convergence tolerance (default 1e-6)
#   --decimals N     .mul value precision (default 2)
#   --fieldw N       .mul column width (default: auto, wide enough for DC)
#   --prefix STR     output base name (default: input with .edf stripped)

use strict;
use warnings;
use PDL;
use PDL::EEG::IO::EDF         qw(read_edf clean_edf_label);
use PDL::EEG::IO::BESA::ASCII qw(write_mul);
use PDL::EEG::ICA             qw(ica_decompose identify_by_reference apply_ica);
use Getopt::Long;

my %o = (nl=>'tanh', seed=>1, maxit=>500, tol=>1e-6, decimals=>2);
GetOptions(\%o, qw(veog=s heog=s eeg=s extra=s eog=s no-eog nl=s seed=i keep=i
                   maxit=i tol=f decimals=i fieldw=i prefix=s)) or die "bad options\n";
my $in = shift @ARGV or die "usage: $0 [options] file.edf\n";

# ---- canonical 10-20 scalp set, with common alias spellings -----------------
my @TEN20 = qw(Fp1 Fp2 F7 F3 Fz F4 F8 T7 C3 Cz C4 T8 P7 P3 Pz P4 P8 O1 O2);
my %ALIAS = (T3=>'T7', T4=>'T8', T5=>'P7', T6=>'P8');
sub canon {
    my $l = shift; $l =~ s/\s+//g;
    for my $c (@TEN20)      { return $c         if lc $l eq lc $c }
    for my $a (keys %ALIAS) { return $ALIAS{$a} if lc $l eq lc $a }
    return undef;
}
# additional non-10-20 electrodes to treat as EEG for the ICA when present
my @EEG_EXTRA = defined $o{extra} ? grep { length } split /,/, $o{extra}
                                  : qw(E A1 A2 X1 BN1 BN2);

sub rms { my $x = shift; my $c = $x - $x->avg; sqrt(($c*$c)->avg) }
sub corr { my ($a,$b) = @_; my $a0=$a-$a->avg; my $b0=$b-$b->avg;
           (($a0*$b0)->avg)/(rms($a)*rms($b)+1e-30) }

# ---- read EDF ---------------------------------------------------------------
my $rec   = read_edf($in);
my $fs    = $rec->{fs};
my $data  = $rec->{data};                          # (n_ch, n_samp) uV
my @clean = map { clean_edf_label($_) } @{ $rec->{labels} };
my %idx;  $idx{$clean[$_]} = $_ for 0 .. $#clean;
printf "read %s : %d ch @ %g Hz, %d samp (%.1f s)\n",
    $in, $data->dim(0), $fs, $data->dim(1), $data->dim(1)/$fs;

# ---- EEG channels for the ICA ----------------------------------------------
my (@eeg_idx, @eeg_lab);
if (defined $o{eeg}) {
    for my $name (split /,/, $o{eeg}) {
        exists $idx{$name} or die "no such EEG channel: $name\n";
        push @eeg_idx, $idx{$name}; push @eeg_lab, $name;
    }
} else {
    my %used;
    for my $want (@TEN20) {                        # keep 10-20 order
        for my $i (0 .. $#clean) {
            next if $used{$i};
            my $c = canon($clean[$i]);
            if (defined $c && $c eq $want) {
                push @eeg_idx, $i; push @eeg_lab, $want; $used{$i} = 1; last;
            }
        }
    }
    for my $x (@EEG_EXTRA) {                        # extra electrodes present
        for my $i (0 .. $#clean) {
            next if $used{$i};
            if (lc $clean[$i] eq lc $x) {
                push @eeg_idx, $i; push @eeg_lab, $clean[$i]; $used{$i} = 1; last;
            }
        }
    }
}
die "found no scalp EEG channels\n" unless @eeg_idx;

# ---- EOG channels to include in the ICA mixture -----------------------------
my (@eog_idx, @eog_lab);
unless ($o{'no-eog'}) {
    if (defined $o{eog}) {
        for my $name (split /,/, $o{eog}) {
            exists $idx{$name} or die "no such EOG channel: $name\n";
            push @eog_idx, $idx{$name}; push @eog_lab, $name;
        }
    } else {
        my %in_eeg = map { $_ => 1 } @eeg_idx;
        for my $i (0 .. $#clean) {
            next if $in_eeg{$i};
            if ($clean[$i] =~ /eog/i) { push @eog_idx, $i; push @eog_lab, $clean[$i] }
        }
    }
}

my @ica_idx = (@eeg_idx, @eog_idx);
my @ica_lab = (@eeg_lab, @eog_lab);
printf "ICA channels (%d): %s\n", scalar @ica_lab, join(' ', @ica_lab);
my %in_ica = map { $_ => 1 } @ica_idx;
my @pass = grep { !$in_ica{$_} } 0 .. $#clean;
printf "passed through unchanged (%d): %s\n",
    scalar @pass, (@pass ? join(' ', @clean[@pass]) : '(none)');

my $Xsub = $data->dice_axis(0, pdl(long, @ica_idx))->sever;  # (n_ica, nt)
my $X    = $Xsub->transpose->sever;                          # (nt, n_ica) for ICA

# ---- locate vEOG and hEOG ---------------------------------------------------
my $find = sub { my $re = shift; for my $i (0..$#clean){ return $i if $clean[$i] =~ $re } undef };
my $veog_i = defined $o{veog} ? ($idx{$o{veog}} // die "no such channel: $o{veog}\n")
           : ($find->(qr/v.?eog/i) // $find->(qr/eog/i));
my $heog_i = defined $o{heog} ? ($idx{$o{heog}} // die "no such channel: $o{heog}\n")
           : $find->(qr/h.?eog/i);
$heog_i = undef if defined $heog_i && defined $veog_i && $heog_i == $veog_i;

my ($vref, $vname);
if (defined $veog_i) { $vref = $data->slice("($veog_i),:")->sever; $vname = $clean[$veog_i]; }
else {
    my @fp = grep { defined } @idx{qw(Fp1 Fp2)};
    die "no vEOG/EOG and no Fp1/Fp2 for a blink reference\n" unless @fp;
    $vref = $data->dice_axis(0, pdl(long,@fp))->average; $vname = 'mean(Fp1,Fp2)';
}
my ($href, $hname);
if (defined $heog_i) { $href = $data->slice("($heog_i),:")->sever; $hname = $clean[$heog_i]; }
printf "blink reference (vertical):   %s\n", $vname;
printf "saccade reference (horizontal): %s\n", (defined $href ? $hname : '(none - vEOG-only fallback)');

# ---- ICA --------------------------------------------------------------------
my %iopt = (nl=>$o{nl}, seed=>$o{seed}, maxit=>$o{maxit}, tol=>$o{tol});
$iopt{keep} = $o{keep} if defined $o{keep};
my $ica = ica_decompose($X, %iopt);
my $S   = $ica->{sources};                          # (nt, k)
printf "ICA: %d comps, iters=%d conv=%.2e (nl=%s seed=%d%s)\n",
    $ica->{k}, $ica->{iters}, $ica->{conv}, $o{nl}, $o{seed},
    (defined $o{keep} ? " keep=$o{keep}" : '');

# correlation of every component with each EOG
my $k = $ica->{k};
my @cv = map { corr($S->slice(":,(".$_.")"), $vref) } 0..$k-1;   # vs vEOG
my @ch = defined $href ? map { corr($S->slice(":,(".$_.")"), $href) } 0..$k-1 : ();

# ---- pick ocular components -------------------------------------------------
my @remove;    # ordered list of [role, ic, cv, ch]
sub argmax_absv { my ($arr,$skip)=@_; my ($bi,$bv)=(-1,-1);
    for my $i (0..$#$arr){ next if $skip->{$i}; my $a=abs $arr->[$i]; if($a>$bv){$bv=$a;$bi=$i} } $bi }

my %used_ic;
# #1 blink = max |corr vEOG|
my $blink = argmax_absv(\@cv, \%used_ic);
$used_ic{$blink} = 1;
push @remove, ['blink', $blink, $cv[$blink], (@ch?$ch[$blink]:undef)];

if (@ch) {
    # #2 horizontal = max |corr hEOG| among the rest
    my $hor = argmax_absv(\@ch, \%used_ic);
    $used_ic{$hor} = 1;
    push @remove, ['horizontal', $hor, $cv[$hor], $ch[$hor]];
    # #3 next ocular = max( |cv|, |ch| ) among the rest
    my @comb = map { my $a=abs $cv[$_]; my $b=abs $ch[$_]; $a>$b?$a:$b } 0..$k-1;
    my $nxt = argmax_absv(\@comb, \%used_ic);
    push @remove, ['next-ocular', $nxt, $cv[$nxt], $ch[$nxt]];
} else {
    # vEOG-only fallback: next two by |corr vEOG|
    for my $role ('vEOG#2','vEOG#3') {
        my $i = argmax_absv(\@cv, \%used_ic); $used_ic{$i}=1;
        push @remove, [$role, $i, $cv[$i], undef];
    }
}

print "selected components to remove:\n";
for my $r (@remove) {
    my ($role,$ic,$v,$h) = @$r;
    printf "  %-11s IC %2d  corr(vEOG)=%+.3f%s\n",
        $role, $ic, $v, (defined $h ? sprintf("  corr(hEOG)=%+.3f",$h) : '');
}

# ---- field width: wide enough for the largest value (DC channels) -----------
my $fieldw;
if (defined $o{fieldw}) { $fieldw = $o{fieldw} }
else {
    my $mx  = $data->abs->max->sclr;
    my $dig = ($mx >= 1) ? int(log($mx)/log(10)) + 1 : 1;
    my $need = $dig + $o{decimals} + 3;              # sign + '.' + one-space margin
    $fieldw = $need > 8 ? $need : 8;
}

# ---- remove #1 / #1-2 / #1-2-3 and write FULL-record .mul --------------------
(my $base = $in) =~ s/\.[Ee][Dd][Ff]$//;
$base = $o{prefix} if defined $o{prefix};
my @tags = ('1', '1-2', '1-2-3');
for my $n (1..3) {
    last if $n > @remove;
    my @bad = map { $remove[$_][1] } 0 .. $n-1;
    my $cleaned = apply_ica($X, $ica, \@bad);        # (nt, n_ica)
    my $cT      = $cleaned->transpose->sever;        # (n_ica, nt)

    my $outdata = $data->copy;                       # full record
    for my $j (0 .. $#ica_idx) {                     # overwrite only ICA columns
        $outdata->slice("(".$ica_idx[$j]."),:") .= $cT->slice("(".$j."),:");
    }
    my $out = { data => $outdata, fs => $fs, labels => [ @clean ],
                t_start => $rec->{t_start} };
    my $file = "$base-ica-$tags[$n-1].mul";
    write_mul($out, $file, decimals => $o{decimals}, fieldw => $fieldw);
    printf "wrote %s : removed IC(s) %s (%s)\n",
        $file, join(',', @bad), join('+', map { $remove[$_][0] } 0..$n-1);
}
