#!/usr/bin/env perl
# examples/heartbeat_resp_quality_pool.pl -- pools the PREFIX_quality.tsv files of heartbeat_resp_quality_vs_belt.pl (several records): how well each ECG-only quality index
# predicts the true error of the ECG respiratory phase, over all windows and within records, and a calibration table for the best one
# (q_agree: disagreement between the H and L phases) -- the true error to expect in a window, given the index, without a belt.
# Windows with more than --maxev (0.3) of their beats inside annotated events are left out (when the stage column is filled).
# usage: perl heartbeat_resp_quality_pool.pl *_quality.tsv [--maxev 0.3]
use strict; use warnings; use Getopt::Long; use List::Util qw(sum);
my %o = (maxev => 0.3); GetOptions(\%o, 'maxev=f') or die;
my (@R, %byrec);
for my $f (@ARGV) { open my $h, '<', $f or die "$f: $!"; my $hd = <$h>; chomp $hd; my @c = split /\t/, $hd;
    while (<$h>) { chomp; my %r; @r{@c} = split /\t/; next if length($r{ev_share} // '') && $r{ev_share} > $o{maxev}; push @R, \%r; push @{ $byrec{ $r{rec} } }, \%r } }
die "no windows\n" unless @R;
sub ranks { my @v = @_; my @i = sort { $v[$a] <=> $v[$b] } 0 .. $#v; my @r; $r[ $i[$_] ] = $_ for 0 .. $#i; @r }
sub spear { my ($x, $y) = @_; my @rx = ranks(@$x); my @ry = ranks(@$y); my ($mx, $my) = (sum(@rx) / @rx, sum(@ry) / @ry);
    my ($s, $sx, $sy) = (0, 0, 0); for (0 .. $#rx) { $s += ($rx[$_] - $mx) * ($ry[$_] - $my); $sx += ($rx[$_] - $mx) ** 2; $sy += ($ry[$_] - $my) ** 2 } $s / sqrt($sx * $sy) }
sub median { my @s = sort { $a <=> $b } @_; @s ? ($s[ $#s / 2 ] + $s[ @s / 2 ]) / 2 : 'nan' }
printf "%d windows from %d records (event share <= %g)\n", scalar @R, scalar keys %byrec, $o{maxev};
printf "%-8s %10s %18s\n", 'index', 'all', 'within records';
for my $q (qw(q_fold q_snr q_band q_agree q_rr)) { my $all = spear([ map { $_->{$q} } @R ], [ map { $_->{true_sd} } @R ]);
    my @w = map { my $rr = $byrec{$_}; @$rr >= 10 ? spear([ map { $_->{$q} } @$rr ], [ map { $_->{true_sd} } @$rr ]) : () } sort keys %byrec;
    printf "%-8s %+10.2f %+10.2f (median of %d)\n", $q, $all, median(@w), scalar @w }
print "\ntrue error [deg] by q_agree (circular SD between the H and L phases in the window):\n";
my @cut = (0, 15, 20, 25, 30, 40, 60, 1e9);
for my $i (0 .. $#cut - 1) { my @s = grep { $_->{q_agree} >= $cut[$i] && $_->{q_agree} < $cut[ $i + 1 ] } @R; next unless @s;
    my @t = map { $_->{true_sd} } @s; my $good = grep { $_ < 45 } @t;
    printf "  q_agree %3s-%-4s %5d windows (%4.1f %%): true error median %5.1f, share < 45 deg %3.0f %%\n", $cut[$i], $cut[ $i + 1 ] > 1e8 ? '' : $cut[ $i + 1 ], scalar @s, 100 * @s / @R, median(@t), 100 * $good / @t }
if (grep { ($_->{stage} // '-') ne '-' } @R) { print "\nby stage: median true error, median q_agree, windows\n";
    for my $st (qw(W N1 N2 N3 R)) { my @s = grep { $_->{stage} eq $st } @R; next unless @s; printf "  %-3s %5.1f %5.1f %6d\n", $st, median(map { $_->{true_sd} } @s), median(map { $_->{q_agree} } @s), scalar @s } }
