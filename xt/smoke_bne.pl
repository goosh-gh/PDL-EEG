#!/usr/bin/env perl
# Real-data smoke test for the --bne re-reference path.
#
#   perl -Ilib xt/smoke_bne.pl /path/to/JJ0090J6.EEG      # nk_to_mul
#   perl -Ilib xt/smoke_bne.pl /path/to/JJ0090J6.edf      # edf_to_mul
#
# Runs the appropriate converter WITH and WITHOUT --bne on a real recording and
# checks structural invariants (no ground-truth values needed):
#   * plain run keeps BN1/BN2; --bne drops them and tags scalp channels -BN
#   * every .mul is self-consistent: label-row token count == Channels + (#Trigger)
#   * data are finite (no NaN/Inf) and re-referenced values differ from raw
# This is an author test (xt/), gated on a real file being passed.

use strict;
use warnings;
use FindBin;
use File::Temp qw(tempdir);
use Test::More;

my $in = shift @ARGV or plan skip_all => "usage: smoke_bne.pl FILE.EEG|FILE.edf";
plan skip_all => "no such file: $in" unless -f $in;

my $is_edf = $in =~ /\.edf$/i;
my $script = "$FindBin::Bin/../examples/" . ($is_edf ? 'edf_to_mul.pl' : 'nk_to_mul.pl');
my $lib    = "$FindBin::Bin/../lib";
plan skip_all => "converter not found: $script" unless -f $script;

my $dir = tempdir(CLEANUP => 1);

sub run_conv {
    my ($out, @opt) = @_;
    my $rc = system($^X, "-I$lib", $script, $in, '--out', $out, @opt);
    return $rc == 0;
}
sub load_mul {
    my $f = shift;
    open my $fh, '<', $f or die "open $f: $!";
    my @L = <$fh>; close $fh;
    my ($ch) = $L[0] =~ /Channels=(\d+)/;
    (my $lab = $L[1]) =~ s/^\s+//; chomp $lab;
    my @tok = split ' ', $lab;
    my @vals = split ' ', ($L[2] // '');
    return { channels => $ch, labels => \@tok, first_row => \@vals, n_lines => scalar @L };
}
sub finite { my $v = shift; $v == $v && $v !~ /inf/i }   # NaN/Inf guard

# ---- plain (no re-reference) ------------------------------------------------
my $plain = "$dir/plain.mul";
ok(run_conv($plain),               'plain conversion exits 0');
my $P = load_mul($plain);
ok(@{$P->{labels}},                'plain .mul has a label row');
my $has_trigger = grep { /^Trigger$/i } @{$P->{labels}};
is(scalar @{$P->{labels}}, $P->{channels} + $has_trigger,
   "plain: label tokens == Channels + trigger ($P->{channels}+$has_trigger)");
ok((grep { /^BN1/ } @{$P->{labels}}), 'plain: BN1 present (not re-referenced)');
ok((!grep { !finite($_) } @{$P->{first_row}}), 'plain: first row all finite');

# ---- --bne (default prop 0.5) ----------------------------------------------
my $bne = "$dir/bne.mul";
ok(run_conv($bne, '--bne'),        '--bne conversion exits 0');
my $B = load_mul($bne);
my $bt = grep { /^Trigger$/i } @{$B->{labels}};
is(scalar @{$B->{labels}}, $B->{channels} + $bt,
   "--bne: label tokens == Channels + trigger ($B->{channels}+$bt)");
ok((!grep { /^BN1(-BN)?$/ } @{$B->{labels}}), '--bne: BN1 dropped');
ok((!grep { /^BN2(-BN)?$/ } @{$B->{labels}}), '--bne: BN2 dropped');
ok((grep { /-BN$/ } @{$B->{labels}}),          '--bne: scalp channels tagged -BN');
ok($B->{channels} < $P->{channels},            '--bne: fewer channels than plain (BN1/BN2 gone)');
ok((!grep { !finite($_) } @{$B->{first_row}}), '--bne: first row all finite');

# re-referenced values should differ from the raw ones (sanity: something changed)
my $changed = 0;
for my $i (0 .. $#{$B->{first_row}}) {
    $changed++ if defined $P->{first_row}[$i]
               && abs(($B->{first_row}[$i] // 0) - ($P->{first_row}[$i] // 0)) > 1e-6;
}
ok($changed > 0, '--bne: re-referenced values differ from raw');

diag("plain:  Channels=$P->{channels}  labels[0..3]=@{$P->{labels}}[0..3]");
diag("--bne:  Channels=$B->{channels}  labels[0..3]=@{$B->{labels}}[0..3]");

done_testing();
