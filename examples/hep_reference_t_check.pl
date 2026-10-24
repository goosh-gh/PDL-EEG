#!/usr/bin/env perl
#
# hep_reference_t_check.pl — run heartbeat_hep_resp_grid.pl on one recording in the five
# conditions used to check a heartbeat-evoked potential against the two weak points of
# the cardiac-field removal, and print one summary:
#
#   1. reference: does the residual A2-A1 (a bipolar difference, the sign of cardiac field
#      left in the ears) stay the same when the head reference changes from the ears to the
#      average of the head channels?  If the cardiac leads take out every cardiac
#      component, it does; if it changes, the reference matters for this recording.
#   2. T wave left in 200..300 ms: the fit in the QRS only (-45,55), in P+QRS (-250,55,
#      nothing after 55 ms used), in the whole beat (-250,600, takes out too much: also
#      brain activity shaped like the leads), and in P+QRS with the ear difference A1-A2 as
#      one more (left-right) cardiac lead.
#
#   perl hep_reference_t_check.pl REC.EEG|REC.edf|REC.vhdr [options] [-- options for the grid script]
#
#     --grid PATH      heartbeat_hep_resp_grid.pl     (the one next to this script)
#     --out PREFIX     prefix of the runs' files and logs           (REC without extension, + _check)
#     --only LIST      run only these conditions (numbers 1..5, e.g. 2,3)          (all)
#     --show LIST      channels in the summary            (Fz,Cz,Pz,C3,C4,T3,T4,T5,T6,F7,F8,O1,O2)
#     after --        passed to every run, e.g.  -- --cardiac ClavL,ClavR   or  -- --exclude F3,F4
#
#   conditions (all with --leads reref --bn-channels):
#     1  qrs       --ref ears --fit=-45,55
#     2  pqrs      --ref ears --fit=-250,55
#     3  pqrs_avg  --ref avg  --fit=-250,55
#     4  wide      --ref ears --fit=-250,600
#     5  pqrs_ear  --ref ears --fit=-250,55 --cardiac A1-A2   (adds A1-A2 to any --cardiac after --)
#
# Each run's full output is kept in PREFIX_<name>.log (and its TSV/PNG files under
# PREFIX_<name>).  The summary takes the residual means without respiration from the
# "residual mean" table of each log; the logs of earlier runs are reused when they exist
# and --rerun is not given, so a long recording need not be run twice.
#
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long qw(:config pass_through);

my ($grid, $out, $only, $show, $rerun) = ("$Bin/heartbeat_hep_resp_grid.pl", undef, '', 'Fz,Cz,Pz,C3,C4,T3,T4,T5,T6,F7,F8,O1,O2', 0);
my @pass;
if (my ($i) = grep { $ARGV[$_] eq '--' } 0 .. $#ARGV) { @pass = @ARGV[$i + 1 .. $#ARGV]; splice @ARGV, $i }
GetOptions('grid=s' => \$grid, 'out=s' => \$out, 'only=s' => \$only, 'show=s' => \$show, 'rerun' => \$rerun)
    or die "bad options\n";
my $file = shift @ARGV or die "usage: $0 REC.EEG|REC.edf|REC.vhdr [options] [-- grid options]\n";
!@ARGV or die "unexpected argument(s): @ARGV (options for the grid script go after --)\n";
-f $file or die "$file: not found\n";
-f $grid or die "$grid: not found (give it with --grid)\n";
($out = $file) =~ s/\.(edf|vhdr|eeg)$//i, $out .= '_check' unless defined $out;

# --cardiac given after -- is merged with A1-A2 for condition 5
my ($pcard) = map { $pass[$_ + 1] } grep { $pass[$_] eq '--cardiac' } 0 .. $#pass - 1;
($pcard) = map { /^--cardiac=(.*)$/ ? $1 : () } @pass unless defined $pcard;
my @pass5 = do { my @k; for (my $j = 0; $j < @pass; $j++) { if ($pass[$j] eq '--cardiac') { $j++; next } next if $pass[$j] =~ /^--cardiac=/; push @k, $pass[$j] } @k };

my @cond = (
    [ 1, 'qrs',      'ears, QRS (-45,55)',           [ '--ref', 'ears', '--fit=-45,55' ],   \@pass ],
    [ 2, 'pqrs',     'ears, P+QRS (-250,55)',        [ '--ref', 'ears', '--fit=-250,55' ],  \@pass ],
    [ 3, 'pqrs_avg', 'average, P+QRS (-250,55)',     [ '--ref', 'avg',  '--fit=-250,55' ],  \@pass ],
    [ 4, 'wide',     'ears, whole beat (-250,600)',  [ '--ref', 'ears', '--fit=-250,600' ], \@pass ],
    [ 5, 'pqrs_ear', 'ears, P+QRS, + A1-A2',         [ '--ref', 'ears', '--fit=-250,55',
                                                       '--cardiac', join(',', grep { defined && length } $pcard, 'A1-A2') ], \@pass5 ],
);
my %want = map { $_ => 1 } grep { length } split /,/, $only;
@cond = grep { !%want || $want{ $_->[0] } } @cond;

my %R;                                             # name -> { ch -> [v1, v2, v3] }, info
for my $c (@cond) {
    my ($n, $name, $desc, $args, $extra) = @$c;
    my $log = "${out}_$name.log";
    if ($rerun || !-s $log || !_complete($log)) {
        my @cmd = ($^X, $grid, $file, '--leads', 'reref', '--bn-channels', @$args, @$extra, '--out', "${out}_$name");
        print "[$n] $desc\n    @cmd\n";
        open my $lh, '>', $log or die "$log: $!";
        open my $ph, '-|', @cmd or die "cannot run $grid: $!";
        while (<$ph>) { print $lh $_ }
        close $ph;
        close $lh;
        $? == 0 or die "[$n] $name failed (exit " . ($? >> 8) . "); see $log\n";
    } else {
        print "[$n] $desc: reusing $log\n";
    }
    $R{$name} = _parse($log);
}

# ---- summary ----------------------------------------------------------------------------
my @w = ('200-300', '300-450', '450-600');
my $f3 = sub { my $v = shift; defined $v ? join(' / ', map { sprintf '%+5.2f', $_ } @$v) : '-' };
print "\n=== $file\n";
for my $c (@cond) {
    my $r = $R{ $c->[1] };
    printf "  [%d] %-30s %s beats, fit share %s\n", $c->[0], $c->[2], $r->{beats} // '?', $r->{share} // '?';
}

print "\nA2-A1 after the fit [uV] ($w[0] / $w[1] / $w[2] ms): left-right cardiac field left in the ears\n";
for my $c (@cond) {
    my $r = $R{ $c->[1] }{ch};
    my $d = ($r->{A1} && $r->{A2}) ? [ map { $r->{A2}[$_] - $r->{A1}[$_] } 0 .. 2 ] : undef;
    printf "  [%d] %-30s %s%s\n", $c->[0], $c->[2], $f3->($d),
        $c->[1] eq 'pqrs_ear' ? '   (0 by construction: A1-A2 is a regressor)' : '';
}
if ($R{pqrs} && $R{pqrs_avg} && $R{pqrs}{ch}{A1} && $R{pqrs_avg}{ch}{A1}) {
    my @a = map { $R{pqrs}{ch}{A2}[$_] - $R{pqrs}{ch}{A1}[$_] } 0 .. 2;
    my @b = map { $R{pqrs_avg}{ch}{A2}[$_] - $R{pqrs_avg}{ch}{A1}[$_] } 0 .. 2;
    my $dmax = 0; for (0 .. 2) { my $x = abs($a[$_] - $b[$_]); $dmax = $x if $x > $dmax }
    printf "  reference check (ears vs average, P+QRS): largest difference %.2f uV -> %s\n", $dmax,
        $dmax <= 0.15 ? 'the same: the cardiac leads cover both references here'
                      : 'DIFFERENT: the reference matters in this recording';
}

my @sh = grep { length } split /,/, $show;
print "\nresidual mean [uV] without respiration ($w[0] / $w[1] / $w[2] ms)\n";
printf "  %-4s %s\n", 'ch', join('', map { sprintf '  %-22s', "[$_->[0]] $_->[1]" } @cond);
for my $ch (@sh) {
    next unless grep { $R{ $_->[1] }{ch}{$ch} } @cond;
    printf "  %-4s %s\n", $ch, join('', map { sprintf '  %-22s', $f3->($R{ $_->[1] }{ch}{$ch}) } @cond);
}
print "\nreading: [2] is the estimate that uses nothing after 55 ms; [4] takes out too much (lower bound);\n",
      "[5] takes out a left-right cardiac field through the ears (it may leave some of the field\n",
      "in the midline: an upper bound there); [3] shows the distribution against the average.\n";

# ---- helpers ----------------------------------------------------------------------------
sub _complete { my $f = shift; open my $h, '<', $f or return 0; local $/; my $t = <$h>; return $t =~ /^residual mean/m && $t =~ /expected by chance/ }
sub _parse {
    my $f = shift;
    open my $h, '<', $f or die "$f: $!";
    my (%ch, $beats, $share, $in);
    while (my $l = <$h>) {
        $beats = $1 if $l =~ /, (\d+) beats \(reject/;
        $share = $1 if $l =~ /explained inside the fit window: without rho median ([\d.]+)/;
        if ($l =~ /^residual mean/) { $in = 1; next }
        if ($in) {
            last if $l =~ /^\s*$/ && %ch;
            # "  Cz     1.30( 3.6)   1.45( 4.0)   ..." : without / with rho for the three windows
            if ($l =~ /^\s+(\S+)\s+(-?[\d.]+)\(/) {
                my $name = $1;
                my @v = $l =~ /(-?\d+\.\d+)\(\s*(?:-?[\d.]+|-)\s*\)/g;
                $ch{$name} = [ @v[0, 2, 4] ] if @v >= 6;
            }
        }
    }
    return { ch => \%ch, beats => $beats, share => $share };
}
