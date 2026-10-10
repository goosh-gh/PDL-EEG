use strict;
use warnings;
use Test::More;
use File::Temp qw(tempfile);
use PDL;
use FindBin;

my $elc = "$FindBin::Bin/data/standard_1020_subset.elc";
plan skip_all => "$elc not found" unless -r $elc;
plan skip_all => 'PDL::Graphics::Cairo with Figure::add_connection / add_text not available'
    unless eval { require PDL::Graphics::Cairo; PDL::Graphics::Cairo::Figure->can('add_connection')
                      && PDL::Graphics::Cairo::Figure->can('add_text') };

use_ok('PDL::EEG::MAP2D', 'plot_joint');

# deterministic average: 19 scalp channels + one channel without a position
my @lab = qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 C3 Cz C4 T4 T5 P3 Pz P4 T6 O1 O2 X1);
my $t   = sequence(221) * 5 - 300;                           # -300 .. 800 ms
my $w   = pdl(map { sin(0.7 * $_) } 0 .. $#lab);
my $avg = $w->dummy(1) * exp(-(($t - 350) / 80) ** 2)->dummy(0) * 2;   # (20, 221)

my ($fh, $png) = tempfile(SUFFIX => '.png', UNLINK => 1);
close $fh;
my $ret = eval { plot_joint(avg => $avg, t => $t, labels => \@lab, montage => $elc,
                            times => [200, 350, 500], nave => 100, outfile => $png) };
ok(!$@, "plot_joint runs ($@)");
is($ret, $png, 'returns the output file name');
ok(-s $png, 'PNG written');

my $fig = eval { plot_joint(avg => $avg, t => $t, labels => \@lab, montage => $elc,
                            times => [200, 350, 500], negative_up => 0, ylim => [-2, 2], clim => 1.5) };
ok(!$@ && ref $fig, "without outfile the figure is returned ($@)");
is(scalar @{ $fig->{_connections} }, 3, 'one connection per requested latency');
# 3 maps + colour bar + waveform panel + electrode key
is(scalar @{ $fig->axes_list }, 6, 'panels: maps, colour bar, waveforms, key');
is($fig->{_connections}[1]{xyB}[0], 350, 'connection ends at the requested latency');

eval { plot_joint(avg => $avg, t => $t, labels => [ @lab[0 .. 5] ], montage => $elc, times => [200]) };
like($@, qr/size mismatch/, 'label count checked');
eval { plot_joint(avg => $avg, t => $t, labels => \@lab, montage => $elc) };
like($@, qr/'times'/, 'times required');

done_testing;
