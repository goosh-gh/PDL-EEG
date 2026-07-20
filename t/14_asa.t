use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use PDL;

use_ok('PDL::EEG::IO::ASA', qw(read_elc parse_ELEC_POS3D_ASA_4AdventCalendar))
    or BAIL_OUT("module will not load");

my $elc = "$Bin/data/standard_1020_subset.elc";
ok(-s $elc, "fixture present: $elc");

my $mon = read_elc($elc);

# ---- shape / metadata -----------------------------------------------------
is($mon->{n}, 28,          'read 28 positions');
is($mon->{unit}, 'mm',     'UnitPosition parsed');
is($mon->{reference}, 'avg','ReferenceLabel parsed');
is($mon->{number_positions}, 28, 'declared NumberPositions matches');
is_deeply([$mon->{coords}->dims], [3, 28], 'coords is (3,28)');
is(scalar @{$mon->{labels}}, 28, 'labels count matches');

# ---- label <-> coordinate association -------------------------------------
is($mon->{labels}[0], 'LPA', 'first label is LPA (fiducial lead)');
is($mon->{labels}[13], 'Cz', 'Cz at expected index');
ok(exists $mon->{pos}{Cz}, 'Cz in name lookup');
is($mon->{pos}{Cz}{index}, 13, 'Cz index consistent');

# ---- anatomy sanity (MNI mm, +x right / +y front / +z up) -----------------
my $cz = $mon->{pos}{Cz};
ok(abs($cz->{x}) < 5,   'Cz near midline (|x|<5mm)');
ok($cz->{z} > 90,       'Cz high on vertex (z>90mm)');
ok($mon->{pos}{T7}{x} < 0, 'T7 on the left  (x<0)');
ok($mon->{pos}{T8}{x} > 0, 'T8 on the right (x>0)');
ok($mon->{pos}{Fpz}{y} > 0, 'Fpz anterior (y>0)');
ok($mon->{pos}{Oz}{y}  < 0, 'Oz posterior (y<0)');

# ---- requested electrode set actually present -----------------------------
my @want = qw(Fp1 Fp2 F7 F3 Fz F4 F8 T7 C3 Cz C4 T8 P7 P3 Pz P4 P8 O1 O2
              Fpz Oz A1 A2 M1 M2);
ok(exists $mon->{pos}{$_}, "electrode $_ present") for @want;

# ---- fiducials ------------------------------------------------------------
is_deeply([sort keys %{$mon->{fiducials}}], [qw(LPA Nz RPA)],
          'three fiducials detected');
ok($mon->{fiducials}{LPA}[0] < 0, 'LPA left of midline');
ok($mon->{fiducials}{RPA}[0] > 0, 'RPA right of midline');
ok($mon->{fiducials}{Nz}[1]  > 0, 'Nz anterior');

# ---- coords row order matches pos lookup ----------------------------------
my $i = $mon->{pos}{Pz}{index};
ok(abs($mon->{coords}->at(0,$i) - $mon->{pos}{Pz}{x}) < 1e-4,
   'coords column agrees with pos lookup');

# ---- advent-compatible shim -----------------------------------------------
my ($h, $epos, $labels, $coords) = parse_ELEC_POS3D_ASA_4AdventCalendar($elc);
is($h->{FileComment}, '# ASA electrode file', 'shim: FileComment');
is($h->{UnitPosition}, 'mm', 'shim: UnitPosition');
is($h->{Cz}{DeviceCh}, 13, 'shim: per-electrode DeviceCh');
ok(abs($h->{Cz}{z} - $cz->{z}) < 1e-4, 'shim: coords agree with read_elc');
is($labels->[13], '  Cz', 'shim: 2-space label padding preserved');
is_deeply([$coords->dims], [3, 28], 'shim: coords (3,28)');
ok(all($coords == $mon->{coords}), 'shim: same coordinate piddle values');

done_testing();
