use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::MAP2D qw(project_positions interpolate_topo);
use PDL::EEG::IO::ASA ();

# Numeric core only (projection + thin-plate spline); no rendering, so
# PDL::Graphics::Cairo is not required for this test.

my $elc = 't/data/standard_1020.elc';
plan skip_all => "$elc not found" unless -r $elc;

my $mon = PDL::EEG::IO::ASA::read_elc($elc);

# 10-20 scalp set (old names; T3/T4/T5/T6 alias to T7/T8/P7/P8 in the montage)
my @chan = qw(Fp1 Fp2 F7 F3 Fz F4 F8 T3 C3 Cz C4 T4 T5 P3 Pz P4 T6 O1 O2);
my %ALIAS = (T3 => 'T7', T4 => 'T8', T5 => 'P7', T6 => 'P8');
my (@x, @y, @z, @keep);
for my $c (@chan) {
    my $k = $mon->{pos}{$c} ? $c : $ALIAS{$c};
    next unless $k && $mon->{pos}{$k};
    my $p = $mon->{pos}{$k};
    push @x, $p->{x}; push @y, $p->{y}; push @z, $p->{z}; push @keep, $c;
}
ok(@keep >= 19, 'matched the 19 standard channels (alias-aware)');
my $coords = cat(pdl(@x), pdl(@y), pdl(@z))->xchg(0, 1);   # (3,N)

my ($px, $py, $ref) = project_positions($coords, margin => 0.2);
my %i; $i{$keep[$_]} = $_ for 0 .. $#keep;

# orientation: nose = +y up, right ear = +x
ok($py->at($i{Fz}) > $py->at($i{Cz}), 'Fz is anterior to Cz (up)');
ok($py->at($i{Cz}) > $py->at($i{O1}), 'Cz is anterior to O1 (down)');
ok($px->at($i{T4}) > 0 && $px->at($i{T3}) < 0, 'T4 right, T3 left');
ok(abs($px->at($i{Cz})) < 0.05, 'Cz near the vertical midline');

# recentring keeps the array roughly centred
ok(abs($py->avg) < 0.06, 'array centroid close to origin after recentering');

# thin-plate spline returns a finite, masked grid
my $v = zeroes(scalar @keep);
$v->set($i{Cz}, 5); $v->set($i{O1}, -3);
my ($field, $inside, $ext, $rclip) = interpolate_topo($px, $py, $v, res => 60, clip => 0.54);
is_deeply([$field->dims], [60, 60], 'field grid is res x res');
ok($field->isfinite->all, 'field is finite everywhere');
ok($inside->sum > 0 && $inside->sum < $field->nelem, 'in-disc mask is a subset');
ok($rclip > 0.5, 'clip radius extends past the head circle');

done_testing;
