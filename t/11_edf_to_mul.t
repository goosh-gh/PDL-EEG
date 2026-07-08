use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::IO::EDF qw(write_edf clean_edf_label);
use File::Temp qw(tempdir);
use FindBin;
# ---------------------------------------------------------------------------
# Integration test for examples/edf_to_mul.pl.
#
# Writes a tiny EDF+ file whose signal labels use the vendor conventions that
# broke the .mul header in P:EEG11 -- an EDF+ type prefix ("EEG "/"POL "), a
# referential "-Ref" marker, a Nihon Kohden "$A1" reference channel, and a
# Trigger channel -- then runs the real CLI and checks the emitted .mul:
#   * EDF+ type prefix stripped        (EEG Fp1-Ref -> Fp1)
#   * referential marker dropped       (-Ref removed)
#   * $-ref channel renamed            ($A1 -> A1_ref)
#   * X1 IS montage-suffixed           (X1 -> X1-BN, matches vendor .mul)
#   * DC / $ref / Trigger un-suffixed
#   * label row is whitespace-free and has exactly n_ch tokens
#   * Trigger written as a column but NOT counted in Channels=
# ---------------------------------------------------------------------------

# --- unit: clean_edf_label (pure Perl, EDF+ label normalisation) ------------
is(clean_edf_label('EEG Fp1-Ref'), 'Fp1',     'type prefix + -Ref stripped');
is(clean_edf_label('EEG A1-Ref'),  'A1',      'scalp A1-Ref -> A1');
is(clean_edf_label('POL DC01'),    'DC01',    'POL prefix stripped');
is(clean_edf_label('POL $A1'),     'A1_ref',  '$A1 -> A1_ref (Perl-safe)');
is(clean_edf_label('POL $A2'),     'A2_ref',  '$A2 -> A2_ref');
is(clean_edf_label('POL X1'),      'X1',      'X1 kept');
is(clean_edf_label('POL E'),       'E',       'E kept');
is(clean_edf_label('Trigger'),     'Trigger', 'no prefix -> unchanged');
is(clean_edf_label('Fp1'),         'Fp1',     'already clean -> idempotent');
is(clean_edf_label('EEG Fp1-REF'), 'Fp1',     '-REF case-insensitive');
is(clean_edf_label('POL foo bar'), 'foo_bar', 'residual whitespace -> _');
is(clean_edf_label(undef),         undef,     'undef passes through');

my $script = "$FindBin::Bin/../examples/edf_to_mul.pl";
my $libdir = "$FindBin::Bin/../lib";

SKIP: {
    skip "edf_to_mul.pl not found", 14 unless -f $script;

    my $dir = tempdir(CLEANUP => 1);

# --- build a 5-channel EDF+ record with the tricky labels -------------------
my $fs  = 1000;
my $ns  = 20;
my $t   = sequence($ns) / $fs;
my $data = pdl(
    (100 * sin($t * 6))->list,      # EEG Fp1-Ref : scalp
    (10 * $t)->list,                # POL X1      : aux (suffixed per vendor)
    (5 * ones($ns))->list,          # POL DC01    : DC trigger input
    (3 * ones($ns))->list,          # POL $A1     : A1 reference value (const)
    (zeroes($ns))->list,            # Trigger     : integer trigger channel
)->reshape($ns, 5)->xchg(0, 1)->sever;

my $rec = {
    data    => $data->float,
    fs      => $fs,
    labels  => ['EEG Fp1-Ref', 'POL X1', 'POL DC01', 'POL $A1', 'Trigger'],
    t_start => '2026-07-05 14:03:19',
};

my $edf = "$dir/tricky.edf";
write_edf($rec, $edf, phys => 'gain');
ok(-s $edf, 'write_edf produced an EDF file');

# --- run the real CLI: whole recording, montage suffix -BN ------------------
my $mul = "$dir/out.mul";
my $rc  = system($^X, "-I$libdir", $script, $edf, '--suffix', '-BN', '--out', $mul);
is($rc, 0, 'edf_to_mul.pl exited 0');
ok(-s $mul, 'edf_to_mul.pl wrote a .mul file');

open my $fh, '<', $mul or die "open $mul: $!";
my @L = <$fh>;
close $fh;

like($L[0], qr/ Channels=4 /,
     'Channels=4 : 5 channels minus the uncounted Trigger');
like($L[0], qr/Time=14:03:19/, 'Time= taken from t_start');

is($L[1], " Fp1-BN X1-BN DC01 A1_ref Trigger\n",
   'label row cleaned: prefix/-Ref stripped, $A1->A1_ref, X1-BN, DC/Trigger bare');

(my $row = $L[1]) =~ s/^ //; chomp $row;
is(scalar(split ' ', $row), 5, 'label row keeps all 5 tokens (incl Trigger)');
unlike($L[1], qr/-Ref/,          'no -Ref remains');
unlike($L[1], qr/\$/,            'no $ remains (Perl-safe names)');
unlike($L[1], qr/(?:EEG|POL)\s/, 'no EDF+ type prefix remains');

# --- --cut in data-coordinate seconds --------------------------------------
my $cutbase = "$dir/cut.mul";
my $rc2 = system($^X, "-I$libdir", $script, $edf,
                 '--suffix', '-BN', '--out', $cutbase,
                 '--cut', '0.005-0.015:seg');
is($rc2, 0, 'edf_to_mul.pl --cut exited 0');

my $cutfile = "$dir/cut_seg.mul";
ok(-s $cutfile, '--cut wrote <base>_seg.mul');
open my $cf, '<', $cutfile or die "open $cutfile: $!";
my @C = <$cf>;
close $cf;
like($C[0], qr/^TimePoints=10 /,
     '--cut "0.005-0.015" -> 10 samples (end-exclusive)');
like($C[0], qr/ Channels=4 /, '--cut output keeps Channels=4');
}

done_testing();
