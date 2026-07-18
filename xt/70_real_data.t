use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::IO::NihonKohden qw(read_nk);

# Real-data regression for .LOG event placement (extblock AND wfmblock).
#
# These recordings are private (clinical data at machine-specific paths), so
# this test is NOT part of `make test`. Run it explicitly, passing the .EEG
# paths after `::` ...
#
#   prove -lv xt/70_real_data.t :: /path/A.EEG /path/B.EEG
#
# ... or via an environment variable (whitespace-separated):
#
#   PDL_EEG_TEST_FILES="/path/A.EEG /path/B.EEG" prove -lv xt/70_real_data.t
#
# With no files it skips, so it is safe to leave in xt/ and in the MANIFEST.
#
# For each file it reads the whole session (all_blocks) and checks that every
# .LOG event received a data-sample position, that every REC START lands exactly
# on a block boundary, and that every event falls inside a real segment in
# non-decreasing order. This is the property that _attach_recstart_samp
# guarantees; a regression (e.g. reverting to wall-clock placement) makes late
# segments drift off their boundaries and fail here.

my @files = @ARGV ? @ARGV
          : defined $ENV{PDL_EEG_TEST_FILES} ? split(' ', $ENV{PDL_EEG_TEST_FILES})
          : ();
my @have = grep { -f $_ } @files;

plan skip_all =>
    'no real .EEG files given (pass paths after :: or set PDL_EEG_TEST_FILES)'
    unless @have;

for my $f (@have) {
    my $r = eval { read_nk($f, all_blocks => 1) };
    ok($r, "read_nk: $f") or do { diag "read failed: $@"; next };

    my $meta = $r->{block_meta} || [];
    my $ev   = $r->{events}     || [];
    my $ns   = eval { $r->{data}->dim(1) } // 0;
    diag sprintf('%s  device=%s  layout=%s  blocks=%d  events=%d  n_samp=%d',
        $f, ($r->{device} // '?'), ($r->{layout} // '?'),
        scalar(@$meta), scalar(@$ev), $ns);

    my $multi = @$meta > 1;
    my @rec   = grep { ($_->{label} // '') =~ /REC\s*START/i } @$ev;

  SKIP: {
        skip "single-block read: events not placed by design ($f)", 3 unless $multi;
        skip "no .LOG events ($f)", 3 unless @$ev;

        my $placed = grep { defined $_->{samp} } @$ev;
        is($placed, scalar(@$ev), "all events placed ($f)");

        my %boundary = map { $_->{start_samp} => 1 } @$meta;
        my $rec_ok = grep { defined $_->{samp} && $boundary{ $_->{samp} } } @rec;
        is($rec_ok, scalar(@rec), "every REC START on a block boundary ($f)");

        my ($inside, $mono, $prev) = (1, 1, -1);
        for my $e (@$ev) {
            next unless defined $e->{samp};
            my $in = 0;
            for my $m (@$meta) {
                $in = 1 if $e->{samp} >= $m->{start_samp}
                        && $e->{samp} <  $m->{start_samp} + $m->{n_samp};
            }
            $inside = 0 unless $in;
            $mono   = 0 if $e->{samp} < $prev;
            $prev   = $e->{samp};
        }
        ok($inside && $mono, "events inside segments, non-decreasing ($f)");
    }
}

done_testing();
