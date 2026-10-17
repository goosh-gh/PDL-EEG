use strict;
use warnings;
use Test::More;
use PDL;
use PDL::EEG::IO::WFDB qw(read_wfdb read_wfdb_annotations);

# Real-data check of PDL::EEG::IO::WFDB on PhysioNet Challenge 2018 records
# (not part of `make test`). Pass record stems (path without extension):
#
#   prove -lv xt/71_wfdb_real.t :: /data/tr03-0005 /data/te03-0024
#   PDL_EEG_WFDB_RECORDS="/data/tr03-0005" prove -lv xt/71_wfdb_real.t
#
# For each record: the header and the .mat agree on the number of signals and
# samples; if PDL::IO::Matlab is installed, its matlab_read of the .mat gives the
# same int16 matrix; if <stem>.arousal exists, its annotation times are
# non-decreasing and lie within the record. With no records it skips.

my @stems = @ARGV ? @ARGV
          : defined $ENV{PDL_EEG_WFDB_RECORDS} ? split(' ', $ENV{PDL_EEG_WFDB_RECORDS})
          : ();
plan skip_all => 'no records given (see the comment at the top)' unless @stems;
my $have_matlab = eval { require PDL::IO::Matlab; 1 };

for my $stem (@stems) {
    my $rec = read_wfdb($stem);
    my ($nsig, $nt) = $rec->{raw}->dims;
    is($nsig, scalar @{ $rec->{labels} }, "$stem: one row per signal");
    is($nt, $rec->{nt}, "$stem: nt");
  SKIP: {
        skip 'PDL::IO::Matlab not installed', 1 unless $have_matlab && -f "$stem.mat";
        my @m = PDL::IO::Matlab::matlab_read("$stem.mat");
        my $M = $m[0];
        $M = $M->xchg(0, 1) if $M->dim(0) != $nsig;
        ok(all($M->short == $rec->{raw}), "$stem: same matrix as PDL::IO::Matlab");
    }
  SKIP: {
        skip 'no .arousal', 2 unless -f "$stem.arousal";
        my $ann = read_wfdb_annotations("$stem.arousal");
        my @t = map { $_->[0] } @$ann;
        ok(!grep({ $t[$_] < $t[$_ - 1] } 1 .. $#t), "$stem: annotation times non-decreasing");
        ok($t[0] >= 0 && $t[-1] < $nt, "$stem: annotations within the record");
    }
}

done_testing();
