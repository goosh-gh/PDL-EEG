use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use Cwd qw(getcwd);

# The synthetic Nihon Kohden fixtures (t/data/test01.eeg, test01.LOG, test02.eeg, test03.eeg)
# are not in the repository: t/mk_synthetic_nk.pl writes them (deterministically, no random
# numbers). Make them here when any is missing, so that make test and prove need no extra step.
my @want = map { "$Bin/data/$_" } qw(test01.eeg test01.LOG test02.eeg test03.eeg);
my @missing = grep { !-f } @want;
if (@missing) {
    my $here = getcwd();
    chdir "$Bin/.." or die "chdir $Bin/..: $!";                  # the generator writes to t/data under the cwd
    my $rc = system($^X, 't/mk_synthetic_nk.pl');
    chdir $here;
    is($rc, 0, 't/mk_synthetic_nk.pl ran (' . scalar(@missing) . ' fixture(s) were missing)');
} else {
    pass('the synthetic fixtures are present');
}
ok(-f $_, "fixture $_") for @want;

done_testing;
