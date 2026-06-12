use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'PDL::EEG',
    AUTHOR           => 'goosh <goosh@example.com>',
    VERSION_FROM     => 'lib/PDL/EEG.pm',
    ABSTRACT         => 'EEG analysis toolkit for PDL -- MNE-Python-inspired',
    LICENSE          => 'perl_5',
    MIN_PERL_VERSION => '5.020',
    PREREQ_PM => {
        'PDL'    => '2.080',
        'Encode' => 0,
        'Carp'   => 0,
    },
    TEST_REQUIRES => {
        'Test::More' => '0.98',
    },
    test => { TESTS => 't/*.t' },
    META_MERGE => {
        'meta-spec' => { version => 2 },
        resources => {
            repository => {
                type => 'git',
                url  => 'https://github.com/goosh-gh/PDL-EEG.git',
                web  => 'https://github.com/goosh-gh/PDL-EEG',
            },
        },
    },
);

package MY;

sub postamble { return <<'END_MAKE'; }
testdata :
	$(PERLRUN) t/mk_synthetic_nk.pl

END_MAKE
