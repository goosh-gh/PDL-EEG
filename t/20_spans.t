use strict; use warnings; use Test::More;
use PDL; use FindBin; use lib "$FindBin::Bin/../lib";
use PDL::EEG::Spans qw(read_spans write_spans spans_to_mask kinds_present);

my @spans = (
  { kind=>'blink',    t0=>1.000, t1=>1.300 },
  { kind=>'hsaccade', t0=>2.500, t1=>2.900 },
  { kind=>'hsaccade', t0=>0.100, t1=>0.200 },
);
my $tmp = "$FindBin::Bin/_spans_tmp.tsv";
is(write_spans($tmp,\@spans), 3, "wrote 3 spans");
my $back = read_spans($tmp);
is(scalar(@$back), 3, "read 3 back");
is($back->[0]{kind}, 'hsaccade', "sorted by t0");
is_deeply([kinds_present($back)], ['blink','hsaccade'], "kinds_present");
my $m_all = spans_to_mask($back, fs=>1000, nt=>4000);
is($m_all->sum, 301+401+101, "any-kind mask count (inclusive)");
my $m_h = spans_to_mask($back, fs=>1000, nt=>4000, kind=>'hsaccade');
is($m_h->sum, 401+101, "hsaccade-only count");
is($m_h->at(2500),1,"start set"); is($m_h->at(2900),1,"end set");
is($m_h->at(2499),0,"before unset"); is($m_h->at(2901),0,"after unset");
my $m_off = spans_to_mask($back, fs=>1000, nt=>4000, abs0=>2000, kind=>'hsaccade');
is($m_off->at(500),1,"abs0 offset maps 2.5s->local 500");
my $m_clip = spans_to_mask([{kind=>'x',t0=>-5,t1=>1.0}], fs=>1000, nt=>4000);
is($m_clip->at(0),1,"neg start clipped"); is($m_clip->sum,1001,"clipped count");
unlink $tmp;
done_testing();
