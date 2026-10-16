use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use PDL;
use PDL::EEG::MAP2D qw(project_positions interpolate_topo);
use PDL::EEG::IO::ASA ();

# Numeric core only (projection + thin-plate spline); no rendering, so
# PDL::Graphics::Cairo is not required for this test.

my $elc = "$Bin/data/standard_1020_subset.elc";
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

# ---- interpolate_topo: cached K against the direct path -------------------------
{
    my $rel = sub {     # max |f1 - f2| / max |f1| over every grid point
        my ($f1, $f2) = @_;
        return (($f1 - $f2)->abs->max / $f1->abs->max)->sclr;
    };
    CORE::srand(29);
    my @layouts = (['montage 19', $px, $py]);
    for my $n (5, 12, 32) {
        push @layouts, ["random $n",
            pdl(map { 0.8 * CORE::rand() - 0.4 } 1 .. $n),
            pdl(map { 0.8 * CORE::rand() - 0.4 } 1 .. $n)];
    }
    PDL::EEG::MAP2D::_topo_cache_clear();
    for my $L (@layouts) {
        my ($name, $lx, $ly) = @$L;
        my $n = $lx->nelem;
        for my $o ([res => 23, clip => 0.5], [res => 60, clip => 0.54],
                   [res => 41, clip => 0.54, anchors => 0], [res => 30, clip => 0.5, extent => 0.7]) {
            my %o = @$o;
            my $tag = "$name, " . join(' ', map { "$_=$o{$_}" } sort keys %o);
            my $vv = pdl(map { 10 * CORE::rand() - 5 } 1 .. $n);
            my ($f1, $i1, $e1, $r1) = interpolate_topo($lx, $ly, $vv, %o, cache => 0);
            my ($f2, $i2, $e2, $r2) = interpolate_topo($lx, $ly, $vv, %o);
            is_deeply([$f2->dims], [$f1->dims], "cached: same dims ($tag)");
            ok($rel->($f1, $f2) < 1e-10, "cached = direct at every grid point ($tag)");
            ok(($i1 == $i2)->all, "cached: same in-disc mask ($tag)");
            is_deeply([$e2, $r2], [$e1, $r1], "cached: same extent and clip ($tag)");
        }
    }

    # same layout, new values: K is not rebuilt
    PDL::EEG::MAP2D::_topo_cache_clear();
    interpolate_topo($px, $py, $v, res => 40, clip => 0.54);
    my $v2 = $v->copy; $v2->set($i{Fz}, -4); $v2->set($i{O2}, 2);
    my ($g2) = interpolate_topo($px, $py, $v2, res => 40, clip => 0.54);
    my ($d2) = interpolate_topo($px, $py, $v2, res => 40, clip => 0.54, cache => 0);
    is(PDL::EEG::MAP2D::_topo_cache_stats()->{builds}, 1, 'new values on the same layout: K built once');
    ok($rel->($d2, $g2) < 1e-10, 'new values on the same layout: still equal to the direct path');

    # another layout, or another resolution: K is built again
    my $qx = $px->copy; $qx->set($i{Cz}, $qx->at($i{Cz}) + 0.05);
    my ($g3) = interpolate_topo($qx, $py, $v2, res => 40, clip => 0.54);
    my ($d3) = interpolate_topo($qx, $py, $v2, res => 40, clip => 0.54, cache => 0);
    is(PDL::EEG::MAP2D::_topo_cache_stats()->{builds}, 2, 'moved sensor: K built again');
    ok($rel->($d3, $g3) < 1e-10, 'moved sensor: equal to the direct path');
    interpolate_topo($px, $py, $v2, res => 41, clip => 0.54);
    is(PDL::EEG::MAP2D::_topo_cache_stats()->{builds}, 3, 'other resolution: K built again');

    # negative control: the K of the first layout does not fit the moved one
    my $Kold = PDL::EEG::MAP2D::_topo_kernel($px, $py, 40, 0.54, 0.64, 16)->{K};
    my $wrong = ($Kold x $v2->dummy(0))->clump(-1)->reshape(40, 40);
    ok($rel->($d3, $wrong) > 1e-3, 'negative control: a stale K gives a different field');

    # the returned mask is the caller's own copy
    my (undef, $m1) = interpolate_topo($px, $py, $v, res => 40, clip => 0.54);
    $m1 .= 0;
    my (undef, $m2) = interpolate_topo($px, $py, $v, res => 40, clip => 0.54);
    ok($m2->sum->sclr > 0, 'changing the returned mask does not change the cached one');

    # the number of layouts kept is bounded
    PDL::EEG::MAP2D::_topo_cache_clear();
    {
        local $PDL::EEG::MAP2D::TOPO_CACHE_MAX = 3;
        interpolate_topo($px, $py, $v, res => $_, clip => 0.54) for 20 .. 25;
        is(PDL::EEG::MAP2D::_topo_cache_stats()->{entries}, 3, 'cache keeps at most TOPO_CACHE_MAX layouts');
        interpolate_topo($px, $py, $v, res => 25, clip => 0.54);
        is(PDL::EEG::MAP2D::_topo_cache_stats()->{builds}, 6, 'the most recent layout is still cached');
        interpolate_topo($px, $py, $v, res => 20, clip => 0.54);
        is(PDL::EEG::MAP2D::_topo_cache_stats()->{builds}, 7, 'the oldest layout was dropped');
    }
    PDL::EEG::MAP2D::_topo_cache_clear();
}

# ---- grid_layout_1020 -----------------------------------------------------------
{
    my $g = PDL::EEG::MAP2D::grid_layout_1020();
    is_deeply([ @$g{qw(nrows ncols)} ], [5, 5], 'grid_layout_1020: 5 x 5');
    my %seen; $seen{ join ',', @{ $g->{cell}{$_} } }++ for @chan;
    is(scalar keys %seen, 19, 'grid_layout_1020: the 19 channels have 19 different cells');
    is_deeply($g->{cell}{Cz},  [2, 2], 'Cz in the centre');
    is_deeply($g->{cell}{Fp1}, [0, 1], 'Fp1 top, left of the middle');
    is_deeply($g->{cell}{O2},  [4, 3], 'O2 bottom, right of the middle');
    is_deeply([ map { $g->{cell}{$_} } qw(T7 T8 P7 P8) ], [ map { $g->{cell}{$_} } qw(T3 T4 T5 T6) ],
              'new names share the cells of the old ones');
    ok(!$seen{ join ',', @{ $g->{map} } }, 'the map cell is not used by a channel');
    $g->{cell}{Cz} = [9, 9];
    is_deeply(PDL::EEG::MAP2D::grid_layout_1020()->{cell}{Cz}, [2, 2], 'each call returns a fresh table');
}

# ---- embedded map with a colour bar of PDL::Graphics::Cairo ---------------------
SKIP: {
    eval { require PDL::Graphics::Cairo; require PDL::Graphics::Cairo::Axes; 1 }
        or skip 'PDL::Graphics::Cairo not installed', 7;
    PDL::Graphics::Cairo::Axes->can('colorbar_reserve')
        or skip 'this PDL::Graphics::Cairo has no colorbar(cmap => ...)', 7;

    my $cm = PDL::EEG::MAP2D::_Cmap->new;
    my @lo = $cm->rgb_at(0); my @mid = $cm->rgb_at(0.5); my @hi = $cm->rgb_at(1);
    ok($lo[2] > $lo[0] && $hi[0] > $hi[2], 'colour scale object: blue at the low end, red at the high end');
    ok((grep { $_ > 0.9 } @mid) == 3, 'colour scale object: near white in the middle');

    my $fig = PDL::Graphics::Cairo::figure(width => 400, height => 300);
    my $gs  = $fig->add_gridspec(1, 1, cell_margins => 0);
    my $ax  = $fig->add_subplot($gs->at(0, 0));
    my @r = PDL::EEG::MAP2D::plot_topomap(values => $v, labels => \@keep, montage => $mon, clim => 3,
        contours => 0, res => 40, _ax => $ax, _fig => $fig, colorbar => { width => 7, reserve => 0 });
    is_deeply(\@r, [-3, 3], 'embedded plot_topomap returns the colour limits');
    is_deeply([ @{ $ax->_colorbar }{qw(vmin vmax)} ], [-3, 3], 'colorbar => { ... }: bar attached with the map limits');
    @r = eval { PDL::EEG::MAP2D::plot_topomap(values => $v, labels => \@keep, montage => $mon, clim => pdl(2.5)->max,
        contours => 0, res => 40, _ax => $fig->add_subplot($gs->at(0, 0)), _fig => $fig) };
    is_deeply(\@r, [-2.5, 2.5], "clim may be a one-element piddle ($@)");
    eval { PDL::EEG::MAP2D::plot_topomap(values => $v, labels => \@keep, montage => $mon, clim => {}, res => 40,
        _ax => $ax, _fig => $fig) };
    like($@, qr/'clim' must be/, 'other references as clim are refused');
    my $bytes = eval { $fig->to_png };
    ok(!$@ && length($bytes) > 0, "embedded map with colour bar renders ($@)");
}

done_testing;
