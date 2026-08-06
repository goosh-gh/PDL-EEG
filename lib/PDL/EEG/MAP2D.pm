package PDL::EEG::MAP2D;

use strict;
use warnings;
use Carp qw(croak carp);
use Exporter 'import';

use PDL;
use PDL::MatrixOps qw(inv);
use PDL::EEG::IO::ASA ();

our $VERSION = '0.01';

our @EXPORT_OK = qw(
    plot_topomap
    project_positions
    interpolate_topo
);

# ---------------------------------------------------------------------------
# PDL::EEG::MAP2D - MNE-style round 2D scalp voltage topography for PDL::EEG
#
#   use PDL::EEG::MAP2D qw(plot_topomap);
#
#   plot_topomap(
#       values  => $avg->slice(":,($t)"),      # (nchan) voltages at one sample
#       labels  => [qw(Fp1 Fp2 F7 F3 Fz ... O2)],
#       montage => '/.../standard_1020.elc',
#       clim    => 5,                           # symmetric +/- uV (default max|v|)
#       contours=> 6,
#       names   => 0,
#       title   => 'SEP 20 ms',
#       outfile => 'topo.png',
#   );
#
# The voltages live in $avg[chan,time]; positions come from an ASA .elc read by
# PDL::EEG::IO::ASA::read_elc.  Old 10-20 names (T3/T4/T5/T6) are aliased to the
# 10-10 names used by standard_1020.elc (T7/T8/P7/P8).
# ---------------------------------------------------------------------------

my $PI = 4 * atan2(1, 1);

# old 10-20  <->  standard_1020.elc (10-10)
my %ALIAS = (
    T3 => 'T7', T4 => 'T8', T5 => 'P7', T6 => 'P8',
    T7 => 'T3', T8 => 'T4', P7 => 'T5', P8 => 'T6',
);

# ColorBrewer RdBu reversed: blue(low) -> white(0) -> red(high), 11 stops.
my @RDBU_R = (
    [0.019, 0.188, 0.380], [0.129, 0.400, 0.674], [0.262, 0.576, 0.764],
    [0.572, 0.772, 0.870], [0.819, 0.898, 0.941], [0.968, 0.968, 0.968],
    [0.992, 0.858, 0.780], [0.956, 0.647, 0.509], [0.839, 0.376, 0.302],
    [0.698, 0.094, 0.168], [0.403, 0.000, 0.121],
);

# ---------------------------------------------------------------------------
# project_positions - 3D electrode coords -> 2D topomap coords
#
#   my ($px, $py, $ref) = project_positions($coords, scalp=>$mask, margin=>0.2);
#
# $coords is (3,N) [x=right, y=anterior, z=superior].  The sphere is fitted, and
# the scale/centre reference, to the *scalp* subset only - and within that, to
# sensors whose polar angle stays above the ears (periph_deg, default 115 deg),
# so deep periocular sensors (EOG, and aliases like X1) can't hijack the scale
# even when their name isn't recognised.  Azimuthal-equidistant projection
# (nose=+y up, right ear=+x); the array is then recentred on the reference
# centroid so it sits symmetrically in the head circle.  Outermost reference
# sensor lands at 0.5/(1+margin); head circle is 0.5.  Third return value is the
# boolean reference mask (true = scalp sensor used for scale/clim).
# ---------------------------------------------------------------------------
sub _sphere_center {
    my ($X, $Y, $Z) = @_;
    my $n  = $X->nelem;
    my $A  = cat(2*$X, 2*$Y, 2*$Z, ones($n))->xchg(0, 1);   # (4,n)
    my $bb = $X**2 + $Y**2 + $Z**2;
    my $At = $A->xchg(0, 1);
    my $c  = (inv($At x $A) x ($At x $bb->dummy(0)))->clump(-1);
    return ($c->at(0), $c->at(1), $c->at(2));
}
sub _azimuthal {
    my ($X, $Y, $Z, $cx, $cy, $cz) = @_;
    my $vx = $X - $cx; my $vy = $Y - $cy; my $vz = $Z - $cz;
    my $rr = sqrt($vx**2 + $vy**2 + $vz**2);
    my $theta = acos(($vz / $rr)->clip(-1, 1));   # 0 at vertex
    my $phi   = atan2($vy, $vx);
    return ($theta, $theta * cos($phi), $theta * sin($phi));
}
sub project_positions {
    my ($coords, %opt) = @_;
    my $margin   = $opt{margin} // 0.20;
    my $recenter = defined $opt{recenter} ? $opt{recenter} : 1;
    my $tcap     = ($opt{periph_deg} // 115) * $PI / 180;
    my $X = $coords->slice("(0),:")->flat;
    my $Y = $coords->slice("(1),:")->flat;
    my $Z = $coords->slice("(2),:")->flat;
    my $N = $X->nelem;
    my $named = defined $opt{scalp} ? $opt{scalp} : ones(long, $N);

    # pass 1: fit on name-scalp, get polar angles
    my $i1 = which($named); $i1 = sequence(long, $N) if $i1->nelem == 0;
    my ($cx, $cy, $cz) = _sphere_center($X->index($i1), $Y->index($i1), $Z->index($i1));
    my ($theta1) = _azimuthal($X, $Y, $Z, $cx, $cy, $cz);

    # reference = name-scalp AND not-deep (above the ear line)
    my $ref = $named & ($theta1 <= $tcap);
    my $ir  = which($ref);
    if ($ir->nelem < 3) { $ref = $named->copy; $ir = $i1; }   # fall back

    # pass 2: refit on the clean reference, project everything
    ($cx, $cy, $cz) = _sphere_center($X->index($ir), $Y->index($ir), $Z->index($ir));
    my (undef, $px, $py) = _azimuthal($X, $Y, $Z, $cx, $cy, $cz);

    # centre the array on the reference centroid (removes front/back bias)
    if ($recenter) {
        $px = $px - $px->index($ir)->avg;
        $py = $py - $py->index($ir)->avg;
    }
    my $ro = sqrt($px**2 + $py**2)->index($ir)->max;
    $ro = 1 if $ro <= 0;
    my $s = 0.5 / ($ro * (1 + $margin));
    return ($px * $s, $py * $s, $ref);
}

# ---------------------------------------------------------------------------
# Thin-plate-spline scattered interpolation (internal)
# ---------------------------------------------------------------------------
sub _tps_solve {
    my ($px, $py, $v) = @_;
    my $n  = $px->nelem;
    my $dx = $px->dummy(0) - $px->dummy(1);
    my $dy = $py->dummy(0) - $py->dummy(1);
    my $r2 = $dx*$dx + $dy*$dy;
    my $K  = 0.5 * $r2 * log($r2 + ($r2 == 0));    # U(r)=r^2 ln r ; diagonal -> 0

    my $L = zeroes($n + 3, $n + 3);
    $L->slice("0:" . ($n-1) . ",0:" . ($n-1)) .= $K;
    for my $i (0 .. $n-1) {
        my ($xi, $yi) = ($px->at($i), $py->at($i));
        $L->set($n+0, $i, 1);   $L->set($i, $n+0, 1);
        $L->set($n+1, $i, $xi); $L->set($i, $n+1, $xi);
        $L->set($n+2, $i, $yi); $L->set($i, $n+2, $yi);
    }
    my $rhs = zeroes($n + 3);
    $rhs->slice("0:" . ($n-1)) .= $v;
    return (inv($L) x $rhs->dummy(0))->clump(-1);   # [w(0..n-1), a0, a1, a2]
}

# ---------------------------------------------------------------------------
# interpolate_topo - scattered sensor values -> dense masked grid
#
#   my ($field, $inside, $extent) = interpolate_topo($px, $py, $vals, %opt);
#
# Returns the interpolated field (res,res) with dim0=x, dim1=y; a boolean
# in-circle mask; and [-EXT,EXT] extent.  Boundary anchors (inverse-distance
# on the head rim) keep the extrapolation to the rim tame.
# ---------------------------------------------------------------------------
sub interpolate_topo {
    my ($px, $py, $vals, %opt) = @_;
    my $res     = $opt{res}     // 220;
    my $R       = $opt{clip} // $opt{head_radius} // 0.5;   # colored-region radius
    my $EXT     = $opt{extent} // ($R + 0.10 > 0.6 ? $R + 0.10 : 0.6);
    my $nanchor = $opt{anchors} // 16;

    # rim anchors, value = inverse-distance-weighted sensor mean
    my ($fx, $fy, $fv) = ($px, $py, $vals);
    if ($nanchor > 0) {
        my $ang = sequence($nanchor) / $nanchor * 2 * $PI;
        my $ax  = $R * cos($ang);
        my $ay  = $R * sin($ang);
        my @av;
        for my $k (0 .. $nanchor - 1) {
            my $d2 = ($px - $ax->at($k))**2 + ($py - $ay->at($k))**2 + 1e-9;
            my $w  = 1 / $d2;
            push @av, ($w * $vals)->sum / $w->sum;
        }
        $fx = $px->append($ax);
        $fy = $py->append($ay);
        $fv = $vals->append(pdl(@av));
    }

    my $c  = _tps_solve($fx, $fy, $fv);
    my $nf = $fx->nelem;
    my $w  = $c->slice("0:" . ($nf-1));
    my ($a0, $a1, $a2) = ($c->at($nf), $c->at($nf+1), $c->at($nf+2));

    my $gl = (sequence($res) / ($res - 1)) * 2 * $EXT - $EXT;   # -EXT..EXT
    my $gx = $gl->dummy(1, $res);     # (res,res) x along dim0
    my $gy = $gl->dummy(0, $res);     # y along dim1
    my $GX = $gx->flat; my $GY = $gy->flat;

    my $DX = $GX->dummy(1, $nf) - $fx->dummy(0, $GX->nelem);
    my $DY = $GY->dummy(1, $nf) - $fy->dummy(0, $GY->nelem);
    my $R2 = $DX*$DX + $DY*$DY;
    my $U  = 0.5 * $R2 * log($R2 + ($R2 == 0));
    my $F  = ($U * $w->dummy(0, $GX->nelem))->xchg(0, 1)->sumover;   # sum over anchors+sensors
    $F += $a0 + $a1 * $GX + $a2 * $GY;

    my $field  = $F->reshape($res, $res);
    my $inside = (sqrt($gx**2 + $gy**2) <= $R);
    return ($field, $inside, $EXT, $R);
}

# vectorized RdBu_r lookup: value grid + limits -> R,G,B grids in [0,1]
sub _colorize {
    my ($field, $lo, $hi) = @_;
    my $nst = scalar @RDBU_R;
    my $sx  = sequence($nst) / ($nst - 1);
    my $sr  = pdl(map { $_->[0] } @RDBU_R);
    my $sg  = pdl(map { $_->[1] } @RDBU_R);
    my $sb  = pdl(map { $_->[2] } @RDBU_R);
    my $t   = (($field - $lo) / (($hi - $lo) || 1))->clip(0, 1)->flat;
    my $R   = $t->interpol($sx, $sr)->reshape($field->dims);
    my $G   = $t->interpol($sx, $sg)->reshape($field->dims);
    my $B   = $t->interpol($sx, $sb)->reshape($field->dims);
    return ($R, $G, $B);
}

# ---------------------------------------------------------------------------
# plot_topomap - the public entry point
# ---------------------------------------------------------------------------
sub plot_topomap {
    my (%a) = @_;

    require PDL::Graphics::Cairo;
    PDL::Graphics::Cairo->import(qw(figure));

    # ---- resolve values (nchan) -------------------------------------------
    my $labels = $a{labels} or croak "plot_topomap: 'labels' arrayref required";
    my $vals;
    if (defined $a{values}) {
        $vals = ref($a{values}) && eval { $a{values}->isa('PDL') }
              ? $a{values}->flat : pdl(@{ $a{values} });
    }
    elsif (defined $a{avg}) {
        my $t = $a{time} // croak "plot_topomap: 'time' index required with 'avg'";
        if (ref($a{avg}) && eval { $a{avg}->isa('PDL') }) {
            $vals = $a{avg}->slice(":,($t)")->flat;
        }
        else {   # array-of-arrays  $avg[chan][time]
            $vals = pdl(map { $_->[$t] } @{ $a{avg} });
        }
    }
    else {
        croak "plot_topomap: supply 'values' (nchan) or 'avg'+'time'";
    }
    croak "plot_topomap: values/labels length mismatch"
        if $vals->nelem != scalar @$labels;

    # ---- montage ----------------------------------------------------------
    my $mon = ref($a{montage}) ? $a{montage}
            : PDL::EEG::IO::ASA::read_elc($a{montage}
                  // croak "plot_topomap: 'montage' (.elc path or read_elc hash) required");

    # ---- match channels to sensor positions (with aliasing) ---------------
    my %extra_periph = map { uc($_) => 1 } @{ $a{periph} || [] };
    my (@x, @y, @z, @val, @keep, @scalp);
    for my $i (0 .. $#$labels) {
        my $name = $labels->[$i];
        next if _is_fiducial($name);
        my $key  = $mon->{pos}{$name}                 ? $name
                 : ($ALIAS{$name} && $mon->{pos}{$ALIAS{$name}}) ? $ALIAS{$name}
                 : undef;
        unless (defined $key) {
            carp "plot_topomap: channel '$name' not in montage - skipped";
            next;
        }
        my $p = $mon->{pos}{$key};
        push @x, $p->{x}; push @y, $p->{y}; push @z, $p->{z};
        push @val, $vals->at($i); push @keep, $name;
        push @scalp, (_is_periph($name) || $extra_periph{uc $name}) ? 0 : 1;
    }
    croak "plot_topomap: no channels matched the montage" unless @keep;
    my $coords = cat(pdl(@x), pdl(@y), pdl(@z))->xchg(0, 1);   # (3,N)
    my $Vv     = pdl(@val);
    my $smask  = pdl(long, @scalp);
    my $N      = scalar @keep;

    # ---- project + interpolate -------------------------------------------
    my ($px, $py, $ref) = project_positions($coords,
        scalp => $smask, margin => ($a{margin} // 0.20),
        recenter => (defined $a{recenter} ? $a{recenter} : 1),
        periph_deg => ($a{periph_deg} // 115));
    my $R_head    = 0.5;                                   # drawn head boundary
    my $overshoot = $a{overshoot} // 0.08;                 # color past the circle
    my $R_clip    = $R_head * (1 + $overshoot);
    # By default only scalp sensors shape the field; periocular (EOG) sensors are
    # plotted but excluded from the interpolation.  eog_interp => 1 lets the eye
    # channels feed the spline too (extends the frontal periphery - useful with a
    # broad montage such as the New York Head).
    my ($ipx, $ipy, $iv) = ($px, $py, $Vv);
    unless ($a{eog_interp}) {
        my $ri = which($ref);
        ($ipx, $ipy, $iv) = ($px->index($ri), $py->index($ri), $Vv->index($ri))
            if $ri->nelem >= 3;
    }
    my ($field, $inside, $EXT) =
        interpolate_topo($ipx, $ipy, $iv, res => ($a{res} // 220), clip => $R_clip);

    # ---- colour limits (from scalp reference, so EOG can't blow the scale) --
    my ($lo, $hi);
    if (defined $a{clim}) {
        ($lo, $hi) = ref($a{clim}) ? @{ $a{clim} } : (-$a{clim}, $a{clim});
    }
    else {
        my $ri  = which($ref);
        my $rv  = $ri->nelem ? $Vv->index($ri) : $Vv;
        my $m = $rv->abs->max; $m = 1e-6 if $m <= 0;
        ($lo, $hi) = (-$m, $m);
    }

    # ---- field -> RGB HWC (white outside the head) ------------------------
    my ($R, $G, $B) = _colorize($field, $lo, $hi);
    my $out = !$inside;
    $R->flat->where($out->flat) .= 1;
    $G->flat->where($out->flat) .= 1;
    $B->flat->where($out->flat) .= 1;
    my $hwc = cat($R, $G, $B)->xchg(0, 1);   # (x,y,3) -> (y,x,3), row0 = y-min

    my $res = ($field->dims)[0];
    my $d2p = sub { my $d = shift; ($d + $EXT) / (2 * $EXT) * $res };

    # ---- figure / axes ----------------------------------------------------
    my $want_cbar = $a{colorbar} // 1;
    my ($fig, $ax, $cax);
    if ($want_cbar) {
        $fig = figure(width => ($a{width} // 640), height => ($a{height} // 600));
        my $gs = $fig->add_gridspec(1, 2, width_ratios => [20, 1],
            left => 0.02, right => 0.82, top => 0.92, bottom => 0.05, wspace => 0.25);
        $ax  = $fig->add_subplot($gs->at(0, 0));
        $cax = $fig->add_subplot($gs->at(0, 1));
    }
    else {
        ($fig, $ax) = PDL::Graphics::Cairo::subplots(1, 1,
            figsize => [ ($a{width} // 560) / 100, ($a{height} // 600) / 100 ]);
    }

    $ax->imshow($hwc, origin => 'lower');
    $ax->set_aspect('equal');
    $ax->axis('off');

    # contour on a NaN-masked copy so lines stop at the head rim
    my $ncont = $a{contours} // 6;
    if ($ncont) {
        # Outside the head: fill with the interior mean (a flat constant region
        # carries no contour crossings), so lines stop at the rim without the
        # BAD-value handling that PDL::Graphics::Cairo's marching squares warns on.
        my $fm = $field->copy;
        my $inmean = $field->flat->where($inside->flat)->avg;
        $fm->flat->where($out->flat) .= $inmean;
        my $gpix = sequence($res) + 0.5;
        eval {
            $ax->contour($gpix, $gpix, $fm->xchg(0, 1),
                levels => $ncont, colors => 'black', lw => 0.6,
                vmin => $lo, vmax => $hi);
            1;
        } or carp "plot_topomap: contour skipped ($@)";
    }

    _draw_outline($ax, $d2p, $R_head, ($a{ear_dy} // 0.03));

    # ---- sensors ----------------------------------------------------------
    my $mark = $a{sensor_color} // '#333333';
    $ax->scatter(pdl(map { $d2p->($px->at($_)) } 0 .. $N-1),
                 pdl(map { $d2p->($py->at($_)) } 0 .. $N-1),
                 color => $mark, s => ($a{sensor_size} // 9));
    if ($a{names}) {
        for my $i (0 .. $N-1) {
            $ax->text($d2p->($px->at($i)) + 3, $d2p->($py->at($i)) + 4,
                $keep[$i], fontsize => ($a{name_size} // 10), color => 'black');
        }
    }
    $fig->suptitle($a{title}, fontsize => ($a{title_size} // 13)) if defined $a{title};

    # ---- colorbar ---------------------------------------------------------
    if ($want_cbar) {
        my $nb = 64;
        my $cf = (sequence($nb) / ($nb - 1)) * ($hi - $lo) + $lo;  # low..high
        my ($cr, $cg, $cb) = _colorize($cf->dummy(0, 1), $lo, $hi); # (1,nb)
        my $bar = cat($cr, $cg, $cb)->xchg(0, 1);                   # (nb,1,3)
        $cax->imshow($bar, origin => 'lower');
        $cax->set_aspect('auto');
        $cax->set_xticks(pdl([]));
        $cax->set_yticks(pdl([]));
        # All three numbers sit at the SAME x (centred), so they stay vertically
        # aligned and one knob slides them together, clear of the bar.  A
        # per-label offset scatters differently across text-engine builds, so we
        # avoid it.  cbar_label_x: larger = further from the bar (default 3.4).
        $cax->xlim(0, 1);
        my $xr = $a{cbar_label_x} // $a{xr} // 3.4;
        my %al = (ha => 'center', halign => 'center', va => 'center', valign => 'center');
        $cax->text($xr, $nb - 1.5,     sprintf('%.1f', $hi), fontsize => 9, %al);
        $cax->text($xr, ($nb - 1) / 2, '0',                  fontsize => 9, %al);
        $cax->text($xr, 1.5,           sprintf('%.1f', $lo), fontsize => 9, %al);
        $cax->set_title($a{unit} // 'uV', fontsize => 9);
    }

    # ---- output -----------------------------------------------------------
    $fig->save($a{outfile}) if $a{outfile};

    my $dev = lc($a{device} // $a{backend} // 'png');
    if ($dev ne 'png') {
        # macOS native display is the giza-server backend now; treat the old
        # osx/aqua/cocoa names as aliases for 'gs' so they don't warn.
        my %map = (osx => 'gs', aqua => 'gs', cocoa => 'gs',
                   gs  => 'gs',  giza => 'gs',
                   gnuplot => 'gnuplot', qt => 'gnuplot', x11 => 'gnuplot');
        $fig->show(backend => ($map{$dev} // $dev));
    }
    return $a{outfile} ? $a{outfile} : (wantarray ? ($fig, $ax) : $fig);
}

# --- classifiers -----------------------------------------------------------
sub _is_fiducial {
    my $n = shift;
    return $n =~ /^(?:LPA|RPA|Nz|Iz|NAS|NASION|CMS|DRL|Ref|Gnd|GND|COMNT|SCALE|
                     Fid(?:Nz|T9|T10)?)$/xi ? 1 : 0;
}
sub _is_periph {
    my $n = shift;   # periocular / EOG sensors sit well below the scalp
    return $n =~ /EOG|^(?:IO|SO|LO|IO[12]|SO[12]|LO[12])$/i ? 1 : 0;
}

# head circle + nose + ears in pixel coords
sub _draw_outline {
    my ($ax, $d2p, $R, $ear_dy) = @_;
    $ear_dy //= 0;
    my $tt = sequence(200) / 199 * 2 * $PI;
    $ax->line(pdl(map { $d2p->($R * cos($_)) } list $tt),
              pdl(map { $d2p->($R * sin($_)) } list $tt),
              color => 'black', lw => 2);

    my $bx = 0.10; my $by = sqrt($R**2 - $bx**2); my $apex = $R + 0.07;
    $ax->line(pdl($d2p->(-$bx), $d2p->(0), $d2p->($bx)),
              pdl($d2p->($by),  $d2p->($apex), $d2p->($by)),
              color => 'black', lw => 2);

    my @ex = (0.497,0.510,0.518,0.5299,0.5419,0.54,0.547,0.532,0.510,0.489);
    my @ey = (0.0555,0.0775,0.0783,0.0746,0.0555,-0.0055,-0.0932,-0.1313,-0.1384,-0.1199);
    my $sc = $R / 0.5;
    # $ear_dy raises the ear graphic so its centre sits on the C-row / midline
    $ax->line(pdl(map { $d2p->( $_ * $sc) } @ex),
              pdl(map { $d2p->($_ * $sc + $ear_dy) } @ey), color => 'black', lw => 2);
    $ax->line(pdl(map { $d2p->(-$_ * $sc) } @ex),
              pdl(map { $d2p->($_ * $sc + $ear_dy) } @ey), color => 'black', lw => 2);
}

1;

__END__

=head1 NAME

PDL::EEG::MAP2D - MNE-style round 2D scalp voltage topography

=head1 SYNOPSIS

    use PDL::EEG::MAP2D qw(plot_topomap);

    # $avg is (nchan, ntime); positions from standard_1020.elc
    plot_topomap(
        avg     => $avg,
        time    => $sample_index,
        labels  => \@channel_names,        # order matches rows of $avg
        montage => '/path/standard_1020.elc',
        clim    => 5,                      # +/- uV; omit for auto max|v|
        contours=> 6,
        title   => 'SEP 20 ms',
        outfile => 'topo.png',
    );

=head1 DESCRIPTION

Draws a circular scalp map with nose and ears from a one-sample voltage vector.
Electrode positions are read from an ASA C<.elc> file via L<PDL::EEG::IO::ASA>.
The sphere is fitted to the scalp sensors, every sensor is projected to 2D by an
azimuthal-equidistant projection (nose = +y up, right ear = +x), the array is
recentred on the scalp centroid, and the field is interpolated with a thin-plate
spline and clipped to the head disc.  Rendering uses L<PDL::Graphics::Cairo>.

Old 10-20 channel names C<T3 T4 T5 T6> are matched to the 10-10 names
C<T7 T8 P7 P8> used by C<standard_1020.elc>.  Fiducials (LPA/RPA/Nz/Iz) are
skipped.  Periocular sensors (names matching C</EOG/>, and any sensor whose
polar angle falls below C<periph_deg>, e.g. an aliased C<X1>) are excluded from
the sphere fit, the scale, and the colour limits; by default they are drawn but
do not feed the interpolation.

=head1 FUNCTIONS

=over 4

=item plot_topomap(%opt)

Renders one topomap.  Data is given as either C<values> (an C<nchan> vector) or
C<avg> + C<time>.  With C<outfile> it writes a PNG and returns the path;
otherwise it returns C<($fig, $ax)>.

Data and montage:

=over 4

=item C<values> => piddle | arrayref

Per-channel voltages at one latency (length C<nchan>).

=item C<avg> => piddle C<(nchan,ntime)> | array-of-arrays C<[chan][time]>

Averaged data; use with C<time>.

=item C<time> => integer

Sample index into C<avg>.

=item C<labels> => arrayref (required)

Channel names in the row order of C<values>/C<avg>.

=item C<montage> => path | hashref (required)

Path to an ASA C<.elc>, or a hash already returned by
C<PDL::EEG::IO::ASA::read_elc>.

=back

Field and layout:

=over 4

=item C<clim> => number | [lo, hi]

Colour limits in µV.  A scalar is symmetric (±). Default: C<max|v|> over the
scalp channels.

=item C<contours> => integer

Number of contour levels (default 6; 0 disables).

=item C<margin> => number

Fraction by which the outermost scalp sensor sits inside the head circle
(default 0.20).

=item C<overshoot> => number

Fraction by which the coloured disc extends outside the head circle
(default 0.08).

=item C<ear_dy> => number

Raise the drawn ears by this much (head-radius units) so their centre sits on
the C-row / midline (default 0.03; use 0 for the plain MNE ear placement).

=item C<recenter> => bool

Centre the array on the scalp centroid (default 1).

=item C<periph_deg> => number

Polar-angle threshold (degrees) above which a sensor is treated as periocular
(default 115).

=item C<eog_interp> => bool

Let periocular sensors feed the interpolation as well (default 0).

=item C<res> => integer

Interpolation grid size (default 220).

=back

Annotation and output:

=over 4

=item C<names> => bool

Draw channel labels (default 0).

=item C<name_size> => number

Channel-label font size (default 10).

=item C<sensor_color>, C<sensor_size>

Electrode marker colour (default C<#333333>) and size (default 9).

=item C<title> => string

Figure title (drawn with C<suptitle>).

=item C<unit> => string

Colour-bar unit label (default C<uV>).

=item C<colorbar> => bool

Draw the colour bar (default 1).

=item C<cbar_label_x> => number

Horizontal position of the three colour-bar numbers (default 3.4, in
bar-widths).  All three share this x and move together; adjust it so the column
clears the bar.  Accepts C<xr> as an alias.

=item C<outfile> => path

Write a PNG and return the path.

=item C<device> => 'png' | 'gs' | 'gnuplot'

Display backend when not (only) writing a file: C<gs> (giza-server window, the
native macOS backend; C<osx>/C<aqua>/C<cocoa> are accepted aliases), C<gnuplot>.
Default C<png>.

=back

=item project_positions($coords, %opt)

C<(3,N)> world coords to two C<(N)> topomap coords, plus a boolean reference
mask (scalp sensors used for scale/limits).  Options: C<scalp>, C<margin>,
C<recenter>, C<periph_deg>.

=item interpolate_topo($px, $py, $vals, %opt)

Scattered sensor values to a dense grid via thin-plate spline; returns
C<($field, $inside, $extent, $clip_radius)>.  Options: C<res>, C<clip>,
C<extent>, C<anchors>.

=back

=head1 SEE ALSO

L<PDL::EEG::IO::ASA>, L<PDL::Graphics::Cairo>.

=cut
