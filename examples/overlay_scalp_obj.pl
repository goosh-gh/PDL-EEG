#!/usr/bin/env perl
# overlay_scalp_obj.pl  —  P:EEG17 / NYHead03
# =============================================================================
# Overlay ASA (.elc) electrodes onto a New York Head SURFACE (scalp / cortex /
# any sa.<x>.vc+tri) and export ONE Wavefront .obj (+ .mtl) for Blender/MeshLab.
#
#   surface   : --surf /sa/head (default)  ->  <path>/vc (mm) + <path>/tri (1-based)
#               swap to /sa/cortex etc. to overlay electrodes on the brain/MRI.
#   electrodes: PDL::EEG::IO::ASA::read_elc  ->  coords [3,N] raw MNI mm
#
# NO ALIGNMENT: P:EEG16 proved .elc and NY Head share the same MNI frame/axes
# (raw residual mean 4.9 mm, no outlier); fiducial align made it worse. So the
# residual you SEE is the real fit, not a mis-registration.
#
# OUTPUT: each electrode is its own named object -> appears by label in the
# Blender outliner. With --labels (default) a readable, outward-facing 3D TEXT
# label is also built as thin ribbon geometry (visible even in Finder preview),
# so a full left/right swap is caught at a glance.
#
#     o scalp     usemtl surface     semi-transparent skin
#     o Fp1       usemtl electrode    octahedron marker
#                 usemtl label        3D text "Fp1"  (omit with --no-labels)
#
# USAGE
#   perl overlay_scalp_obj.pl --mat sa_nyhead.mat --elc standard_1020.elc \
#        --out nyhead_overlay.obj [--surf /sa/head] [--radius 4] \
#        [--scalp-alpha 0.35] [--no-labels] [--label-size 6] [--no-stats]
#   perl overlay_scalp_obj.pl --selftest        # no PDL / no data needed
# =============================================================================
use strict;
use warnings;
use Getopt::Long;

# ---- CLI --------------------------------------------------------------------
my $mat         = 'sa_nyhead.mat';
my $elc         = undef;
my $out         = 'nyhead_overlay.obj';
my $surf        = '/sa/head';   # HDF5 group holding vc + tri
my $radius      = 4.0;          # electrode marker radius (mm)
my $scalp_alpha = 0.35;         # surface opacity: 1 opaque, 0 invisible
my $labels      = 1;            # build visible 3D text labels
my $label_size  = 6.0;          # cap height of label text (mm)
my $stats       = 1;            # report electrode -> nearest-vertex distance
my $selftest    = 0;

GetOptions(
    'mat=s'         => \$mat,
    'elc=s'         => \$elc,
    'out=s'         => \$out,
    'surf=s'        => \$surf,
    'radius=f'      => \$radius,
    'scalp-alpha=f' => \$scalp_alpha,
    'labels!'       => \$labels,
    'label-size=f'  => \$label_size,
    'stats!'        => \$stats,
    'selftest'      => \$selftest,
) or die "bad options\n";

# =============================================================================
# GEOMETRY + WRITER  — pure Perl (so --selftest runs without PDL)
# =============================================================================

# tiny 3-vector helpers
sub _sub  { [ $_[0][0]-$_[1][0], $_[0][1]-$_[1][1], $_[0][2]-$_[1][2] ] }
sub _add  { [ $_[0][0]+$_[1][0], $_[0][1]+$_[1][1], $_[0][2]+$_[1][2] ] }
sub _scal { my($v,$s)=@_; [ $v->[0]*$s, $v->[1]*$s, $v->[2]*$s ] }
sub _cross{ my($a,$b)=@_; [ $a->[1]*$b->[2]-$a->[2]*$b->[1],
                            $a->[2]*$b->[0]-$a->[0]*$b->[2],
                            $a->[0]*$b->[1]-$a->[1]*$b->[0] ] }
sub _len  { sqrt($_[0][0]**2 + $_[0][1]**2 + $_[0][2]**2) }
sub _norm { my $v=shift; my $l=_len($v)||1e-12; [ $v->[0]/$l,$v->[1]/$l,$v->[2]/$l ] }

# Octahedron radius $r at ($cx,$cy,$cz): 6 verts, 8 outward-wound faces (1-based local).
sub octahedron {
    my ( $cx, $cy, $cz, $r ) = @_;
    my @v = (
        [ $cx+$r,$cy,   $cz    ], [ $cx,   $cy+$r,$cz    ], [ $cx-$r,$cy,   $cz    ],
        [ $cx,   $cy-$r,$cz    ], [ $cx,   $cy,   $cz+$r ], [ $cx,   $cy,   $cz-$r ],
    );
    my @f = (
        [5,1,2],[5,2,3],[5,3,4],[5,4,1],
        [6,2,1],[6,3,2],[6,4,3],[6,1,4],
    );
    return ( \@v, \@f );
}

# ---- stroke font (unit em: baseline y=0, cap y=1) ---------------------------
# each glyph: { w => advance, s => [ [ [x,y],... polyline ], ... ] }
my %GLYPH = (
    '0' => { w=>0.6, s=>[ [[0,0],[0.6,0],[0.6,1],[0,1],[0,0]], [[0,0],[0.6,1]] ] },
    '1' => { w=>0.5, s=>[ [[0.15,0.8],[0.32,1],[0.32,0]], [[0.1,0],[0.55,0]] ] },
    '2' => { w=>0.6, s=>[ [[0,0.78],[0.15,0.97],[0.45,0.97],[0.6,0.7],[0,0],[0.6,0]] ] },
    '3' => { w=>0.6, s=>[ [[0,0.9],[0.5,1],[0.6,0.75],[0.3,0.55],[0.6,0.32],[0.4,0],[0,0.1]] ] },
    '4' => { w=>0.6, s=>[ [[0.45,1],[0,0.32],[0.6,0.32]], [[0.45,1],[0.45,0]] ] },
    '5' => { w=>0.6, s=>[ [[0.6,1],[0,1],[0,0.55],[0.4,0.62],[0.6,0.35],[0.4,0],[0,0.06]] ] },
    '6' => { w=>0.6, s=>[ [[0.55,0.9],[0.28,1],[0.05,0.68],[0,0.3],[0.15,0.05],[0.45,0],[0.6,0.25],[0.45,0.46],[0.1,0.44]] ] },
    '7' => { w=>0.6, s=>[ [[0,1],[0.6,1],[0.25,0]] ] },
    '8' => { w=>0.6, s=>[ [[0.3,0.55],[0.06,0.7],[0.3,1],[0.54,0.7],[0.3,0.55],[0.06,0.3],[0.3,0],[0.54,0.3],[0.3,0.55]] ] },
    '9' => { w=>0.6, s=>[ [[0.5,0.55],[0.15,0.56],[0,0.75],[0.15,1],[0.45,0.96],[0.6,0.7],[0.55,0.3],[0.3,0],[0.05,0.1]] ] },
    'A' => { w=>0.6, s=>[ [[0,0],[0.3,1],[0.6,0]], [[0.12,0.4],[0.48,0.4]] ] },
    'C' => { w=>0.6, s=>[ [[0.6,0.85],[0.35,1],[0.1,0.8],[0,0.5],[0.1,0.2],[0.35,0],[0.6,0.15]] ] },
    'F' => { w=>0.55,s=>[ [[0,0],[0,1],[0.55,1]], [[0,0.55],[0.45,0.55]] ] },
    'I' => { w=>0.5, s=>[ [[0.25,0],[0.25,1]], [[0.05,1],[0.45,1]], [[0.05,0],[0.45,0]] ] },
    'M' => { w=>0.65,s=>[ [[0,0],[0,1],[0.32,0.4],[0.65,1],[0.65,0]] ] },
    'N' => { w=>0.6, s=>[ [[0,0],[0,1],[0.6,0],[0.6,1]] ] },
    'O' => { w=>0.6, s=>[ [[0,0],[0.6,0],[0.6,1],[0,1],[0,0]] ] },
    'P' => { w=>0.6, s=>[ [[0,0],[0,1],[0.5,1],[0.6,0.75],[0.5,0.5],[0,0.5]] ] },
    'T' => { w=>0.6, s=>[ [[0.3,0],[0.3,1]], [[0,1],[0.6,1]] ] },
    'Z' => { w=>0.6, s=>[ [[0,1],[0.6,1],[0,0],[0.6,0]] ] },
    # lowercase (x-height 0.6, ascender/descender used where needed)
    'h' => { w=>0.55,s=>[ [[0.05,1],[0.05,0]], [[0.05,0.5],[0.28,0.62],[0.5,0.48],[0.5,0]] ] },
    'p' => { w=>0.55,s=>[ [[0.05,0.6],[0.05,-0.35]], [[0.05,0.55],[0.35,0.6],[0.5,0.4],[0.35,0.2],[0.05,0.2]] ] },
    'z' => { w=>0.5, s=>[ [[0,0.6],[0.48,0.6],[0,0],[0.48,0]] ] },
    '?' => { w=>0.5, s=>[ [[0,0],[0.5,0],[0.5,1],[0,1],[0,0]] ] },   # fallback box
);

# string -> (\@polylines2d, total_width_em); undefined chars -> box, warned once
sub label_strokes {
    my ($str) = @_;
    my @out;
    my @missing;
    my $cx = 0;
    for my $ch ( split //, $str ) {
        if ( $ch eq ' ' ) { $cx += 0.4; next; }
        my $g = $GLYPH{$ch};
        if ( !$g ) { push @missing, $ch; $g = $GLYPH{'?'}; }
        for my $stroke ( @{ $g->{s} } ) {
            push @out, [ map { [ $_->[0] + $cx, $_->[1] ] } @$stroke ];
        }
        $cx += $g->{w} + 0.18;
    }
    warn "label: no glyph for '" . join( '', @missing ) . "'\n" if @missing;
    return ( \@out, $cx > 0 ? $cx - 0.18 : 0 );
}

# Build readable 3D text for $str at $pos, plane facing outward from $ctr.
# Returns (\@verts,\@faces) (faces 1-based LOCAL). Text is centred on the
# outward point and floated clear of the marker.
sub label_geometry {
    my ( $pos, $ctr, $str, $size, $radius ) = @_;
    my $n   = _norm( _sub( $pos, $ctr ) );
    my $ref = [ 0, 0, 1 ];
    my $x   = _cross( $ref, $n );
    if ( _len($x) < 1e-6 ) { $x = _cross( [ 0, 1, 0 ], $n ); }
    $x = _norm($x);
    my $y = _norm( _cross( $n, $x ) );                 # in-plane up
    my $o = _add( $pos, _scal( $n, 2.0 * $radius + 1 ) );

    my ( $strokes, $W ) = label_strokes($str);
    my $hw = 0.06;                                     # ribbon half-width (em)
    my ( @V, @F );
    my $to3d = sub {
        my ( $u, $v ) = @_;
        my $uu = ( $u - $W / 2 ) * $size;
        my $vv = ( $v - 0.5 ) * $size;
        return [ $o->[0] + $uu*$x->[0] + $vv*$y->[0],
                 $o->[1] + $uu*$x->[1] + $vv*$y->[1],
                 $o->[2] + $uu*$x->[2] + $vv*$y->[2] ];
    };
    for my $poly (@$strokes) {
        for my $i ( 0 .. $#$poly - 1 ) {
            my ( $p0, $p1 ) = ( $poly->[$i], $poly->[ $i + 1 ] );
            my ( $dx, $dy ) = ( $p1->[0]-$p0->[0], $p1->[1]-$p0->[1] );
            my $L  = sqrt( $dx*$dx + $dy*$dy ) || 1e-9;
            my ( $px, $py ) = ( -$dy/$L*$hw, $dx/$L*$hw );
            my $b = scalar @V;
            push @V,
                $to3d->( $p0->[0]+$px, $p0->[1]+$py ),
                $to3d->( $p1->[0]+$px, $p1->[1]+$py ),
                $to3d->( $p1->[0]-$px, $p1->[1]-$py ),
                $to3d->( $p0->[0]-$px, $p0->[1]-$py );
            push @F, [ $b+1, $b+2, $b+3 ], [ $b+1, $b+3, $b+4 ];
        }
    }
    return ( \@V, \@F );
}

sub write_obj {
    my (%a)   = @_;
    my $file  = $a{out};
    my $sv    = $a{scalp_v};
    my $sf    = $a{scalp_f};
    my $elecs = $a{electrodes};
    my $r     = $a{radius};
    my $alpha = $a{scalp_alpha};
    my $ctr   = $a{centroid} // [ 0, 0, 0 ];
    my $do_lab= $a{labels};
    my $lsize = $a{label_size} // 6;

    ( my $mtl = $file ) =~ s/\.obj$//i; $mtl .= '.mtl';
    ( my $mtl_base = $mtl ) =~ s{.*/}{};

    open my $M, '>', $mtl or die "write $mtl: $!";
    printf $M "newmtl surface\nKa 0.20 0.18 0.16\nKd 0.86 0.74 0.66\nKs 0.10 0.10 0.10\n"
            . "Ns 8\nillum 2\nd %.3f\nTr %.3f\n\n", $alpha, 1 - $alpha;
    print  $M "newmtl electrode\nKa 0.10 0.00 0.00\nKd 0.90 0.12 0.12\nKs 0.60 0.60 0.60\n"
            . "Ns 40\nillum 2\nd 1.000\nTr 0.000\n\n";
    print  $M "newmtl label\nKa 0.02 0.02 0.02\nKd 0.05 0.05 0.05\nKs 0.0 0.0 0.0\n"
            . "Ns 1\nillum 1\nd 1.000\nTr 0.000\n";
    close $M;

    open my $O, '>', $file or die "write $file: $!";
    print  $O "# New York Head surface + ASA electrode overlay\n";
    print  $O "# P:EEG17 / NYHead03  (overlay_scalp_obj.pl)\n";
    printf $O "mtllib %s\n\n", $mtl_base;

    printf $O "o scalp\nusemtl surface\n";
    printf $O "v %.6f %.6f %.6f\n", @$_ for @$sv;
    print  $O "\n";
    printf $O "f %d %d %d\n",       @$_ for @$sf;
    print  $O "\n";

    my $base = scalar @$sv;
    for my $e (@$elecs) {
        ( my $name = defined $e->{label} ? $e->{label} : '' ) =~ s/\s+/_/g;
        $name = 'elec' unless length $name;
        printf $O "o %s\n", $name;

        my ( $lv, $lf ) = octahedron( $e->{x}, $e->{y}, $e->{z}, $r );
        printf $O "usemtl electrode\n";
        printf $O "v %.6f %.6f %.6f\n", @$_ for @$lv;
        printf $O "f %d %d %d\n", $_->[0]+$base, $_->[1]+$base, $_->[2]+$base for @$lf;
        $base += scalar @$lv;

        if ($do_lab) {
            my ( $tv, $tf ) =
                label_geometry( [ $e->{x}, $e->{y}, $e->{z} ], $ctr,
                                $e->{label} // $name, $lsize, $r );
            if (@$tv) {
                printf $O "usemtl label\n";
                printf $O "v %.6f %.6f %.6f\n", @$_ for @$tv;
                printf $O "f %d %d %d\n", $_->[0]+$base, $_->[1]+$base, $_->[2]+$base for @$tf;
                $base += scalar @$tv;
            }
        }
        print $O "\n";
    }
    close $O;
    return $base;
}

sub centroid {
    my ($V) = @_;
    return [ 0, 0, 0 ] unless @$V;
    my @c = ( 0, 0, 0 );
    for my $v (@$V) { $c[0]+=$v->[0]; $c[1]+=$v->[1]; $c[2]+=$v->[2]; }
    return [ map { $_ / @$V } @c ];
}

sub scalp_stats {
    my ( $V, $elecs ) = @_;
    return unless @$elecs && @$V;
    my @dist;
    my ( $wi, $wd ) = ( 0, -1 );
    for my $ei ( 0 .. $#$elecs ) {
        my $e = $elecs->[$ei];
        my $best = 1e30;
        for my $v (@$V) {
            my ( $dx,$dy,$dz ) = ( $e->{x}-$v->[0], $e->{y}-$v->[1], $e->{z}-$v->[2] );
            my $d2 = $dx*$dx + $dy*$dy + $dz*$dz;
            $best = $d2 if $d2 < $best;
        }
        my $d = sqrt $best;
        push @dist, $d;
        ( $wi, $wd ) = ( $ei, $d ) if $d > $wd;
    }
    my @s = sort { $a <=> $b } @dist;
    my $sum = 0; $sum += $_ for @s;
    printf STDERR "electrode -> nearest surface vertex: mean %.2f  median %.2f  max %.2f mm\n",
        $sum/@s, $s[ int(@s/2) ], $s[-1];
    printf STDERR "  worst: %s (%.2f mm)\n", $elecs->[$wi]{label}, $wd;
}

# =============================================================================
# PDL LOADERS  (lazy; only needed for real runs)
# =============================================================================
sub _as_Nx3 {
    my ( $p, $name ) = @_;
    my @d = $p->dims;
    die "$name: expected 2-D, got ${\ scalar @d}-D (@d)\n" unless @d == 2;
    return $p if $d[1] == 3;
    if ( $d[0] == 3 ) { warn "$name stored as (3,N); transposing\n"; return $p->transpose->copy; }
    die "$name: no axis of length 3 in (@d)\n";
}

sub _face_base {
    my ( $tri, $nv ) = @_;
    my ( $min, $max ) = ( $tri->min, $tri->max );
    if ( $min < 0.5 ) {
        die sprintf "tri max %d out of range for %d verts (0-based)\n", $max, $nv
            if $max > $nv - 1 + 0.5;
        return 0;
    }
    die sprintf "tri max %d out of range for %d verts (1-based)\n", $max, $nv
        if $max > $nv + 0.5;
    return 1;
}

sub load_surface {
    my ( $file, $path ) = @_;
    $path =~ s{/$}{};
    require PDL;           PDL->import;
    require PDL::IO::HDF5; PDL::IO::HDF5->import;

    my $h = PDL::IO::HDF5->new($file);
    my ( $vc, $tri );
    eval {
        $vc  = _as_Nx3( $h->dataset("$path/vc")->get,  "$path/vc"  );
        $tri = _as_Nx3( $h->dataset("$path/tri")->get, "$path/tri" );
        1;
    } or die
        "could not read $path/vc + $path/tri ($@)"
      . "  -> is '$path' a surface (has vc/tri)?\n"
      . "     list candidates with:  h5ls -r $file | grep -iE 'vc|tri'\n"
      . "     a volumetric MRI (no tri) needs an isosurface first.\n";

    my $nv = $vc->dim(0);
    my $nf = $tri->dim(0);
    printf STDERR "surface %s: %d vertices, %d faces\n", $path, $nv, $nf;

    my @x = $vc->slice(':,(0)')->list;
    my @y = $vc->slice(':,(1)')->list;
    my @z = $vc->slice(':,(2)')->list;
    my @V = map { [ $x[$_], $y[$_], $z[$_] ] } 0 .. $#x;

    my $add = _face_base( $tri, $nv ) ? 0 : 1;
    my @a = $tri->slice(':,(0)')->list;
    my @b = $tri->slice(':,(1)')->list;
    my @c = $tri->slice(':,(2)')->list;
    my @F = map {
        [ int($a[$_]+0.5)+$add, int($b[$_]+0.5)+$add, int($c[$_]+0.5)+$add ]
    } 0 .. $#a;

    return ( \@V, \@F );
}

sub load_asa_electrodes {
    my ($file) = @_;
    require PDL;               PDL->import;
    require PDL::EEG::IO::ASA; PDL::EEG::IO::ASA->import('read_elc');

    my $elc = read_elc($file);
    # --- read_elc interface (P:EEG15): adjust these 2 lines if yours differs
    my $coords = $elc->{coords};    # PDL [3,N] raw MNI mm
    my $labels = $elc->{labels};    # arrayref of N labels
    # -----------------------------------------------------------------------
    my @d = $coords->dims;
    $coords = $coords->transpose->copy if @d == 2 && $d[0] != 3 && $d[1] == 3;
    my $N = $coords->dim(1);

    my @ex = $coords->slice('(0),')->list;
    my @ey = $coords->slice('(1),')->list;
    my @ez = $coords->slice('(2),')->list;

    my @E;
    for my $i ( 0 .. $N - 1 ) {
        my $lab = ( ref $labels eq 'ARRAY' && defined $labels->[$i] )
            ? $labels->[$i] : sprintf( 'E%03d', $i + 1 );
        push @E, { label => $lab, x => $ex[$i], y => $ey[$i], z => $ez[$i] };
    }
    printf STDERR "electrodes: %d from %s\n", $N, $file;
    return \@E;
}

# =============================================================================
# SELF-TEST
# =============================================================================
sub _parse_obj {
    my ($file) = @_;
    open my $F, '<', $file or die "read $file: $!";
    my ( @v, @groups, $cur );
    while (<$F>) {
        if    (/^o\s+(\S+)/)                 { $cur = { name=>$1, faces=>[], mats=>{} }; push @groups, $cur; }
        elsif (/^usemtl\s+(\S+)/)            { $cur->{mats}{$1}++ if $cur; }
        elsif (/^v\s+(\S+)\s+(\S+)\s+(\S+)/) { push @v, [ $1, $2, $3 ]; }
        elsif (/^f\s+(\d+)\s+(\d+)\s+(\d+)/) { push @{ $cur->{faces} }, [ $1, $2, $3 ]; }
    }
    close $F;
    return ( \@v, \@groups );
}

sub run_selftest {
    my ( $pass, $fail ) = ( 0, 0 );
    my $ok = sub { my ($c,$m)=@_; $c ? $pass++ : do { $fail++; print "  NOT OK: $m\n" } };

    # octahedron outward winding
    {
        my ( $v, $f ) = octahedron( 0, 0, 0, 1 );
        my $out = 1;
        for my $face (@$f) {
            my @p = map { $v->[ $_-1 ] } @$face;
            my $n = _cross( _sub($p[1],$p[0]), _sub($p[2],$p[0]) );
            my $ct = [ map { ($p[0][$_]+$p[1][$_]+$p[2][$_])/3 } 0..2 ];
            $out = 0 if $n->[0]*$ct->[0]+$n->[1]*$ct->[1]+$n->[2]*$ct->[2] <= 0;
        }
        $ok->( $out, "octahedron faces wound outward" );
    }

    # font covers a representative 10-05 montage
    {
        my @chars = split //, 'Fp1Fp2F7F3FzF4F8T7C3CzC4T8P7P3PzP4P8O1OzO2IzNzM1M2A1A2AFhFTCPPOTP';
        my %seen; $seen{$_}++ for @chars;
        my @miss = grep { $_ ne ' ' && !$GLYPH{$_} } keys %seen;
        $ok->( !@miss, "font covers montage chars (missing: @miss)" );
    }

    my ( $sv, $sf ) = octahedron( 0, 0, 0, 90 );
    my $ctr = centroid($sv);
    my @elecs = (
        { label=>'Fp1', x=> 40, y=> 70, z=>30 },
        { label=>'Fp2', x=>-40, y=> 70, z=>30 },
        { label=>'Cz',  x=>  0, y=>  0, z=>90 },
        { label=>'O 1', x=> 30, y=>-70, z=>20 },
    );
    my $Nsv = @$sv; my $Ne = @elecs;

    # labels ON
    my $tmp = "/tmp/st_lab_$$.obj";
    my $tot = write_obj( out=>$tmp, scalp_v=>$sv, scalp_f=>$sf, electrodes=>\@elecs,
        radius=>5, scalp_alpha=>0.4, centroid=>$ctr, labels=>1, label_size=>6 );
    my ( $V, $G ) = _parse_obj($tmp);

    $ok->( @$V == $tot,          "vertex count matches return ($tot)" );
    $ok->( @$G == 1 + $Ne,       "object count = surface + N electrodes" );
    $ok->( $G->[0]{name} eq 'scalp', "first object is surface" );
    $ok->( @$V > $Nsv + 6*$Ne,   "labels ON adds text geometry beyond markers" );

    my $inrange = 1; my $finite = 1;
    for my $vv (@$V) { $finite = 0 if grep { !/^-?\d+(?:\.\d+)?$/ } @$vv; }
    for my $g (@$G) { for my $f (@{$g->{faces}}) { $inrange = 0 if grep { $_<1 || $_>@$V } @$f; } }
    $ok->( $inrange, "all face indices within [1,Nverts]" );
    $ok->( $finite,  "all label/marker coords finite" );

    my $mats_ok = 1;
    for my $ei ( 0..$Ne-1 ) {
        my $m = $G->[$ei+1]{mats};
        $mats_ok = 0 unless $m->{electrode} && $m->{label};
    }
    $ok->( $mats_ok, "each electrode uses electrode + label materials" );
    $ok->( ( grep { $_->{name} eq 'O_1' } @$G ), "whitespace label sanitised to '_'" );
    unlink $tmp, ($tmp =~ s/\.obj$/.mtl/r);

    # labels OFF
    my $tmp2 = "/tmp/st_nolab_$$.obj";
    my $tot2 = write_obj( out=>$tmp2, scalp_v=>$sv, scalp_f=>$sf, electrodes=>\@elecs,
        radius=>5, scalp_alpha=>0.4, centroid=>$ctr, labels=>0 );
    my ( $V2, $G2 ) = _parse_obj($tmp2);
    $ok->( @$V2 == $Nsv + 6*$Ne, "labels OFF -> exactly surface + 6/electrode verts" );
    my $no_lab_mat = 1;
    for my $g ( @$G2 ) { $no_lab_mat = 0 if $g->{mats}{label}; }
    $ok->( $no_lab_mat, "labels OFF -> no label material used" );
    unlink $tmp2, ($tmp2 =~ s/\.obj$/.mtl/r);

    printf "\nselftest: %d passed, %d failed\n", $pass, $fail;
    exit( $fail ? 1 : 0 );
}

# =============================================================================
# MAIN
# =============================================================================
run_selftest() if $selftest;

die "need --elc FILE (ASA .elc electrodes)\n" unless defined $elc;
die "no such .mat: $mat\n"                    unless -f $mat;
die "no such .elc: $elc\n"                    unless -f $elc;

my ( $sv, $sf ) = load_surface( $mat, $surf );
my $elecs       = load_asa_electrodes($elc);
my $ctr         = centroid($sv);
scalp_stats( $sv, $elecs ) if $stats;

my $total = write_obj(
    out => $out, scalp_v => $sv, scalp_f => $sf, electrodes => $elecs,
    radius => $radius, scalp_alpha => $scalp_alpha, centroid => $ctr,
    labels => $labels, label_size => $label_size,
);

( my $mtl = $out ) =~ s/\.obj$/.mtl/;
printf "wrote %s (+ %s): %d surf verts, %d surf faces, %d electrodes%s, %d total verts\n",
    $out, $mtl, scalar(@$sv), scalar(@$sf), scalar(@$elecs),
    ( $labels ? " +labels" : "" ), $total;
