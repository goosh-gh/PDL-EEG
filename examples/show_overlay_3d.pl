#!/usr/bin/env perl
#
# show_overlay_3d.pl — electrodes_overlay.xyz を GS3D で 3D 重畳表示 / .obj 書き出し。
#
#   overlay_nyhead.pl が書く "name x y z SET" (SET = NY | ELC) を読み、
#   NY Head 19ch と standard_1020.elc 19ch を色分けで重ねる。両者は同一 MNI 枠
#   なので生座標のまま重ねれば ~5mm の近接ペアが 19 組見える。
#
#   NY = 赤   ELC = シアン   変位線 = 黄
#
#   ラベル(--labels, 既定 ON): 前後左右の入れ替わりが無いことを目視で保証するため、
#     左半球(x<0)の名前は .elc 側から、右半球・正中(x>=0)の名前は NY 側から出す。
#     左に "T7"(.elc)・右に "T8"(NY) が正しい位置に出れば L/R も A/P も健全。
#
#   perl -I<P:G:C>/lib show_overlay_3d.pl [--xyz electrodes_overlay.xyz]
#       [--no-labels] [--no-lines] [--obj overlay.obj] [--no-show]
#
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PDL;

my %opt = (xyz => 'electrodes_overlay.xyz', labels => 1, lines => 1,
           obj => undef, show => 1, radius => 3.0);
while (@ARGV) {
    my $a = shift @ARGV;
    if    ($a eq '--xyz')      { $opt{xyz}    = shift @ARGV }
    elsif ($a eq '--labels')   { $opt{labels} = 1 }
    elsif ($a eq '--no-labels'){ $opt{labels} = 0 }
    elsif ($a eq '--no-lines') { $opt{lines}  = 0 }
    elsif ($a eq '--obj')      { $opt{obj}    = shift @ARGV }
    elsif ($a eq '--radius')   { $opt{radius} = shift @ARGV }
    elsif ($a eq '--no-show')  { $opt{show}   = 0 }
    else  { die "unknown option: $a\n" }
}

# ---- read the 2-set point cloud ------------------------------------------
my (%NY, %ELC, @order);
open my $fh, '<', $opt{xyz} or die "open $opt{xyz}: $!";
while (<$fh>) {
    next if /^\s*#/ || /^\s*$/;
    my ($nm, $x, $y, $z, $set) = split;
    next unless defined $set;
    if    ($set eq 'NY')  { $NY{$nm}  = [$x, $y, $z]; push @order, $nm unless grep { $_ eq $nm } @order }
    elsif ($set eq 'ELC') { $ELC{$nm} = [$x, $y, $z] }
}
close $fh;
die "no NY/ELC points in $opt{xyz}\n" unless %NY && %ELC;
printf "overlay: %d NY + %d ELC electrodes from %s\n",
    scalar(keys %NY), scalar(keys %ELC), $opt{xyz};

# ---- build points / colors / labels, tracking 1-based vertex indices -----
my $C_NY  = [0.95, 0.35, 0.30];   # aka
my $C_ELC = [0.30, 0.75, 0.95];   # cyan
my $C_SEG = [0.85, 0.85, 0.35];   # yellow

my (@X, @Y, @Z, @R, @G, @B, @lab);
my (%nyv, %elcv);                 # name -> 1-based vertex index (for .obj lines)
my $vi = 0;
my $add = sub {
    my ($p, $c, $name) = @_;
    push @X, $p->[0]; push @Y, $p->[1]; push @Z, $p->[2];
    push @R, $c->[0]; push @G, $c->[1]; push @B, $c->[2];
    push @lab, $name;
    return ++$vi;
};
my $hemi = sub { $NY{$_[0]}[0] < 0 ? 'L' : 'R' };   # hemisphere from NY x (elc agrees)

# NY points: label only right/midline names (from NY)
for my $nm (@order) {
    my $l = ($opt{labels} && $hemi->($nm) eq 'R') ? $nm : '';
    $nyv{$nm} = $add->($NY{$nm}, $C_NY, $l);
}
# ELC points: label only left names (from .elc)
for my $nm (grep { $ELC{$_} } @order) {
    my $l = ($opt{labels} && $hemi->($nm) eq 'L') ? $nm : '';
    $elcv{$nm} = $add->($ELC{$nm}, $C_ELC, $l);
}

my $points = pdl([ \@X, \@Y, \@Z ]);   # (N,3)
my $colors = pdl([ \@R, \@G, \@B ]);   # (N,3)
my %scene = (points => $points, colors => $colors);
$scene{labels} = \@lab if $opt{labels};

# ---- displacement segments (NY -- ELC) -----------------------------------
my @seg = grep { $elcv{$_} } @order;    # names present in both sets
if ($opt{lines} && @seg) {
    my @c = map { [] } 0 .. 8;          # 9 columns: x0 y0 z0 x1 y1 z1 r g b
    for my $nm (@seg) {
        my ($a, $b) = ($NY{$nm}, $ELC{$nm});
        push @{$c[$_]},   $a->[$_] for 0 .. 2;
        push @{$c[3+$_]}, $b->[$_] for 0 .. 2;
        push @{$c[6]}, $C_SEG->[0]; push @{$c[7]}, $C_SEG->[1]; push @{$c[8]}, $C_SEG->[2];
    }
    $scene{lines} = pdl(\@c);           # (M,9)
}

# ---- optional .obj export (octahedron markers + .mtl, Blender-friendly) ----
# giza-server/P:G:C have no built-in .obj export (giza-server saves the current
# 2D view as PDF/SVG; GS3D is a point/line viewer). Blender's OBJ importer also
# ignores per-vertex 'v x y z r g b' colour and shows bare points as loose
# vertices — so instead each electrode is a small octahedron marker with a
# material, grouped into objects (NY / ELC / links). Colour survives via .mtl.
# When the scalp surface is added later, append its verts and 'f i j k' faces.
if ($opt{obj}) {
    (my $mtlpath = $opt{obj}) =~ s/\.obj$//i; $mtlpath .= '.mtl';
    (my $mtlref  = $mtlpath) =~ s{.*/}{};
    my $r = $opt{radius};
    open my $o, '>', $opt{obj} or die "write $opt{obj}: $!";
    print $o "# NY(red) + ELC(cyan) electrodes as octahedron markers; links = displacement\n";
    print $o "mtllib $mtlref\n";
    my $vc = 0;                              # running 1-based vertex count
    my @D = ([1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]);   # +x -x +y -y +z -z
    my @F = ([1,3,5],[3,2,5],[2,4,5],[4,1,5],[3,1,6],[2,3,6],[4,2,6],[1,4,6]);
    my $oct = sub {
        my ($name, $mat, $c) = @_;
        print $o "o $name\nusemtl $mat\n";
        printf $o "v %.3f %.3f %.3f\n", $c->[0]+$r*$_->[0], $c->[1]+$r*$_->[1], $c->[2]+$r*$_->[2] for @D;
        printf $o "f %d %d %d\n", $vc+$_->[0], $vc+$_->[1], $vc+$_->[2] for @F;
        $vc += 6;
    };
    $oct->("${_}_NY",  'mNY',  $NY{$_} ) for @order;
    $oct->("${_}_ELC", 'mELC', $ELC{$_}) for grep { $ELC{$_} } @order;
    if ($opt{lines} && @seg) {
        print $o "o links\nusemtl mSeg\n";
        for my $nm (@seg) {
            printf $o "v %.3f %.3f %.3f\n", @{$NY{$nm}};
            printf $o "v %.3f %.3f %.3f\n", @{$ELC{$nm}};
            printf $o "l %d %d\n", $vc+1, $vc+2;
            $vc += 2;
        }
    }
    close $o;
    open my $m, '>', $mtlpath or die "write $mtlpath: $!";
    print $m "newmtl mNY\nKd 0.95 0.35 0.30\n\n";
    print $m "newmtl mELC\nKd 0.30 0.75 0.95\n\n";
    print $m "newmtl mSeg\nKd 0.85 0.85 0.35\n";
    close $m;
    printf "wrote %s + %s  (%d markers r=%.1fmm%s)\n",
        $opt{obj}, $mtlref, scalar(@order) + scalar(grep { $ELC{$_} } @order),
        $r, ($opt{lines} && @seg) ? " + links" : "";
}

# ---- render ---------------------------------------------------------------
exit 0 unless $opt{show};
unless (eval { require PDL::Graphics::Cairo::Driver::GS3D; 1 }) {
    die "needs PDL::Graphics::Cairo (Driver::GS3D) on \@INC:\n$@";
}
print "GS3D:  NY = red,  ELC = cyan,  displacement = yellow\n";
print "       labels: left = .elc, right/mid = NY (L/R & A/P sanity)\n";
print "       drag to rotate, r = canonical, p = perspective/ortho, close to quit\n";
PDL::Graphics::Cairo::Driver::GS3D->new(
    scene => \%scene, width => 750, height => 750,
)->run;
