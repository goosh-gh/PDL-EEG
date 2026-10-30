#!/usr/bin/env perl
#
# nyhead_atlas_sagittal.pl — the New York Head laid over a side view of an anatomy atlas (Open
# Anatomy Atlas / BodyParts3D, the face to the left) in the midsagittal plane, scaled to the atlas
# head, to see where the NYHead neck electrodes and the neck model end fall on a whole body
# (neck, spine, chest), and where own neck sites placed by tape-measure lengths fall.
#
#   perl examples/nyhead_atlas_sagittal.pl SKIN.png BONE.png --mat ~/src/NYHead/sa_nyhead.mat
#   perl examples/nyhead_atlas_sagittal.pl SKIN.png BONE.png --points neck_points_atlas.tsv --mark 'Cv12=85,Cv45=125,BN1=165,Th12=200'
#   perl examples/nyhead_atlas_sagittal.pl --selftest
#
#   SKIN BONE        the two atlas views as for atlas_neck_levels.pl (same camera, skin / bones)
#   --mat FILE       sa_nyhead.mat                                         (~/src/NYHead/sa_nyhead.mat)
#   --points TSV     NYHead points (name, x, y, z, kind) from PDL-IO-NYHead nyhead_neck_electrodes.pl
#                    --points-out: kind point (own sites), profile_mesh / profile_extended (the back
#                    midline), level (atlas levels put on NYHead); drawn after the same transform
#   --mark LIST      NAME=MM,...  own sites by their length from the inion along the atlas outline
#                    (the subject's tape measure; the atlas scaled by --c7-mm)
#                                                                            (Cv12=85,Cv45=125,BN1=165,Th12=200)
#   --inion-row, --tips, --box, --smooth, --c7-mm    as for atlas_neck_levels.pl
#   --fit-above MM   NYHead head outline points used in the fit: z above Iz + MM  (0)
#   --fit-rows PX    atlas outline pixels used in the fit: rows above the atlas inion + PX (80)
#   --neck-pivot X,Y the centre of the neck turn, atlas SKIN pixels (with --neck angle and arc: the
#                    upper end of the NYHead chord --ny-neck as placed by the head fit; otherwise 25
#                    mm forward and 10 mm up of the C1 tip, about the joint between skull and C1)
#   --neck-band MM   the turn grows from 0 at the height of the pivot to full at MM below it (10)
#   --neck MODE      arc: each NYHead neck site by its length along the NYHead back midline,
#                    put at the same length along the atlas back outline (B); the neck is drawn
#                    turned as with angle;
#                    angle: turn the neck so that the direction of the back of the NYHead neck
#                    (--ny-neck) is that of the atlas neck (--atlas-neck), and read the sites where
#                    they then fall;
#                    fit: turn it to bring the NYHead back of the neck closest to the atlas back
#                    outline (positions, which also take up differences of neck length);
#                    none: the head similarity only. arc and angle need the profile rows of --points
#                                                                            (arc)
#   --ny-neck A,B    the NYHead back-of-the-neck chord, each end: a length along its back midline
#                    profile from Iz (mm), an electrode name (its projection on the profile),
#                    'end' (the last mesh vertex, where /sa/head already turns forward into its
#                    cap) or @name (the electrode itself; Nk2 lies on the surface of the FEM
#                    model below that turn)                                  (Exxz,@Nk2)
#   --atlas-neck A,B the atlas back-of-the-neck chord: lengths from the atlas inion along its back
#                    outline, mm at the subject scale (--c7-mm)              (85,125: Cv12, Cv45 = Nk2)
#   --nasion X,Y     the atlas nasion, SKIN pixels (the inion is that of atlas_levels.pl)
#                                                                            (170,212: the
#                    Open Anatomy Atlas view of 2026-10-28, as used for the nasion-inion arc)
#   --head-fit F     all: fit the whole head above Iz (its outline); back: only its back half
#                    (behind the middle of Iz and Nz, or Fpz; the scale held at that of all: the
#                    back alone is nearly one arc and turns freely along it); IN: Iz and Nz (or Fpz)
#                    put on the atlas inion and nasion (the head tilt by the I-N line); all three
#                    are printed, F is used (IN: with Iz and Nz on the atlas landmarks the head
#                    outline RMS was 6.2 mm against 5.5 for all, and Exxz came 2.4 mm from the
#                    atlas back outline against 22 mm: the outline fit does not take up the tilt of
#                    the atlas head, chin up 14 deg against NYHead)          (IN)
#   --no-labels      the picture without the text labels (the marks, lines, axes and title stay),
#                    to place labels by hand
#   --out PREFIX     PREFIX_sagittal.png and PREFIX_atlas_in_nyhead.tsv     (nyhead_atlas)
#   --strings S      auto | h5dump | lowlevel (how the NYHead names are read) (auto)
#
# The head fit (A): the cut of the NYHead scalp mesh (/sa/head) by the plane x = 0 above the
# height of Iz (forehead, vertex, back of the head, face down to that height), or its back half,
# is matched to the outline of the atlas body above its inion (+ --fit-rows) by a similarity (one
# scale, a turn, a shift) found by iterated closest points both ways (point-to-line steps; pairs
# farther than 2.5 x the median left out), starting from Iz on the atlas inion and the NYHead
# vertex height on the atlas head height.  Printed for both: the scale, turn, RMS, and how far the
# moved Iz, Nz (or Fpz) and Exxz lie from the atlas inion, nasion and back outline.  The atlas
# outline is a silhouette (the outermost skin), the NYHead cut is the midline.
#
# The angles: in each body on its own, the I-N line (inion -> nasion; NYHead Iz -> Nz or Fpz) and
# the back-of-the-neck line (NYHead: --ny-neck, from the height of the Exx ring to Nk2, below the
# forward turn of the end of /sa/head; atlas: --atlas-neck), and the angle between them; their difference is the difference of the
# head-neck posture.  After the head fit, the turn of the neck that makes the two neck lines
# parallel (the head kept as fitted) is the neck turn of --neck angle and arc.  (A direction taken
# below the deepest point of the atlas neck follows the slope into the upper back, not the neck:
# tried, -24 deg.)
#
# The neck sites (B, --neck arc): the neck length relative to the head differs between people (the
# atlas neck is short against NYHead and the subject YJ0394VB: inion -> C7 at the head scale 118
# against 165 mm), and no turn takes that up, so the sites are placed by their length along the
# back midline.  The level is read twice: with the NYHead neck taken as the subject's (Nk2 at
# 124.6 mm and the subject's Cv45 at 125 mm lie together), and at the head scale (the atlas body).
# A fit of positions (--neck fit) mixes posture and length and swings a longer neck into the body;
# a free scale of the neck as well collapses it onto a short stretch of the outline (tried: 0.58,
# -26 deg).
#
# Written: PREFIX_atlas_in_nyhead.tsv, the atlas levels (skin points and tips) and the atlas back
# outline taken into NYHead MNI mm by the head fit and the neck turn (x = 0; kinds level, tip,
# profile_atlas), and the picture PREFIX_sagittal.png.  The atlas is one person in one posture.

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use PDL;
use PDL::EEG::ElectrodeLocation::Polyline qw(arc_length point_at_length project_point);
use PDL::EEG::ElectrodeLocation::AtlasView qw(read_png foreground body_outline atlas_levels level_of);

my $TIPS = 'C1=415,331;C2=422,388;C3=419,414;C4=415,437;C5=422,458;C6=435,474;C7=449,502;'
         . 'T1=456,536;T2=463,571;T3=482,610';
my %o = (mat => "$ENV{HOME}/src/NYHead/sa_nyhead.mat", points => undef, mark => 'Cv12=85,Cv45=125,BN1=165,Th12=200',
         inion_row => 266, tips => $TIPS, box => '15,8', smooth => 4, c7_mm => 165, mm_per_px => undef,
         fit_above => 0, fit_rows => 80, out => 'nyhead_atlas', strings => 'auto', selftest => 0,
         neck_pivot => undef, neck_band => 10, neck => 'arc', ny_neck => 'Exxz,@Nk2', atlas_neck => '85,125',
         nasion => '170,212', head_fit => 'IN', labels => 1);
GetOptions('mat=s' => \$o{mat}, 'points=s' => \$o{points}, 'mark=s' => \$o{mark}, 'inion-row=i' => \$o{inion_row},
           'tips=s' => \$o{tips}, 'box=s' => \$o{box}, 'smooth=i' => \$o{smooth}, 'c7-mm=f' => \$o{c7_mm},
           'fit-above=f' => \$o{fit_above}, 'fit-rows=i' => \$o{fit_rows}, 'out=s' => \$o{out},
           'strings=s' => \$o{strings}, 'selftest' => \$o{selftest},
           'neck-pivot=s' => \$o{neck_pivot}, 'neck-band=f' => \$o{neck_band}, 'neck=s' => \$o{neck}, 'ny-neck=s' => \$o{ny_neck}, 'atlas-neck=s' => \$o{atlas_neck},
           'nasion=s' => \$o{nasion}, 'head-fit=s' => \$o{head_fit}, 'labels!' => \$o{labels})
    or die "bad options (see the header of $0)\n";
if ($o{selftest}) { exit(selftest() ? 0 : 1) }
@ARGV == 2 or die "usage: $0 SKIN.png BONE.png [options]  (see the header)\n";
my ($fskin, $fbone) = @ARGV;
-f $o{mat} or die "not found: $o{mat} (give it with --mat)\n";
require PDL::IO::NYHead;

# ---- the atlas ------------------------------------------------------------------------
my $skin = read_png($fskin);
my $bone = read_png($fbone);
my $at   = atlas_levels(\%o, $skin, $bone);                    # inion, levels, back outline, mm/px (C7)
my $B    = body_outline(foreground($skin));                    # (m,2) px
my $Bfit = $B->dice_axis(0, which($B->slice(':,(1)') <= $at->{inion_skin}->at(1) + $o{fit_rows}));

# ---- NYHead ---------------------------------------------------------------------------
my $ny  = PDL::IO::NYHead->new($o{mat}, strings => $o{strings});
my @lab = @{ $ny->electrode_labels };
my $pos = $ny->electrode_pos->double;                           # (ne,3)
my %ix; @ix{@lab} = 0 .. $#lab;
exists $ix{Iz} or die "no electrode Iz in $o{mat}\n";
my $hd  = $ny->surface('head');
my $cut = mesh_cut_x($hd->{vc}->double, $hd->{tri0}, 0);        # (n,3)
my $iz  = $pos->slice("($ix{Iz}),");
my $C2  = $cut->slice(':,1:2');                                 # (n,2): y, z
my $Pfit = $C2->dice_axis(0, which($cut->slice(':,(2)') >= $iz->at(2) + $o{fit_above}));
printf "NYHead: /sa/head cut by x = 0: %d points, %d above Iz%+g mm used; atlas outline: %d pixels, %d above the inion + %d rows used\n",
    $cut->dim(0), $Pfit->dim(0), $o{fit_above}, $B->dim(0), $Bfit->dim(0), $o{fit_rows};

# ---- the head fit (A): the whole head, and the back half of the head only ---------------
# the atlas inion is that of atlas_levels (inion_skin); the nasion is --nasion (SKIN px)
my $nas = pdl(split /,/, $o{nasion});
my $nyN = exists $ix{Nz} ? 'Nz' : exists $ix{Fpz} ? 'Fpz' : die "no electrode Nz or Fpz in $o{mat}\n";
my $nz  = $pos->slice("($ix{$nyN}),");
my $top_px = $Bfit->slice(':,(1)')->min;
my $s0 = ($at->{inion_skin}->at(1) - $top_px) / ($Pfit->slice(':,(1)')->max - $iz->at(2));   # px per mm
my $T0 = { s => $s0, th => 0, t => to_uv($at->{inion_skin}) - $s0 * $iz->slice('1:2') };
my $ymid = 0.5 * ($iz->at(1) + $nz->at(1));
my %region = (all => $Pfit, back => $Pfit->dice_axis(0, which($Pfit->slice(':,(0)') < $ymid)));
my %fit;
print "\nhead fit (NYHead mm per atlas px, turn, RMS; then how far the moved NYHead landmarks lie from the atlas ones, NYHead mm)\n";
# IN: the two landmarks alone (Iz on the atlas inion, Nz or Fpz on the atlas nasion): the head
# tilt by the I-N line, the scale by the I-N distance
my $in_fit = do {
    my ($pI, $pN) = ($iz->slice('1:2'), $nz->slice('1:2'));
    my ($qI, $qN) = (to_uv($at->{inion_skin}), to_uv($nas));
    my ($dp, $dq) = ($pN - $pI, $qN - $qI);
    my $sc = sqrt(($dq**2)->sum / ($dp**2)->sum)->sclr;
    my $th = (atan2($dq->at(1), $dq->at(0)) - atan2($dp->at(1), $dp->at(0)));
    my $F = { s => $sc, th => $th, t => pdl(0, 0) };
    $F->{t} = $qI - apply_sim($F, $pI->dummy(0))->slice('(0),');
    my (undef, $dd) = nearest(apply_sim($F, $Pfit), to_uv($Bfit));   # the head outline against it
    $F->{rms} = sqrt(($dd**2)->avg)->sclr;
    $F;
};
for my $r ('all', 'back', 'IN') {
    # the back half alone is nearly an arc of one curvature: its scale and shift trade off, so the
    # scale is held at that of the whole head (tried free: 769 mm/px) and the turn and shift fitted
    my $Tr = $r eq 'all' ? fit_similarity_icp($region{$r}, to_uv($Bfit), $T0)
           : $r eq 'back' ? fit_similarity_icp($region{$r}, to_uv($Bfit), $fit{all}, fix_s => $fit{all}{s}, cut_front => 1)
           : $in_fit;
    my $mp = sub { from_uv(apply_sim($Tr, $_[0]->dummy(0)))->slice('(0),') };
    my $dI = sqrt((($mp->($iz->slice('1:2')) - $at->{inion_skin})**2)->sum)->sclr / $Tr->{s};
    my $dN = sqrt((($mp->($nz->slice('1:2')) - $nas)**2)->sum)->sclr / $Tr->{s};
    my $dX = exists $ix{Exxz} ? (project_point($at->{P}, $mp->($pos->slice("($ix{Exxz}),1:2")), arc_length($at->{P})))[1] / $Tr->{s} : -1;
    my $what = $r eq 'all' ? sprintf('%3d points above Iz', $region{all}->dim(0))
             : $r eq 'back' ? sprintf('%3d points above Iz, behind y %.0f; scale held', $region{back}->dim(0), $ymid)
             : "Iz, $nyN on the atlas inion, nasion; head outline RMS against it";
    printf "  %-4s %s: %.4f mm/px, turned %+5.1f deg, RMS %.1f mm; Iz -> atlas inion %.1f, %s -> atlas nasion %.1f, Exxz -> atlas back outline %.1f\n",
        $r, $what, 1 / $Tr->{s}, 57.29578 * $Tr->{th}, $Tr->{rms} / $Tr->{s}, $dI, $nyN, $dN, $dX;
    $fit{$r} = $Tr;
}
exists $fit{ $o{head_fit} } or die "--head-fit all | back | IN, not '$o{head_fit}'\n";
our $T = $fit{ $o{head_fit} };
printf "used: --head-fit %s; atlas scale: head %.4f mm/px (NYHead head), subject inion -> C7 = %g mm: %.4f mm/px (ratio %.2f)\n",
    $o{head_fit}, 1 / $T->{s}, $o{c7_mm}, $at->{mm}, $at->{mm} * $T->{s};

# ---- the --points file ------------------------------------------------------------------
my (@ptin, @prof);                                              # [name, (y,z)]; [mesh|extended, (y,z)]
if (defined $o{points}) {
    open my $fh, '<', $o{points} or die "$o{points}: $!";
    while (<$fh>) {
        chomp; next if /^\s*(#|$)/;
        my ($n, $x, $y, $z, $k) = split /\t/;
        next if $n eq 'name';
        $k //= 'point';
        if    ($k eq 'point')                     { push @ptin, [ $n, pdl($y, $z) ] }
        elsif ($k =~ /^profile_(mesh|extended)$/) { push @prof, [ $1, pdl($y, $z) ] }
    }
}

# ---- the neck angle ----------------------------------------------------------------------
# NYHead below the pivot is turned about it (the weight of the turn grows from 0 at the height of
# the pivot to 1 at --neck-band mm below it), by the angle that brings the NYHead back of the
# neck (the /sa/head back midline below Iz, and Nk2 on the surface of the finite-element model)
# closest to the atlas back outline.
my $P = $at->{P}; my $cum = arc_length($P); my $sI = $at->{sI};
my %lv = map { $_->{name} => $_ } @{ $at->{levels} };
my $pv = defined $o{neck_pivot} ? pdl(split /,/, $o{neck_pivot})
       : exists $lv{C1} ? $lv{C1}{tip} + pdl(-25 * $T->{s}, -10 * $T->{s})   # 25 mm forward, 10 mm up of the C1 tip
       : die "--neck-pivot needed (no level C1)\n";
my $NM = { pivot => to_uv($pv), phi => 0, k => 1, band => $o{neck_band} };
$NM->{z0} = inv_sim($T, to_uv($pv)->dummy(0))->at(0, 1);       # NYHead height of the pivot
my @nk = map { $_->[1] } grep { $_->[0] eq 'mesh' } @prof;
my $Nset = @nk >= 2 ? do {                                     # the mesh profile (its vertices), every 2 mm along it
    my $pl = cat(@nk)->transpose; my $pc = arc_length($pl);
    point_at_length($pl, sequence(int($pc->at(-1) / 2) + 1) * 2, $pc) } : do {                     # without --points: the cut, back half, Iz .. 95 mm below
    my $w = which(($C2->slice(':,(0)') < $iz->at(1) + 40) & ($C2->slice(':,(1)') < $iz->at(2)) & ($C2->slice(':,(1)') > $iz->at(2) - 95));
    $C2->dice_axis(0, $w) };
$Nset = $Nset->dice_axis(0, which($Nset->slice(':,(1)') < $NM->{z0}));
$Nset = $Nset->glue(0, $pos->slice("($ix{Nk2}),1:2")->dummy(0)) if exists $ix{Nk2};
my $rms0 = neck_rms($NM, 0, 1, $Nset, $P, $cum);
my $ang  = sub { my ($u, $v) = @_; 57.29578 * atan2($u->at(0) * $v->at(1) - $u->at(1) * $v->at(0), ($u * $v)->sum)->sclr };
my $unit = sub { $_[0] / sqrt(($_[0]**2)->sum) };
my $horiz = sub { 57.29578 * atan2($_[0]->at(1), $_[0]->at(0)) };              # from forward, + up
my $vert  = sub { 57.29578 * atan2($_[0]->at(0), -$_[0]->at(1)) };             # from straight down, + forward
my ($pl, $pc);                                                  # the NYHead back midline profile (mesh, then extended)
if (@prof >= 2) { $pl = cat(map { $_->[1] } @prof)->transpose; $pc = arc_length($pl) }
# ---- the angle between the I-N line and the back-of-the-neck line, in each body ----------
if (@nk >= 2) {
    my $plm = cat(@nk)->transpose; my $pcm = arc_length($plm);
    # the two ends of the chord: a length along the mesh profile (mm), 'end' (the last vertex of the
    # mesh, where /sa/head already turns forward into its cap), an electrode name (its projection on
    # the profile) or '@name' (the electrode itself, e.g. @Nk2 on the surface of the FEM model)
    my $endpt = sub {
        my ($e) = @_;
        return $pos->slice("($ix{$1}),1:2")->copy if $e =~ /^@(.+)$/ && exists $ix{$1};
        my $L = $e eq 'end' ? $pcm->at(-1) : exists $ix{$e} ? (project_point($plm, $pos->slice("($ix{$e}),1:2"), $pcm))[0] : $e;
        $L =~ /^[\d.]+$/ or die "--ny-neck A,B: a length (mm), 'end', an electrode name or \@name; got '$e'\n";
        return point_at_length($plm, pdl($L), $pcm)->slice('(0),')->copy;
    };
    my ($na, $nb) = split /,/, $o{ny_neck};
    my $ny2 = cat($endpt->($na), $endpt->($nb))->transpose->copy;                       # (2,2): (y,z)
    sqrt((($ny2->slice('(1),') - $ny2->slice('(0),'))**2)->sum) > 10 or die "--ny-neck $o{ny_neck}: the two ends are less than 10 mm apart\n";
    my $dny = $unit->($ny2->slice('(1),') - $ny2->slice('(0),'));
    my $iny = $unit->($nz->slice('1:2') - $iz->slice('1:2'));
    my ($aa, $ab) = split /,/, $o{atlas_neck};
    my $at2 = point_at_length($P, pdl($aa, $ab) / $at->{mm} + $sI, $cum);               # (2,2) px
    my $dat = $unit->(to_uv($at2->slice('(1),') - $at2->slice('(0),')));
    my $iat = $unit->(to_uv($nas - $at->{inion_skin}));
    my ($hn, $ha) = ($ang->($iny, $dny), $ang->($iat, $dat));
    print "\nhead and neck: the I-N line (inion -> nasion) and the back-of-the-neck line, each body on its own\n";
    printf "  NYHead (MNI y,z): Iz -> %s %+.1f deg from horizontal; neck %s (%.1f,%.1f -> %.1f,%.1f) %+.1f deg from vertical; neck from I-N %.1f deg\n",
        $nyN, $horiz->($iny), $o{ny_neck}, $ny2->slice('(0),')->list, $ny2->slice('(1),')->list, $vert->($dny), $hn;
    # how much the lower end matters: the same upper end to the last mesh vertex, the vertex before
    # it, and Nk2
    my @alt = (['end', $plm->slice('(-1),')], ['the vertex before the end', $plm->slice('(-2),')],
               exists $ix{Nk2} ? (['@Nk2', $pos->slice("($ix{Nk2}),1:2")]) : ());
    printf "    lower end: %s\n", join('; ', map { sprintf '%s %+.1f deg', $_->[0], $vert->($_->[1] - $ny2->slice('(0),')) } @alt);
    printf "  atlas (image):    inion -> nasion %+.1f deg from horizontal; neck %g -> %g mm from the inion (subject scale) %+.1f deg from vertical; neck from I-N %.1f deg\n",
        $horiz->($iat), $aa, $ab, $vert->($dat), $ha;
    printf "  the atlas neck is turned %+.1f deg against its head compared with NYHead (+ = towards the face)\n", $ha - $hn;
    my $duv = apply_sim({ s => 1, th => $T->{th}, t => pdl(0, 0) }, $dny->dummy(0))->slice('(0),');
    my $phi = $ang->($duv, $dat) / 57.29578;
    printf "  after the head fit (turn %+.1f deg; the I-N lines alone would turn %+.1f deg) the neck turn is %+.1f deg\n",
        57.29578 * $T->{th}, $ang->($iny, $iat), 57.29578 * $phi;
    if ($o{neck} eq 'angle' || $o{neck} eq 'arc') {
        $NM->{phi} = $phi;
        # the turn about the upper end of the NYHead chord (without --neck-pivot): that point stays
        # where the head fit put it (with IN, Exxz 2.4 mm from the atlas outline), and the chord
        # below it is made parallel to the atlas one.  (About a point near the skull-C1 joint the
        # chord turns parallel but moves off the skin: tried, 10-30 mm.)
        if (!defined $o{neck_pivot}) {
            $NM->{pivot} = apply_sim($T, $ny2->slice('(0),')->dummy(0))->slice('(0),')->copy;
            $NM->{z0} = $ny2->at(0, 1);
            $pv = from_uv($NM->{pivot});
        }
        my $m2 = from_uv(ny_to_uv($NM, $ny2));
        my @off = map { (project_point($P, $m2->slice("($_),"), $cum))[1] / $T->{s} } 0, 1;
        my @len = map { ((project_point($P, $m2->slice("($_),"), $cum))[0] - $sI) * $at->{mm} } 0, 1;
        printf "  turned about (%.0f, %.0f) px (NYHead z %.1f mm, band %g mm): the NYHead chord ends lie %.1f and %.1f mm (NYHead) off the atlas outline, at %.0f and %.0f mm from the atlas inion (subject scale)\n",
            $pv->list, $NM->{z0}, $NM->{band}, @off, @len;
    }
}
elsif ($o{neck} eq 'angle' || $o{neck} eq 'arc') { die "--neck $o{neck} needs the back midline profile in --points (profile_mesh rows)\n" }
if ($o{neck} eq 'fit') { @$NM{qw(phi k)} = fit_neck($NM, $Nset, $P, $cum, 0) if $Nset->dim(0) >= 3 }
elsif ($o{neck} !~ /^(angle|arc|none)$/) { die "--neck arc | angle | fit | none, not '$o{neck}'\n" }
($o{neck} ne 'arc' || defined $pl) or die "--neck arc needs the back midline profile in --points\n";

# ---- NYHead electrodes and points on the atlas ----------------------------------------
my (@rows, @pt);                                                # [name, px (2)]
if ($o{neck} eq 'arc') {
    # B: each site by its length along the NYHead back midline profile (from Iz; an electrode off
    # the profile is taken at its projection), put at the same length along the atlas back outline
    # (NYHead neck = subject neck, the subject scale); the head-scale reading beside it
    my $on_arc = sub {
        my ($yz) = @_;
        my ($sn, $off) = project_point($pl, $yz, $pc);
        my $px = point_at_length($P, pdl($sn / $at->{mm} + $sI), $cum)->slice('(0),');
        return ($px, $sn, $off);
    };
    print "\nNYHead sites by their length along the NYHead back midline (B): atlas level with the NYHead neck = the subject's neck (subject scale), and at the head scale (the atlas body)\n";
    printf "  %-8s %9s %8s  %-24s %s\n", 'site', 'length', 'off', 'level (subject scale)', 'level (head scale)';
    my $line = sub {
        my ($n, $yz, $tag) = @_;
        my ($px, $sn, $off) = $on_arc->($yz);
        printf "  %-8s %6.1f mm %5.1f mm  %-24s %s%s\n", $n, $sn, $off, level_of($at, $sn),
            level_of($at, $sn * $T->{s} * $at->{mm}), $tag;
        return $px;
    };
    for my $n (grep { /^Nk\d$|^Iz$|^Exxz$/ } @lab) {
        my $off = (project_point($pl, $pos->slice("($ix{$n}),1:2"), $pc))[1];
        my $tag = abs($pos->at($ix{$n}, 0)) > 5 ? sprintf('   (x %+.0f, off the midline: projected, no meaning on the back)', $pos->at($ix{$n}, 0))
                : $off > 30 ? '   (front of the neck: projected, no meaning on the back)' : '';
        push @rows, [ $n, $line->($n, $pos->slice("($ix{$n}),1:2"), $tag) ];
    }
    push @pt, [ $_->[0], $line->($_->[0], $_->[1], '   (--points)') ] for @ptin;
}
else {
my $on_atlas = sub {                                            # NYHead (y,z) -> px, and the place on the atlas back outline
    my ($yz) = @_;
    my $px = from_uv(ny_to_uv($NM, $yz->dummy(0)))->slice('(0),');
    my ($s, $d) = project_point($P, $px, $cum);
    my $L = $s - $sI;
    return ($px, $L / $T->{s}, $L * $at->{mm}, $d);
};
print "\nNYHead sites on the atlas back outline (length from the atlas inion; NYHead mm by the head fit, subject mm by --c7-mm)\n";
printf "  %-8s %10s %10s %8s  %s\n", 'site', 'NYHead mm', 'subject mm', 'off mm', 'atlas level (subject mm)';
for my $n (grep { /^Nk\d$|^Iz$|^Exxz$/ } @lab) {
    my ($px, $Lny, $Lsub, $d) = $on_atlas->($pos->slice("($ix{$n}),1:2"));
    my $lat = abs($pos->at($ix{$n}, 0)) > 5 ? sprintf(' (x %+.0f, projected)', $pos->at($ix{$n}, 0)) : '';
    printf "  %-8s %10.1f %10.1f %8.1f  %s%s\n", $n, $Lny, $Lsub, $d / $T->{s}, ($n =~ /^Nk[2-4]$|^Exxz$/ ? level_of($at, $Lsub) : '-'), $lat;
    push @rows, [ $n, $px ];
}
for my $p (@ptin) {
    my ($px, $Lny, $Lsub, $d) = $on_atlas->($p->[1]);
    printf "  %-8s %10.1f %10.1f %8.1f  %s   (--points)\n", $p->[0], $Lny, $Lsub, $d / $T->{s}, level_of($at, $Lsub);
    push @pt, [ $p->[0], $px ];
}
}

# ---- the atlas taken into NYHead ------------------------------------------------------
my $tsv = "$o{out}_atlas_in_nyhead.tsv";
open my $fo, '>', $tsv or die "$tsv: $!";
printf $fo "# nyhead_atlas_sagittal.pl %s %s: atlas -> NYHead MNI mm by the head fit (%.4f mm/px, turn %+.2f deg) and the neck turn (%+.2f deg about NYHead z %.1f mm, --neck $o{neck}), x = 0\n",
    $fskin, $fbone, 1 / $T->{s}, 57.29578 * $T->{th}, 57.29578 * $NM->{phi}, $NM->{z0};
print $fo join("\t", qw(name x y z kind)), "\n";
for my $l (@{ $at->{levels} }) {
    my $q  = point_at_length($P, pdl($l->{s} + $sI), $cum)->slice('(0),');
    for my $w ([ "atlas_$l->{name}", $q, 'level' ], [ "tip_$l->{name}", $l->{tip}, 'tip' ]) {
        my $yz = uv_to_ny($NM, to_uv($w->[1]->dummy(0)))->slice('(0),');
        printf $fo "%s\t0.00\t%.2f\t%.2f\t%s\n", $w->[0], $yz->list, $w->[2];
    }
}
my $keep = which(($cum - $sI >= 0) & (($cum - $sI) / $T->{s} <= 300));
my $Pyz  = uv_to_ny($NM, to_uv($P->dice_axis(0, $keep)));
printf $fo "atlasback%04d\t0.00\t%.2f\t%.2f\tprofile_atlas\n", $_, $Pyz->slice("($_),")->list for 0 .. $Pyz->dim(0) - 1;
close $fo;
print "\nwrote $tsv\n";

plot_overlay();
exit 0;

# =====================================================================================
# atlas pixels (x right = back, y down) <-> a right-handed frame like NYHead (y forward, z up)
sub to_uv   { my ($p) = @_; return -$p }                        # (u, v) = (-x, -y)
sub from_uv { my ($q) = @_; return -$q }
sub apply_sim {                                                 # (n,2) NYHead (y,z) -> (n,2) atlas (u,v)
    my ($T, $p) = @_;
    my ($c, $s) = (cos($T->{th}), sin($T->{th}));
    my ($a, $b) = ($p->slice(':,(0)'), $p->slice(':,(1)'));
    return cat($T->{s} * ($c * $a - $s * $b), $T->{s} * ($s * $a + $c * $b)) + $T->{t}->dummy(0);
}
sub inv_sim {                                                   # (n,2) atlas (u,v) -> NYHead (y,z)
    my ($T, $q) = @_;
    my $d = ($q - $T->{t}->dummy(0)) / $T->{s};
    my ($c, $s) = (cos($T->{th}), sin($T->{th}));
    my ($a, $b) = ($d->slice(':,(0)'), $d->slice(':,(1)'));
    return cat($c * $a + $s * $b, -$s * $a + $c * $b);
}

# neck_rms($NM, $phi, $k, $N (n,2) NYHead, $P (m,2) atlas back outline px, $cum): RMS px from the
# turned and scaled neck points to the outline
sub neck_rms {
    my ($NM, $phi, $k, $N, $P, $cum) = @_;
    my $m = from_uv(ny_to_uv({ %$NM, phi => $phi, k => $k }, $N));
    my $s2 = 0;
    for my $i (0 .. $m->dim(0) - 1) { my (undef, $d) = project_point($P, $m->slice("($i),"), $cum); $s2 += $d * $d }
    return sqrt($s2 / $m->dim(0));
}
# principal_dir($pts (n,2)[, $px]): the unit direction of least-squares line through the points,
# pointing down (z decreasing for NYHead (y,z); with $px, image y increasing)
sub principal_dir {
    my ($p, $px) = @_;
    my $c = $p - $p->average->dummy(0);
    my ($sxx, $syy, $sxy) = (($c->slice(':,(0)')**2)->sum, ($c->slice(':,(1)')**2)->sum, ($c->slice(':,(0)') * $c->slice(':,(1)'))->sum);
    my $phi = 0.5 * atan2(2 * $sxy, $sxx - $syy);
    my $d = pdl(cos($phi), sin($phi));
    $d = -$d if ($px ? $d->at(1) < 0 : $d->at(1) > 0);
    return $d;
}

# fit_neck(..., $stretch): the angle (and with $stretch the scale) of least neck_rms: a grid of
# 2 deg over -44 .. 44 (and 0.05 over 0.50 .. 1.50), then 0.2 deg (and 0.005) around the best
sub fit_neck {
    my ($NM, $N, $P, $cum, $stretch) = @_;
    my @best = (0, 1); my $bc = neck_rms($NM, 0, 1, $N, $P, $cum);
    my $try = sub { my ($g, $k) = @_; my $c = neck_rms($NM, $g, $k, $N, $P, $cum); (@best = ($g, $k), $bc = $c) if $c < $bc };
    my @ks = $stretch ? map { 0.5 + 0.05 * $_ } 0 .. 20 : (1);
    for my $k (@ks) { $try->($_ / 57.29578, $k) for map { 2 * $_ } -22 .. 22 }
    my ($g0, $k0) = @best;
    my @kf = $stretch ? map { $k0 + 0.005 * $_ } -10 .. 10 : (1);
    for my $k (@kf) { $try->($g0 + $_ / 5 / 57.29578, $k) for -10 .. 10 }
    return @best;
}

# ny_to_uv($NM, $yz (n,2)): NYHead (y,z) -> atlas (u,v) by the head similarity $T, then the neck:
# about $NM->{pivot} turned by phi x w and scaled by 1 + (k - 1) w, w = (z0 - z) / band clipped
# to 0 .. 1 (the head is not moved, the neck below the band turned and scaled fully)
sub ny_to_uv {
    my ($NM, $yz) = @_;
    my $q = apply_sim($T, $yz);
    my $k = $NM->{k} // 1;
    return $q if $NM->{phi} == 0 && $k == 1;
    my $w = (($NM->{z0} - $yz->slice(':,(1)')) / $NM->{band})->clip(0, 1);
    return turn_about($q, $NM->{pivot}, $NM->{phi} * $w, 1 + ($k - 1) * $w);
}
# uv_to_ny($NM, $q (n,2)): the inverse.  The weight w of a point depends on the NYHead height being
# sought: for a trial w the point turned back by phi x w gives a height and so a weight; the
# difference h(w) = weight - w is >= 0 at w = 0 and <= 0 at w = 1, so w is found by bisection
# (all points at once, 50 halvings)
sub uv_to_ny {
    my ($NM, $q) = @_;
    my $k = $NM->{k} // 1;
    return inv_sim($T, $q) if $NM->{phi} == 0 && $k == 1;
    my $back = sub { my ($w) = @_; inv_sim($T, turn_about($q, $NM->{pivot}, -$NM->{phi} * $w, 1 / (1 + ($k - 1) * $w))) };
    my ($lo, $hi) = (zeroes($q->dim(0)), ones($q->dim(0)));
    for (1 .. 50) {
        my $mid = 0.5 * ($lo + $hi);
        my $h = (($NM->{z0} - $back->($mid)->slice(':,(1)')) / $NM->{band})->clip(0, 1) - $mid;
        my $up = ($h > 0);
        $lo = $up * $mid + (1 - $up) * $lo;
        $hi = $up * $hi + (1 - $up) * $mid;
    }
    return $back->(0.5 * ($lo + $hi));
}
# turn_about($q (n,2), $c (2), $a (n) or scalar[, $sc (n) or scalar]): each point turned about $c
# by its angle (and its distance from $c scaled by $sc)
sub turn_about {
    my ($q, $c, $a, $sc) = @_;
    my $d = ($q - $c->dummy(0)) * (defined $sc ? (ref $sc ? $sc->dummy(1) : $sc) : 1);
    my ($ca, $sa) = (cos($a), sin($a));
    my ($x, $y) = ($d->slice(':,(0)'), $d->slice(':,(1)'));
    return cat($ca * $x - $sa * $y, $sa * $x + $ca * $y) + $c->dummy(0);
}

# similarity_to_lines($p (k,2), $q (k,2), $nq (k,2)): the similarity q ~ s R p + t that makes
# the sum of squared distances from the moved p to the lines through q with normals nq least.
# With a = s cos th, b = s sin th the distance n . (sRp + t - q) is linear in (a, b, tx, ty):
# one 4 x 4 least-squares solve.
sub similarity_to_lines {
    my ($p, $q, $nq) = @_;
    my ($px, $py) = ($p->slice(':,(0)'), $p->slice(':,(1)'));
    my ($nx, $ny_) = ($nq->slice(':,(0)'), $nq->slice(':,(1)'));
    my $A = cat($nx * $px + $ny_ * $py, $ny_ * $px - $nx * $py, $nx, $ny_);   # (k,4)
    my $r = ($nq * $q)->mv(1, 0)->sumover;                                     # (k): n . q
    my $AtA = ($A->dummy(1) * $A->dummy(2))->sumover;                          # (k,4,4) -> (4,4): sum over k (dim 0)
    my $Atr = ($A * $r->dummy(1))->sumover;                                    # (k,4) -> (4)
    my $x = (PDL::MatrixOps::inv($AtA) x $Atr->dummy(0))->flat;                # (4): a, b, tx, ty
    my ($a, $b) = ($x->at(0), $x->at(1));
    return { s => sqrt($a**2 + $b**2), th => atan2($b, $a), t => $x->slice('2:3')->copy };
}

# rigid_to_lines($p, $q, $nq, $s, $T): as similarity_to_lines with the scale held at $s: the turn
# and shift by Gauss-Newton steps from $T (residual n . (s R(th) p + t - q), linear in t)
sub rigid_to_lines {
    my ($p, $q, $nq, $s, $T) = @_;
    my ($th, $t) = ($T->{th}, $T->{t}->copy);
    my ($px, $py) = ($p->slice(':,(0)'), $p->slice(':,(1)'));
    my ($nx, $ny_) = ($nq->slice(':,(0)'), $nq->slice(':,(1)'));
    for (1 .. 10) {
        my ($c, $sn) = (cos($th), sin($th));
        my $mx = $s * ($c * $px - $sn * $py) + $t->at(0);
        my $my = $s * ($sn * $px + $c * $py) + $t->at(1);
        my $r  = $nx * ($mx - $q->slice(':,(0)')) + $ny_ * ($my - $q->slice(':,(1)'));          # (k)
        my $dth = $nx * $s * (-$sn * $px - $c * $py) + $ny_ * $s * ($c * $px - $sn * $py);       # dr / dth
        my $J = cat($dth, $nx, $ny_);                                                            # (k,3)
        my $JtJ = ($J->dummy(1) * $J->dummy(2))->sumover;                                        # (3,3)
        my $Jtr = ($J * $r->dummy(1))->sumover;                                                  # (3)
        my $d = (PDL::MatrixOps::inv($JtJ) x $Jtr->dummy(0))->flat;
        $th -= $d->at(0); $t -= $d->slice('1:2');
        last if $d->abs->max < 1e-10;
    }
    return { s => $s, th => $th, t => $t };
}

# outline_normals($q (m,2), $r): unit normals of an unordered outline (pixels of a boundary):
# for every point, the direction of least spread of its neighbours within $r
sub outline_normals {
    my ($q, $r) = @_;
    my $m = $q->dim(0);
    my $N = zeroes($m, 2);
    for my $i (0 .. $m - 1) {
        my $d  = $q - $q->slice("($i),")->dummy(0);                            # (m,2)
        my $w  = which((($d**2)->mv(1, 0)->sumover) <= $r * $r);
        my $dd = $d->dice_axis(0, $w); $dd = $dd - $dd->average->dummy(0);
        my ($cxx, $cyy, $cxy) = (($dd->slice(':,(0)')**2)->sum, ($dd->slice(':,(1)')**2)->sum,
                                 ($dd->slice(':,(0)') * $dd->slice(':,(1)'))->sum);
        my $phi = 0.5 * atan2(2 * $cxy, $cxx - $cyy);                          # the tangent
        $N->slice("($i),") .= pdl(-sin($phi), cos($phi));
    }
    return $N;
}

# nearest($a (n,2), $b (m,2)): for every point of $a the index of the nearest point of $b and
# the distance (one point of $a at a time against all of $b)
sub nearest {
    my ($A, $Bq) = @_;
    my $n = $A->dim(0);
    my $d2 = zeroes($n); my $jj = zeroes(long, $n);
    for my $i (0 .. $n - 1) {
        my $e = (($Bq - $A->slice("($i),")->dummy(0))**2)->mv(1, 0)->sumover;   # (m): over the coordinates (dim 1)
        my $j = $e->minimum_ind->sclr;
        $jj->set($i, $j); $d2->set($i, $e->at($j));
    }
    return ($jj, sqrt($d2));
}

# fit_similarity_icp($p (n,2) NYHead, $q (m,2) atlas uv, $T0): iterated closest points, both
# ways: every NYHead point is paired with its nearest atlas point, and every atlas point (within
# the box of the moved NYHead points) with its nearest moved NYHead point, so that neither set
# can shrink onto a part of the other; in each direction the pairs farther than 2.5 x the median
# distance are left out.  Each step solves the point-to-line similarity (similarity_to_lines,
# the atlas normals), which does not creep along the outline as point-to-point steps do.  The
# RMS reported is the point-to-point distance of the NYHead -> atlas pairs kept.
sub fit_similarity_icp {
    my ($p, $q, $T, %opt) = @_;                                 # fix_s => scale (turn and shift only), cut_front => 1
    my $nq = outline_normals($q, 4);
    my $last = 9e9;
    my ($rms, $nu);
    for my $it (1 .. 100) {
        my $m = apply_sim($T, $p);                              # (n,2)
        my ($jf, $df) = nearest($m, $q);                        # NYHead -> atlas
        # atlas -> NYHead only for the atlas points on the sides where the NYHead set is not cut: it
        # is cut at a height (no atlas point below the lowest moved NYHead point, 3 px allowed) and,
        # with cut_front, halfway to the face (none in front of the most forward one).  Limiting the
        # other sides too lets NYHead shrink and slide (tried: 6.3 deg off and Iz 26 mm off).
        my ($mu, $mv) = ($m->slice(':,(0)'), $m->slice(':,(1)'));
        my ($qu, $qv) = ($q->slice(':,(0)'), $q->slice(':,(1)'));
        my $keep = ($qv >= $mv->min - 3);
        $keep &= ($qu <= $mu->max + 3) if $opt{cut_front};
        my $qa = which($keep);
        my ($jr, $dr) = nearest($q->dice_axis(0, $qa), $m);
        my $uf = which($df <= 2.5 * $df->median + 1e-9);
        my $ur = which($dr <= 2.5 * $dr->median + 1e-9);
        $jr = $jr->index($ur); $ur = $qa->index($ur);
        my $iq = $jf->index($uf)->append($ur);                  # the atlas point of every pair
        my $P2 = $p->dice_axis(0, $uf)->glue(0, $p->dice_axis(0, $jr));
        $T   = defined $opt{fix_s} ? rigid_to_lines($P2, $q->dice_axis(0, $iq), $nq->dice_axis(0, $iq), $opt{fix_s}, $T)
             : similarity_to_lines($P2, $q->dice_axis(0, $iq), $nq->dice_axis(0, $iq));
        $rms = sqrt(($df->index($uf)**2)->avg)->sclr;
        $nu  = $uf->nelem;
        $T->{iter} = $it;
        last if abs($last - $rms) < 1e-5;
        $last = $rms;
    }
    @$T{qw(rms n_used)} = ($rms, $nu);
    return $T;
}

# mesh_cut_x($vc (nv,3), $tri0 (nf,3) 0-based, $x0): the points where the mesh edges cross the
# plane x = x0 (and vertices on it), (n,3), not ordered
sub mesh_cut_x {
    my ($vc, $tri, $x0) = @_;
    my @pts;
    for my $e ([0, 1], [1, 2], [2, 0]) {
        my $va = $vc->dice_axis(0, $tri->slice(":,($e->[0])"));
        my $vb = $vc->dice_axis(0, $tri->slice(":,($e->[1])"));
        my ($xa, $xb) = ($va->slice(':,(0)') - $x0, $vb->slice(':,(0)') - $x0);
        my $w = which($xa * $xb < 0);
        next if $w->isempty;
        my $t = ($xa->index($w) / ($xa->index($w) - $xb->index($w)))->dummy(1);
        push @pts, $va->dice_axis(0, $w) + $t * ($vb->dice_axis(0, $w) - $va->dice_axis(0, $w));
    }
    my $on = which(abs($vc->slice(':,(0)') - $x0) < 1e-9);
    push @pts, $vc->dice_axis(0, $on) unless $on->isempty;
    @pts or die "mesh_cut_x: the plane x = $x0 does not cut the mesh\n";
    my $P = shift @pts; $P = $P->glue(0, $_) for @pts;
    # an edge shared by two triangles gives the same point twice: keep one (to 1e-6 mm)
    my $k = ($P * 1e6)->rint;
    my $ord = $k->slice(':,(2)')->qsorti;                         # sort by z, then drop rows equal to the previous one
    ($P, $k) = ($P->dice_axis(0, $ord), $k->dice_axis(0, $ord));
    my $ord2 = ($k->slice(':,(1)') + 1e9 * $k->slice(':,(2)'))->qsorti;
    ($P, $k) = ($P->dice_axis(0, $ord2), $k->dice_axis(0, $ord2));
    my $same = (($k->slice('1:-1,:') == $k->slice('0:-2,:'))->mv(1, 0)->andover) if $P->dim(0) > 1;
    return $P if !defined $same;
    return $P->dice_axis(0, pdl(long, 0)->append(which(!$same) + 1));
}

sub plot_overlay {
    if (!eval { require PDL::Graphics::Cairo; 1 }) { print "PDL::Graphics::Cairo not found: no PNG\n"; return }
    my ($nx, $ny) = ($skin->dim(1), $skin->dim(2));
    my $fig = PDL::Graphics::Cairo::figure(width => 2 * $nx + 160, height => 2 * $ny + 100);
    my $ax  = $fig->axes;
    $ax->imshow($skin->reorder(2, 1, 0));
    my $Y  = sub { $ny - $_[0] };
    my $px = sub { from_uv(ny_to_uv($NM, $_[0])) };              # (n,2) NYHead (y,z) -> (n,2) px, neck turned
    # NYHead /sa/head at x = 0: with a neck turn only the head (above the pivot height): the turn by
    # height would also swing the face and jaw below it, which belong to the head
    my $hp = $NM->{phi} == 0 ? $C2 : $C2->dice_axis(0, which($C2->slice(':,(1)') >= $NM->{z0}));
    my $cp = $px->($hp);
    $ax->scatter($cp->slice(':,(0)'), $Y->($cp->slice(':,(1)')), color => 'lime', s => 3);
    # the back midline: the mesh part (green); the extension of --points (orange) is not drawn with
    # --neck arc, whose sites go by length along the atlas outline
    for my $k (['mesh', 'green'], $o{neck} eq 'arc' ? () : (['extended', 'orange'])) {
        my @q = map { $_->[1] } grep { $_->[0] eq $k->[0] } @prof;
        next if @q < 2;
        my $qq = $px->(cat(@q)->transpose);                      # cat of (2)-vectors is (2,n): points to dim 0
        $ax->line($qq->slice(':,(0)'), $Y->($qq->slice(':,(1)')), color => $k->[1], lw => 2);
    }
    # labels: placed after all marks; a label that would overlap one already placed goes one line
    # down (then further).  Sizes in image pixels: a line is 9 px, a character 3.4 px (11 pt on a
    # figure of twice the image size)
    my @placed;                                                  # [x0, x1, ytop, ybot] image px
    my $label = sub {
        my ($x, $y, $txt, $col, $fs, %a) = @_;                   # $y image px (down); ha right: text ends at x
        my $w  = 3.4 * length($txt) * $fs / 11; my $hgt = 9 * $fs / 11;
        my ($x0, $x1) = ($a{ha} // 'left') eq 'right' ? ($x - $w, $x) : ($x, $x + $w);
        my $dy = 0;
        for my $try (0 .. 8) {
            my ($t, $b) = ($y + $dy - $hgt, $y + $dy);
            last unless grep { $_->[0] < $x1 && $_->[1] > $x0 && $_->[2] < $b && $_->[3] > $t } @placed;
            $dy += $hgt;
        }
        push @placed, [ $x0, $x1, $y + $dy - $hgt, $y + $dy ];
        $ax->text($x, $Y->($y + $dy), $txt, color => $col, fontsize => $fs, (defined $a{ha} ? (ha => $a{ha}) : ()));
    };
    my $mid = which(abs($pos->slice(':,(0)')) < 3);              # NYHead midline electrodes
    my $ep  = $px->($pos->dice_axis(0, $mid)->slice(':,1:2'));
    $ax->scatter($ep->slice(':,(0)'), $Y->($ep->slice(':,(1)')), color => 'white', s => 2);
    my @labels;                                                  # [x, y, text, colour, size, ha]
    my %side = (Fpz => ['right', 0], Fz => ['right', -5]);       # left of the point (Fz also up): not over the brain
    for my $n (grep { exists $ix{$_} } qw(Fpz Fz Cz Pz Oz)) {
        my $q = $px->($pos->slice("($ix{$n}),1:2")->dummy(0))->slice('(0),');
        my ($ha, $dy) = @{ $side{$n} // ['left', 0] };
        push @labels, [ $q->at(0), $q->at(1) + $dy, $ha eq 'right' ? "$n " : " $n", 'white', 10, $ha ];
    }
    my %rp = map { $_->[0] => $_->[1] } @rows;
    my @rn = map { $_->[0] } @rows;
    if (exists $rp{Nk3} && exists $rp{Nk4} && sqrt((($rp{Nk3} - $rp{Nk4})**2)->sum) < 2) {   # symmetric: the same place
        @rn = grep { $_ ne 'Nk4' } @rn;
    }
    for my $n (@rn) {
        my $q = $rp{$n};
        $ax->scatter(pdl($q->at(0)), pdl($Y->($q->at(1))), color => 'yellow', s => 5);
        push @labels, [ $q->at(0), $q->at(1), $n eq 'Nk3' && !grep({ $_ eq 'Nk4' } @rn) && exists $rp{Nk4} ? ' Nk3/Nk4' : " $n", 'yellow', 11, 'left' ];
    }
    for my $l (@{ $at->{levels} }) {
        my ($x, $y) = $l->{tip}->list;
        $ax->scatter(pdl($x), pdl($Y->($y)), color => 'red', s => 4);
        push @labels, [ $x, $y, "$l->{name} ", 'red', 10, 'right' ];
    }
    for my $kv (grep { $_ ne '' } split /,/, $o{mark}) {         # own sites by length on the atlas
        my ($n, $L) = split /=/, $kv;
        my $q = point_at_length($P, pdl($L / $at->{mm} + $sI), $cum)->slice('(0),');
        $ax->scatter(pdl($q->at(0)), pdl($Y->($q->at(1))), color => 'gold', s => 5);
        push @labels, [ $q->at(0), $q->at(1), " $n $L", 'gold', 11, 'left' ];
    }
    for my $p (@pt) {                                            # own sites as placed on NYHead
        my ($n, $q) = @$p;
        $ax->scatter(pdl($q->at(0)), pdl($Y->($q->at(1))), color => 'magenta', s => 5);
        push @labels, [ $q->at(0), $q->at(1), " $n (NYHead)", 'magenta', 10, 'left' ];
    }
    my @gq = map { $_->[1] } grep { $_->[0] eq 'mesh' } @prof;   # what the green line is: right of it, 2-3 lines below Iz
    if (@gq >= 2) {
        my $e = $px->($gq[0]->dummy(0))->slice('(0),');
        push @labels, [ $e->at(0) + 10, $e->at(1) + 22, ' NYHead /sa/head back', 'green', 9, 'left' ];
    }
    $ax->scatter(pdl($pv->at(0)), pdl($Y->($pv->at(1))), color => 'cyan', s => 6);
    push @labels, [ $pv->at(0), $pv->at(1), sprintf(' neck pivot %+.1f deg', 57.29578 * $NM->{phi}), 'cyan', 10, 'left' ];
    if ($o{labels}) { $label->($_->[0], $_->[1], $_->[2], $_->[3], $_->[4], ha => $_->[5]) for @labels }
    $ax->title(sprintf 'NYHead x = 0 (lime) on the atlas, %.3f mm/px by the head, neck %+.1f deg (--neck %s); subject C7 = %g mm: %.3f mm/px',
        1 / $T->{s}, 57.29578 * $NM->{phi}, $o{neck}, $o{c7_mm}, $at->{mm});
    my $f = "$o{out}_sagittal.png";
    $fig->save($f);
    print "wrote $f\n";
}

# =====================================================================================
# a head-and-neck outline (an ellipse on a neck) as the "atlas" boundary, and as "NYHead" the
# same curve sampled in mm and moved by a known similarity; the fit must find it again
sub selftest {
    my $ok = 1;
    my $chk = sub { my ($c, $m) = @_; printf "  [%s] %s\n", $c ? 'PASS' : 'FAIL', $m; $ok &&= $c };
    my $PI = 3.14159265358979;
    # NYHead-like outline in mm: an egg (longer at the back), y forward, z up, 720 points
    my $a  = sequence(720) * 2 * $PI / 720;
    my $ry = 95 + 10 * (cos($a) < 0) * (-cos($a))                 # the back half a little longer
           + 18 * exp(-(($a - 0.15)**2) / (2 * 0.06**2));         # and a nose, which fixes the turn as on a real head
    my $Pn = cat($ry * cos($a), 75 * sin($a) + 10 * sin(2 * $a) * (sin($a) > 0));   # (720,2)
    my $Ttrue = { s => 1.65, th => 7 / 57.29578, t => pdl(-320, -180) };
    my $Q  = apply_sim($Ttrue, $Pn);                              # atlas (u,v)
    # the atlas outline as integer pixels (what body_outline gives), duplicates removed
    my $Bp = from_uv($Q)->rint;
    my $key = $Bp->slice(':,(0)') * 10000 + $Bp->slice(':,(1)');
    $Bp = $Bp->dice_axis(0, $key->uniqind);
    my $T0 = { s => 1.4, th => 0, t => apply_sim($Ttrue, pdl(0, 0)->dummy(0))->slice('(0),') + pdl(12, -9) };
    my $use = which($Pn->slice(':,(1)') > -20);                   # the part above "Iz"
    # the atlas pixels above the "inion" + 30 rows, as --fit-rows does
    my $yI  = from_uv(apply_sim($Ttrue, pdl(0, -20)->dummy(0)))->at(0, 1);
    my $Bf  = $Bp->dice_axis(0, which($Bp->slice(':,(1)') <= $yI + 30));
    my $T  = fit_similarity_icp($Pn->dice_axis(0, $use), to_uv($Bf), $T0);
    $chk->(abs($T->{s} / 1.65 - 1) < 0.01, sprintf 'scale %.4f px/mm (true 1.65), from a start of 1.4', $T->{s});
    $chk->(abs(57.29578 * $T->{th} - 7) < 0.5, sprintf 'turn %+.2f deg (true +7), from 0', 57.29578 * $T->{th});
    $chk->($T->{rms} < 0.6, sprintf 'RMS %.2f px (pixel rounding only)', $T->{rms});
    # the back half only, the scale held at the true one, from a start 4 deg and 10 px off
    my $ub = which(($Pn->slice(':,(1)') > -20) & ($Pn->slice(':,(0)') < 0));
    my $Tb = fit_similarity_icp($Pn->dice_axis(0, $ub), to_uv($Bf), { s => 1.65, th => 3 / 57.29578, t => $Ttrue->{t} + pdl(8, -6) },
                                fix_s => 1.65, cut_front => 1);
    my $eb = sqrt(((apply_sim($Tb, $Pn) - $Q)**2)->mv(1, 0)->sumover)->max;
    $chk->(abs(57.29578 * $Tb->{th} - 7) < 0.5 && $eb < 2, sprintf 'back half, scale held: turn %+.2f deg (true +7), whole outline within %.2f px', 57.29578 * $Tb->{th}, $eb);
    my $back = inv_sim($T, apply_sim($T, $Pn));
    $chk->(($back - $Pn)->abs->max < 1e-9, 'inv_sim undoes apply_sim');
    my $e = sqrt(((apply_sim($T, $Pn) - $Q)**2)->mv(1, 0)->sumover)->max;
    $chk->($e < 1.5, sprintf 'the whole outline (also the part not fitted) within %.2f px of the truth', $e);
    # similarity_to_lines on exact pairs: any normals give the similarity back
    my $T2 = similarity_to_lines($Pn, $Q, outline_normals($Q, 4));
    $chk->(abs($T2->{s} - 1.65) < 1e-9 && abs($T2->{th} - $Ttrue->{th}) < 1e-9 && ($T2->{t} - $Ttrue->{t})->abs->max < 1e-7,
        'similarity_to_lines: exact pairs give the similarity back');
    my $Nq = outline_normals($Q, 4);                              # normals against the exact ones of the curve
    my $tg = $Q->slice('2:-1,:') - $Q->slice('0:-3,:');           # central differences along the ordered curve
    my $cosang = abs(($Nq->slice('1:-2,:') * $tg)->mv(1, 0)->sumover / sqrt(($tg**2)->mv(1, 0)->sumover));
    $chk->($cosang->pct(0.95) < 0.05, sprintf 'outline_normals: perpendicular to the curve (|cos| 95 %% %.3f; max %.3f at the sharpest bend, the nose)',
        $cosang->pct(0.95), $cosang->max);
    # the neck turn: inverse, and a known angle found again
    local $main::T = { s => 1.2, th => 0.05, t => pdl(-300, -200) };
    my $NMt = { pivot => apply_sim($main::T, pdl(0, -60)->dummy(0))->slice('(0),'), z0 => -60, band => 40, phi => 12 / 57.29578, k => 0.75 };
    my $G = cat(map { $_->flat } (xvals(21, 21) * 8 - 80), (yvals(21, 21) * 8 - 140))->copy;   # (441,2) y, z
    my $rt = uv_to_ny($NMt, ny_to_uv($NMt, $G));
    $chk->(($rt - $G)->abs->max < 1e-6, sprintf 'uv_to_ny undoes ny_to_uv (max %.1e mm)', ($rt - $G)->abs->max);
    my $nl = cat(-100 - 0.2 * (sequence(60) * 2), -40 - sequence(60) * 2)->copy;        # a neck back line, z -40 .. -158
    my $Pa = from_uv(ny_to_uv($NMt, cat(-100 - 0.2 * (sequence(300) * 0.5), -40 - sequence(300) * 0.5)));  # the "atlas" outline
    my $Nn = $nl->dice_axis(0, which($nl->slice(':,(1)') < -60));
    my ($ph, $kk) = fit_neck({ %$NMt, phi => 0, k => 1 }, $Nn, $Pa, arc_length($Pa), 1);
    $chk->(abs(57.29578 * $ph - 12) < 0.3 && abs($kk - 0.75) < 0.01, sprintf 'neck angle and scale found again: %+.2f deg x %.3f (true +12 x 0.75)', 57.29578 * $ph, $kk);
    # principal_dir: a noisy line 10 deg forward of straight down, both frames
    my $tt = sequence(50);
    my $ln = cat(10 + $tt * sin(10 / 57.29578) + 0.3 * sin($tt * 1.7), -$tt * cos(10 / 57.29578));   # (y, z), going down
    my $pd = principal_dir($ln);
    $chk->(abs(57.29578 * atan2($pd->at(0), -$pd->at(1)) - 10) < 0.5 && $pd->at(1) < 0,
        sprintf 'principal_dir: %+.2f deg from straight down (true +10), pointing down', 57.29578 * atan2($pd->at(0), -$pd->at(1)));
    my $pp = principal_dir(-$ln->slice('-1:0,:'), 1);              # as image pixels (y down), listed bottom first
    $chk->($pp->at(1) > 0, 'principal_dir with $px: points down the image');
    # mesh_cut_x on a tetrahedron-like strip: a square from (-1,0,0)..(1,0,0)..(1,0,1)..(-1,0,1)
    my $vc  = pdl([[-1, 0, 0], [1, 0, 0], [1, 0, 1], [-1, 0, 1]])->transpose;   # (4,3): point, coordinate
    my $tri = pdl(long, [[0, 1, 2], [0, 2, 3]])->transpose;                      # (2,3)
    my $c   = mesh_cut_x($vc, $tri, 0);
    my $zs  = $c->slice(':,(2)')->qsort;
    $chk->($c->dim(0) == 3 && abs($zs->at(0)) < 1e-12 && abs($zs->at(1) - 0.5) < 1e-12 && abs($zs->at(2) - 1) < 1e-12,
        sprintf 'mesh_cut_x: a square cut at x = 0 gives z = %s', join(' ', map { sprintf '%.2f', $_ } $zs->list));
    print $ok ? "selftest passed\n" : "selftest FAILED\n";
    return $ok;
}
