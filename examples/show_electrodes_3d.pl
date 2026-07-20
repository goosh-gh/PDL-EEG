#!/usr/bin/env perl
#
# show_electrodes_3d.pl — 頭皮電極位置(ASA .elc)を 3D 表示する。
#
#   perl examples/show_electrodes_3d.pl                       # 同梱 fixture, GS3D
#   perl examples/show_electrodes_3d.pl --elc standard_1020.elc
#   perl examples/show_electrodes_3d.pl --backend trid        # Advent 原典の TriD
#   perl examples/show_electrodes_3d.pl --backend gs3d        # P:G:C Driver::GS3D
#   perl examples/show_electrodes_3d.pl --fiducials           # LPA/RPA/Nz も描く
#
# 座標は PDL::EEG::IO::ASA::read_elc で読む。coords は (3,N)。
#   TriD : points3d([$x,$y,$z]) + Labels     (Advent 2024 Day 12 と同じ道)
#   GS3D : Driver::GS3D->new(scene=>{points=>[N,3],...})->run
#          正準 seed (0,0,0,1) / R=diag(-1,-1,1) はモジュール側が持つ。
#
# 色は左右で非対称に付ける(対称配色だと鏡映バグが原理的に見えないため):
#   左半球=オレンジ系 / 右半球=シアン系 / 正中=白。
#
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PDL;
use PDL::EEG::IO::ASA qw(read_elc);

my %opt = (backend => 'gs3d', labels => undef, fiducials => 0,
           elc => "$Bin/../t/data/standard_1020_subset.elc");
while (@ARGV) {
    my $a = shift @ARGV;
    if    ($a eq '--elc')         { $opt{elc}       = shift @ARGV }
    elsif ($a eq '--backend')     { $opt{backend}   = lc shift @ARGV }
    elsif ($a eq '--labels')      { $opt{labels}    = 1 }
    elsif ($a eq '--no-labels')   { $opt{labels}    = 0 }
    elsif ($a eq '--fiducials')   { $opt{fiducials} = 1 }
    elsif ($a eq '-h' or $a eq '--help') { exec 'perldoc', $0 if -t STDOUT; die "see header\n" }
    else  { die "unknown option: $a\n" }
}

my $mon = read_elc($opt{elc});
printf "read %d positions from %s (unit=%s, ref=%s)\n",
    $mon->{n}, $mon->{file}, $mon->{unit}, $mon->{reference};

# ---- select electrodes (drop fiducials unless asked) ----------------------
my %is_fid = map { $_ => 1 } keys %{ $mon->{fiducials} };
my @idx = grep { $opt{fiducials} || !$is_fid{ $mon->{labels}[$_] } } 0 .. $mon->{n} - 1;
my @name = map { $mon->{labels}[$_] } @idx;
my $sel  = $mon->{coords}->dice_axis(1, pdl(long, \@idx));   # (3, M)
my $M    = scalar @idx;
die "no electrodes to plot\n" unless $M;

# ---- left/right/midline colouring -----------------------------------------
my @rgb;
for my $j (0 .. $M - 1) {
    my $x = $sel->at(0, $j);
    if    (abs($x) < 5) { push @rgb, [1.0, 1.0,  1.0 ] }   # midline: white
    elsif ($x < 0)      { push @rgb, [0.9, 0.55, 0.25] }   # left:  orange
    else                { push @rgb, [0.3, 0.7,  0.85] }   # right: cyan
}

my $backend = $opt{backend};
if ($backend eq 'gs3d') { show_gs3d($sel, \@name, \@rgb, \%opt) }
elsif ($backend eq 'trid') { show_trid($sel, \@name, \@rgb, \%opt) }
else { die "unknown backend '$backend' (use gs3d|trid)\n" }

# ===========================================================================
# GS3D — PDL::Graphics::Cairo::Driver::GS3D (giza-server 3D viewer)
# ===========================================================================
sub show_gs3d {
    my ($xyz, $names, $rgb, $o) = @_;
    unless (eval { require PDL::Graphics::Cairo::Driver::GS3D; 1 }) {
        die "GS3D backend needs PDL::Graphics::Cairo (Driver::GS3D):\n$@";
    }
    # GS3D wants points/colors as (N,3); read_elc gives (3,N) -> transpose.
    my $points = $xyz->xchg(0, 1)->sever;                 # (M,3)
    my $colors = pdl([ [ map { $_->[0] } @$rgb ],
                       [ map { $_->[1] } @$rgb ],
                       [ map { $_->[2] } @$rgb ] ]);       # (M,3)
    my $show_labels = defined $o->{labels} ? $o->{labels} : 1;   # GS3D: on
    my $drv = PDL::Graphics::Cairo::Driver::GS3D->new(
        scene  => { points => $points, colors => $colors,
                    labels => ($show_labels ? $names : undef) },
        width  => 700, height => 700,
    );
    print "GS3D: drag to rotate, r=canonical, p=perspective/ortho, close to quit\n";
    $drv->run;
}

# ===========================================================================
# TriD — Advent Calendar 2024 Day 12 path (points3d + Labels)
# ===========================================================================
sub show_trid {
    my ($xyz, $names, $rgb, $o) = @_;
    # TriD core: OpenGL::Modern backend on a current (GitHub-master) install.
    eval "use PDL::Graphics::TriD; 1"
        or die "TriD backend needs PDL::Graphics::TriD:\n$@";

    # Labels are a SEPARATE, optional module — and the sole hazard here. Current
    # TriD master dropped Labels.pm entirely; if an old classic-OpenGL Labels.pm
    # lingers on disk, merely loading it drags classic OpenGL into a process that
    # already has OpenGL::Modern, and the two GL XS layers collide (bogus
    # GL_LIGHTING AUTOLOAD kills rendering). So only load it when labels are
    # actually requested — --no-labels must never touch it.
    # TriD labels default OFF: current TriD master ships no Labels.pm, and a
    # stale classic-OpenGL Labels.pm collides with OpenGL::Modern and kills
    # rendering. Only load it if labels are EXPLICITLY requested (--labels).
    my $want_labels = defined $o->{labels} ? $o->{labels} : 0;
    if ($want_labels && !eval "use PDL::Graphics::TriD::Labels; 1") {
        warn "TriD labels unavailable, rendering without them:\n$@";
        $want_labels = 0;
    }

    # coords is (3,M): take the three ROWS as length-M piddles. NB dog() would
    # split along the highest dim (M) and hand back M little 3-vectors instead.
    my $x = $xyz->slice("(0),")->sever;                   # x of all electrodes
    my $y = $xyz->slice("(1),")->sever;                   # y
    my $z = $xyz->slice("(2),")->sever;                   # z
    # Colours as three length-M channel piddles [R,G,B] (Advent's form).
    my $r = pdl(map { $_->[0] } @$rgb);
    my $g = pdl(map { $_->[1] } @$rgb);
    my $b = pdl(map { $_->[2] } @$rgb);

    print "TriD: drag to rotate; press 'q' in the window to close\n";
    PDL::Graphics::TriD::points3d([ $x, $y, $z ], [ $r, $g, $b ],
                                  { PointSize => 8 });
    # points3d displays the window and waits for 'q' (TriD's implicit twiddle),
    # exactly as the Advent script does — no explicit twiddle3d() call needed.
    # (The public loop function is twiddle3d(), not twiddle(); calling the
    # latter is what raised "Undefined subroutine ...::twiddle".)
    if ($want_labels) {
        PDL::Graphics::TriD::hold3d();
        my @padded = map { "  $_" } @$names;              # Advent 2-space pad
        PDL::Graphics::TriD::graph_object(
            PDL::Graphics::TriD::Labels->new([ $x, $y, $z ],
                                             { Strings => \@padded }));
        PDL::Graphics::TriD::release3d();
    }
}

__END__

=head1 NAME

show_electrodes_3d.pl - display ASA .elc scalp electrodes in 3D

=head1 SYNOPSIS

  perl examples/show_electrodes_3d.pl [--elc FILE] [--backend gs3d|trid]
                                      [--labels|--no-labels] [--fiducials]

Reads electrode coordinates with L<PDL::EEG::IO::ASA/read_elc> and renders
them, either through the giza-server 3D viewer (L<PDL::Graphics::Cairo>'s
C<Driver::GS3D>) or through L<PDL::Graphics::TriD> as in the PDL Advent
Calendar 2024, Day 12. Fiducials (LPA/RPA/Nz) are hidden unless C<--fiducials>.

=cut
