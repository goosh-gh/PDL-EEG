#!/usr/bin/env perl
# examples/read_nihonkohden.pl
#
# Usage:
#   perl -Ilib examples/read_nihonkohden.pl subject.eeg [--plot] [--block N]
#              [--sec S] [--nch N] [--uv U] [--chans LIST] [--aux MODE]
#              [--cut a-b[:name]] [--cut-clock HH:MM:SS-HH:MM:SS[:name]] [--no-lazy]
#
# Reads either Nihon Kohden layout (wfmblock / extblock) -- dispatch is entirely
# inside read_nk(), so this viewer needs no format knowledge.
#
# READ STRATEGY
#   default        : whole recording (all blocks concatenated). read_nk() reports
#                    event times in whole-recording coordinates, so loading only
#                    block 0 would leave the data and the events in different
#                    coordinate systems.
#   --cut/-clock   : LAZY. Only the blocks needed to cover the requested range are
#                    read; a range that straddles a block boundary pulls in every
#                    block it touches and they are glued together. Blocks past the
#                    end of the range are never read. Peak memory is "blocks in
#                    range + 1 transient", not the whole recording.
#   --no-lazy      : read the whole recording even with --cut (old behaviour).
#   --block N      : read ONE waveform block (wfmblock/EEG-1100C). --cut is then
#                    interpreted in THAT BLOCK's coordinates, like nk_to_mul.pl.
#
# Display options:
#   --sec S     : seconds per screen (default 10, also under --cut; the horizontal
#                 slider scrolls this window. Pass --sec <cut width> to see the
#                 whole cut at once -- but then the slider has nothing to scroll.)
#   --nch N     : number of channels to show (default 8; ignored if --chans)
#   --uv U      : initial µV per channel division for EEG traces (default 100;
#                 sets the gain slider's starting position, 10..500)
#   --chans L   : comma-separated channel NAMES to show, in order
#                 e.g. --chans Fp1,Cz,Pz,DC01,DC02,STIM
#   --aux MODE  : how to scale "aux" channels (DC*/STIM/PAD/COM/$*/BN*/Pulse/CO2...)
#                   same  : draw as-is, same µV scale as EEG (default)
#                   auto  : auto-scale each aux channel to fit its own slot.
#                           The scale then DIFFERS PER CHANNEL, so each aux row's
#                           label gets "(x)" appended: the real excursion that one
#                           calibration-bar height represents on that channel.
#                           Read it. Auto-scaling blows a dead channel's noise up
#                           to full slot height, which looks identical to a live
#                           trigger until you see that its bar means 50 uV and the
#                           neighbouring one means 3 V.
#                   <N>   : fixed N µV/div for aux channels (e.g. --aux 2000)
#
# RANGE OPTIONS -- SAME SPEC AS nk_to_mul.pl / nk_to_edf.pl, so a range you eyeball
# here can be pasted straight into the converters:
#   --cut       "a-b[:name]"                  DATA-COORDINATE seconds
#   --cut-clock "HH:MM:SS-HH:MM:SS[:name]"    WALL-CLOCK (vendor viewer times)
# ':name' has no output file here, so it is used as a tag in the figure title.
# A comma-separated list is accepted for spec compatibility but only the FIRST
# range is shown (one window, one range).
#
# The time axis is in whole-recording data-coordinate seconds (NOT reset to 0 by
# --cut), so the numbers you read off the plot are exactly what you feed back to
# --cut / nk_to_mul.pl.
#
# EEG channels always use the slider-controlled µV/div. Only aux channels are
# affected by --aux.

use strict;
use warnings;
use PDL;
use PDL::EEG::IO::NihonKohden qw(read_nk block_extents);
use Getopt::Long;
use POSIX ();

binmode(STDOUT, ':encoding(UTF-8)');   # .LOG annotations may be Japanese
binmode(STDERR, ':encoding(UTF-8)');

my ($plot, $allblk, $blk, $cut, $cutclock);
my $lazy     = 1;          # --no-lazy to disable
my $gap_samples = 0;       # butt-join blocks; no synthetic samples in the data
my $nsec     = undef;      # seconds per screen
my $nch      = 8;
my $uv_init  = 100;        # initial µV/div (EEG) -> initial gain-slider position
my $chans    = '';
my $aux_mode = 'same';     # same | auto | <µV/div>

GetOptions(
    'plot'                 => \$plot,
    'block=i'              => \$blk,
    'allblocks!'           => \$allblk,
    'lazy!'                => \$lazy,
    'gap-samples=i'        => \$gap_samples,
    'sec=f'                => \$nsec,
    'nch=i'                => \$nch,
    'uv=f'                 => \$uv_init,
    'chans=s'              => \$chans,
    'aux=s'                => \$aux_mode,
    'cut=s'                => \$cut,
    'cutclock|cut-clock=s' => \$cutclock,
) or die "bad options\n";

my $file = shift @ARGV
    or die "Usage: $0 subject.eeg [--plot] [--block N] [--sec S] "
         . "[--nch N] [--uv U] [--chans LIST] [--aux MODE]\n"
         . "          [--cut a-b[:name]] [--cut-clock HH:MM:SS-HH:MM:SS[:name]] [--no-lazy]\n"
         . "  (default: whole recording; with --cut only the needed blocks are read)\n";

die "--block and --allblocks are mutually exclusive\n"
    if defined $blk && $allblk;
die "--cut and --cut-clock are mutually exclusive\n"
    if defined $cut && defined $cutclock;

# ---------------------------------------------------------------------------
# Parse the range SPEC now (pure syntax; fs / t_start are not needed yet).
# ---------------------------------------------------------------------------
my $T = qr/\d{1,2}:\d{2}:\d{2}|\d+(?:\.\d+)?/;   # HH:MM:SS or (fractional) seconds
my ($cut_a, $cut_b, $cut_name, $is_clock);

if (defined $cut || defined $cutclock) {
    $is_clock = defined $cutclock ? 1 : 0;
    my $spec  = $is_clock ? $cutclock : $cut;
    my @parts = grep { length } split /\s*,\s*/, $spec;
    warn "only the first range is displayed (" . scalar(@parts) . " given)\n"
        if @parts > 1;
    my $part = $parts[0] // '';
    $part =~ /^\s*($T)\s*-\s*($T)(?::(.*))?$/
        or die "bad range '$part' (expected a-b[:name]"
             . ($is_clock ? " with HH:MM:SS" : " in seconds") . ")\n";
    ($cut_a, $cut_b, $cut_name) = ($1, $2, $3);
    undef $cut_name unless defined $cut_name && length $cut_name;
}
my $have_cut = defined $cut_a;

# "YYYY-MM-DD HH:MM:SS" -> epoch seconds (undef if unparseable)
sub ts_epoch {
    my $ts = shift // '';
    my ($Y,$M,$D,$h,$m,$sec) = $ts =~ /^(\d{4})-(\d{2})-(\d{2})[ T](\d{2}):(\d{2}):(\d{2})/
        or return undef;
    my $e = POSIX::mktime($sec, $m, $h, $D, $M - 1, $Y - 1900, 0, 0, -1);
    return (defined $e && $e >= 0) ? $e : undef;
}

# Build the break list from block extents (header-only, cheap).
sub make_bounds {
    my ($ext, $fs) = @_;
    my @b;
    for my $i (0 .. $#$ext - 1) {
        my ($a, $z) = ($ext->[$i], $ext->[$i + 1]);
        my ($ea, $ez) = (ts_epoch($a->{t_start}), ts_epoch($z->{t_start}));
        # dt = wall-clock time that elapsed with NOTHING recorded, i.e. how much
        # real time the vertical line skips over:
        #     (start of block b) - (start of block a) - (duration of block a)
        my $dt = (defined $ea && defined $ez)
               ? ($ez - $ea) - $a->{n_samp} / $fs
               : undef;
        push @b, { samp => $z->{start_samp},           # == $a->{end_samp}
                   prev => $a->{index}, next => $z->{index}, dt => $dt,
                   kind => 'block' };
    }
    return @b;
}
# extblock (EEG-1200A) files are ONE physical data block, so block_extents finds
# no boundaries -- yet the recording is still cut into segments, and the data are
# gap-REMOVED continuous. The only record of the breaks is the .LOG: each segment
# opens with a "REC START" event carrying both a wall-clock t and an epoch-derived
# data-sample position (samp). Inside a segment wall-clock and data advance 1:1,
# so the dead time accumulated before a segment is
#     dead = t - samp/fs
# and the gap AT a break is the jump in that quantity. (This mirrors what
# block_ranges(source => 'log') does, plus the elapsed time it does not report.)
sub log_bounds {
    my ($events, $fs) = @_;
    my @a = grep { defined $_->{samp} && defined $_->{t} }
            grep { ($_->{label} // '') =~ /REC\s*START/i } @{ $events // [] };
    return () unless @a;

    my %seen;
    @a = grep { !$seen{ $_->{samp} }++ }
         sort { $a->{samp} <=> $b->{samp} || $a->{t} <=> $b->{t} } @a;

    my (@out, $seg);
    $seg = 0;
    for my $e (@a) {
        push @out, { samp => $e->{samp}, t => $e->{t},
                     prev => $seg - 1, next => $seg, kind => 'seg' }
            if $e->{samp} > 0 && $seg > 0;
        $seg++;
    }
    return @out;
}

# The .LOG gives only an ESTIMATE of where a break sits in the data.
# _attach_epoch_samp() maps epoch -> sample with  samp = (epoch-1) * n_total/max_epoch,
# i.e. it assumes every segment is the same length. They are not (this recording has
# segments of 296 s and of 30 s), so the estimate drifts by seconds -- enough to draw
# the line in the wrong place, and enough to make a computed gap come out NEGATIVE.
#
# The break itself is unambiguous in the data: gap-removed concatenation makes every
# channel jump at once. So search a window around each estimate for the largest
# simultaneous jump and snap the line to it. Only non-aux (EEG) channels are used --
# DC/trigger channels have large legitimate steps (TTL edges) that would win.
sub snap_bounds {
    my ($data, $abs0, $bounds, $fs, $labels, $n_ch_valid, $radius_s) = @_;
    my $n = $data->dim(1);
    my @eeg = grep { $_ < $n_ch_valid && !is_aux($labels->[$_]) } 0 .. $data->dim(0) - 1;
    return unless @eeg;
    my $r = int(($radius_s // 8) * $fs);

    for my $b (@$bounds) {
        my $loc = $b->{samp} - $abs0;                     # estimate, local samples
        my $w0  = $loc - $r;  $w0 = 0      if $w0 < 0;
        my $w1  = $loc + $r;  $w1 = $n - 1 if $w1 > $n - 1;
        next if $w1 - $w0 < 3;                            # not in the loaded window
        next if $loc < -$r || $loc > $n + $r;

        my $seg  = $data->slice("," . $w0 . ":" . $w1)->slice(pdl(\@eeg), 'X');
        my $diff = ($seg->slice(",1:-1") - $seg->slice(",0:-2"))->abs->sumover;
        my $med  = $diff->medover->sclr;
        my $mx   = $diff->max->sclr;
        my $i    = $diff->maximum_ind->sclr;

        # A real break dwarfs ordinary sample-to-sample movement. If it does not,
        # the estimate is all we have -- say so rather than pretending.
        if ($med > 0 && $mx > 20 * $med) {
            $b->{samp}    = $abs0 + $w0 + $i + 1;         # first sample of new segment
            $b->{snapped} = 1;
        }
    }

    # Recompute the elapsed gap from the (now exact) sample positions:
    #   dt = wall-clock elapsed - data elapsed
    # .LOG times are whole seconds, so dt carries about +/-1 s of quantisation.
    for my $i (0 .. $#$bounds) {
        my $b = $bounds->[$i];
        my ($t0, $s0) = $i == 0 ? (0, 0)
                      : ($bounds->[$i-1]{t}, $bounds->[$i-1]{samp});
        next unless defined $b->{t} && defined $t0;
        $b->{dt} = ($b->{t} - $t0) - ($b->{samp} - $s0) / $fs;
    }
    return;
}

# {data} from read_nk is uniformly MICROVOLTS -- every channel, DC included.
#
# Do NOT scale by $rec->{units} here. {units} is the dimension a channel should be
# EXPORTED in (write_edf writes DC as mV because +/-12002900 uV does not fit EDF's
# 8-character physical_min field); it is NOT the unit the data is in. Treating it
# as the latter multiplies every DC row by 1000 and reports a 2 V trigger as
# "2.06e+03 V", which is exactly what this code did until it was caught.
#
# So: microvolts in, a readable unit out.
sub fmt_scale {
    my $v = shift;
    return '?' unless defined $v && $v > 0;
    # NB: double quotes -- '\x{B5}' in single quotes is the literal 7 characters
    return sprintf("%.3g \x{B5}V", $v)        if $v < 1_000;
    return sprintf("%.3g mV", $v / 1_000)     if $v < 1_000_000;
    return sprintf("%.3g V",  $v / 1_000_000);
}

sub fmt_dt {
    my $dt = shift;
    return '?' unless defined $dt;
    return sprintf('%.1fs', $dt)              if $dt <  60;
    return sprintf('%dm%02ds', int($dt / 60), int($dt) % 60) if $dt < 3600;
    return sprintf('%dh%02dm', int($dt / 3600), int($dt / 60) % 60);
}

sub start_clock_sec {
    my $ts = shift // '';
    return $ts =~ /(\d{2}):(\d{2}):(\d{2})/ ? $1 * 3600 + $2 * 60 + $3 : 0;
}
# One token -> sample index. Identical arithmetic to nk_to_mul.pl's $to_samp.
sub tok_to_samp {
    my ($tok, $fs, $start_clk) = @_;
    if ($tok =~ /^(\d+):(\d+):(\d+)$/) {                 # HH:MM:SS (wall-clock)
        return int((($1 * 3600 + $2 * 60 + $3) - $start_clk) * $fs + 0.5);
    }
    return int((0 + $tok) * $fs + 0.5);                  # bare seconds
}

# ---------------------------------------------------------------------------
# Read.
#
# $abs0 = absolute sample index (whole-recording concatenated coordinates) of
#         local sample 0 of $data. The plot's time axis is ($abs0 + i) / $fs.
# ---------------------------------------------------------------------------
my ($rec, $data, $abs0, $read_note);
# Recording breaks. Blocks are butt-joined (gap_samples => 0), so a break is a
# zero-width discontinuity at an absolute sample index -- not a stretch of fake
# samples. Each entry also carries the block indices either side and the REAL
# elapsed wall-clock gap, recovered from the per-block t_start:
#     dt = epoch(t_start[b+1]) - epoch(t_start[b]) - n_samp[b]/fs
my @bounds;      # [ { samp, prev, next, dt }, ... ]  ABSOLUTE samples
my @blocks;      # block extents, for "which block am I looking at?"

my $lazy_cut = 0;                               # 1 = the buffer IS already the cut

if (defined $blk) {
    # ---- single block: --cut is in THIS block's coordinates (nk_to_mul style) -
    $rec       = read_nk($file, block => $blk);
    $data      = $rec->{data};
    $abs0      = 0;
    $read_note = sprintf("block %d/%d only", $blk, $rec->{n_blocks} - 1);
    @blocks    = ({ index => $blk, kind => 'block', start_samp => 0,
                    end_samp => $data->dim(1), n_samp => $data->dim(1),
                    t_start => $rec->{t_start} });
    # no breaks inside a single block
}
elsif (!$have_cut || !$lazy) {
    # ---- whole recording ----------------------------------------------------
    $rec       = read_nk($file, all_blocks => 1, gap_samples => $gap_samples);
    $data      = $rec->{data};
    $abs0      = 0;
    @blocks    = map { { %$_, kind => 'block' } }
                 @{ block_extents($file, gap_samples => $gap_samples) };
    @bounds    = make_bounds(\@blocks, $rec->{fs});
    $read_note = sprintf("all %d blocks concat, %d break%s",
                         $rec->{n_blocks}, scalar @bounds,
                         (@bounds == 1 ? '' : 's'));
}
else {
    # ---- LAZY: read only the blocks the range touches ------------------------
    # block_extents() reads the control-block address table and the per-block
    # headers ONLY (no samples), so we can decide which blocks to read without
    # touching the waveform data. Its start_samp/end_samp are in exactly the
    # coordinates read_nk(all_blocks => 1) produces, INCLUDING the gap_samples
    # zero padding between blocks -- so the same padding has to be re-inserted
    # here, or --cut coordinates would drift 100 samples per crossed boundary
    # and would no longer match nk_to_mul.pl.
    my $ext = block_extents($file, gap_samples => $gap_samples);
    my $nb  = scalar @$ext;
    my $fsz = $ext->[0]{fs};
    my $end = $ext->[-1]{end_samp};

    my $sc = 0;
    if ($is_clock) {
        $sc = start_clock_sec($ext->[0]{t_start});
    }
    my ($lo, $hi) = ( tok_to_samp($cut_a, $fsz, $sc), tok_to_samp($cut_b, $fsz, $sc) );
    ($lo, $hi) = ($hi, $lo) if $lo > $hi;
    $lo = 0 if $lo < 0;
    die "cut range is empty\n" if $hi - $lo < 1;
    die sprintf("cut range %s-%s is outside the recording (0.000-%.3f s)\n",
                $cut_a, $cut_b, $end / $fsz)
        if $lo >= $end;
    warn sprintf("cut clipped to the recording: end %.3f s requested, "
               . "%.3f s available\n", $hi / $fsz, $end / $fsz)
        if $hi > $end;
    $hi = $end if $hi > $end;

    my $nch  = $ext->[0]{n_ch};
    my @need = grep { $_->{end_samp} > $lo && $_->{start_samp} < $hi } @$ext;

    # A cut can start or end INSIDE a gap (or even lie entirely within one). The
    # all_blocks buffer has zeros there, so pad with zeros to keep the lazy read
    # sample-for-sample identical to read_nk(all_blocks => 1) sliced to [lo,hi).
    my (@piece, $first_abs);
    if (!@need) {
        $first_abs = $lo;
        @piece     = (zeroes(float, $nch, $hi - $lo));
    }
    else {
        for my $i (0 .. $#need) {
            if ($i > 0) {                                # inter-block gap padding
                my $gap = $need[$i]{start_samp} - $need[$i - 1]{end_samp};
                push @piece, zeroes(float, $nch, $gap) if $gap > 0;
            }
            my $r = read_nk($file, block => $need[$i]{index});
            $rec //= { %$r };                            # metadata from first block
            push @piece, $r->{data};
            undef $r;
        }
        $first_abs = $need[0]{start_samp};
        if ($lo < $first_abs) {                          # starts inside a gap
            unshift @piece, zeroes(float, $nch, $first_abs - $lo);
            $first_abs = $lo;
        }
        if ($hi > $need[-1]{end_samp}) {                 # ends inside a gap
            push @piece, zeroes(float, $nch, $hi - $need[-1]{end_samp});
        }
    }
    $rec //= read_nk($file, block => 0);                 # gap-only cut: need labels
    delete $rec->{data};

    my $nkeep = scalar @need;
    my $d = shift @piece;
    $d = $d->glue(1, @piece) if @piece;                  # straddling blocks joined here
    @piece = ();

    my $lo_rel = $lo - $first_abs;  $lo_rel = 0          if $lo_rel < 0;
    my $hi_rel = $hi - $first_abs;  $hi_rel = $d->dim(1) if $hi_rel > $d->dim(1);
    $d = $d->slice(":," . $lo_rel . ":" . ($hi_rel - 1))->sever;

    $data      = $d;
    $abs0      = $first_abs + $lo_rel;
    $lazy_cut  = 1;
    @blocks    = map { { %$_, kind => 'block' } } @$ext;
    @bounds    = make_bounds($ext, $fsz);
    my $n_in   = grep { $_->{samp} > $abs0 && $_->{samp} < $abs0 + $data->dim(1) } @bounds;
    $read_note = sprintf("lazy: %d of %d blocks read (%d break%s in view)",
                         $nkeep, $nb, $n_in, ($n_in == 1 ? '' : 's'));
}

my $fs         = $rec->{fs};
my @labels     = @{ $rec->{labels} };
my @events     = @{ $rec->{events} // [] };
my $n_ch_valid = $rec->{n_ch_valid};
my $n_total    = $data->dim(1);                  # LOCAL length
my $duration   = $n_total / $fs;

# $rec->{t_start} is the start of the first LOADED segment, not of the recording
# (read_nk(block => N) rebases it, matching select_block). Report both.
my $t_rec_start = @blocks ? $blocks[0]{t_start} : $rec->{t_start};
printf "Device : %s\n", $rec->{device};
printf "Start  : %s  (recording)%s\n", $t_rec_start,
    ($rec->{t_start} ne $t_rec_start ? "   loaded from $rec->{t_start}" : '');
# One physical block (extblock / a single-block wfmblock) means block_extents
# found no boundaries -- but the recording can still be segmented. Fall back to
# the .LOG "REC START" markers, which is the only record of the breaks there.
if (!@bounds && @events) {
    @bounds = log_bounds(\@events, $fs);
    if (@bounds) {
        # Rebuild the "which region am I in?" list from the .LOG segments, so the
        # bottom-left label still works when no break is on screen.
        my @st = (0, map { $_->{samp} } @bounds);
        @blocks = map { { index => $_, kind => 'seg', start_samp => $st[$_],
                          end_samp => ($_ < $#st ? $st[$_ + 1] : 9**18) } } 0 .. $#st;
        # The .LOG positions are estimates (uniform-epoch assumption); snap them
        # to the real discontinuity in the loaded data.
        snap_bounds($data, $abs0, \@bounds, $fs, \@labels, $n_ch_valid, 8);
        my $ns = grep { $_->{snapped} } @bounds;
        $read_note .= sprintf(", %d .LOG segment break%s (%d snapped to the data)",
                              scalar @bounds, (@bounds == 1 ? '' : 's'), $ns);
    }
}

printf "Loaded : %d ch (%d valid), %d samples = %.1f s @ %g Hz  [%s]\n",
    $data->dim(0), $n_ch_valid, $n_total, $duration, $fs, $read_note;
printf "Events : %s\n",
    join(', ', map { "$_->{label}(\@$_->{t}s)" } @events) || 'none';
printf "Breaks : %s\n",
    @bounds
      ? join(', ', map { sprintf("%s %d @ %.1fs%s (%s skipped)",
                                 ($_->{kind} eq 'seg' ? 'seg' : 'block'),
                                 $_->{next}, $_->{samp} / $fs,
                                 ($_->{kind} eq 'seg'
                                    ? ($_->{snapped} ? '' : '~') : ''),
                                 fmt_dt($_->{dt})) }
                   @bounds)
      : 'none';

# ---------------------------------------------------------------------------
# Scrollable region (LOCAL, end-exclusive). Under lazy loading the buffer already
# IS the cut. Under --no-lazy / --block the cut still has to be applied here.
# ---------------------------------------------------------------------------
my ($r_lo, $r_hi) = (0, $n_total);

if ($have_cut && !$lazy_cut) {
    my $sc = start_clock_sec($rec->{t_start});
    my ($lo, $hi) = ( tok_to_samp($cut_a, $fs, $sc), tok_to_samp($cut_b, $fs, $sc) );
    ($lo, $hi) = ($hi, $lo) if $lo > $hi;
    my ($lo_req, $hi_req) = ($lo, $hi);
    $lo = 0        if $lo < 0;
    $hi = $n_total if $hi > $n_total;
    die sprintf("cut range %s-%s is outside the loaded data (0.000-%.3f s)%s\n",
                $cut_a, $cut_b, $duration,
                defined $blk ? " -- --cut is in block ${blk}'s coordinates" : '')
        if $lo >= $n_total || $hi - $lo < 1;
    warn sprintf("cut clipped: %.3f-%.3f s requested, %.3f-%.3f s available\n",
                 $lo_req / $fs, $hi_req / $fs, $lo / $fs, $hi / $fs)
        if $lo_req < 0 || $hi_req > $n_total;
    ($r_lo, $r_hi) = ($lo, $hi);
}

printf "Window : %.3f - %.3f s (data coords)%s\n",
    ($r_lo + $abs0) / $fs, ($r_hi + $abs0) / $fs,
    (defined $cut_name ? "  '$cut_name'" : '');

my $n_region   = $r_hi - $r_lo;
my $sec_region = $n_region / $fs;

# Window width is ALWAYS 10 s by default -- also under --cut. Defaulting to the
# whole cut width would make n_show == n_region, i.e. max_offset == 0, and the
# time slider would be dead. Pass --sec <cut width> to see it all at once.
$nsec //= 10;
$nsec = $sec_region if $nsec > $sec_region;

# --- aux channel classification (by label) ---
sub is_aux {
    my $l = shift // '';
    return 1 if $l =~ /^DC/i;
    return 1 if $l eq 'STIM' || $l eq 'PAD' || $l eq 'COM';
    return 1 if $l =~ /^\$/;                  # reference channels $A1, $Cz, ...
    return 1 if $l =~ /^BN/i;
    return 1 if $l =~ /^(Pulse|CO2|SpO2|EtCO2)/i;
    return 0;
}

# --- resolve which channel indices to display ---
my @sel;
if (length $chans) {
    my %by_name;
    for my $i (0 .. $#labels) { $by_name{lc $labels[$i]} //= $i }
    for my $nm (split /\s*,\s*/, $chans) {
        $nm =~ s/^\s+|\s+$//g;                # tolerate "--chans Fp1, Pz, DC03"
        next unless length $nm;
        my $i = $by_name{lc $nm};
        if (defined $i) { push @sel, $i }
        else {
            # comma-separated, so the list can be pasted straight back into --chans
            warn "channel '$nm' not found. Available:\n  --chans \""
               . join(',', @labels) . "\"\n";
        }
    }
    die "no valid channels selected\n" unless @sel;
} else {
    my $np = ($n_ch_valid < $nch) ? $n_ch_valid : $nch;
    @sel = (0 .. $np - 1);
}
my $n_plot     = scalar @sel;
my @sel_labels = @labels[@sel];

exit unless $plot;

# --- Load Cairo once, outside render loop ---
require PDL::Graphics::Cairo;
PDL::Graphics::Cairo->import(qw(subplots));
require PDL::Graphics::Cairo::Driver::GS;

my $gs = PDL::Graphics::Cairo::Driver::GS->new(
    width  => 1100,
    height => 80 + $n_plot * 70,
);

# gain slider: uv = 10 * 50**(1-sv)  =>  sv = 1 - ln(uv/10)/ln(50)
my $sv_init = 1.0 - log(($uv_init > 0 ? $uv_init : 100) / 10) / log(50);
$sv_init = 0 if $sv_init < 0;
$sv_init = 1 if $sv_init > 1;

$gs->show_interactive(
    init   => { 0 => 0.0, 1 => $sv_init },
    render => sub {
        my ($state, $w, $h) = @_;
        # warn sprintf("SLIDER state0=%.4f state1=%.4f\n", $state->{0}//-1, $state->{1}//-1); ### debug for slider position version 0.02
        # --- Slider 0: time offset within [r_lo, r_hi) (LOCAL samples) ---
        my $n_show     = int($nsec * $fs);
        $n_show        = 1         if $n_show < 1;
        $n_show        = $n_region if $n_show > $n_region;
        my $max_offset = ($n_region - $n_show) / $fs;
        my $t_off      = ($state->{0} // 0.0) * $max_offset;
        my $s0         = $r_lo + int($t_off * $fs);
        $s0            = $r_hi - $n_show if $s0 + $n_show > $r_hi;
        $s0            = $r_lo if $s0 < $r_lo;
        my $s1         = $s0 + $n_show - 1;
        $s1            = $r_hi - 1 if $s1 >= $r_hi;
        $s0            = $s1 if $s0 > $s1;

        # --- Slider 1: gain (µV/div) for EEG, logarithmic 10..500 ---
        my $sv         = $state->{1} // 0.5;
        my $uv_per_div = 10 * (50 ** (1.0 - $sv));
        my $spacing    = $uv_per_div * 2.0;

        my ($fig, $ax) = subplots(1, 1, width => $w, height => $h);

        my $plot_w   = $w - 120; $plot_w = 200 if $plot_w < 200;
        my $n_raw    = $s1 - $s0 + 1;
        my $step_max = 5;
        my $step     = int($n_raw / $plot_w) || 1;
        $step        = $step_max if $step > $step_max;
        my $t;
        my @cal;                            # per-row: what one cal-bar height means

        for my $i (0 .. $n_plot - 1) {
            my $ci  = $sel[$i];
            my $raw = $data->slice("($ci),$s0:$s1");
            my $sig = ($step > 1) ? $raw->slice("0:-1:$step") : $raw;
            # ABSOLUTE data-coordinate seconds, so the axis matches --cut
            $t //= (sequence($sig->dim(0)) * $step + $s0 + $abs0) / $fs;

            my $offset = ($n_plot - 1 - $i) * $spacing;
            my $aux    = is_aux($sel_labels[$i]);

            # $cal = the REAL signal excursion that one calibration-bar height
            # (uv_per_div) represents on THIS channel. For EEG that is just
            # uv_per_div, but --aux auto rescales every aux channel to fit its own
            # slot, so the bar means something different on each one -- and a dead
            # channel's noise gets blown up to look exactly like a live signal.
            # Without the number next to the label, the plot is unreadable.
            my ($y, $cal);
            if ($aux && $aux_mode ne 'same') {
                my $mean = $sig->avg;
                if ($aux_mode eq 'auto') {
                    my $amp = ($sig - $mean)->abs->max; $amp = 1 if $amp <= 0;
                    $y    = ($sig - $mean) / $amp * ($spacing * 0.4) + $offset;
                    # spacing = 2*uv_per_div, so spacing*0.4 = 0.8*uv_per_div
                    $cal  = $amp / 0.8;
                } else {
                    my $val = $aux_mode + 0; $val = 1 if $val <= 0;
                    $y    = ($sig - $mean) * ($uv_per_div / $val) + $offset;
                    $cal  = $val;
                }
            } else {
                $y = $sig + $offset;
            }
            $cal[$i] = defined $cal ? $cal : $uv_per_div;   # EEG: the slider gain

            $ax->line($t, $y, color => ($aux ? '#CC0000' : '#1565C0'), lw => 0.7);
            $ax->axhline($offset, color => '#DDDDDD', lw => 0.3);
        }

        my $x_lo = ($s0 + $abs0) / $fs;
        my $x_hi = $t->max;

        # --- recording breaks -------------------------------------------------
        # A break is a zero-width discontinuity, not a stretch of samples. Draw a
        # DASHED orange line (dashed so it cannot be mistaken for an EMG spike or
        # a trace), label the block on EITHER SIDE at the bottom (the trigger/DC
        # rows are the least cluttered place, and a label per side removes the
        # "is N the block that ends or the one that starts?" ambiguity), and put
        # the REAL elapsed wall-clock gap on top -- the thing the old 100-sample
        # zero pad was standing in for.
        # The break annotation lives BELOW the traces, in two rows: the block
        # names either side of the line, and under them the time it skips. It used
        # to sit above the top trace, where it landed on the data.
        my $y_lo   = -$spacing * 0.92;                # bottom of the dashed line
        my $y_hi   = ($n_plot - 0.2) * $spacing;
        my $y_lab  = -$spacing * 0.52;                # row 1: block N-1 | block N
        my $y_dt   = -$spacing * 0.80;                # row 2: "46.0s skipped"
        my $span   = $x_hi - $x_lo;
        my $BRK    = '#FF9800';                       # not used by traces/grid

        my $dashed = sub {                            # no ls=> in P:G:C: draw it
            my ($x, $y0, $y1) = @_;
            my $nseg = 26;
            my $h    = ($y1 - $y0) / $nseg;
            for my $k (0 .. $nseg - 1) {
                next if $k % 2;
                $ax->line(PDL->new([$x, $x]),
                          PDL->new([$y0 + $k * $h, $y0 + ($k + 0.62) * $h]),
                          color => $BRK, lw => 1.4);
            }
        };

        my $n_seen = 0;
        for my $b (@bounds) {
            my $bx = $b->{samp} / $fs;
            next if $bx <= $x_lo || $bx >= $x_hi;
            $n_seen++;
            $dashed->($bx, $y_lo, $y_hi);

            # P:G:C text() now honours halign (P:G:C14): anchor by right/left/
            # center instead of hand-subtracting an estimated text width.
            my $kw   = $b->{kind} eq 'seg' ? 'seg' : 'block';
            my $lft  = sprintf('%s %d', $kw, $b->{prev});
            my $rgt  = sprintf('%s %d', $kw, $b->{next});
            my $dts  = sprintf("\x{25B2} %s skipped", fmt_dt($b->{dt}));
            my $u    = $span / $plot_w;               # data units per pixel
            my $pad  = 10 * $u;                        # clear of the line itself (data coords)

            # left label ends just left of the line; right label starts just right of it.
            $ax->text($bx - $pad, $y_lab, $lft,
                      fontsize => 8, color => $BRK, halign => 'right', valign => 'bottom');
            $ax->text($bx + $pad, $y_lab, $rgt,
                      fontsize => 8, color => $BRK, halign => 'left', valign => 'bottom');
            # centred under the line, on the row below the block names
            $ax->text($bx, $y_dt, $dts,
                      fontsize => 8, color => $BRK, halign => 'center', valign => 'bottom');
        }

        # No break in view? Then the block number is not shown anywhere above, so
        # say which block we are actually in, at the bottom-left.
        if (!$n_seen && @blocks > 1) {
            my $mid = ($s0 + $s1) / 2 + $abs0;
            my ($cur) = grep { $mid >= $_->{start_samp} && $mid < $_->{end_samp} } @blocks;
            $cur //= $blocks[0];
            $ax->text($x_lo + $span * 0.004, $y_lab,
                      sprintf('%s %d', $cur->{kind}, $cur->{index}),
                      fontsize => 8, color => '#9E9E9E',
                      halign => 'left', valign => 'bottom');
        }

        # --- EEG calibration bar in a right gutter (off the data) ---
        my $gutter = ($x_hi - $x_lo) * 0.05 + 1e-9;
        my $bar_t  = $x_hi + $gutter * 0.5;
        $ax->line(PDL->new([$bar_t, $bar_t]), PDL->new([0, $uv_per_div]),
                  color => '#CC0000', lw => 2.0);

        # ONE tick per channel, with a two-line label: the name, and under it the
        # scale that one calibration-bar height actually represents on THAT row.
        # (Two separate ticks would put a tick mark at the scale row, where there
        # is nothing.) Every row gets a number because --aux auto normalises each
        # aux row to its own slot: a dead channel's noise is stretched to full
        # height and looks exactly like a live trigger, until you see that its bar
        # means 16 mV and its neighbour's means 4 V.
        #
        # But two lines only fit if there is room. At 34 channels each row is
        # ~15 px and the second line lands on the next channel's name. Rather than
        # shrinking the font until nothing is legible, work out what fits:
        #
        #   >= 30 px/row : two lines everywhere
        #   >= 18 px/row : two lines on AUX rows only -- the EEG gain is already
        #                  in the suptitle, so nothing is lost by dropping it
        #    < 18 px/row : names only
        my $row_px = ($h - 120) / ($n_plot || 1);       # ~120 px of margins
        my $lines  = $row_px >= 30 ? 'all' : $row_px >= 18 ? 'aux' : 'none';

        my @tick_vals = map { ($n_plot - 1 - $_) * $spacing } 0 .. $n_plot - 1;
        my @tick_labels = map {
            my $want = $lines eq 'all' ? 1
                     : $lines eq 'aux' ? is_aux($sel_labels[$_])
                     : 0;
            $want ? $sel_labels[$_] . "\n(" . fmt_scale($cal[$_]) . ')'
                  : $sel_labels[$_];
        } 0 .. $n_plot - 1;
        $ax->yticks(PDL->new(\@tick_vals), [@tick_labels]);
        $ax->xlim($x_lo, $x_hi + $gutter);
        $ax->ylim(-$spacing * 1.0, ($n_plot - 0.2) * $spacing);
        $ax->xlabel(sprintf("Time (s, data coords)  [%.1f-%.1fs  step=%d]",
            $x_lo, $x_hi, $step));

        $fig->suptitle(
            sprintf("%s  %s  %s  %gHz  | EEG %.0fµV/div  aux=%s  %gs/screen%s",
                $rec->{device},
                (defined $blk ? "blk=$blk" : 'all blocks'),
                $t_rec_start, $fs, $uv_per_div, $aux_mode, $nsec,
                ($have_cut
                    ? sprintf("  cut=%.1f-%.1fs%s",
                              ($r_lo + $abs0) / $fs, ($r_hi + $abs0) / $fs,
                              (defined $cut_name ? " '$cut_name'" : ''))
                    : '')),
            fontsize => 9
        );
        $fig->tight_layout;
        return $fig;
    },
);
