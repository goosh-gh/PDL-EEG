#!/usr/bin/env perl
# ---------------------------------------------------------------------------
# sep_hfo_tfa.pl - SEP high-frequency-oscillation time-frequency map, in PDL
#
# Pure-PDL counterpart of an MNE-Python compute_tfr(method="morlet") pipeline:
# read EEGLAB SEP data, epoch around the stimulus, run a complex Morlet CWT
# (PDL::EEG::TFA), z-score against the pre-stimulus baseline, and render the
# time-frequency heatmap with PDL::Graphics::Cairo.
#
#   # runs anywhere, no data needed (synthetic HFO bursts):
#   perl sep_hfo_tfa.pl --demo -o hfo_demo.png
#
#   # real BIDS/EEGLAB data (one or more runs, contralateral channel auto-picked):
#   perl sep_hfo_tfa.pl \
#       /path/sub-001_task-median_run-03_eeg.set \
#       /path/sub-001_task-median_run-05_eeg.set \
#       -o hfo_C4.png
#
#   # superlet instead of plain Morlet (sharper HFO bursts):
#   perl sep_hfo_tfa.pl --demo --superlet -o hfo_slt.png
#
# Data loading assumes the common BIDS/EEGLAB layout: a continuous <base>.set
# with a companion <base>.fdt (float32, nbchan x nsamples) and BIDS sidecars
# (_eeg.json for SamplingFrequency, _channels.tsv, _events.tsv).  A v7.3
# (HDF5) .set is read via PDL::IO::HDF5 when no .fdt is present.
# ---------------------------------------------------------------------------
use strict;
use warnings;
use PDL;
use PDL::EEG::TFA qw(tfr_morlet tfr_superlet tfr_stat apply_baseline);
use Getopt::Long;
use POSIX qw(ceil);

# ---- options ----
my $out       = 'hfo_tfa.png';
my $demo      = 0;
my $superlet  = 0;
my $sfreq_cli = 0;
my $ch_cli    = '';
my $fmin      = 70;
my $fmax      = 900;
my $fstep     = 5;
my $n_cycles  = 7;
my $tmin      = -0.050;    # epoch start (s)
my $tmax      =  0.120;    # epoch end   (s)
my $artifact  =  4.0;      # stimulus-artifact blank half-width (ms)
my $event_type = '';       # BIDS trial_type filter (empty = all events)
my $ymin;                  # display frequency floor (default = fmin)
my $ymax;                  # display frequency ceiling (default = fmax)
my $ytick     = 50;        # y-axis tick step (Hz)
my $xtick     = 5;         # x-axis tick step (ms)
my $cmap      = '';        # colormap override (default jet; decomp default coolwarm)
my $itc       = 0;         # render inter-trial coherence instead of power
my $stat      = 0;         # render across-trial reliability z instead of power
my $decomp    = 0;         # three-panel total / evoked / induced decomposition
my $sig       = '';        # significance contour levels, e.g. "3,5" (empty=off)
my $excise    = 0;         # in-data artifact interpolation half-width (ms; 0=off)
my $bl_min;                # z-score baseline start (ms; default = tmin)
my $bl_max;                # z-score baseline end   (ms; default = -artifact)

GetOptions(
    'o|out=s'      => \$out,
    'demo'         => \$demo,
    'superlet'     => \$superlet,
    'sfreq=f'      => \$sfreq_cli,
    'channel=s'    => \$ch_cli,
    'fmin=f'       => \$fmin,
    'fmax=f'       => \$fmax,
    'fstep=f'      => \$fstep,
    'ncycles=f'    => \$n_cycles,
    'tmin=f'       => \$tmin,
    'tmax=f'       => \$tmax,
    'artifact=f'   => \$artifact,
    'event-type=s' => \$event_type,
    'ymin=f'       => \$ymin,
    'ymax=f'       => \$ymax,
    'ytick=f'      => \$ytick,
    'xtick=f'      => \$xtick,
    'cmap=s'       => \$cmap,
    'itc'          => \$itc,
    'stat'         => \$stat,
    'decomp'       => \$decomp,
    'sig=s'        => \$sig,
    'excise=f'     => \$excise,
    'bl-min=f'     => \$bl_min,
    'bl-max=f'     => \$bl_max,
) or die "bad options\n";

my @sets = @ARGV;

die "  [error] use only one of --stat / --itc (they are different maps)\n"
    if $stat && $itc;

# ---- assemble target-channel epochs (nepoch, ntime) and time axis ----
my ($epochs, $times, $sfreq, $ch_name);

if ($demo || !@sets) {
    ($epochs, $times, $sfreq, $ch_name) = synth_epochs($tmin, $tmax);
    print "DEMO: synthetic SEP-HFO, $epochs->[0] ... using $ch_name\n" if 0;
    printf "DEMO synthetic: %d epochs, sfreq=%g Hz, channel=%s\n",
           $epochs->dim(0), $sfreq, $ch_name;
} else {
    my @all;
    for my $set (@sets) {
        my $run = load_run($set, sfreq => $sfreq_cli);
        $sfreq = $run->{sfreq};
        my $ci = pick_channel($run, $ch_cli);
        $ch_name = $run->{ch_names}[$ci];
        my $trace = $run->{data}->slice("($ci),:");   # (nsamp)
        my @ev = filter_events($run, $event_type);
        my ($ep, $tt) = epoch_trace($trace, \@ev, $sfreq, $tmin, $tmax, $artifact, $excise);
        push @all, $ep if defined $ep && $ep->dim(0) > 0;
        $times = $tt;
        printf " %-40s -> %d epochs (ch %s)\n",
               (split m{/}, $set)[-1], (defined $ep ? $ep->dim(0) : 0), $ch_name;
    }
    die "no epochs collected\n" unless @all;
    $epochs = _vstack_epochs(@all);   # concat runs along the epoch axis
    printf "TOTAL: %d epochs, sfreq=%g Hz, channel=%s\n",
           $epochs->dim(0), $sfreq, $ch_name;
}

# ---- time-frequency transform ----
if ($fmax >= 0.45 * $sfreq) {                 # keep the grid safely below Nyquist
    my $newmax = int(0.45 * $sfreq / $fstep) * $fstep;
    warn sprintf("  [warn] fmax=%g Hz too close to Nyquist (%g Hz); clamping to %g Hz\n",
                 $fmax, $sfreq / 2, $newmax);
    $fmax = $newmax;
}
my $freqs = $fmin + $fstep * sequence(int(($fmax - $fmin) / $fstep) + 1);
my $bl_lo = defined($bl_min) ? $bl_min / 1000.0 : $tmin;
my $bl_hi = defined($bl_max) ? $bl_max / 1000.0 : -($artifact / 1000.0);
my @sig_levels = length($sig) ? (map { $_ + 0 } split /[,\s]+/, $sig) : ();
my $t0 = time;

# ---- three-panel decomposition: total / evoked / induced ----
if ($decomp) {
    my %r = tfr_morlet($epochs, $sfreq, $freqs, n_cycles => $n_cycles,
                       output => 'all');       # one CWT pass yields all three
    printf "decomposition (total/evoked/induced): %d epochs x %d freqs x %d samples in %.1fs\n",
           $epochs->dim(0), $freqs->nelem, $epochs->dim(1), time - $t0;
    # Normalise all three panels by the SAME, well-conditioned reference: the
    # total-power baseline SD (per frequency).  z-scoring evoked/induced against
    # their OWN baseline is wrong -- the evoked baseline is ~0 with ~0 variance
    # (no pre-stimulus phase locking), so dividing by its SD explodes to
    # meaningless values.  Using the total-power baseline SD keeps every panel
    # in the same, interpretable units (SDs of total baseline fluctuation).
    my $bidx = which(($times >= $bl_lo) & ($times <= $bl_hi));
    my $mu0  = $r{power}->dice_axis(1, $bidx)->xchg(0,1)->avgover;          # (nf)
    my $sd0  = (($r{power}->dice_axis(1,$bidx) - $mu0->dummy(1))**2)
                   ->xchg(0,1)->avgover->sqrt;
    $sd0->where($sd0 <= 0) .= 1;
    my $me   = $r{evoked} ->dice_axis(1, $bidx)->xchg(0,1)->avgover;
    my $mi   = $r{induced}->dice_axis(1, $bidx)->xchg(0,1)->avgover;
    my $zt = ($r{power}   - $mu0->dummy(1)) / $sd0->dummy(1);
    my $ze = ($r{evoked}  - $me ->dummy(1)) / $sd0->dummy(1);
    my $zi = ($r{induced} - $mi ->dummy(1)) / $sd0->dummy(1);
    render_triptych(
        [ ['Total power',              $zt],
          ['Evoked (phase-locked)',    $ze],
          ['Induced (non-phase-locked)', $zi] ],
        $times, $freqs, $ch_name, $artifact, $out,
        n => $epochs->dim(0), sig => \@sig_levels, cmap => ($cmap || undef),
        ymin => $ymin, ymax => $ymax, ytick => $ytick, xtick => $xtick);
    print "wrote $out\n";
    exit 0;
}

my ($disp, $clabel);

if ($stat) {
    my %s = tfr_stat($epochs, $sfreq, $freqs, $times, [$bl_lo, $bl_hi],
                     n_cycles => $n_cycles);
    $disp   = $s{z};
    $clabel = 'Across-trial reliability (z)';
    printf "reliability stat: %d epochs x %d freqs x %d samples in %.1fs\n",
           $epochs->dim(0), $freqs->nelem, $epochs->dim(1), time - $t0;
} elsif ($itc) {
    warn "  [warn] --itc uses the Morlet transform (superlet has no phase/ITC)\n"
        if $superlet;
    $disp   = tfr_morlet($epochs, $sfreq, $freqs, n_cycles => $n_cycles,
                         output => 'itc');
    $clabel = 'Inter-trial coherence';
    printf "morlet ITC: %d epochs x %d freqs x %d samples in %.1fs\n",
           $epochs->dim(0), $freqs->nelem, $epochs->dim(1), time - $t0;
} else {
    my $power = $superlet
        ? tfr_superlet($epochs, $sfreq, $freqs, base_cycles => 3,
                       order_min => 1, order_max => 5)
        : tfr_morlet($epochs, $sfreq, $freqs, n_cycles => $n_cycles);
    printf "%s TFA: %d epochs x %d freqs x %d samples in %.1fs\n",
           ($superlet ? 'superlet' : 'morlet'),
           $epochs->dim(0), $freqs->nelem, $epochs->dim(1), time - $t0;
    $disp   = apply_baseline($power, $times, [$bl_lo, $bl_hi], mode => 'zscore');
    $clabel = 'Power (z vs baseline)';
}

# ---- render ----
render_heatmap($disp, $times, $freqs, $ch_name, $artifact, $out,
               superlet => $superlet, n => $epochs->dim(0), itc => $itc,
               stat => $stat, sig => \@sig_levels, cmap => ($cmap || undef),
               ymin => $ymin, ymax => $ymax, ytick => $ytick, xtick => $xtick,
               clabel => $clabel);
print "wrote $out\n";

# ===========================================================================
# Rendering
# ===========================================================================
sub render_heatmap {
    my ($z, $times, $freqs, $ch, $aw, $outfile, %opt) = @_;
    require PDL::Graphics::Cairo;
    PDL::Graphics::Cairo->import(qw(subplots));

    my $tms = $times * 1000;

    # colour scaling from the true analysis region only: (+aw .. +70 ms),
    # excluding the stimulus artifact and the right edge transient.
    my $vmask = ($tms > $aw) & ($tms < 70);
    my $dv    = $z->dice_axis(1, which($vmask));
    my ($vmin, $vmax) = ($dv->min->sclr, $dv->max->sclr);
    $vmin = 0 if $opt{itc};                       # ITC is in [0,1]
    my $cmap = $opt{cmap} || 'jet';
    if ($cmap =~ /coolwarm|rdbu/i && !$opt{itc}) { # diverging -> symmetric about 0
        my $m = $dv->abs->max->sclr; $m = 1 if $m <= 0;
        ($vmin, $vmax) = (-$m, $m);
    }

    my $xe = _centers_to_edges($tms);
    my $ye = _centers_to_edges($freqs);

    my ($fig, $ax) = subplots(1, 1, figsize => [10, 6]);
    $ax->set_facecolor('#e0e0e0');
    $ax->pcolormesh($xe, $ye, $z, cmap => $cmap, vmin => $vmin, vmax => $vmax);
    # significance contours (|z| >= level): positive solid black, negative dashed white
    if ($opt{sig} && @{ $opt{sig} }) {
        my @pos = grep { $_ > 0 } @{ $opt{sig} };
        $ax->contour($tms, $freqs, $z, levels => [@pos],
                     colors => 'black', lw => 1.1) if @pos;
        $ax->contour($tms, $freqs, $z, levels => [map { -$_ } @pos],
                     colors => 'white', lw => 0.9, linestyles => 'dashed') if @pos;
    }
    # opaque grey blank over the stimulus-artifact window.  (P:G:C pcolormesh
    # has no NaN transparency, so we cover rather than blank the data.)
    $ax->axvspan(-$aw, $aw, color => '#c9c9c9', alpha => 1.0);
    $ax->axvline(0, color => 'black', lw => 1.2);
    my $ylo = defined($opt{ymin}) ? $opt{ymin} : $freqs->min->sclr;
    my $yhi = defined($opt{ymax}) ? $opt{ymax} : $freqs->max->sclr;
    my $ys  = $opt{ytick} || 50;
    my $xs  = $opt{xtick} || 5;
    $ax->xlim(-20, 80);
    $ax->ylim($ylo, $yhi);
    my @yt; for (my $v = $ys * ceil($ylo / $ys); $v <= $yhi + 1e-6; $v += $ys) { push @yt, $v }
    my @xt; for (my $v = -20; $v <= 80 + 1e-6; $v += $xs) { push @xt, $v }
    $ax->xticks(\@xt);
    $ax->yticks(\@yt);
    my $method = $opt{superlet} ? 'superlet' : 'Morlet';
    $method = 'Morlet ITC'         if $opt{itc};
    $method = 'Morlet reliability' if $opt{stat};
    $ax->set_title(sprintf("SEP-HFO %s TFA (N=%d) - %s", $method, $opt{n}, $ch));
    $ax->xlabel('Time (ms)');
    $ax->ylabel('Frequency (Hz)');
    $ax->colorbar(label => ($opt{clabel} // 'Power (z vs baseline)'));
    $fig->tight_layout;
    $fig->save($outfile);
}

sub _centers_to_edges {
    my ($c) = @_;
    my $n = $c->nelem;
    my $e = zeroes($n + 1);
    $e->slice("1:-2") .= 0.5 * ($c->slice("0:-2") + $c->slice("1:-1"));
    $e->set(0,  2 * $c->at(0)      - $e->at(1));
    $e->set($n, 2 * $c->at($n - 1) - $e->at($n - 1));
    return $e;
}

# Three side-by-side heatmaps (total / evoked / induced) in one figure,
# sharing a colour scale so the panels are directly comparable.
sub render_triptych {
    my ($panels, $times, $freqs, $ch, $aw, $outfile, %opt) = @_;
    require PDL::Graphics::Cairo;
    PDL::Graphics::Cairo->import(qw(subplots));

    my $tms = $times * 1000;
    my $xe  = _centers_to_edges($tms);
    my $ye  = _centers_to_edges($freqs);
    my $ylo = defined($opt{ymin}) ? $opt{ymin} : $freqs->min->sclr;
    my $yhi = defined($opt{ymax}) ? $opt{ymax} : $freqs->max->sclr;
    my $ys  = $opt{ytick} || 50;
    my $xs  = $opt{xtick} || 5;
    my @yt; for (my $v = $ys * ceil($ylo / $ys); $v <= $yhi + 1e-6; $v += $ys) { push @yt, $v }
    my @xt; for (my $v = -20; $v <= 80 + 1e-6; $v += $xs) { push @xt, $v }

    # signed z: use a diverging colormap centred at 0 (0 = white/neutral,
    # +red / -blue), NOT jet (whose brightest colour sits at the middle).
    # Scale each panel independently and symmetrically about 0 -- evoked power
    # is much smaller than total/induced, so a shared scale hides it.
    my $cmap  = $opt{cmap} || 'coolwarm';
    my $vmask = ($tms > $aw) & ($tms < 70);
    my $fidx  = which(($freqs >= $ylo) & ($freqs <= $yhi));   # displayed band only

    my ($fig, @ax) = subplots(1, 3, figsize => [16, 4.8]);
    for my $i (0 .. 2) {
        my ($title, $z) = @{ $panels->[$i] };
        my $ax = $ax[$i];
        my $d  = $z->dice_axis(0, $fidx)->dice_axis(1, which($vmask));
        my $m  = $d->abs->max->sclr; $m = 1 if $m <= 0;   # symmetric limit
        $ax->set_facecolor('#e0e0e0');
        $ax->pcolormesh($xe, $ye, $z, cmap => $cmap, vmin => -$m, vmax => $m);
        if ($opt{sig} && @{ $opt{sig} }) {
            my @pos = grep { $_ > 0 } @{ $opt{sig} };
            $ax->contour($tms, $freqs, $z, levels => [@pos],
                         colors => 'black', lw => 1.0) if @pos;
        }
        $ax->axvspan(-$aw, $aw, color => '#c9c9c9', alpha => 1.0);
        $ax->axvline(0, color => 'black', lw => 1.0);
        $ax->xlim(-20, 80);
        $ax->ylim($ylo, $yhi);
        $ax->xticks(\@xt);
        $ax->yticks(\@yt);
        $ax->set_title($title);
        $ax->xlabel('Time (ms)');
        $ax->ylabel('Frequency (Hz)') if $i == 0;
        $ax->colorbar(label => 'z vs total baseline');
    }
    $fig->suptitle(sprintf("SEP-HFO decomposition (N=%d) - %s", $opt{n}, $ch))
        if $fig->can('suptitle');
    $fig->tight_layout;
    $fig->save($outfile);
}

# ===========================================================================
# Epoching
# ===========================================================================
sub epoch_trace {
    my ($trace, $ev, $sfreq, $tmin, $tmax, $aw, $excise) = @_;
    my $s0 = int($tmin * $sfreq + ($tmin < 0 ? -0.5 : 0.5));
    my $s1 = int($tmax * $sfreq + ($tmax < 0 ? -0.5 : 0.5));
    my $len = $s1 - $s0 + 1;
    my $nsamp = $trace->dim(0);
    my $times = ($s0 + sequence($len)) / $sfreq;
    my $bl_hi = -($aw / 1000.0);
    my $bl_mask = ($times < $bl_hi);              # pre-stimulus baseline samples

    my @keep;
    for my $s (@$ev) {
        my $a = $s + $s0;
        my $b = $s + $s1;
        next if $a < 0 || $b >= $nsamp;
        my $seg = $trace->slice("$a:$b")->sever;
        _excise_artifact($seg, $times, $excise) if $excise && $excise > 0;
        push @keep, $seg;
    }
    return (undef, $times) unless @keep;

    my $ep = _vstack(map { $_->dummy(0) } @keep);  # (nepoch, len)
    # per-epoch DC removal on the pre-stimulus window
    if (which($bl_mask)->nelem) {
        my $bl = $ep->dice_axis(1, which($bl_mask))->xchg(0,1)->avgover; # (nepoch)
        $ep = $ep - $bl->dummy(1);
    }
    return ($ep, $times);
}

# Replace the stimulus-artifact window [-ex,+ex] ms with a straight line
# between the samples just outside it (removes the transient from the DATA,
# so the wavelet no longer smears it into neighbouring latencies).
sub _excise_artifact {
    my ($seg, $times, $ex_ms) = @_;
    my $ex = $ex_ms / 1000.0;
    my $lo = which($times <= -$ex);
    my $hi = which($times >=  $ex);
    return unless $lo->nelem && $hi->nelem;
    my $ia = $lo->max->sclr;
    my $ib = $hi->min->sclr;
    return unless $ib > $ia + 1;
    my ($va, $vb) = ($seg->at($ia), $seg->at($ib));
    my $n = $ib - $ia;
    for my $k (1 .. $n - 1) {
        $seg->set($ia + $k, $va + ($vb - $va) * $k / $n);
    }
}

sub _vstack_epochs {
    my @blocks = @_;                   # each (nepoch_i, len)
    my $len   = $blocks[0]->dim(1);
    my $total = 0; $total += $_->dim(0) for @blocks;
    my $out = zeroes(double, $total, $len);
    my $off = 0;
    for my $b (@blocks) {
        my $n = $b->dim(0);
        $out->slice("$off:@{[$off+$n-1]},:") .= $b;
        $off += $n;
    }
    return $out;
}

sub _vstack {
    my @rows = @_;                     # each (1, len) or (len)
    @rows = map { $_->ndims == 1 ? $_->dummy(0) : $_ } @rows;
    my $len = $rows[0]->dim(1);
    my $out = zeroes(double, scalar(@rows), $len);
    $out->slice("($_),:") .= $rows[$_]->slice("(0),:") for 0 .. $#rows;
    return $out;
}

# ===========================================================================
# EEGLAB / BIDS loading
# ===========================================================================
sub load_run {
    my ($set, %opt) = @_;
    -e $set or die "not found: $set\n";
    (my $base = $set) =~ s/\.set$//;
    my $fdt = "$base.fdt";

    my %r;
    # --- metadata from BIDS sidecars (fall back to CLI / HDF5) ---
    (my $stem = $base) =~ s/_eeg$//;
    $r{sfreq} = $opt{sfreq} || _read_json_sfreq("${stem}_eeg.json") || 0;
    my ($names, $types) = _read_channels_tsv("${stem}_channels.tsv");

    if (-e $fdt) {
        die "sfreq unknown (no _eeg.json); pass --sfreq\n" unless $r{sfreq};
        die "no _channels.tsv beside $set\n" unless @$names;
        my $nch   = scalar @$names;
        my $flat  = _read_float32($fdt);
        my $total = $flat->nelem;
        my $nsamp = int($total / $nch);
        warn sprintf("  [warn] %s: %d floats not divisible by %d channels "
                   . "(channels.tsv vs .fdt mismatch?)\n",
                   (split m{/}, $fdt)[-1], $total, $nch) if $total % $nch;
        $r{data}     = $flat->slice("0:@{[$nch*$nsamp-1]}")->reshape($nch, $nsamp);
        $r{ch_names} = $names;
        $r{ch_types} = $types;
        $r{events}   = _read_events_tsv("${stem}_events.tsv", $r{sfreq});
        printf "  loaded %s: %d ch x %d samp = %.1f s @ %g Hz, %d events\n",
               (split m{/}, $fdt)[-1], $nch, $nsamp, $nsamp / $r{sfreq},
               $r{sfreq}, scalar @{ $r{events} };
    } else {
        ($r{data}, $r{ch_names}, $r{ch_types}, $r{sfreq}, $r{events}) =
            _read_set_hdf5($set, $r{sfreq}, $names, $types);
    }
    die "no channels loaded from $set\n" unless @{ $r{ch_names} };
    return \%r;
}

sub _read_float32 {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "open $path: $!\n";
    local $/; my $bytes = <$fh>; close $fh;
    my $n = int(length($bytes) / 4);
    my $p = zeroes(float, $n);
    my $ref = $p->get_dataref;
    $$ref = substr($bytes, 0, $n * 4);
    $p->upd_data;
    return $p->double;
}

sub _read_json_sfreq {
    my ($json) = @_;
    return 0 unless -e $json;
    open my $fh, '<', $json or return 0;
    local $/; my $t = <$fh>; close $fh;
    return $1 if $t =~ /"SamplingFrequency"\s*:\s*([0-9.eE+-]+)/;
    return 0;
}

sub _read_channels_tsv {
    my ($tsv) = @_;
    my (@names, @types);
    return (\@names, \@types) unless -e $tsv;
    open my $fh, '<', $tsv or return (\@names, \@types);
    my $hdr = <$fh>; chomp $hdr; $hdr =~ s/\r$//;
    my @cols = split /\t/, $hdr;
    my %idx  = map { lc($cols[$_]) => $_ } 0 .. $#cols;
    my $ni = $idx{name} // 0;
    my $ti = $idx{type};
    while (my $line = <$fh>) {
        chomp $line; $line =~ s/\r$//;
        my @f = split /\t/, $line;
        push @names, $f[$ni];
        push @types, defined($ti) ? uc($f[$ti] // '') : 'EEG';
    }
    close $fh;
    return (\@names, \@types);
}

sub _read_events_tsv {
    my ($tsv, $sfreq) = @_;
    my @ev;
    return \@ev unless -e $tsv;
    open my $fh, '<', $tsv or return \@ev;
    my $hdr = <$fh>; chomp $hdr; $hdr =~ s/\r$//;
    my @cols = split /\t/, $hdr;
    my %idx  = map { lc($cols[$_]) => $_ } 0 .. $#cols;
    my $oi = $idx{onset} // 0;
    my $tti = $idx{trial_type};
    while (my $line = <$fh>) {
        chomp $line; $line =~ s/\r$//;
        my @f = split /\t/, $line;
        my $s = int($f[$oi] * $sfreq + 0.5);
        push @ev, { sample => $s, type => (defined $tti ? $f[$tti] : '') };
    }
    close $fh;
    return \@ev;
}

# v7.3 (HDF5) .set fallback: read /EEG/data and /EEG/srate
sub _read_set_hdf5 {
    my ($set, $sfreq, $names, $types) = @_;
    eval { require PDL::IO::HDF5; 1 }
        or die "$set has no companion .fdt and PDL::IO::HDF5 is not available.\n"
             . "Re-export from EEGLAB with a separate .fdt, or install PDL::IO::HDF5.\n";
    my $h5 = PDL::IO::HDF5->new($set)
        or die "cannot open $set as HDF5 (is it a v7.3 .set?)\n";
    my $data = $h5->dataset('/EEG/data')->get;        # HDF5 {pnts,nbchan}->PDL(nbchan,pnts) after xchg
    $data = $data->xchg(0,1) if $data->dim(0) > $data->dim(1);  # heuristics; keep (nch,nsamp)
    unless ($sfreq) {
        my $sr = eval { $h5->dataset('/EEG/srate')->get };
        $sfreq = $sr->at(0) if defined $sr && $sr->nelem;
    }
    die "sfreq unknown; pass --sfreq\n" unless $sfreq;
    # events from BIDS sidecar still preferred; empty here
    (my $stem = $set) =~ s/\.set$//; $stem =~ s/_eeg$//;
    my $events = _read_events_tsv("${stem}_events.tsv", $sfreq);
    $names = [ map { "ch$_" } 1 .. $data->dim(0) ] unless @$names;
    return ($data->double, $names, $types, $sfreq, $events);
}

sub filter_events {
    my ($run, $type) = @_;
    my @ev = @{ $run->{events} };
    @ev = grep { $_->{type} eq $type } @ev if length $type;
    return map { $_->{sample} } @ev;
}

sub pick_channel {
    my ($run, $want) = @_;
    my @names = @{ $run->{ch_names} };
    my %pos   = map { uc($names[$_]) => $_ } 0 .. $#names;
    if (length $want) {
        return $pos{uc $want} if exists $pos{uc $want};
        die "channel '$want' not found\n";
    }
    # contralateral median-nerve SEP: prefer right-hemisphere central sites
    for my $c (qw(C4 CP4 C4P C3 CP3 C3P CZ)) {
        return $pos{$c} if exists $pos{$c};
    }
    # first EEG-typed channel
    my @types = @{ $run->{ch_types} };
    for my $i (0 .. $#names) { return $i if ($types[$i] // 'EEG') eq 'EEG' }
    return 0;
}

# ===========================================================================
# Synthetic demo data (phase-locked HFO bursts riding an N20-like deflection)
# ===========================================================================
sub synth_epochs {
    my ($tmin, $tmax) = @_;
    my $PI = 4 * atan2(1, 1);
    my $sfreq = 4000;
    my $ne    = 200;
    my $nt    = int(($tmax - $tmin) * $sfreq) + 1;
    my $t     = $tmin + sequence($nt) / $sfreq;
    srand(7);
    my $ep = zeroes(double, $ne, $nt);
    for my $e (0 .. $ne - 1) {
        my $n20 = 6   * exp(-(($t-0.020)**2)/(2*0.004**2))  * sin(2*$PI*45*($t-0.020));
        my $b1  =       exp(-(($t-0.022)**2)/(2*0.005**2))  * sin(2*$PI*120*$t);
        my $b2  = 0.8 * exp(-(($t-0.020)**2)/(2*0.0025**2)) * sin(2*$PI*275*$t);
        # small phase-locked sigma-burst (~600 Hz on N20): tiny power, strong phase locking
        my $sig = 0.35* exp(-(($t-0.020)**2)/(2*0.0018**2)) * sin(2*$PI*600*$t);
        # non-phase-locked (induced) burst: random phase each trial
        my $ind = 0.9 * exp(-(($t-0.030)**2)/(2*0.006**2))  * sin(2*$PI*200*$t + rand()*2*$PI);
        $ep->slice("($e),:") .= $n20 + $b1 + $b2 + $sig + $ind + grandom($nt) * 2.0;
    }
    return ($ep, $t, $sfreq, 'C4 (synthetic)');
}
