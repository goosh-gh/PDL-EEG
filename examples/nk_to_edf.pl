#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use PDL;
use PDL::EEG::IO::NihonKohden qw(read_nk block_ranges select_block select_range clock_to_samp);
use PDL::EEG::IO::EDF         qw(write_edf);

# nk_to_edf.pl - convert a Nihon Kohden .EEG file to EDF / EDF+
#
#   perl nk_to_edf.pl INPUT.EEG [OUTPUT.edf] [options]
#
#   --plain            write plain EDF (no annotation channel; drops events)
#   --phys auto|gain|N physical scaling (default auto)
#   --allblocks        concatenate all waveform blocks into one EDF
#   --blocks SPEC      write ONE EDF per selected block (0-based). SPEC is a
#                      comma/range list, e.g. "0", "2,3", "1-4", "0,2-3", "all".
#                      Output: <base>_bNN.edf per block. Blocks come from real
#                      waveform blocks (wfmblock) or, for a single-block file
#                      (extblock/EEG-1200A), from .LOG "REC START" segments
#                      (placed by epoch, so gaps are handled correctly).
#   --cut SPEC         write ONE EDF per range in DATA-coordinate seconds:
#                      "a-b[:name],c-d[:name]" e.g. "0-300,300-590:task2".
#   --cut-clock SPEC   same, but ranges are WALL-CLOCK "HH:MM:SS" (as shown in
#                      the vendor viewer), e.g. "14:06:14-14:07:15:task2".
#                      Converted to data samples via epoch anchors.
#   --subject  STR     local subject identification (alias: --patient)
#   --recording STR    local recording identification
#   --equipment STR    acquisition equipment string (default: vendor + device,
#                      e.g. Nihon_Kohden_EEG-1200A_V01.00)
#   --rename SPEC      rename channels before writing, for recordings where an input
#                      named after a 10-20 site carried an electrode elsewhere.
#                      SPEC is "OLD=NEW,OLD=NEW" or a file: a TSV whose first two
#                      columns are OLD and NEW (more columns, e.g. the site, are
#                      ignored; '#' starts a comment, a header line "recorded ..." is
#                      skipped). Every OLD must be in the file, and no NEW may equal
#                      a label that stays. Applies to every mode.
#                      e.g. --rename YJ0394VB_rename.tsv

my %opt = (phys => 'auto', plain => 0, allblocks => 0, subject => '', recording => '',
           equipment => undef, gapsamples => 100, rename => undef);
GetOptions(\%opt, 'plain', 'phys=s', 'allblocks', 'blocks=s', 'cut=s', 'cutclock|cut-clock=s',
                  'subject|patient=s', 'recording=s', 'equipment=s', 'gapsamples=i', 'rename=s')
    or die "bad options\n";
my %RENAME = defined $opt{rename} ? parse_rename($opt{rename}) : ();

my $in  = shift @ARGV or die "usage: $0 INPUT.EEG [OUTPUT.edf] [--options]\n";
my $out = shift @ARGV;
unless (defined $out) { ($out = $in) =~ s/\.[^.]+$//; $out .= '.edf'; }

# subject default: input file basename without extension (e.g. subject)
unless (length $opt{subject}) {
    (my $stem = $in) =~ s{.*/}{};      # strip directory
    $stem =~ s/\.[^.]+$//;             # strip extension
    $opt{subject} = $stem;
}

my $phys = ($opt{phys} =~ /^-?\d/) ? 0 + $opt{phys} : $opt{phys};
sub write_one {
    my ($rec, $path) = @_;
    write_edf($rec, $path,
        plus => $opt{plain} ? 0 : 1, phys => $phys,
        subject => $opt{subject}, recording => $opt{recording},
        equipment => $opt{equipment});
    printf "wrote %s (%s): %d ch x %d samp\n",
        $path, ($opt{plain} ? 'EDF' : 'EDF+C'), $rec->{data}->dim(0), $rec->{data}->dim(1);
}

# --rename: "OLD=NEW,..." or a TSV file (first two columns OLD, NEW)
sub parse_rename {
    my $spec = shift;
    my @pairs;
    if (-f $spec) {
        open my $fh, '<', $spec or die "$spec: $!\n";
        while (my $l = <$fh>) {
            $l =~ s/#.*//; $l =~ s/\s+$//;
            next unless $l =~ /\S/;
            my ($old, $new) = map { s/^\s+|\s+$//gr } (split /\t/, $l)[0, 1];
            next if lc $old eq 'recorded';                     # header line
            die "$spec: line '$l' has no second column (columns are tab-separated)\n"
                unless defined $new && length $new;
            push @pairs, [ $old, $new ];
        }
    } else {
        for my $p (split /\s*,\s*/, $spec) {
            my ($old, $new) = split /=/, $p, 2;
            die "--rename: '$p' is not OLD=NEW\n" unless defined $new && length $old && length $new;
            push @pairs, [ $old, $new ];
        }
    }
    my (%m, %to);
    for (@pairs) {
        my ($old, $new) = @$_;
        die "--rename: $old given twice\n" if exists $m{$old};
        die "--rename: two channels renamed to $new\n" if $to{$new}++;
        length $new <= 16 or die "--rename: '$new' is longer than the 16 characters of an EDF label\n";
        $m{$old} = $new;
    }
    return %m;
}
sub rename_channels {
    my $rec = shift;
    return unless %RENAME;
    my @lab = @{ $rec->{labels} };
    my %have = map { $_ => 1 } @lab;
    for my $old (sort keys %RENAME) {
        die "--rename: no channel '$old' in $in (channels: @lab)\n" unless $have{$old};
    }
    for my $new (values %RENAME) {
        die "--rename: '$new' is already a channel that is not renamed\n"
            if $have{$new} && !exists $RENAME{$new};
    }
    my @newlab = map { exists $RENAME{$_} ? $RENAME{$_} : $_ } @lab;
    $rec->{labels} = \@newlab;
    print "  renamed: ", join(', ', map { "$_ -> $RENAME{$_}" } grep { exists $RENAME{$_} } @lab), "\n";
}

# block spec parser: "2,3" | "1-4" | "0,2-3" | "all" (0-based)
sub parse_blocks {
    my ($spec, $nb) = @_;
    return (0 .. $nb - 1) if lc($spec) eq 'all';
    my %pick;
    for my $tok (split /[,\s]+/, $spec) {
        next unless length $tok;
        if ($tok =~ /^(\d+)-(\d+)$/) {
            my ($a, $b) = ($1, $2); ($a, $b) = ($b, $a) if $a > $b;
            $pick{$_} = 1 for $a .. $b;
        } elsif ($tok =~ /^(\d+)$/) {
            $pick{$1} = 1;
        } else {
            warn "ignoring bad block token '$tok'\n";
        }
    }
    return grep { $_ >= 0 && $_ < $nb } sort { $a <=> $b } keys %pick;
}

# ----- cut mode: arbitrary ranges (data-coordinate or wall-clock) ------------
if (defined $opt{cut} || defined $opt{cutclock}) {
    my $rec = read_nk($in, all_blocks => 1, gap_samples => 0);
    rename_channels($rec);
    my $fs  = $rec->{fs};
    my $n   = $rec->{data}->dim(1);
    printf "read %s : %s [%s], %d ch @ %g Hz, %.1f s data\n",
        $in, ($rec->{device}//'?'), ($rec->{layout}//'?'), $rec->{data}->dim(0), $fs, $n/$fs;

    my $is_clock = defined $opt{cutclock};
    my $spec     = $is_clock ? $opt{cutclock} : $opt{cut};

    # start-of-recording wall-clock seconds (for absolute HH:MM:SS input)
    my $start_clk = 0;
    if ($is_clock && ($rec->{t_start} // '') =~ /(\d{2}):(\d{2}):(\d{2})/) {
        $start_clk = $1 * 3600 + $2 * 60 + $3;
    }
    my $to_samp = sub {
        my $tok = shift;
        if ($tok =~ /^(\d+):(\d+):(\d+)$/) {            # HH:MM:SS (wall-clock)
            my $wall = ($1*3600 + $2*60 + $3) - $start_clk;
            return clock_to_samp($rec, $wall);
        }
        return $is_clock ? clock_to_samp($rec, 0 + $tok)   # bare seconds = wall-clock
                         : int((0 + $tok) * $fs + 0.5);    # bare seconds = data coord
    };

    (my $base = $out) =~ s/\.edf$//i;
    my $k = 0;
    for my $part (split /\s*,\s*/, $spec) {
        next unless length $part;
        my ($range, $name) = split /:/, $part, 2;
        my ($a, $b) = split /-/, $range, 2;
        unless (defined $a && defined $b) { warn "skip bad range '$part'\n"; next; }
        my ($lo, $hi) = ($to_samp->($a), $to_samp->($b));
        ($lo, $hi) = ($hi, $lo) if $lo > $hi;
        my $sub  = select_range($rec, $lo, $hi);
        my $file = defined $name && length $name
                 ? sprintf('%s_%s.edf', $base, $name)
                 : sprintf('%s_cut%02d.edf', $base, $k);
        printf "  cut %s: samples %d..%d  (%.3f-%.3f s data)\n",
            ($name // "#$k"), $lo, $hi, $lo/$fs, $hi/$fs;
        write_one($sub, $file);
        $k++;
    }
    exit 0;
}

# ----- block-per-file mode ---------------------------------------------------
if (defined $opt{blocks}) {
    my $rec = read_nk($in, all_blocks => 1, gap_samples => 0);  # clean boundaries
    rename_channels($rec);
    printf "read %s : %s [%s], %d ch @ %g Hz\n",
        $in, ($rec->{device}//'?'), ($rec->{layout}//'?'), $rec->{data}->dim(0), $rec->{fs};

    my $ranges = block_ranges($rec);
    my $nb = @$ranges;
    my $fs = $rec->{fs};
    printf "detected %d block(s):\n", $nb;
    for my $r (@$ranges) {
        printf "  [%d] %8.3f-%8.3f s  (%6d samp)  start=%s\n",
            $r->{index}, $r->{start}/$fs, $r->{end}/$fs, $r->{end}-$r->{start},
            ($r->{t_start}//'?');
    }

    my @idx = parse_blocks($opt{blocks}, $nb);
    die "no valid block indices in '$opt{blocks}' (have 0..".($nb-1).")\n" unless @idx;

    (my $base = $out) =~ s/\.edf$//i;
    for my $i (@idx) {
        my $sub = select_block($rec, $i);
        write_one($sub, sprintf('%s_b%02d.edf', $base, $i));
    }
    exit 0;
}

# ----- single-EDF mode (default / --allblocks) -------------------------------
my $rec = read_nk($in,
    $opt{allblocks} ? (all_blocks => 1, gap_samples => $opt{gapsamples}) : ());
rename_channels($rec);

printf "read %s\n  format : %s [%s]\n  signal : %d ch x %d samp @ %g Hz\n  events : %d\n",
    $in, ($rec->{device} // '?'), ($rec->{layout} // '?'),
    $rec->{data}->dim(0), $rec->{data}->dim(1), $rec->{fs},
    ($rec->{events} ? scalar @{ $rec->{events} } : 0);

if ($opt{allblocks} && $rec->{t_block_starts} && @{ $rec->{t_block_starts} } > 1) {
    my $fs = $rec->{fs};
    printf "  blocks : %d, starts(s) = %s\n",
        scalar @{ $rec->{t_block_starts} },
        join(', ', map { sprintf('%.3f', $_ / $fs) } @{ $rec->{t_block_starts} });
}

write_one($rec, $out);
