#!/usr/bin/env perl
# blink_ged_clean.pl - GED eye-blink removal for a continuous EEG recording.
#
#   perl blink_ged_clean.pl TASK.edf [VOLUNTARY_BLINKS.edf] [options]
#     --veog LABEL     vEOG channel label            (default: vEOG)
#     --trim-sec S     drop the first S seconds of TASK (default: 0)
#     --n N            remove N blink components       (default: 1)
#     --reg G          shrinkage on R                  (default: 0.05)
#     --out FILE       write cleaned EEG (subtracted); extension picks the format:
#                      .edf -> EDF (PDL::EEG::IO::EDF::write_edf), .mul -> BESA
#                      ASCII (PDL::EEG::IO::BESA::ASCII::write_mul)
#     --out-removed F  write the removed artifact (same extension rules)
#     --exclude L,L..  channels to pass through unchanged (not filtered).
#     --montage-only   write only the filtered channels (subset EDF) instead of
#                      the whole recording.
#
# By DEFAULT every signal channel is filtered -- EOG included -- and the written
# EDF is the WHOLE recording: filtered channels are replaced by their cleaned
# versions and every other channel (DC trigger channels, Event/Marker/annotation
# channels, and any --exclude channel) passes through unchanged. This keeps the
# DC trigger channels, which are needed to epoch/average the task.
#
# DC trigger channels and Event/Marker/annotation channels are NEVER filtered.
# --exclude adds channels (e.g. a noisy aux or ECG) to the pass-through set.
# vEOG is always used for blink detection; --exclude vEOG keeps that role but
# writes its original (unfiltered) signal. --montage-only writes just the
# filtered channels as a subset EDF.
#
# TASK.edf is the recording to clean. If VOLUNTARY_BLINKS.edf is given, the blink
# *pattern* (S) is estimated from it (a dedicated voluntary-blink session, same
# montage); otherwise the pattern is estimated from TASK's own blinks. The
# background R is always taken from TASK's blink-free segments.
#
# Core call:
#   ($subtracted, $removed, $op) = remove_blinks($task, veog=>$i, fit=>$voluntary);
#   $subtracted (nch,nt) = cleaned data;  $removed = $task - $subtracted.
# Both inputs are (nch, nt) piddles with identical channel order/count.

use strict;
use warnings;
use PDL;
use FindBin;
use lib "$FindBin::Bin/../lib";
use PDL::EEG::IO::EDF qw(read_edf);
# clean_edf_label is optional: use your EDF.pm's if it exports one, else a small
# built-in fallback (strips EDF+ type codes like "EEG Fp1-Ref" -> "Fp1").
my $clean_label = PDL::EEG::IO::EDF->can('clean_edf_label') || sub {
    my $l = shift; return $l unless defined $l;
    $l =~ s/^\s+//; $l =~ s/\s+$//;
    $l =~ s/^(?:EEG|EOG|ECG|EKG|EMG|MEG|POL|DC)\s+//;
    $l =~ s/-Ref$//i; $l =~ s/\s+/_/g; return $l;
};

# load_edf: accept either read_edf return style and normalise to ($data,$labels,$fs)
#   * record hashref  { data=>(nch,nt), labels=>[...], fs=>... }  (current EDF.pm)
#   * list ($data, \@labels, $fs)                                (older readers)
sub load_edf {
    my ($f) = @_;
    my @r = read_edf($f);
    if (@r == 1 && ref $r[0] eq 'HASH') {
        my $h = $r[0];
        return ($h->{data}, $h->{labels}, $h->{fs}, $h->{units}, $h->{events}, $h->{t_start});
    }
    return (@r[0,1,2], undef, undef, undef);   # ($data, \@labels, $fs)
}


use PDL::EEG::GED     qw(remove_blinks blink_evoked detect_blinks blink_free);
use Getopt::Long;

my ($veog_label, $trim_sec, $n_remove, $reg, $out_file, $out_removed, $montage_only, $exclude_opt) =
    ('vEOG', 0, 1, 0.05, undef, undef, 0, undef);
GetOptions('veog=s'=>\$veog_label, 'trim-sec=f'=>\$trim_sec,
           'n=i'=>\$n_remove, 'reg=f'=>\$reg,
           'out=s'=>\$out_file, 'out-removed=s'=>\$out_removed,
           'montage-only'=>\$montage_only, 'exclude=s'=>\$exclude_opt) or die;
my $task_file = shift @ARGV; defined $task_file or die
    "usage: $0 TASK.edf [VOLUNTARY_BLINKS.edf] [--veog L --trim-sec S --n N]\n".
    "               [--out CLEAN.edf] [--out-removed BLINK.edf]\n".
    "               [--exclude LABEL,LABEL,...] [--montage-only]\n";
my $vol_file  = shift @ARGV;

# Channels never filtered (hard pass-through): DC trigger channels and any
# Event / Marker / annotation / status / trigger channel.
sub hard_pass { my $l = shift;
    $l =~ /^DC\d/i || $l =~ /(?:mark|event|annot|status|trig)/i }

# ---- load ---- (read_edf returns a record hashref; labels normalised)
my ($task_raw, $task_labels, $sf, $task_units, $task_events, $task_start) = load_edf($task_file);
my $task  = $task_raw->double;                     # (nch,nt) uV
my @tlab  = map { $clean_label->($_) } @$task_labels;
printf "TASK %s: %d ch, %.0f Hz, %.1f s
", $task_file, scalar(@tlab), $sf, $task->dim(1)/$sf;
if ($trim_sec > 0) {
    my $s = int($trim_sec*$sf);
    $task = $task->slice(":,$s:-1")->sever;
    printf "  trimmed first %.1f s -> %.1f s
", $trim_sec, $task->dim(1)/$sf;
}
my %thave = map { $_ => 1 } @tlab;
my %tix;  $tix{$tlab[$_]} = $_ for 0..$#tlab;

my ($vol, $vlab, $vhave, $vix);
if (defined $vol_file) {
    my ($vol_raw, $vol_labels, $vsf) = load_edf($vol_file);
    $vol  = $vol_raw->double;
    my @vl = map { $clean_label->($_) } @$vol_labels;
    $vlab = \@vl; $vhave = { map {$_=>1} @vl }; $vix = {}; $vix->{$vl[$_]} = $_ for 0..$#vl;
    printf "FIT  %s: %d ch, %.0f Hz, %.1f s (voluntary blinks)
",
           $vol_file, scalar(@vl), $vsf, $vol->dim(1)/$vsf;
}

# ---- montage: every signal channel is filtered by default (EOG included) ----
die "vEOG label '$veog_label' not in TASK\n" unless $thave{$veog_label};
my %excl = map { (my $x=$_) =~ s/^\s+|\s+$//g; ($x => 1) }
           grep { length } split /,/, ($exclude_opt // '');
# vEOG drives blink detection, so it is always kept in the fit even if excluded
# from the output (its output is then the original signal).
my (@labels_used, @tsel, @vsel);
for my $i (0 .. $#tlab) {
    my $l = $tlab[$i];
    next if hard_pass($l);
    next if $excl{$l} && $l ne $veog_label;
    if (defined $vol) { next unless $vhave->{$l}; }
    push @labels_used, $l;
    push @tsel, $i;
    push @vsel, $vix->{$l} if defined $vol;
}
my $Xt = $task->dice_axis(0, pdl(long,@tsel))->sever;
my $Xv = defined $vol ? $vol->dice_axis(0, pdl(long,@vsel))->sever : undef;
my %u; $u{$labels_used[$_]} = $_ for 0..$#labels_used;
my $vi = $u{$veog_label};
{
    my @cleaned = grep { !$excl{$_} } @labels_used;
    my @passed  = grep {  $excl{$_} } @labels_used;
    printf "montage: %d channels filtered%s -- fit %s\n",
        scalar(@cleaned),
        (@passed ? " (passthrough: ".join(',',@passed).")" : ""),
        (defined $vol ? "from voluntary session" : "within task");
}

# ---- clean ----
my ($sub, $removed, $op) = remove_blinks($Xt, veog=>$vi, fit=>$Xv, sfreq=>$sf,
                                         n=>$n_remove, reg=>$reg);
printf "\nfit on %d blinks; blink component #%d (|r(vEOG)|=%.3f)\n",
       $op->{n_blinks}, $op->{blink_comp}, $op->{veog_corr}->at($op->{blink_comp});
{ my @c = $op->{veog_corr}->qsort->slice("-3:-1")->list;
  printf "top-3 component |r(vEOG)| = %.3f %.3f %.3f  (only the blink should be high)\n",
         reverse @c; }

# ---- diagnostics: blink-evoked reduction + blink-free preservation ----
sub p2p { my($e,$i)=@_; my $r=$e->slice("($i),"); $r->max-$r->min }
my $pk   = detect_blinks($Xt->slice("($vi),")->sever, sfreq=>$sf);
my $free = blink_free($Xt->slice("($vi),")->sever, sfreq=>$sf);
my $eo   = blink_evoked($Xt,  $pk, veog=>$vi, sfreq=>$sf, normalize=>0);
my $ec   = blink_evoked($sub, $pk, veog=>$vi, sfreq=>$sf, normalize=>0);
printf "\ntask has %d blinks. blink-evoked peak-to-peak (uV) raw -> clean [%% removed]\n", $pk->nelem;
for my $ch (qw/Fp1 Fp2 Fz F7 F8/) { next unless exists $u{$ch}; my $i=$u{$ch};
    my ($r,$c)=(p2p($eo,$i),p2p($ec,$i));
    printf "  %-4s %6.1f -> %6.1f  [%3.0f%%]\n", $ch,$r,$c,100*(1-$c/$r); }
print "posterior (should stay small)\n";
for my $ch (qw/O1 O2 Pz Cz/) { next unless exists $u{$ch}; my $i=$u{$ch};
    printf "  %-4s %6.1f -> %6.1f\n", $ch,p2p($eo,$i),p2p($ec,$i); }
print "\nblink-free preservation (var kept, ~1.00 = untouched)\n";
for my $ch (qw/Fp1 Fp2 Cz Pz O1 O2/) { next unless exists $u{$ch}; my $i=$u{$ch};
    my $o=$Xt->slice("($i),")->index($free);
    my $c=$sub->slice("($i),")->index($free);
    printf "  %-4s %.3f\n", $ch, sum(($c-$c->avg)**2)/sum(($o-$o->avg)**2); }

# ---- optional output (EDF via your EDF.pm writer) ----
if (defined $out_file || defined $out_removed) {
    die "write_edf not available in your PDL::EEG::IO::EDF; cannot save.\n"
        unless PDL::EEG::IO::EDF->can('write_edf');

    # output dispatch on file extension: .mul -> BESA ASCII (canonical writer),
    # anything else -> EDF.
    my $w = sub {
        my ($rec, $path) = @_;
        if ($path =~ /\.mul$/i) {
            eval { require PDL::EEG::IO::BESA::ASCII; 1 }
                or die "'.mul' output needs PDL::EEG::IO::BESA::ASCII: $@";
            PDL::EEG::IO::BESA::ASCII::write_mul($rec, $path);
        } else {
            PDL::EEG::IO::EDF::write_edf($rec, $path);
        }
    };

    # EDF+ annotation onsets shift with the trim (DC trigger channels are data
    # columns and are already trimmed consistently). Start time is passed through
    # so the writer keeps the recording date/time instead of a placeholder.
    my $ev = [ map { { %$_, onset => $_->{onset} - $trim_sec } }
               grep { !defined $_->{onset} || $_->{onset} >= $trim_sec }
               @{ $task_events || [] } ];

    if (!$montage_only) {
        # DEFAULT: whole recording. Every input channel is kept; the filtered
        # channels are replaced by their cleaned versions, everything else
        # (DC/ECG/markers + --exclude) passes through unchanged.
        die "full-record output needs per-channel units from read_edf (your reader\n".
            "returned none; use the record-hashref read_edf, or pass --montage-only)\n"
            unless ref $task_units eq 'ARRAY';
        my $orig  = $task->copy;                       # full (nch,nt) uV, trimmed
        my $clean = $task->copy;
        for my $i (0 .. $#labels_used) {
            next if $excl{$labels_used[$i]};            # --exclude channels stay original
            $clean->slice("(".$tsel[$i]."),") .= $sub->slice("($i),");
        }
        if (defined $out_file) {
            $w->({ data=>$clean->float, fs=>$sf, labels=>$task_labels,
                   units=>$task_units, events=>$ev, t_start=>$task_start }, $out_file);
            printf "wrote cleaned recording -> %s  (%d ch; %d filtered, rest passed through)\n",
                   $out_file, scalar(@$task_labels), scalar(grep { !$excl{$_} } @labels_used);
        }
        if (defined $out_removed) {
            my $art = $orig - $clean;                   # artifact on filtered ch, 0 elsewhere
            $w->({ data=>$art->float, fs=>$sf, labels=>$task_labels,
                   units=>$task_units, events=>$ev, t_start=>$task_start }, $out_removed);
            printf "wrote removed artifact  -> %s\n", $out_removed;
        }
    } else {
        # --montage-only: subset EDF of just the filtered channels.
        my @keep = grep { !$excl{$labels_used[$_]} } 0 .. $#labels_used;  # drop --exclude (e.g. vEOG)
        my $kidx = pdl(long, @keep);
        my @klab = @labels_used[@keep];
        my @u    = ('uV') x scalar(@klab);
        if (defined $out_file) {
            $w->({ data=>$sub->dice_axis(0,$kidx)->sever->float, fs=>$sf,
                   labels=>\@klab, units=>\@u, t_start=>$task_start }, $out_file);
            printf "wrote cleaned montage  -> %s  (%d ch)\n", $out_file, scalar(@klab);
        }
        if (defined $out_removed) {
            $w->({ data=>$removed->dice_axis(0,$kidx)->sever->float, fs=>$sf,
                   labels=>\@klab, units=>\@u, t_start=>$task_start }, $out_removed);
            printf "wrote removed montage  -> %s\n", $out_removed;
        }
    }
} else {
    print "(nothing saved; add --out CLEAN.edf [--out-removed BLINK.edf] [--montage-only])\n";
}

printf "returned in-memory: subtracted %s, removed %s (subtracted + removed = input montage)\n",
       join('x',$sub->dims), join('x',$removed->dims);
