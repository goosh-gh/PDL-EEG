#!/usr/bin/env perl
# examples/artifact_remove.pl
#
# Remove ocular / stereotyped artifacts from an EEG recording with GED spatial
# filters (PDL::EEG::GED) and write the cleaned recording back out.
#
# Artifact classes are declared on the command line, one --ged per class, in the
# order they should be removed:
#
#   --ged NAME[:CHANS[:SOURCE]]
#
#   NAME    a kind marked in the spans file by raw_plot.pl (e.g. blink,
#           hsaccade, slowA, vert_EOG, sweat_like). Any name is allowed.
#   CHANS   OPTIONAL. Signed channel list used ONLY to give the estimated axis a
#           meaning (a correlation report, and the source for --regress). It does
#           NOT build S or R. e.g. F7,Fp1,vEOG  or  Fp1,Fp2,-A1,-A2. Use 'none'
#           (or omit) for an axis with no natural reference -- then it is judged
#           by its eigenvalue gap and by inspecting its topography in the review.
#   SOURCE  OPTIONAL. How S (the artifact covariance) is built:
#             span  (default) covariance of the NAME segments;
#             blink           peak-normalised blink-evoked covariance, from
#                             --blink-from (a blink-only recording, same
#                             montage) or the target's own 'blink' segments.
#
# With no --ged at all, the default is a single blink class:
#   --ged blink::blink   (needs --blink-from or 'blink' spans)
#
# For each class, S comes from that artifact and R (the reference covariance)
# comes from the segments that are NOT that artifact and not 'bad', KEEPING every
# other declared artifact in R. A GED filter is orthogonal (in the R metric)
# only to what lives in R, so keeping the other artifacts in R is what stops one
# class's removal from disturbing another. Each S also excludes any samples that
# overlap another declared class, so classes do not contaminate each other's S.
#
# Removal is either spatial (default: PDL::EEG::GED::apply_ged, one operator per
# class, applied in the given order -- each is an exact, bounded oblique
# projection) or temporal (--regress: PDL::EEG::Regress::regress_out on the class
# reference time courses -- stronger, but also removes brain that co-varies in
# time). 'bad' segments are left exactly as recorded (restored after removal),
# to be excluded later at averaging.
#
# Before writing, a review image (waveform overlay of original vs cleaned, plus
# one 2D scalp topography per class showing the removed pattern) is rendered and
# the write is confirmed on the terminal. --no-review skips it; --yes writes
# without asking.
#
# Data layout is (nch, nt) throughout, matching PDL::EEG::GED and read_edf.
#
# Usage:
#   perl -Ilib examples/artifact_remove.pl task.edf \
#        --ged blink:vEOG:blink --blink-from voluntary_blinks.edf \
#        --ged sweat_like:Fp1,Fp2,-A1,-A2 \
#        --ged hsaccade:hEOG \
#        --spans task.edf.spans.tsv --montage standard_1020.elc \
#        --output edf --out task-clean.edf
#
# Options:
#   --ged SPEC        artifact class (repeatable, ordered); see above
#   --spans FILE      span annotations (default: <infile>.spans.tsv)
#   --blink-from FILE blink-only recording for a blink-SOURCE class
#   --regress         temporal regression instead of spatial projection
#   --montage FILE    .elc electrode file for the review topographies
#   --no-review       do not render the review image or ask
#   --yes             render the review but write without asking
#   --review-out PATH review image path (default <stem>-review.png)
#   --output mul|edf / --out PATH
#   --eeg/--extra/--eog/--veog/--heog   channel selection
#   --blink-thresh K / --blink-refr MS / --blink-half MS
#   --shrink G        R shrinkage passed to ged_operator reg (0.05)
#   --gap-min G       warn if an operator's eigenvalue gap is below G (3)
#   --decimals N / --fieldw N   .mul precision / column width (auto)

use strict;
use warnings;
use PDL;
use PDL::EEG::IO::NihonKohden qw(read_nk);
use PDL::EEG::IO::EDF         qw(read_edf write_edf clean_edf_label);
use PDL::EEG::IO::BESA::ASCII qw(write_mul);
use PDL::EEG::Spans           qw(read_spans spans_to_mask);
use PDL::EEG::GED             qw(ged_cov ged_operator apply_ged
                                 detect_blinks blink_evoked);
use PDL::EEG::Regress         qw(regress_out);
use Getopt::Long;

my %o = (output=>'mul', shrink=>0.05,
         'blink-thresh'=>2.5, 'blink-refr'=>200, 'blink-half'=>150,
         'gap-min'=>3.0, decimals=>2, 'eye-corr'=>0.4, 'lf-cut'=>5);
GetOptions(\%o, qw(ged=s@ spans=s blink-from=s regress
                   montage=s no-review yes review-out=s
                   rank eye-corr=f rank-top=i dump-components=s rank-out=s lf-cut=f
                   output=s out=s eeg=s extra=s eog=s veog=s heog=s
                   blink-thresh=f blink-refr=f blink-half=f shrink=f gap-min=f
                   decimals=i fieldw=i)) or die "bad options\n";
my $in = shift @ARGV or die "usage: $0 [options] file.edf|file.eeg\n";
$o{output} =~ /^(mul|edf)$/ or die "--output must be mul or edf\n";
(my $stem = $in) =~ s/\.[^.]+$//;   # output basename
!defined $o{'dump-components'} || $o{'dump-components'} =~ /^(mul|edf)$/
    or die "--dump-components must be mul or edf\n";

# ---- channel model ----------------------------------------------------------
my @TEN20 = qw(Fp1 Fp2 F7 F3 Fz F4 F8 T7 C3 Cz C4 T8 P7 P3 Pz P4 P8 O1 O2);
my %ALIAS = (T3=>'T7', T4=>'T8', T5=>'P7', T6=>'P8');
sub canon { my $l=shift; $l=~s/\s+//g;
    for my $c (@TEN20){return $c if lc$l eq lc$c} for my $a (keys %ALIAS){return $ALIAS{$a} if lc$l eq lc$a} undef }
my @EEG_EXTRA = defined $o{extra} ? grep {length} split /,/,$o{extra} : qw(E A1 A2 X1 BN1 BN2);
sub is_eog  { my $l=shift//''; $l=~/eog/i || $l=~/^X1$/i }
sub eog_kind{ my $l=shift//''; return 'h' if $l=~/h\s*eog/i; return 'v' if $l=~/v\s*eog/i||$l=~/^X1$/i; '?' }

# select GED channels (10-20 + extras + EOG); returns (\@idx, \@lab)
sub select_ged_channels {
    my ($labels) = @_;
    my @lab = @$labels; my %idx; $idx{$lab[$_]}=$_ for 0..$#lab;
    my (@i,@l,%used);
    if (defined $o{eeg}) { for (split /,/,$o{eeg}){ exists $idx{$_} or die "no channel $_\n"; push @i,$idx{$_};push @l,$_;$used{$idx{$_}}=1 } }
    else {
        for my $w (@TEN20){ for my $k (0..$#lab){ next if $used{$k}; if((canon($lab[$k])//'') eq $w){push @i,$k;push @l,$w;$used{$k}=1;last} } }
        for my $x (@EEG_EXTRA){ for my $k (0..$#lab){ next if $used{$k}; if(lc$lab[$k] eq lc$x){push @i,$k;push @l,$lab[$k];$used{$k}=1;last} } }
    }
    my @eog = defined $o{eog} ? (split /,/,$o{eog})
            : (grep { !$used{$idx{$_}} && is_eog($_) } @lab);
    for my $name (@eog){ next unless exists $idx{$name}; next if $used{$idx{$name}}; push @i,$idx{$name};push @l,$name;$used{$idx{$name}}=1 }
    return (\@i, \@l);
}

# read a recording (edf or nk) -> {data (nch,nt), fs, labels(clean), rec, device}
sub load_any {
    my ($file) = @_;
    if ($file =~ /\.edf$/i) {
        my $r = read_edf($file);
        return { data=>$r->{data}, fs=>$r->{fs},
                 labels=>[ map { clean_edf_label($_) } @{$r->{labels}} ],
                 rec=>$r, device=>'EDF'.(defined $r->{edf_type}?" ($r->{edf_type})":'') };
    }
    my $r = read_nk($file, all_blocks=>1);
    return { data=>$r->{data}, fs=>$r->{fs}, labels=>$r->{labels},
             rec=>$r, device=>$r->{device}//'NK' };
}

sub _corr { my($a,$b)=@_; my $x=$a->flat-$a->avg; my $y=$b->flat-$b->avg;
    sum($x*$y)/(sqrt(sum($x*$x)*sum($y*$y))+1e-30) }
sub _std  { my $x=shift; my $c=$x-$x->avg; sqrt(($c*$c)->avg) }

# component-k activation time course of an operator on X (nch,nt) -> (nt).
# Matches apply_ged's own Y = W' X, component index along dim1.
sub comp_tc {
    my ($op, $X, $k) = @_;
    my $Y = $op->{W}->transpose x $X->transpose;       # (nt, nch)
    return $Y->slice(":,($k)")->flat->sever;
}

# excess kurtosis of a time course (blink/spike components run high)
sub _kurtosis {
    my ($y) = @_; my $c = $y->flat - $y->flat->avg;
    my $v = ($c*$c)->avg; return 0 if $v <= 0;
    return ($c**4)->avg / ($v*$v) - 3;
}
# low-frequency power ratio via a boxcar low-pass proxy (no FFT dependency):
# var(smoothed)/var(raw), smoothing window ~ fs/(2*cut). Ocular activity is
# low-frequency, so eye components score high; muscle scores low.
sub _lf_ratio {
    my ($y, $fs, $cut) = @_;
    my $win = int($fs/(2*($cut||5))); $win = 1 if $win < 1;
    my $n = $y->dim(0); my $h = int($win/2);
    my $csp = append(pdl($y->type,0), $y->flat->cumusumover);   # (n+1) prefix sums
    my $lo  = (sequence($n)-$h)->clip(0,$n);
    my $hi  = (sequence($n)+$h+1)->clip(0,$n);
    my $sm  = ($csp->index($hi) - $csp->index($lo)) / ($hi-$lo);
    my $vy = (($y->flat-$y->flat->avg)**2)->avg;  return 0 if $vy<=0;
    return (($sm-$sm->avg)**2)->avg / $vy;
}

# ---- target -----------------------------------------------------------------
my $T   = load_any($in);
my $fs  = $T->{fs};
my ($ti,$tl) = select_ged_channels($T->{labels});
my $Xt  = $T->{data}->dice_axis(0, pdl(long,@$ti))->sever;   # (nsub, nt)
my $nt  = $Xt->dim(1);
printf "target: %s  %d samp @ %g Hz  | GED channels (%d): %s\n",
    $T->{device}, $nt, $fs, scalar @$tl, join(' ',@$tl);

# EOG references, both as global label index and as position within the subset
my %sub_pos; $sub_pos{$tl->[$_]} = $_ for 0..$#$tl;
my ($veog_g) = grep { eog_kind($T->{labels}[$_]) eq 'v' } 0..$#{$T->{labels}};
my ($heog_g) = grep { eog_kind($T->{labels}[$_]) eq 'h' } 0..$#{$T->{labels}};
$veog_g = (grep { $T->{labels}[$_] eq $o{veog} } 0..$#{$T->{labels}})[0] if defined $o{veog};
$heog_g = (grep { $T->{labels}[$_] eq $o{heog} } 0..$#{$T->{labels}})[0] if defined $o{heog};
my $veog = defined $veog_g ? $T->{data}->slice("($veog_g),:")->flat : undef;
my $heog = defined $heog_g ? $T->{data}->slice("($heog_g),:")->flat : undef;
my $veog_sub = defined $veog_g ? $sub_pos{ $T->{labels}[$veog_g] } : undef;   # row in $Xt

# signed z-scored reference time course from a CHANS spec (meaning only) -> (nt)
sub build_ref {
    my ($spec) = @_;
    return undef unless defined $spec && length $spec && lc($spec) ne 'none';
    my $ref = zeroes($nt); my $used = 0;
    for my $tok (split /,/, $spec) {
        $tok =~ s/^\s+|\s+$//g; next unless length $tok;
        my $sgn = 1; $sgn = -1 if $tok =~ s/^-//; $tok =~ s/^\+//;
        my ($gi) = grep { $T->{labels}[$_] eq $tok } 0..$#{$T->{labels}};
        unless (defined $gi) { warn "  (ref channel '$tok' not found - skipped)\n"; next }
        my $x = $T->{data}->slice("($gi),:")->flat;
        $ref = $ref + $sgn * ($x - $x->avg) / (_std($x) || 1);
        $used++;
    }
    return $used ? $ref : undef;
}

# ---- spans / masks ----------------------------------------------------------
my $spans_path = $o{spans} // "$in.spans.tsv";
my $spans = (-e $spans_path) ? read_spans($spans_path) : [];
printf "spans:  %d from %s\n", scalar @$spans, (-e $spans_path ? $spans_path : "(none)");
my $m_bad = spans_to_mask($spans, fs=>$fs, nt=>$nt, kind=>'bad');

# ---- parse --ged classes (ordered) ------------------------------------------
my @specs = @{ $o{ged} || [] };
@specs = ('blink::blink') unless @specs;             # default: blink only
my @classes;
for my $spec (@specs) {
    my ($name,$chans,$source) = split /:/, $spec, 3;
    defined $name && length $name or die "--ged: empty NAME in '$spec'\n";
    $source = (defined $source && length $source) ? $source : 'span';
    $source =~ /^(span|blink)$/ or die "--ged: SOURCE must be 'span' or 'blink' (got '$source')\n";
    push @classes, { name=>$name, chans=>(defined $chans?$chans:''), source=>$source };
}
# masks for every declared span-name (blink-source names may be absent -> zero)
my %mask; for my $c (@classes) {
    $mask{$c->{name}} //= spans_to_mask($spans, fs=>$fs, nt=>$nt, kind=>$c->{name});
}
sub union_other {                                    # OR of all OTHER classes' masks
    my ($self) = @_; my $u = zeroes(byte,$nt);
    for my $c (@classes) { next if $c->{name} eq $self; $u = $u | $mask{$c->{name}} }
    return $u;
}
sub gapmsg { my ($g)=@_; sprintf "gap=%.2f%s",$g,($g<$o{'gap-min'}?"  <-- WARNING < --gap-min=$o{'gap-min'}":'') }

# ---- build each class: S, R, operator, chosen component, tc, pattern --------
my @arts;
for my $c (@classes) {
    my ($name,$chans,$source) = @{$c}{qw(name chans source)};
    my $ref     = build_ref($chans);
    my $refname = (defined $ref) ? $chans : ($source eq 'blink' ? 'vEOG' : undef);
    my $refuse  = defined $ref ? $ref : ($source eq 'blink' ? $veog : undef);

    if ($source eq 'blink') {
        my ($src,$vi_src,$srclabel);
        if (defined $o{'blink-from'}) {
            my $B = load_any($o{'blink-from'});
            my ($bi,$bl) = select_ged_channels($B->{labels});
            @$bl==@$tl && join(',',@$bl) eq join(',',@$tl)
                or die "blink-from channels differ from target (need same montage):\n  target: @$tl\n  blink : @$bl\n";
            my ($bvg) = grep { eog_kind($B->{labels}[$_]) eq 'v' } 0..$#{$B->{labels}};
            defined $bvg or die "blink-from has no vertical EOG\n";
            $src=$B->{data}->dice_axis(0,pdl(long,@$bi))->sever;
            my @m=grep { $bl->[$_] eq $B->{labels}[$bvg] } 0..$#$bl; $vi_src=$m[0];
            $srclabel=$o{'blink-from'};
        } elsif ($mask{$name}->sum >= 50 && defined $veog_sub) {
            $src=$Xt; $vi_src=$veog_sub; $srclabel="target '$name' spans";
        } else { warn "  $name: no --blink-from and no usable '$name' spans - skipping\n"; next; }

        my $veog_src = $src->slice("($vi_src),:")->flat;
        my $peaks = detect_blinks($veog_src, sfreq=>$fs,
                        thresh_k=>$o{'blink-thresh'}, refractory_ms=>$o{'blink-refr'});
        if ($peaks->nelem < 8) { warn "  $name: too few blinks (".$peaks->nelem.") - skipping\n"; next; }
        my $evk = blink_evoked($src,$peaks,veog=>$vi_src,sfreq=>$fs,half_ms=>$o{'blink-half'},normalize=>1);
        my $S   = ged_cov($evk);
        # R on the target: exclude this class's spans, bad, and detected-blink windows
        my $m_bw = zeroes(byte,$nt);
        if (defined $veog) {
            my $tp = detect_blinks($veog, sfreq=>$fs, thresh_k=>$o{'blink-thresh'}, refractory_ms=>$o{'blink-refr'});
            my $half = int($o{'blink-half'}*$fs/1000);
            for my $p ($tp->list) { my $a=$p-$half; $a=0 if $a<0; my $b=$p+$half; $b=$nt-1 if $b>=$nt; $m_bw->slice("$a:$b").=1; }
        }
        my $Rmask = ($mask{$name} | $m_bad | $m_bw) == 0;
        $Rmask->sum >= 100 or die "  $name: too few reference samples for R\n";
        my $R  = ged_cov($Xt->dice_axis(1, which($Rmask)));
        my $op = ged_operator($S,$R,reg=>$o{shrink});
        # identify blink component by |corr| of comp activation on the raw evoked
        # with the vEOG-evoked time course
        my $evk_raw = blink_evoked($src,$peaks,veog=>$vi_src,sfreq=>$fs,half_ms=>$o{'blink-half'},normalize=>0);
        my $vy = $evk_raw->slice("($vi_src),:")->flat; $vy=$vy-$vy->avg;
        my $Yevk = $op->{W}->transpose x $evk_raw->transpose;
        my ($bc,$best)=(0,-1);
        for my $k (0..$op->{nch}-1){ my $cx=$Yevk->slice(":,($k)")->flat; $cx=$cx-$cx->avg;
            my $cc=abs(sum($cx*$vy))/sqrt(sum($cx**2)*sum($vy**2)+1e-30); ($bc,$best)=($k,$cc) if $cc>$best; }
        my $tc = comp_tc($op,$Xt,$bc);
        my $ev = $op->{evals};
        my $gap = abs($ev->at(-2))>0 ? $ev->at(-1)/abs($ev->at(-2)) : 1e30;
        printf "  %-10s %s  %s  (%d blinks from %s)\n", $name, gapmsg($gap),
            (defined $refuse ? sprintf("corr(%s)=%.3f",$refname,_corr($tc,$refuse)) : "corr=n/a"),
            $peaks->nelem, $srclabel;
        push @arts, { name=>$name, op=>$op, comp=>$bc, tc=>$tc,
                      pattern=>$op->{A}->slice("($bc),:")->flat->sever };
    }
    else {   # span source
        my $m_self  = $mask{$name};
        my $m_other = union_other($name);
        my $Smask = $m_self & (($m_other | $m_bad) == 0);
        my $Rmask = ($m_self | $m_bad) == 0;
        if ($Smask->sum < 50) { warn "  $name: too few '$name' samples (".$Smask->sum.") - skipping\n"; next; }
        $Rmask->sum >= 100 or die "  $name: too few reference samples for R\n";
        my $S  = ged_cov($Xt->dice_axis(1, which($Smask)));
        my $R  = ged_cov($Xt->dice_axis(1, which($Rmask)));
        my $op = ged_operator($S,$R,reg=>$o{shrink});
        my $comp = $op->{nch}-1;                        # top eigenvalue = last (ascending)
        my $tc = comp_tc($op,$Xt,$comp);
        my $ev = $op->{evals};
        my $gap = abs($ev->at(-2))>0 ? $ev->at(-1)/abs($ev->at(-2)) : 1e30;
        printf "  %-10s %s  %s\n", $name, gapmsg($gap),
            (defined $refuse ? sprintf("corr(%s)=%.3f",$refname,_corr($tc,$refuse))
                             : "corr=n/a (inspect topography)");
        push @arts, { name=>$name, op=>$op, comp=>$comp, tc=>$tc,
                      pattern=>$op->{A}->slice("($comp),:")->flat->sever };
    }
}
@arts or die "nothing to remove (no usable artifact axes)\n";

# ---- diagnostic mode: rank every component, do NOT remove/write -------------
if ($o{rank}) {
    my @dump;        # selected components to draw / write: {class,k,y,pattern,lam}
    for my $a (@arts) {
        my $op = $a->{op}; my $nch = $op->{nch}; my $ev = $op->{evals};
        printf "\n=== class '%s' : per-component diagnostics (removal candidate = #%d) ===\n",
            $a->{name}, $a->{comp};
        printf "  %-4s %10s %8s %8s %8s %7s %7s  %s\n",
            'comp','lambda','r(vEOG)','r(hEOG)','kurt','lf','|a|max','flag';
        # rank by eigenvalue, largest first (evals are ascending)
        my @order = reverse(0..$nch-1);
        my $eye_c = $o{'eye-corr'};
        for my $k (@order) {
            my $y  = comp_tc($op,$Xt,$k);
            my $rv = defined $veog ? _corr($y,$veog) : 0;
            my $rh = defined $heog ? _corr($y,$heog) : 0;
            my $ku = _kurtosis($y);
            my $lf = _lf_ratio($y,$fs,$o{'lf-cut'});
            my $pat= $op->{A}->slice("($k),:")->flat;
            my $amax = $pat->abs->max->sclr;
            my @fl;
            push @fl,'eye-v' if defined $veog && abs($rv)>=$eye_c;
            push @fl,'eye-h' if defined $heog && abs($rh)>=$eye_c;
            push @fl,'CAND'  if $k==$a->{comp};
            my $flag = join(',',@fl);
            printf "  %-4d %10.3g %8.3f %8.3f %8.2f %7.3f %7.2g  %s\n",
                $k,$ev->at($k),$rv,$rh,$ku,$lf,$amax,$flag;
            # select for png/dump: flagged (eye-*) OR the removal candidate
            if (@fl && $flag ne '') {
                my $eye = grep { /^eye-/ } @fl;
                push @dump, { class=>$a->{name}, k=>$k, y=>$y, pattern=>$pat->sever,
                              lam=>$ev->at($k), rv=>$rv, rh=>$rh }
                    if $eye || $k==$a->{comp};
            }
        }
        # optionally include top-K by eigenvalue even if unflagged
        if ($o{'rank-top'}) {
            for my $k (@order[0 .. ($o{'rank-top'}-1 < $#order ? $o{'rank-top'}-1 : $#order)]) {
                next if grep { $_->{class} eq $a->{name} && $_->{k}==$k } @dump;
                my $y=comp_tc($op,$Xt,$k);
                push @dump, { class=>$a->{name}, k=>$k, y=>$y,
                              pattern=>$op->{A}->slice("($k),:")->flat->sever,
                              lam=>$ev->at($k),
                              rv=>(defined $veog?_corr($y,$veog):0),
                              rh=>(defined $heog?_corr($y,$heog):0) };
            }
        }
    }
    printf "\nselected %d component(s) for figure/dump\n", scalar @dump;

    # --- diagnostic PNG: one row per selected component -----------------------
    if (defined $o{montage} && @dump) {
        my $rout = $o{'rank-out'} // "$stem-rank.png";
        my $ok = eval {
            require PDL::Graphics::Cairo; PDL::Graphics::Cairo->import(qw(figure));
            require PDL::EEG::MAP2D;      PDL::EEG::MAP2D->import(qw(plot_topomap));
            my $R = scalar @dump;
            my $fig = figure(width=>1200, height=>90+150*$R);
            my $W=2000; my $step=int($nt/$W)||1;
            my $t = (sequence(int(($nt+$step-1)/$step))*$step)/$fs;
            for my $r (0..$R-1) {
                my $d = $dump[$r];
                my $top = 1 - 0.02 - $r*(1/$R);
                my $bot = 1 - 0.02 - ($r+0.9)*(1/$R);
                # topo (left)
                my $gtl = $fig->add_gridspec(1,1, left=>0.02, right=>0.16, top=>$top, bottom=>$bot);
                my $pax = $fig->add_subplot($gtl->at(0,0));
                plot_topomap(values=>$d->{pattern}, labels=>$tl, montage=>$o{montage},
                             _ax=>$pax, _fig=>$fig, colorbar=>0, names=>0,
                             title=>sprintf("%s #%d",$d->{class},$d->{k}));
                # waveforms (right): top-3 representative channels (|a| largest) + component y
                my $gwr = $fig->add_gridspec(1,1, left=>0.20, right=>0.99, top=>$top, bottom=>$bot);
                my $ax  = $fig->add_subplot($gwr->at(0,0));
                my $ord = $d->{pattern}->abs->qsorti->slice("-1:0");     # channels by |a| desc
                my @rep = $ord->slice("0:2")->list;
                my $scale = _std($Xt->flat)||1; my $sp = 4*$scale;
                for my $j (0..$#rep) {
                    my $ci = $rep[$j];
                    my $w1 = $Xt->slice("($ci),0:-1:$step")->flat;
                    my $tt = $t->slice("0:".($w1->dim(0)-1));
                    my $off = (3-$j)*$sp;
                    $ax->axhline($off, color=>'#EEEEEE', lw=>0.3);
                    $ax->line($tt, $w1+$off, color=>'#1565C0', lw=>0.5);
                    $ax->text($tt->at(0), $off+0.3*$sp, $tl->[$ci], fontsize=>8);
                }
                # component y (scaled) at the bottom slot
                my $y1 = $d->{y}->slice("0:-1:$step");
                my $ty = $t->slice("0:".($y1->dim(0)-1));
                my $yz = ($y1-$y1->avg)/((_std($y1)||1)) * (0.8*$sp);
                $ax->axhline(0, color=>'#EEEEEE', lw=>0.3);
                $ax->line($ty, $yz, color=>'#C62828', lw=>0.6);
                $ax->text($ty->at(0), 0.3*$sp, "y (comp)", fontsize=>8, color=>'#C62828');
                $ax->set_yticks(pdl([]));
                $ax->set_title(sprintf("%s #%d  lambda=%.3g  r(vEOG)=%.2f r(hEOG)=%.2f  [rep: %s]",
                    $d->{class},$d->{k},$d->{lam},$d->{rv},$d->{rh}, join(' ',map{$tl->[$_]}@rep)),
                    fontsize=>9);
            }
            $fig->save($rout); 1;
        };
        if ($ok) { printf "rank figure: %s\n", $rout; }
        else     { warn "rank figure failed ($@)\n"; }
    } elsif (@dump && !defined $o{montage}) {
        warn "(--montage not given; skipping rank topographies)\n";
    }

    # --- optional component dump (time courses y as channels) -----------------
    if (defined $o{'dump-components'} && @dump) {
        my $fmt = $o{'dump-components'};
        my $D = zeroes($Xt->type, scalar(@dump), $nt);
        my @lab; my @un;
        for my $r (0..$#dump) {
            $D->slice("($r),:") .= $dump[$r]{y};
            push @lab, "GED_$dump[$r]{class}_$dump[$r]{k}";
            push @un, 'a.u.';                       # component scores are NOT uV
        }
        my $dout = "$stem-components.$fmt";
        if ($fmt eq 'edf') {
            write_edf({ data=>$D, fs=>$fs, labels=>\@lab, units=>\@un,
                        t_start=>$T->{rec}{t_start} }, $dout);
        } else {
            write_mul({ data=>$D, fs=>$fs, labels=>\@lab,
                        t_start=>$T->{rec}{t_start} }, $dout, decimals=>4, fieldw=>12);
        }
        printf "component dump: %s  (%d component scores, arbitrary units, NOT uV)\n",
            $dout, scalar @dump;
    }
    exit 0;
}

# ---- removal ----------------------------------------------------------------
my $cleaned_sub;
if ($o{regress}) {
    ($cleaned_sub) = regress_out($Xt, map { $_->{tc} } @arts);          # temporal
    printf "removal: temporal regression (PDL::EEG::Regress) on %d source(s)\n", scalar @arts;
}
elsif (@arts == 1) {
    $cleaned_sub = apply_ged($arts[0]{op}, $Xt, idx=>[ $arts[0]{comp} ]);
    printf "removal: spatial projection, 1 component (%s)\n", $arts[0]{name};
}
else {
    # sequential single-operator projection, in the --ged order. Each apply_ged
    # is one operator with W'A=I, so each step is an exact, bounded oblique
    # projection (no cross-operator (W'A)^{-1} to ill-condition on overlap).
    my $cur = $Xt;
    for my $a (@arts) { $cur = apply_ged($a->{op}, $cur, idx=>[ $a->{comp} ]); }
    $cleaned_sub = $cur;
    printf "removal: spatial projection, %d axes in order (%s)\n",
        scalar @arts, join('+', map { $_->{name} } @arts);
}

# 'bad' segments are left EXACTLY as recorded (removal ran over the whole record)
if ($m_bad->sum) {
    my $bm = $m_bad->dummy(0);
    $cleaned_sub = $cleaned_sub*(1-$bm) + $Xt*$bm;
}

# ---- review image (waveform overlay + per-class topographies) ---------------
sub render_review {
    my ($png) = @_;
    require PDL::Graphics::Cairo; PDL::Graphics::Cairo->import(qw(figure));
    my $have_topo = defined $o{montage};
    my $fig = figure(width=>1200, height=>($have_topo?840:600));

    # waveform band (top)
    my $gw = $fig->add_gridspec(1,1, left=>0.05, right=>0.99,
                 top=>0.96, bottom=>($have_topo?0.42:0.06));
    my $ax = $fig->add_subplot($gw->at(0,0));
    my $W=2000; my $step=int($nt/$W)||1;
    my $t = (sequence(int(($nt+$step-1)/$step))*$step)/$fs;
    my $scale = _std($Xt->flat) || 1; my $spacing = 6*$scale;
    my $n = scalar @$tl;
    for my $i (0..$n-1) {
        my $o1=$Xt->slice("($i),0:-1:$step")->flat;
        my $c1=$cleaned_sub->slice("($i),0:-1:$step")->flat;
        my $tt=$t->slice("0:".($o1->dim(0)-1));
        my $off=($n-1-$i)*$spacing;
        $ax->axhline($off, color=>'#EEEEEE', lw=>0.3);
        $ax->line($tt,$o1+$off, color=>'#BDBDBD', lw=>0.5);
        $ax->line($tt,$c1+$off, color=>'#1565C0', lw=>0.5);
    }
    $ax->set_title("original (grey) vs cleaned (blue)  |  removed: ".join('+',map{$_->{name}}@arts));
    $ax->set_yticks(pdl([])); $ax->set_xlabel('time (s)');

    if ($have_topo) {
        require PDL::EEG::MAP2D; PDL::EEG::MAP2D->import(qw(plot_topomap));
        my $k = scalar @arts;
        my $gt = $fig->add_gridspec(1,$k, left=>0.04, right=>0.99,
                     top=>0.34, bottom=>0.04, wspace=>0.20);
        for my $c (0..$k-1) {
            my $pax = $fig->add_subplot($gt->at(0,$c));
            plot_topomap(values=>$arts[$c]{pattern}, labels=>$tl,
                         montage=>$o{montage}, _ax=>$pax, _fig=>$fig,
                         colorbar=>0, names=>0, title=>$arts[$c]{name});
        }
    }
    $fig->save($png);
    return $png;
}

unless ($o{'no-review'}) {
    my $rout = $o{'review-out'} // "$stem-review.png";
    my $ok = eval { render_review($rout); 1 };
    if ($ok) { printf "review image: %s\n", $rout; }
    else     { warn "review image failed ($@); continuing\n"; }
    unless ($o{yes}) {
        print "write $o{output} output now? [y/N] ";
        my $ans = <STDIN>;
        unless (defined $ans && $ans =~ /^\s*y/i) { print "aborted (no file written)\n"; exit 0; }
    }
}

# ---- assemble full record (only GED channels changed) + write ---------------
my $outdata = $T->{data}->copy;
for my $j (0 .. $#$ti) { $outdata->slice("(".$ti->[$j]."),:") .= $cleaned_sub->slice("($j),:"); }

my $out = $o{out} // "$stem-clean.$o{output}";
if ($o{output} eq 'edf') {
    my %rec = %{ $T->{rec} };
    $rec{data} = $outdata;
    write_edf(\%rec, $out);
} else {
    my $fieldw = $o{fieldw};
    unless (defined $fieldw) {
        my $mx=$T->{data}->abs->max->sclr; my $dig=($mx>=1)?int(log($mx)/log(10))+1:1;
        my $need=$dig+$o{decimals}+3; $fieldw=$need>8?$need:8;
    }
    write_mul({ data=>$outdata, fs=>$fs, labels=>$T->{labels},
                t_start=>$T->{rec}{t_start} }, $out,
              decimals=>$o{decimals}, fieldw=>$fieldw);
}
printf "wrote %s  (removed: %s;  bad segments left in place)\n",
    $out, join('+', map { $_->{name} } @arts);
