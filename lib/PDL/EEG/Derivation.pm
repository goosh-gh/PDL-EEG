package PDL::EEG::Derivation;

use strict;
use warnings;
use Carp qw(croak carp);
use PDL;
use Exporter 'import';

our @EXPORT_OK = qw(derive bne rereference);
our $VERSION   = '0.01';

=head1 NAME

PDL::EEG::Derivation - Linear channel derivations (re-reference, bipolar, ...)

=head1 SYNOPSIS

  use PDL::EEG::IO::NihonKohden qw(read_nk);
  use PDL::EEG::Derivation      qw(bne derive);

  my $rec = read_nk('JJ0090J6.EEG', all_blocks => 1);

  # Balanced non-cephalic (BNE) re-reference: y = x - (prop*BN1 + (1-prop)*BN2)
  my $bn = bne($rec, prop => 0.5, suffix => '-BN');   # BN1=V, BN2=S

  # ...or a general linear derivation y = M x  (M given as rows of weights):
  my $bip = derive($rec, [[1,-1, 0], [0,1,-1]], ['Fp1-Fp2','Fp2-F3']);

=head1 DESCRIPTION

An EEG re-reference, bipolar montage, common-average reference (CAR) and, given
electrode geometry elsewhere, a surface Laplacian / interpolation are all the
same object: a B<linear derivation> C<y = M x>, where C<x> is the recorded data
C<[n_ch, n_samp]> and C<M> is an C<[n_out, n_in]> weight matrix. This module
provides that primitive (L</derive>) plus constructors for common schemes
(L</bne>, L</rereference>). Geometry-based operators (CSD, spline interpolation)
belong in their own modules and feed their matrix to L</derive>.

=head2 Why BNE cancels the acquisition reference

Nihon Kohden acquires against a system reference (here C<Avr(C3,C4)>), so the
recorded channel is C<x_i = s_i - s_ref>. Re-referencing to
C<r = prop*V + (1-prop)*S> gives C<y_i = s_i - r>. Because the recorded BN
channels are also referenced to C<s_ref>,

  prop*x_V + (1-prop)*x_S = (prop*V + (1-prop)*S) - (prop+(1-prop))*s_ref
                         = r - s_ref            (since prop + (1-prop) = 1)

so C<y_i = x_i - (prop*x_V + (1-prop)*x_S) = s_i - r> exactly: the acquisition
reference C<s_ref> cancels and need not be known. This relies only on the
weights summing to 1, which L</bne> enforces by construction.

=head1 FUNCTIONS

=head2 derive($rec, \@rows, \@labels, %opt)

Apply a linear derivation. C<@rows> is C<n_out> arrayrefs each holding C<n_in>
weights, so C<data_new[o] = sum_i rows[o][i] * data[i]>. C<@labels> gives the
C<n_out> output channel names. Returns a new record hashref (data replaced,
C<fs>/C<t_start> carried over).

=cut

sub derive {
    my ($rec, $rows, $labels, %opt) = @_;
    croak "derive: record hashref required" unless ref $rec eq 'HASH';
    my $data = $rec->{data};
    croak "derive: \$rec->{data} must be a 2-D PDL"
        unless eval { $data->isa('PDL') } && $data->ndims == 2;

    my $n_in  = $data->dim(0);
    my $n_samp = $data->dim(1);
    my $n_out = scalar @$rows;
    croak "derive: got $n_out label(s) for $n_out output row(s)"
        unless @$labels == $n_out;

    for my $o (0 .. $n_out - 1) {
        croak "derive: row $o has ${\ scalar @{$rows->[$o]}} weights, need $n_in"
            unless @{ $rows->[$o] } == $n_in;
    }

    my $out = zeroes(float, $n_out, $n_samp);
    for my $o (0 .. $n_out - 1) {
        my $acc = zeroes(float, $n_samp);
        my $w   = $rows->[$o];
        for my $i (0 .. $n_in - 1) {
            next unless $w->[$i];                       # skip zero weights
            $acc = $acc + $w->[$i] * $data->slice("($i),:");
        }
        my $orow = $out->slice("($o),:");
        $orow .= $acc;
    }

    return {
        %$rec,
        data   => $out,
        labels => [ @$labels ],
        ($opt{meta} ? %{ $opt{meta} } : ()),
    };
}

=head2 bne($rec, %opt)

Balanced non-cephalic (BNE) re-reference. Options:

  prop     => 0.5        weight on V (BN1); S (BN2) gets 1-prop  (must sum to 1)
  v        => 'BN1'      label of the V (vertebral) reference electrode
  s        => 'BN2'      label of the S (sternal) reference electrode
  exclude  => qr/.../    channels passed through unchanged (not re-referenced);
                         default qr/^(?:DC\d+|Trigger)$/i (non-EEG channels)
  drop_ref => 1          drop the V/S channels from the output (default 1)
  suffix   => undef      appended to re-referenced channel labels (e.g. '-BN')

=cut

sub bne {
    my ($rec, %opt) = @_;
    my $prop = defined $opt{prop} ? $opt{prop} : 0.5;
    my $vlab = defined $opt{v} ? $opt{v} : 'BN1';
    my $slab = defined $opt{s} ? $opt{s} : 'BN2';
    my $excl = defined $opt{exclude} ? $opt{exclude} : qr/^(?:DC\d+|Trigger)$/i;
    my $drop = exists $opt{drop_ref} ? $opt{drop_ref} : 1;
    my $suf  = $opt{suffix};

    my @labels = @{ $rec->{labels} // croak "bne: \$rec->{labels} required" };
    my %idx; $idx{ $labels[$_] } = $_ for 0 .. $#labels;
    my $vi = $idx{$vlab}; my $si = $idx{$slab};
    croak "bne: reference electrode '$vlab' not found in labels" unless defined $vi;
    croak "bne: reference electrode '$slab' not found in labels" unless defined $si;

    my $wv = $prop;
    my $ws = 1 - $prop;
    # weights sum to 1 by construction -> acquisition reference cancels.

    my $n_in = scalar @labels;
    my (@rows, @out_labels);
    for my $i (0 .. $n_in - 1) {
        next if $drop && ($i == $vi || $i == $si);      # reference source chans

        my @w = (0) x $n_in;
        if ($labels[$i] =~ $excl) {                     # DC/Trigger: pass through
            $w[$i] = 1;
            push @out_labels, $labels[$i];
        } else {                                        # EEG: subtract BNE ref
            $w[$i]  = 1;
            $w[$vi] -= $wv;
            $w[$si] -= $ws;
            push @out_labels, (defined $suf ? "$labels[$i]$suf" : $labels[$i]);
        }
        push @rows, \@w;
    }

    return derive($rec, \@rows, \@out_labels,
                  meta => { reference => 'BNE', bne_prop => $prop });
}

=head2 rereference($rec, \@ref_labels, %opt)

Re-reference every EEG channel to the (equally weighted) average of the named
channels: C<y_i = x_i - mean(ref)>. C<\@ref_labels> may be a single label
(single-electrode reference) or several (e.g. linked-ears, or the whole set for
a common average reference). Accepts the same C<exclude>/C<suffix> options as
L</bne>.

=cut

sub rereference {
    my ($rec, $ref_labels, %opt) = @_;
    my $excl = defined $opt{exclude} ? $opt{exclude} : qr/^(?:DC\d+|Trigger)$/i;
    my $suf  = $opt{suffix};

    my @labels = @{ $rec->{labels} // croak "rereference: labels required" };
    my %idx; $idx{ $labels[$_] } = $_ for 0 .. $#labels;

    my @rlist = ref $ref_labels ? @$ref_labels : ($ref_labels);
    my @ri;
    for my $r (@rlist) {
        croak "rereference: reference channel '$r' not found" unless defined $idx{$r};
        push @ri, $idx{$r};
    }
    my $rw = 1 / @ri;                                   # equal weights, sum to 1

    my $n_in = scalar @labels;
    my (@rows, @out_labels);
    for my $i (0 .. $n_in - 1) {
        my @w = (0) x $n_in;
        if ($labels[$i] =~ $excl) {
            $w[$i] = 1;
            push @out_labels, $labels[$i];
        } else {
            $w[$i] = 1;
            $w[$_] -= $rw for @ri;
            push @out_labels, (defined $suf ? "$labels[$i]$suf" : $labels[$i]);
        }
        push @rows, \@w;
    }
    return derive($rec, \@rows, \@out_labels,
                  meta => { reference => join('+', @rlist) });
}

1;

=head1 SEE ALSO

L<PDL::EEG::IO::NihonKohden>, L<PDL::EEG::IO::BESA::ASCII>

=head1 AUTHOR

goosh

=cut
