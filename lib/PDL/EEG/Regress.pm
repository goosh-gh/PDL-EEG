package PDL::EEG::Regress;

# Temporal nuisance regression for EEG: remove, from every channel, the
# least-squares linear fit on a set of source time courses.
#
# The sources are method-agnostic -- GED component activations, ICA component
# time courses, or raw EOG channels. With EOG channels as sources this is the
# classic Gratton-Coles ocular regression; with GED/ICA source time courses it
# is the temporal counterpart to a spatial-projection removal.
#
# It removes ALL variance time-correlated with the sources, so it is more
# aggressive than a spatial projection (it can also remove brain signal that
# happens to co-vary in time). Hand it genuine artifact sources only, and
# prefer clean source estimates (e.g. GED/ICA activations) over raw EOG.
#
# The sources are orthonormalised (modified Gram-Schmidt) before removal, so
# correlated sources are partialled out jointly rather than double-subtracted.
# Per-channel and per-source means are removed in the fit; the channel DC level
# is preserved in the output.
#
# Data layout: X is (nch, nt), matching PDL::EEG::GED and read_edf.

use strict;
use warnings;
use PDL;
use Exporter 'import';
our @EXPORT_OK = qw(regress_out);
our $VERSION = '0.01';

# regress_out($X, @sources) -> $cleaned  (or ($cleaned,$removed) in list ctx)
#   $X       : (nch, nt)
#   @sources : one or more (nt) time courses to regress out
sub regress_out {
    my ($X, @sources) = @_;
    die "regress_out: X must be 2-D (nch,nt)\n" unless $X->ndims == 2;
    die "regress_out: need at least one source\n" unless @sources;
    my $nt = $X->dim(1);

    my @u;                                         # orthonormal basis of sources
    for my $s (@sources) {
        my $v = $s->flat;
        die "regress_out: source length != nt\n" unless $v->nelem == $nt;
        $v = $v - $v->avg;
        for my $w (@u) { $v = $v - sum($v * $w) * $w; }   # w already unit-norm
        my $nrm = sqrt(sum($v * $v));
        next if $nrm < 1e-12;                       # drop a source with no new info
        push @u, $v / $nrm;
    }
    die "regress_out: sources are degenerate (no usable direction)\n" unless @u;

    my $removed = zeroes($X->type, $X->dims);       # (nch,nt)
    for my $w (@u) {
        my $c = ($X * $w->dummy(0))->mv(1,0)->sumover;      # (nch) = X . w per channel
        $removed = $removed + $c->dummy(1) * $w->dummy(0);  # (nch,nt)
    }
    my $cleaned = $X - $removed;
    return wantarray ? ($cleaned, $removed) : $cleaned;
}

1;

__END__

=head1 NAME

PDL::EEG::Regress - temporal nuisance / EOG regression for EEG

=head1 SYNOPSIS

    use PDL::EEG::Regress qw(regress_out);

    # $X : (nch, nt).  Sources are (nt) time courses.
    # Classic EOG regression -- feed the EOG channels themselves:
    my $clean = regress_out($X, $veog, $heog);

    # Or remove GED / ICA source activations (label-free: you need not know
    # which component is "vertical" or "horizontal", only that they are ocular):
    my ($clean, $removed) = regress_out($X, $blink_tc, $sacc_tc);

=head1 DESCRIPTION

Removes, from every channel of C<$X> (nch,nt), the least-squares linear fit on
the supplied source time courses, returning the residual. The sources are
orthonormalised first, so correlated sources are partialled out jointly. This
is temporal regression: it removes all variance time-correlated with the
sources, which is more aggressive than a spatial projection and can remove
genuine brain signal that co-varies with the sources. Prefer clean source
estimates and hand it only true artifact sources.

Data layout is C<(nch, nt)>, matching L<PDL::EEG::GED> and C<read_edf>.

=head1 SEE ALSO

L<PDL::EEG::GED> (spatial-projection removal), L<PDL::EEG::ICA>.

=cut
