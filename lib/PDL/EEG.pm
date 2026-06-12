package PDL::EEG;

use strict;
use warnings;

our $VERSION = '0.01';

=head1 NAME

PDL::EEG - EEG analysis toolkit for PDL (MNE-Python inspired)

=head1 SYNOPSIS

  use PDL::EEG::IO::NihonKohden qw(read_nk);

  my $rec = read_nk('patient.eeg');
  # $rec->{data}   [n_ch, n_samples] float32 µV
  # $rec->{fs}     sampling rate Hz
  # $rec->{labels} channel labels
  # $rec->{events} [{t=>$sec, label=>$str}, ...]

=head1 DESCRIPTION

PDL::EEG is a growing EEG analysis framework for Perl/PDL,
aiming to provide functionality comparable to MNE-Python.

Current modules:

  PDL::EEG::IO::NihonKohden   Read Nihon Kohden *.eeg binary files
  PDL::EEG::IO::EDF            (planned) Read EDF/EDF+ files
  PDL::EEG::Epochs             (planned) Epoch extraction
  PDL::EEG::Evoked             (planned) Averaged ERP
  PDL::EEG::Viewer             (planned) Interactive raw viewer (raw.show())
  PDL::EEG::Preprocessing::ICA (planned) Independent component analysis
  PDL::EEG::TimeFreq::STFT     (planned) Time-frequency analysis

=head1 SEE ALSO

L<PDL>, L<PDL::EEG::IO::NihonKohden>

=head1 AUTHOR

goosh

=head1 LICENSE

Same terms as Perl itself.

=cut

1;
