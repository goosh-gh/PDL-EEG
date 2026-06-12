# PDL::EEG

EEG analysis toolkit for PDL — MNE-Python inspired.

## Status

Version 0.01 — initial release. Currently implements:

- **`PDL::EEG::IO::NihonKohden`** — read Nihon Kohden EEG-1100C binary files (`.EEG`)

Planned: `PDL::EEG::IO::EDF`, `PDL::EEG::Epochs`, `PDL::EEG::Evoked`, `PDL::EEG::Viewer`

## Requirements

- Perl 5.20+
- PDL 2.080+
- PDL::Graphics::Cairo (for `--plot`, from [goosh-gh/PDL-Graphics-Cairo](https://github.com/goosh-gh/PDL-Graphics-Cairo))
- giza_server (for interactive viewer, from [goosh-gh/giza-server](https://github.com/goosh-gh/giza-server))

## Installation

```bash
perl Makefile.PL
make
make testdata   # generate synthetic test data
make test
make install
```

## Usage

```perl
use PDL::EEG::IO::NihonKohden qw(read_nk);

my $rec = read_nk('patient.EEG');
# $rec->{data}      [n_ch, n_samples] float32 µV
# $rec->{fs}        sampling rate Hz (auto-detected)
# $rec->{labels}    channel names (from .21E file)
# $rec->{t_start}   "YYYY-MM-DD HH:MM:SS"
# $rec->{events}    [{t => $sec, label => $str}, ...]
# $rec->{device}    "EEG-1100C V01.00"
# $rec->{n_blocks}  number of waveform blocks in file

# Read a specific block
my $rec2 = read_nk('patient.EEG', block => 2);
```

### Interactive viewer

```bash
perl -Ilib examples/read_nihonkohden.pl patient.EEG --plot
perl -Ilib examples/read_nihonkohden.pl patient.EEG --plot --block 2 --nch 16 --sec 10
```

Horizontal slider: time scroll  
Vertical slider: gain (µV/div)

### Verify a file

```bash
perl -Ilib examples/verify_read.pl patient.EEG
```

## Format notes

Confirmed from real EEG-1100C hardware (JE-921A amplifier):

- ADC: uint16 offset binary, `0.09765625 µV/bit` (fixed, not stored in file)
- Sampling rate: lower 14 bits of u16 at `wfmblock+0x1A`
- Channel table: 10-byte entries at `wfmblock+0x2F`
- Data start: `wfmblock+0x171`
- Events: 6-digit ASCII seconds in `.LOG` file

## License

Same terms as Perl itself (Artistic License 2.0 or GPL-1+).  
Format knowledge derived from EDFbrowser (GPL-2) used as reference only.

## See also

- [goosh-gh/PDL-Graphics-Cairo](https://github.com/goosh-gh/PDL-Graphics-Cairo)
- [goosh-gh/giza-server](https://github.com/goosh-gh/giza-server)
- [goosh-gh/PDL-IO-PNG](https://github.com/goosh-gh/PDL-IO-PNG)
