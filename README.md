# PDL::EEG

EEG analysis toolkit for PDL — MNE-Python inspired.

## Status

Version 0.02. Currently implements:

- **`PDL::EEG::IO::NihonKohden`** — read Nihon Kohden EEG binary files (`.EEG`).
  Both on-disk layouts are supported and auto-selected from the file's format
  signature:
  - **`wfmblock`** (legacy) — `EEG-1100A/B/C`, `EEG-2100`, `QI-403A`, `DAE-2100D`
  - **`extblock`** (newer) — `EEG-1200A` (e.g. the EEG-1290 recorder / JE-92NX headbox)

Planned: `PDL::EEG::IO::EDF`, `PDL::EEG::Evoked`, `PDL::EEG::Epochs`, `PDL::EEG::Viewer`.

## Supported Nihon Kohden formats

The on-disk layout is chosen from the **format signature** at offset `0x0000`,
never from the physical recorder model (which is not stored in the file):

| Layout | Signature examples | Channel info location |
|--------|--------------------|-----------------------|
| `wfmblock` | `EEG-1100C V01.00`, `EEG-2100 V01.00`, `QI-403A V01.00`, `DAE-2100D V01.30` | inside the wfmblock |
| `extblock` | `EEG-1200A V01.00` | extended-block chain via `ext_address` (`0x03EE`) |

`read_nk()` dispatches automatically. For unseen signatures it falls back to a
structural check (`ext_address` non-zero ⇒ `extblock`), so future firmware
variants (e.g. `EEG-1200B/C`, or an `EEG-1214` migrating to `EEG-1200A`) read
without code changes; signatures that are not Nihon Kohden are rejected.

Two helpers are also exported:

```perl
use PDL::EEG::IO::NihonKohden qw(read_nk nk_layout nk_format_hint);

my ($sig, $layout, $how) = nk_layout('subjct.EEG');   # authoritative (reads the file)
# e.g. ("EEG-1200A V01.00", "extblock", "table")

my ($guess_sig, $guess_layout, $note) = nk_format_hint('EEG-1290');
# NON-authoritative hint from a recorder model name; always confirm via nk_layout()
```

> **Model vs. format.** A recorder can change the format it emits after a
> software update, so `nk_format_hint()` is only a pre-flight guess. The
> signature in the file always wins.

## Requirements

- Perl 5.20+
- PDL 2.080+
- [PDL::Graphics::Cairo](https://github.com/goosh-gh/PDL-Graphics-Cairo) — for interactive viewer and notebook plots
- [giza-server](https://github.com/goosh-gh/giza-server) — for native-window interactive viewer (optional)
- [App-PDL-Notebook](https://github.com/goosh-gh/App-PDL-Notebook) — for browser-based notebook viewer (optional)

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
my $rec = read_nk('subjct.EEG');
# $rec->{data}      [n_ch, n_samples] float32 µV
# $rec->{fs}        sampling rate Hz (auto-detected)
# $rec->{labels}    channel names (from .21E; last row is PAD/STIM marker)
# $rec->{t_start}   "YYYY-MM-DD HH:MM:SS"
# $rec->{events}    [{t => $sec, label => $str}, ...]  (.LOG session annotations)
# $rec->{gains}     [n_ch] µV/bit
# $rec->{device}    "EEG-1100C V01.00" / "EEG-1200A V01.00"
# $rec->{n_blocks}  number of waveform blocks
# extblock also returns: layout, ch_hw_idx, stim_index, t_block_starts, n_samp_per_block

# Read a specific block (wfmblock files may hold several)
my $rec2 = read_nk('subjct.EEG', block => 2);
```

### Channel labels

Labels come from the `.21E` file. Names absent from `[ELECTRODE]` but present in
`[REFERENCE]` (the amplifier's reference derivations) are filled in as a
fallback, and Nihon Kohden's `$`-prefixed reference names are normalized to a
Perl/filename-safe, collision-free `_ref` suffix:

| in `.21E` | reported label | meaning |
|-----------|----------------|---------|
| `A1` (electrode) | `A1` | A1 ear electrode |
| `$A1` (reference) | `A1_ref` | A1 used as reference input |
| `$AV` | `AV_ref` | average reference |

The last channel is the hardware marker: `PAD` (wfmblock, zero-filled) or `STIM`
(extblock, raw marker codes).

### Stimulus triggers

Per-trial stimulus triggers are recorded as **TTL levels on the DC channels**,
not in `.LOG` (which holds only session annotations) and not on the `STIM`
marker channel. The DC jack numbering differs by format:

| Format | Trigger DC channels |
|--------|---------------------|
| `EEG-1100C` | `DC03`–`DC06` |
| `EEG-1200A` | `DC01`–`DC04` |

`read_nk()` passes the `.21E` names through unchanged, so each file reports the
correct DC labels.

### Interactive viewer (giza-server, native window)

```bash
perl -Ilib examples/read_nihonkohden.pl subject.EEG --plot
perl -Ilib examples/read_nihonkohden.pl subject.EEG --plot --block 2 --sec 10
```

Options:

- `--sec S` — seconds per screen (default 10)
- `--nch N` — number of channels (default 8; ignored if `--chans`)
- `--uv U` — µV per division for EEG traces (default 100)
- `--chans LIST` — comma-separated channel **names**, in order
  (e.g. `--chans Fp1,Cz,Pz,DC01,DC02,STIM`; any label, incl. DC/STIM/`*_ref`)
- `--aux MODE` — scaling for aux channels (`DC*` / `STIM` / `PAD` / `COM` /
  `*_ref` / `BN*` / `Pulse` / `CO2` …):
  - `same` (default) — draw as-is at the EEG scale; TTL may overlap neighbours
    (often useful for lining a trigger up against the EEG)
  - `auto` — auto-scale each aux channel to its own slot (TTL squares stay
    visible, no overlap)
  - `<N>` — give aux channels a fixed `N` µV/div (e.g. `--aux 2000`)

Horizontal slider: time scroll · Vertical slider: EEG gain (µV/div).
EEG traces are blue, aux channels red.

### Interactive viewer (App-PDL-Notebook, browser)

Full MNE `raw.plot()`-style viewer running inside a notebook cell.
No giza-server required. ~7 ms/frame with LTTB downsampling.

```perl
# In a notebook cell:
use lib '/path/to/PDL-EEG/lib';
use lib '/path/to/PDL-Graphics-Cairo/lib';
use PDL;
local @ARGV = ('subjct.EEG');          # omit for synthetic demo
do '/path/to/App-PDL-Notebook/examples/notebook_eeg_raw.pl';
```

Browser controls: Position (time scroll), Window (ms), Gain (10–1000 µV/div),
Ch offset (channel scroll), Neg-up toggle.

### Verify a file

```bash
perl -Ilib examples/verify_read.pl subjct.EEG
```

## Format notes

Confirmed from real hardware: EEG-1100C (JE-921A amplifier) and EEG-1290
(JE-92NX headbox, `EEG-1200A` format).

Common to both layouts:

- ADC: uint16 offset binary, center `0x8000`; `µV = (raw − 32768) × gain`
- Sampling rate: lower 14 bits of the u16 at `data_block+0x1A`
- Signature `0x0000`; control-block list at `0x0091`/`0x0092`; `ext_address` at `0x03EE`

`wfmblock` (EEG-1100C):

- Channel table: 10-byte entries at `wfmblock+0x2F`; data start `wfmblock+0x171`
- Gain: fixed `0.09765625 µV/bit`

`extblock` (EEG-1200A):

- Channels via chain `ext → ext+18 → +20`; count at `eb3+68`, indices at
  `eb3+72+i*10`; data start `eb3+72+(n_ch−1)*10`
- Sample-interleaved; `n_samples` computed from file size
- Per-channel gain: EEG/micro `0.09765624 µV/bit`, DC/other `0.36629984` (mV range);
  the trailing `STIM` marker channel is raw

Events: 6-digit ASCII seconds in `.LOG` (session annotations; Shift-JIS labels
supported). Per-trial triggers are on the DC channels (see above).

## License

Same terms as Perl itself (Artistic License 2.0 or GPL-1+).
Format knowledge derived from EDFbrowser (GPL-2) and Brainstorm's `in_fopen_nk.m`
used as reference only; this is a clean-room Perl implementation.

## See also

- [goosh-gh/PDL-Graphics-Cairo](https://github.com/goosh-gh/PDL-Graphics-Cairo)
- [goosh-gh/App-PDL-Notebook](https://github.com/goosh-gh/App-PDL-Notebook)
- [goosh-gh/giza-server](https://github.com/goosh-gh/giza-server)
- [goosh-gh/PDL-IO-PNG](https://github.com/goosh-gh/PDL-IO-PNG)
