# PDL::EEG — Nihon Kohden / EDF / BESA EEG toolkit

Read Nihon Kohden Neurofax recordings in PDL, resolve headbox-independent
trigger/channel labels, re-reference (incl. balanced non-cephalic), and export
to EDF/EDF+ or BESA ASCII multiplexed (`.mul`).

## Requirements

- Perl ≥ 5.36 and [PDL](https://pdl.perl.org/) (tested against PDL 2.085).
- On macOS/MacPorts, build Cocoa-dependent extras with
  `./configure CC=clang OBJC=clang PKG_CONFIG=/opt/local/bin/pkg-config`.
- The readers assume a **little-endian** host (Apple Silicon, x86-64, ARM64 all
  qualify); binary buffers are interpreted directly as native `ushort`/`short`.

## Modules

| Package | Role |
|---------|------|
| `PDL::EEG::IO::NihonKohden` | Reader for `.EEG` (EEG-1100C `wfmblock` + EEG-1200A `extblock`). Options: `all_blocks`, `block`, `label_map`. Returns `data [n_ch,n_samp]` µV, `fs`, `labels`, `t_start`, `events`, `gains`, `t_block_starts`, `gap_bounds`, `device`, `layout`, `system_reference`, `last_pattern`. |
| `PDL::EEG::IO::NihonKohden::PTN` | Parse Neurofax `.PTN` montage files (1100C + 1200A). |
| `PDL::EEG::IO::NihonKohden::Montage` | `.LOG` montage name + `.PTN` + signal → `label_map`; `resolve_labels`. |
| `PDL::EEG::IO::EDF` | `write_edf` (EDF / EDF+C) and `read_edf` (round-trips the `read_nk` contract); `clean_edf_label` normalises EDF+ signal labels. |
| `PDL::EEG::IO::BESA::ASCII` | `write_mul` — BESA ASCII multiplexed (`.mul`) export. |
| `PDL::EEG::Derivation` | `derive` (general linear derivation `y = M·x`), `bne` (balanced non-cephalic re-reference), `rereference` (single/linked/average). |
| `PDL::EEG::Signal` | Device-independent square-pulse / TTL detector. |

## Command-line tools

| Tool | Role |
|------|------|
| `examples/read_nihonkohden.pl` | Interactive viewer (`--block/--sec/--nch/--chans/--aux`, optional Cairo plot); dispatch is inside `read_nk`, so it needs no format knowledge |
| `examples/nk_to_edf.pl` | NK `.EEG` → EDF/EDF+ (`--subject`, `--equipment`, `--allblocks`) |
| `examples/nk_to_mul.pl` | NK `.EEG` → BESA `.mul` (`--cut`, `--cut-clock`, `--suffix`, `--bne`) |
| `examples/edf_to_mul.pl` | EDF → BESA `.mul` (`--chans`, `--cut`, `--cut-clock`, `--suffix`, `--bne`) |
| `xt/verify_read.pl` | Real-data (or synthetic) `read_nk` sanity check, independent of `make test` |
| `xt/smoke_bne.pl` | Author smoke test: `--bne` on a real `.EEG`/`.edf` |

## Quick start

```perl
use PDL::EEG::IO::NihonKohden qw(read_nk);
use PDL::EEG::IO::EDF         qw(write_edf);
use PDL::EEG::IO::BESA::ASCII qw(write_mul);
use PDL::EEG::Derivation      qw(bne);

my $rec = read_nk('JJ0090J6.EEG', all_blocks => 1);   # data[n_ch,n_samp] µV
write_edf($rec, 'out.edf');                            # EDF+C, events → annotations
write_mul($rec, 'out.mul');                            # BESA ASCII multiplexed

# balanced non-cephalic re-reference, then export
my $bn = bne($rec, prop => 0.5, suffix => '-BN');      # y = x − (p·BN1 + (1−p)·BN2)
write_mul($bn, 'out_bne.mul');
```

### BESA `.mul` export (CLI)

```
perl -Ilib examples/nk_to_mul.pl  JJ0090J6.EEG
perl -Ilib examples/edf_to_mul.pl JJ0090J6.edf --suffix -BN
perl -Ilib examples/nk_to_mul.pl  JJ0090J6.EEG --cut "21-376:b0b1_21_376"
perl -Ilib examples/nk_to_mul.pl  JJ0090J6.EEG --bne          # re-reference to BNE
```

- `--cut a-b[:name],…` writes one `.mul` per range in data-coordinate seconds;
  `--cut-clock HH:MM:SS-HH:MM:SS[:name]` uses wall-clock (continuous recordings).
- EDF+ labels are cleaned on the way in (`EEG Fp1-Ref` → `Fp1`, `POL DC01` →
  `DC01`, `$A1` → `A1_ref`), so the `.mul` label row is whitespace-free and its
  token count matches `Channels=`.
- The dedicated **Trigger** channel is written as a column but **not counted in
  `Channels=`** (matching the vendor export; pass `count_trigger => 1` to
  include it).

### Re-referencing / balanced non-cephalic (BNE)

Nihon Kohden acquires against a system reference (`Avr(C3,C4)`; see
`$rec->{system_reference}`), so a recorded channel is `x_i = s_i − s_ref`.
Re-referencing to `r = p·BN1 + (1−p)·BN2` gives `y_i = x_i − (p·BN1 + (1−p)·BN2)`;
because the weights sum to 1, the acquisition reference **cancels exactly** and
need not be known. `bne()` auto-detects BN1/BN2, drops them from the output,
passes DC/Trigger through unchanged, and tags re-referenced channels `-BN`.

`--bne` on the CLIs is **off by default** (data written as recorded). When used,
provenance is recorded in the `.mul` header as `SegmentName=BNE_prop<value>`
(a standard BESA field). Default `prop = 0.5` — the BN balance is applied in
analog hardware at electrode application, so 0.5 means "no extra digital
re-balance"; other values apply one.

## Trigger / channel-label resolution (headbox-independent)

Trigger/DC channel names are **not** derivable from the recording format alone:

- The same trigger line is `DC03–06` on the 1100C headbox and `DC01–04` on the
  1200A; a fixed-name search is a landmine.
- The authoritative display names live in the **montage** (`.PTN`), which labels
  the four TTL lines `TrigBit0/2/4/8`; the electrode table calls them `DCxx`.
- **Which recorded `ch_idx` carries a trigger is only visible in the signal** —
  the `.PTN` gives the count/names but stores `G1=0`, not the channel index.

`resolve_labels` combines all three:

```
.LOG  ──montage_from_log──▶ "IIA"
.PTN  ──parse_ptn────────▶ trigger names [TrigBit0,2,4,8] (count = 4)
.EEG  ──detect_square_pulses(n=4)──▶ ch_idx that actually pulse (needs all_blocks)
        zip names(montage order) ⟷ triggers(ch_idx order)
              → label_map { ch_idx => name } → read_nk(label_map => …)
```

```perl
use PDL::EEG::IO::NihonKohden::Montage qw(resolve_labels);
my $r = resolve_labels($rec, ptn_dir => 'YJ0394VB.PTN');
# $r->{montage} "IIA"; $r->{label_map} { 45=>'TrigBit0', … }
my $rec2 = read_nk($f, all_blocks=>1, label_map => $r->{label_map});
```

`resolve_labels` is an API in `PDL::EEG::IO::NihonKohden::Montage` (there is no
dedicated CLI). Pass `names => [qw(DC03 DC04 DC05 DC06)]` to use physical box
labels instead of the montage's `TrigBit*` names, or pin `label_map` by hand.

## File-format reference

`docs/nihon_kohden_files.md` documents every file in a Neurofax recording
bundle (`.EEG/.21E/.LOG/.CN3/.PTN/.bam/…`) and what each carries, including
where the system reference and per-segment display montage live.

## Tests

```
make test        # 12 files (t/06 reserved/skipped), ~185 subtests
```

`t/01_nihonkohden` `t/02_edf` (write + read_edf round-trip) `t/03_ptn`
`t/04_signal` `t/05_montage` `t/06_reserved` (skip: reserved) `t/07_blocks`
`t/08_epoch` `t/09_i18n` `t/10_besa_ascii` `t/11_edf_to_mul` `t/12_derivation`.

`xt/smoke_bne.pl FILE.EEG|FILE.edf` is an author test that runs a converter with
and without `--bne` on a real recording and checks structural invariants.

## Performance

Binary I/O interprets buffers directly through PDL's data pointer rather than
building multi-million-element Perl lists: the `wfmblock`/`extblock` readers and
`read_edf`/`write_edf` all use `get_dataref` byte copies, and `write_mul`
formats each flushed block with a single `sprintf`. These keep large-recording
conversion memory-bounded and several times faster than the naïve approach.

## Honest caveats

- **Full recording required for trigger detection.** Triggers must fire to the
  rail to separate from EEG; run `all_blocks=>1`. A line that never toggles in
  the window won't be detected.
- **Name↔channel order is an assumption.** Montage trigger names are zipped onto
  detected triggers sorted by ascending `ch_idx`; verify once per headbox.
  `label_map` overrides.
- **The BNE balance ratio is not stored in the recording.** It is applied in
  analog hardware (a potentiometer at electrode application), so no file in the
  bundle carries a digital ratio; `bne()` therefore takes `prop` as a parameter
  (default 0.5). Verified exhaustively across `.21E/.LOG/.CN3/.11D/.bam/.EEG`
  header and all 36 `.PTN` montages.
- **Per-segment display montage lives in `.CN3`**; the *recording* montage name
  is in `.LOG`/`.21E [LASTPATTERN]`. The *export-time* review montage is not
  recoverable from the files.
- **Absolute DC µV is not calibrated.** Trigger *edges* (hence ERP epoching) are
  gain-independent, so this does not affect event extraction.
- **`read_edf` assumes one sample rate.** All non-annotation signals must share a
  single rate; the EDF+ annotation channel is parsed into `events` and excluded
  from `data`. A mixed-rate EDF is read with a `carp` warning, using the first
  signal's rate. EDF permits per-signal rates, but multi-rate reads are not
  implemented. (This affects only third-party EDFs; files written by `write_edf`
  are always single-rate.)
