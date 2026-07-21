# PDL::EEG — Nihon Kohden / EDF / BESA EEG toolkit

Read Nihon Kohden Neurofax recordings in PDL, resolve headbox-independent
trigger/channel labels, re-reference (incl. balanced non-cephalic), and export
to EDF/EDF+ or BESA ASCII multiplexed (`.mul`).

## Requirements

- Perl ≥ 5.36 and [PDL](https://pdl.perl.org/) (tested against PDL 2.085+).
- On macOS/MacPorts, build Cocoa-dependent extras with
  `./configure CC=clang OBJC=clang PKG_CONFIG=/opt/local/bin/pkg-config`.
- The readers assume a **little-endian** host (Apple Silicon, x86-64, ARM64 all
  qualify); binary buffers are interpreted directly as native `ushort`/`short`.

## Modules

| Package | Role |
|---------|------|
| `PDL::EEG::IO::NihonKohden` | Reader for `.EEG` (EEG-1100 `wfmblock` + EEG-1200 `extblock`, incl. multi-segment recordings). Options: `all_blocks`, `block`, `label_map`, `dc_base`. Returns `data [n_ch,n_samp]` **µV (all channels, DC included)**, `fs`, `labels`, `units` (per-channel export dimension `uV`/`mV`/`code`), `t_start`, `events`, `gains` (µV/bit), `n_samp_per_block`, `block_meta`, `t_block_starts`, `gap_bounds`, `device`, `layout`, `system_reference`, `last_pattern`. |
| `PDL::EEG::IO::NihonKohden::PTN` | Parse Neurofax `.PTN` montage files (1100C + 1200A). |
| `PDL::EEG::IO::NihonKohden::Montage` | `.LOG` montage name + `.PTN` + signal → `label_map`; `resolve_labels`. |
| `PDL::EEG::IO::EDF` | `write_edf` (EDF / EDF+C) and `read_edf` (round-trips the `read_nk` contract); `clean_edf_label` normalises EDF+ signal labels. |
| `PDL::EEG::IO::BESA::ASCII` | `write_mul` — BESA ASCII multiplexed (`.mul`) export. |
| `PDL::EEG::IO::ASA` | `read_elc` — read ASA electrode-position files (`.elc`). Returns `coords [3,N]` (native unit, MNI mm), parallel `labels`, name→xyz `pos`, `unit`/`reference`, and auto-detected `fiducials` (LPA/RPA/Nz). Robust to indented blocks/CRLF; coordinates parsed vectorised. `parse_ELEC_POS3D_ASA_4AdventCalendar` is a drop-in shim for the PDL Advent Calendar 2024 (Day 12) parser. |
| `PDL::EEG::Derivation` | `derive` (general linear derivation `y = M·x`), `bne` (balanced non-cephalic re-reference), `rereference` (single/linked/average). |
| `PDL::EEG::Signal` | Device-independent square-pulse / TTL detector. |

## Command-line tools

| Tool | Role |
|------|------|
| `examples/read_nihonkohden.pl` | Interactive viewer (`--block/--sec/--nch/--chans/--aux`, optional Cairo plot); dispatch is inside `read_nk`, so it needs no format knowledge |
| `examples/nk_to_edf.pl` | NK `.EEG` → EDF/EDF+ (`--subject`, `--equipment`, `--allblocks`) |
| `examples/nk_to_mul.pl` | NK `.EEG` → BESA `.mul` (`--cut`, `--cut-clock`, `--suffix`, `--bne`) |
| `examples/edf_to_mul.pl` | EDF → BESA `.mul` (`--chans`, `--cut`, `--cut-clock`, `--suffix`, `--bne`) |
| `examples/mul_to_nk.pl` | Diff a vendor `.mul` against `read_nk` (round-trip check); `--solve-bne` recovers the BN balance from the vendor's own export |
| `examples/find_bn_balance.pl`, `examples/find_bn_diff.pl` | Search NK header files for where a known BN balance is stored (investigative; see caveats) |
| `xt/verify_read.pl` | Real-data (or synthetic) `read_nk` sanity check, independent of `make test` |
| `xt/smoke_bne.pl` | Author smoke test: `--bne` on a real `.EEG`/`.edf` |
| `examples/show_electrodes_3d.pl` | 3D scalp-electrode viewer over `read_elc` (GS3D or TriD backend); ships a 28-point fixture. Needs `PDL::Graphics::Cairo` (GS3D) or `PDL::Graphics::TriD`. |
| `examples/dump_nyhead19.pl` | Extract New York Head 19ch + fiducials from `sa_nyhead.mat` (`/sa/locs_3D_orig`) into `nyhead19.txt`; built-in fiducial sanity check. Needs `PDL::IO::HDF5` + the NY Head `.mat`. |
| `examples/overlay_nyhead.pl` | Overlay `standard_1020.elc` onto NY Head 19ch via `read_elc`: raw residual + fiducial-frame-aligned residual (mm) + worst-channel, writes `electrodes_overlay.xyz`. `--selftest` validates the alignment math. |
| `examples/show_overlay_3d.pl` | GS3D 3D overlay of the two electrode sets with per-electrode displacement segments (left labels from `.elc`, right/mid from NY, an L/R & A/P sanity check); `--obj` exports a Blender-ready `.obj`+`.mtl` (octahedron markers + materials). |
| `examples/overlay_scalp_obj.pl` | Overlay `.elc` electrodes onto a NY Head **surface** and export one Blender/MeshLab `.obj`+`.mtl`. `--surf` selects the mesh (`/sa/head` scalp, `/sa/cortex75K` cortex); electrodes drop on unaligned (same MNI frame). Each electrode is its own named object with an optional outward-facing 3D **text label** (`--labels`/`--no-labels`) — an L/R check readable even in Finder preview. `--stats` reports electrode→nearest-vertex distance; `--selftest` needs no PDL or data. Needs `PDL::IO::HDF5` + the NY Head `.mat`. |
| `xt/70_real_data.t` | Real-data event-placement regression (`extblock` + `wfmblock`); pass `.EEG` paths after `::` |

## Quick start

```perl
use PDL::EEG::IO::NihonKohden qw(read_nk);
use PDL::EEG::IO::EDF         qw(write_edf);
use PDL::EEG::IO::BESA::ASCII qw(write_mul);
use PDL::EEG::Derivation      qw(bne);

my $rec = read_nk('subject.EEG', all_blocks => 1);   # data[n_ch,n_samp] µV
write_edf($rec, 'out.edf');                            # EDF+C, events → annotations
write_mul($rec, 'out.mul');                            # BESA ASCII multiplexed

# balanced non-cephalic re-reference, then export.
# prop is REQUIRED: the BN balance is a hardware setting, not stored in the file.
# Measure it once with examples/mul_to_nk.pl --solve-bne, or read it off the amp.
my $bn = bne($rec, prop => 0.71, suffix => '-BN');     # y = x − (p·BN1 + (1−p)·BN2)
write_mul($bn, 'out_bne.mul');
```

### BESA `.mul` export (CLI)

```
perl -Ilib examples/nk_to_mul.pl  subject.EEG
perl -Ilib examples/edf_to_mul.pl subject.edf --suffix -BN
perl -Ilib examples/nk_to_mul.pl  subject.EEG --cut "21-376:b0b1_21_376"
perl -Ilib examples/nk_to_mul.pl  subject.EEG --bne          # re-reference to BNE
```

- `--cut a-b[:name],…` writes one `.mul` per range in data-coordinate seconds;
  `--cut-clock HH:MM:SS-HH:MM:SS[:name]` uses wall-clock, mapped to samples
  through `block_meta` (piecewise, break-aware): a time that lands in a
  recording gap clamps to the last real sample before it, and a range never
  leaks the next block's data across a break. For `wfmblock` files add
  `--allblocks` so `block_meta` spans every segment.
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
(a standard BESA field).

**`prop` is required — there is no safe default.** The BN balance is set on the
amplifier (a front-panel value the operator dials in at recording time), and it
is **not written to any file in the bundle** (see caveats). Two machines here
measured **0.71** and **0.64** (logged as 0.65), confirming it is per-machine /
per-session. An earlier version of this toolkit defaulted to `0.5`; that value
was never correct for a real recording and only looked harmless because
`BN1 ≈ BN2` in calibration segments. If you do not know the balance, recover it
from a vendor `.mul` export (next section).

#### Recovering the balance from a vendor `.mul` (`--solve-bne`)

The Nihon Kohden viewer's own `.mul` export is already BN re-referenced. If you
have one, `examples/mul_to_nk.pl` measures the balance the recorder actually
used:

```
perl -Ilib examples/mul_to_nk.pl vendor.m01 --eeg subject.EEG --solve-bne
```

It aligns the `.mul` against `read_nk(all_blocks=>1)` (the `.mul` is a
hand-selected range, so the offset is found by search, not assumed), then
regresses `raw − mul` onto `BN1`/`BN2`. Because that residual is one common
signal on every scalp channel — a reference difference — the fit is exact:
weights that sum to 1 (confirming the model) with a residual at the ADC step.
The recovered `prop` cross-checks against `|p−0.5|·rms(BN1−BN2)` to sub-percent.

The tool is also a general **round-trip check**: matching the vendor export
channel-for-channel is independent confirmation that block boundaries, channel
order and gains are correct — including across recording breaks, which nothing
in this toolkit could otherwise self-verify.

## Trigger / channel-label resolution (headbox-independent)

Trigger/DC channel names are **not** derivable from the recording format alone:

- The same trigger line is `DC03–06` on the EEG-1100 family and `DC01–04` on the
  EEG-1200 family; a fixed-name search is a landmine. `read_nk` keys the default
  DC numbering on the **format signature** at offset 0 (not the on-disk layout,
  and not the enclosing directory name — `NKT/EEG2100/` is a folder, the
  signature is `EEG-1200A V01.00`). Signatures outside the 1100/1200 families
  have no assumed numbering: `read_nk` **croaks** rather than mislabel a trigger,
  unless a `.21e` names the channels or you pass `dc_base => 1|3`.
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
my $r = resolve_labels($rec, ptn_dir => 'subject.PTN');
# $r->{montage} "IIA"; $r->{label_map} { 45=>'TrigBit0', … }
my $rec2 = read_nk($f, all_blocks=>1, label_map => $r->{label_map});
```

`resolve_labels` is an API in `PDL::EEG::IO::NihonKohden::Montage` (there is no
dedicated CLI). Pass `names => [qw(DC03 DC04 DC05 DC06)]` to use physical box
labels instead of the montage's `TrigBit*` names, or pin `label_map` by hand.

## Multi-segment recordings & recording breaks

EEG-1200 `extblock` recordings are **not one continuous stream**. At every
recording break the recorder re-emits a 442-byte channel-info block
(`72 + (n_ch−1)·10` bytes) into the sample stream, and the gaps between segments
are real. `read_nk` detects these embedded headers, treats each span as its own
block, and reports per-segment geometry:

```perl
my $rec = read_nk($f, all_blocks => 1);
$rec->{n_samp_per_block};   # [205000, 176000, 30000, …]
$rec->{block_meta};         # [{ index, start_samp, n_samp, t_start }, …]
```

The viewer marks each break with the **real elapsed time skipped**
(`epoch(t_start[b+1]) − epoch(t_start[b]) − n_samp[b]/fs`), e.g. `▲ 46.0s
skipped`, and loads only the segments a `--cut` range touches.

`.LOG` events are placed at their true data-sample position (`{samp}`/`{t_data}`)
against these exact segment boundaries. The `.LOG` elapsed-seconds field also
counts the paused time between recording blocks, so it cannot be treated as a
wall-clock offset into the data; events are instead anchored to the `REC START`
markers that open each segment (every `REC START` lands exactly on a block
boundary) and offset within the segment, where wall-clock and data advance 1:1.
This holds for both `extblock` and, for `all_blocks` reads, `wfmblock`.

`block_extents($file)` returns the same per-segment table without reading sample
data, for quick inspection:

```
perl -Ilib -MPDL -MPDL::EEG::IO::NihonKohden=block_extents -e '
  my $e = block_extents($ARGV[0]);
  printf "%d segments: %s\n", scalar @$e,
    join(", ", map { sprintf("%.1fs", $_->{n_samp}/$_->{fs}) } @$e);
' subject.EEG
```

> Data written from a multi-segment 1200-family file by **any earlier version**
> of this toolkit is wrong past the first break and must be regenerated.

## File-format reference

`docs/nihon_kohden_files.md` documents every file in a Neurofax recording
bundle (`.EEG/.21E/.LOG/.CN3/.PTN/.bam/…`) and what each carries, including
where the system reference and per-segment display montage live.

## Electrode positions & 3D (ASA `.elc`)

`PDL::EEG::IO::ASA::read_elc` reads ASA electrode files (e.g. mne-python's
`standard_1020.elc`) into a `(3,N)` coordinate piddle plus a name→xyz lookup and
detected fiducials. Five optional examples build on it, covering single-montage
3D display and coregistration against the New York Head forward model:

```
# view one montage in 3D (colour by hemisphere)
perl -Ilib examples/show_electrodes_3d.pl --elc standard_1020.elc --backend gs3d

# NY Head coregistration: dump 19ch from the .mat, overlay, then view / export
perl -Ilib examples/dump_nyhead19.pl   sa_nyhead.mat nyhead19.txt
perl -Ilib examples/overlay_nyhead.pl  --elc standard_1020.elc --ny nyhead19.txt
perl -I<P:G:C>/lib examples/show_overlay_3d.pl --obj overlay.obj   # 3D + Blender .obj

# overlay electrodes on a NY Head surface (scalp or cortex) → Blender/MeshLab .obj
perl -Ilib examples/overlay_scalp_obj.pl --elc standard_1020.elc \
    --mat sa_nyhead.mat --surf /sa/head      --out nyhead_scalp.obj
perl -Ilib examples/overlay_scalp_obj.pl --elc standard_1020.elc \
    --mat sa_nyhead.mat --surf /sa/cortex75K --out nyhead_cortex.obj --no-stats
```

`standard_1020.elc` and the NY Head 19ch are both MNI mm on the same axis
convention, so their **raw (unaligned) residual is already small** (~5 mm mean,
no channel above ~11 mm) — direct confirmation that the electrode correspondence
is correct, with no fitting. `overlay_nyhead.pl --selftest` validates the
fiducial-frame alignment independently (a known transform recovers to 0 mm). The
3D tools need `PDL::Graphics::Cairo` (GS3D), and the NY Head dump needs
`PDL::IO::HDF5`; the core `PDL::EEG::IO::ASA` reader needs only PDL.

`overlay_scalp_obj.pl` drops the electrodes straight onto a surface mesh
(`<group>/vc`+`/tri`) with no alignment, so the mm offset you see is the true
electrode-to-surface fit. In `sa_nyhead.mat` only `/sa/head` (1082 verts) and
`/sa/cortex75K` (74 382 verts) carry vertices — the lower-resolution `cortexNK`
groups are faces-only. Each electrode becomes its own named object carrying an
outward-facing 3D text label (`--no-labels` to omit), so a left/right swap is
obvious in any viewer, Finder Quick Look included. (h5ls reports these datasets
transposed: `{3,N}` on disk is `(N,3)` in PDL.)

## Tests

```
make test        # 14 files (t/06 reserved/skipped), 340 subtests
```

`t/01_nihonkohden` `t/02_edf` `t/03_ptn` `t/04_signal` `t/05_montage`
`t/06_reserved` (skip: reserved) `t/07_blocks` (block_extents + multi-segment
extblock regression) `t/08_epoch` (event placement + wall-clock→sample mapping)
`t/09_i18n` `t/10_besa_ascii` `t/11_edf_to_mul` `t/12_derivation`
`t/13_edf_roundtrip` (µV round-trip incl. DC, per-signal EDF dimension)
`t/14_asa` (ASA `.elc` reader: parse, fiducials, name lookup, Advent shim,
against a 28-point real-coordinate fixture).

`xt/70_real_data.t` is a real-data regression (not part of `make test`; needs
private recordings). Pass `.EEG` paths and it checks event placement on real
`extblock` and `wfmblock` files — every `REC START` on a block boundary, every
event inside a segment:

```
prove -lv xt/70_real_data.t :: /path/A.EEG /path/B.EEG
```

Test fixtures are generated by `perl t/mk_synthetic_nk.pl`; `--long[=SEC]` also
writes larger scrollable files (`t/data/*_long.eeg`, git-ignored — not fixtures).

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
- **The BN balance is not stored in any recording file.** It is a value the
  operator sets on the amplifier, and the bundle records only the *choice* to
  reference to BN (`.21E [REFERENCE] = $BN`), never the ratio. Searched
  exhaustively — value scan, BCD, integer-percent, per-mil, and raw byte-diff of
  two recordings with known-different balances — across `.EEG/.21E/.PNT/.LOG`
  headers, with no field found. `bne()` therefore **requires** `prop`. The one
  reliable way to recover it after the fact is `mul_to_nk.pl --solve-bne`
  against a vendor `.mul`; failing that, read it off the amplifier or your notes.
  (`examples/find_bn_balance.pl` / `find_bn_diff.pl` are the search tools, kept
  for when a controlled two-recording diff — same machine, balance changed —
  becomes available.)
- **Per-segment display montage lives in `.CN3`**; the *recording* montage name
  is in `.LOG`/`.21E [LASTPATTERN]`. The *export-time* review montage is not
  recoverable from the files.
- **DC channels are calibrated in µV** (366.30 µV/bit, i.e. the ±12 V input
  range; confirmed against the vendor `.mul`, whose DC columns are integer
  multiples of 366.22 µV). `read_nk` returns **every** channel in µV, DC
  included. Because a ±12 V DC line is ±12 002 913 µV and EDF's `physical_min`
  field is only 8 characters, `write_edf` gives each signal its own physical
  dimension — EEG in `uV`, DC in `mV` — and `read_edf` normalises back to µV.
  BESA `.mul` has a single `Bins/uV`, so DC there is written in µV at full
  magnitude; pass `exclude => [grep /^DC/]` if you only want the EEG scaled
  sensibly.
- **`read_edf` assumes one sample rate.** All non-annotation signals must share a
  single rate; the EDF+ annotation channel is parsed into `events` and excluded
  from `data`. A mixed-rate EDF is read with a `carp` warning, using the first
  signal's rate. EDF permits per-signal rates, but multi-rate reads are not
  implemented. (This affects only third-party EDFs; files written by `write_edf`
  are always single-rate.)
