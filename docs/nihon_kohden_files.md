# Nihon Kohden EEG-1200A recording bundle — file reference

Files written by a Neurofax EEG-1200A (headbox JE-208A / JE-92NX) for one
recording, e.g. `subject.*`. Basename is shared; each extension carries a
different part of the session. Verified against the `subject` MMN recording
(2026-07-02). "★" marks files that carry reference / montage information.

## Core signal & metadata

| File | Kind | Contents |
|------|------|----------|
| `.EEG` | binary | The waveform. EEG-1100C uses `wfmblock`, EEG-1200A uses `extblock`. Header holds the format signature (`EEG-1200A V01.00`), the linked `.PNT`, the acquisition start timestamp (`YYYYMMDDhhmmssn`), the headbox model (`JE-208A`), and per-segment start times. Read by `read_nk`. |
| `.pnt` / `.PNT` | binary | Patient / study identification (patient ID, study name e.g. `MMN Trigger`, dates, protocol comments). |
| `.21E` | text (Shift_JIS, CRLF) | ★ Electrode & reference table. Sections: `[ELECTRODE]` (`0000=Fp1` … including `0020=BN1`, `0021=BN2`, `0037=BN`, `0038=AV`), `[REFERENCE]` (`$`-prefixed reference pseudo-electrodes `$BN`, `$AV`, `$Cz`, `$A1`, `$A2` …), `[SD_DEF]` (standard-derivation weights — **empty** here), `[SYSTEM_SETUP]` (**`SystemReference=C3,C4`**, `DeviceName=<JE-92NX>`), `[LASTPATTERN]` (**`PATTERN=36`** = recording-time display montage, `REFERENCE=-1`). Parsed by `_read_21e`; `read_nk` now also exposes `system_reference` and `last_pattern`. |

## Events, log & montage timeline

| File | Kind | Contents |
|------|------|----------|
| `.LOG` | binary | ★ Event log. Per-segment records (45-byte entries): task markers (`task1`…`task5`, `安静開眼`/`安静閉眼`, CP932), `REC START MMN EEG/CAL`, `Recording Gap`, and a reference-state note **`A1+A2 OFF`** re-stated at every segment start. Times are 6-digit ASCII seconds from recording start; epoch number at byte 42. Read by `_read_log`. |
| `.LO2` | XML | `KohdenNeurology` secondary log: acquisition start/end times, (empty) `EventArea`. |
| `.EVT` | text | Trigger/event table; header only (`Tmu Code TriNo`) — empty in this session. |
| `.CN3` | binary | ★ Per-segment **display-montage** records (one block per recording segment; scalp channels stored as electrode indices into `.21E`, special channels — `vEOG`, `hEOG`, `TrigBit0/2/4/8`, `rVEOG`, `nose` — inline) plus segment start times. This is where a **mid-recording montage change would be captured**. (Older EEG-1100C used `.CN2`; EEG-1200A uses `.CN3`.) |
| `.CMT` | binary | Comment file (references `.CN3` and the device signature). |

## Montage library & display configuration

| File | Kind | Contents |
|------|------|----------|
| `.PTN/Pattern_0NN.PTN` (×36) | binary | ★ Montage-definition library. Header `EEG-1000/9000 Pattern Info File`; montage name at offset 0x80 (`IA`, `IIA`, … `VIIID`, `FREE A`–`FREE D`, `MMN`, `PSG-D`, `PSG-T`, `MSLT`); reference groups (`AV1`…`AV16`, equal-weight 0/1 inclusion masks); 80-byte channel records from 0x410. `[LASTPATTERN] PATTERN=36` selects `Pattern_036.PTN` = `FREE D`. |
| `.11D` | text | ★ DC-channel calibration (`DCxx_Coefficient/BaseLine/Offset/Pitch/Unit/Enable`) and `[DisplayMontageComent]` — all 36 montages (`PatternA1`…`PatternDF`) as per-electrode `M`/`C` maps (`C` = electrodes forming that montage's common/reference set; **0/1 flags, no weights**). |
| `.11S` | binary | Waveform-review (`TRACE`) data. |
| `.bam` | binary | ★ Per-channel applied reference — lists **`Avr(C3,C4)`** for every channel, confirming acquisition against the average of C3 and C4. |

## Trend, QC & misc

| File | Kind | Contents |
|------|------|----------|
| `.BFT` | binary | Band/trend index: per-segment start/end timestamps + signature. |
| `.TRD/Trend.config` | XML | Trend & DSA (density spectral array) **display** config: heart-rate trend items, DSA groups (`DSA C3,C4`, `DSA All`, …) with channel masks, `VoltMax`, colour tables. No reference definition. |
| `.PATCHINFO` | XML | Network / data-loss recovery info (`SenderInfo`, GUID, IP) — despite the name, **not** electrode patching. |
| `.MISC/AccuracyCheck/*/…` | XML + PDF | Periodic accuracy / impedance / noise / sensitivity / filter QC report (per channel). Restates `<SystemReference>C3-C4</SystemReference>`; no BN/ratio. |
| `DskUUID.vol` | binary | Disk volume identifier. |

## Reference & the "-BN" suffix

- **System (acquisition) reference is `Avr(C3,C4)`** — stated in `.21E [SYSTEM_SETUP]` and confirmed per-channel in `.bam`. Recorded samples are `x_i = s_i − (s_C3+s_C4)/2`.
- **`$BN`** (balanced non-cephalic) is a *named* reference (`.21E [REFERENCE]`) built from the body electrodes **BN1** (vertebral, "V") and **BN2** (sternal, "S"). A display/export referenced to `$BN` renders each channel as `«electrode»-BN` (hence `Fp1-BN` in a BESA `.mul` export).
- **The BN1/BN2 combination ratio is not stored anywhere in the recording bundle** — checked exhaustively: `.21E` (incl. empty `[SD_DEF]`), `.LOG`, `.CN3`, `.11D`, `.bam`, `.EEG` header, and all 36 `.PTN` montages. No montage references BN1+BN2, and reference masks are equal-weight (0/1). The balance is applied in the **analog domain** (a potentiometer at electrode application, per the classic Stephenson–Gibbs method); digitally `$BN` is a simple combination of the already-balanced BN1/BN2, so there is no digital ratio to save.

## What `read_nk` exposes (relevant subset)

`data` `[n_ch, n_samp]` µV · `fs` · `labels` · `t_start` · `events` · `gains` ·
`gap_bounds` · `t_block_starts` · `device` · `layout` ·
**`system_reference`** (e.g. `"C3,C4"`) · **`last_pattern`** (e.g. `36`).

## Downstream re-reference (`PDL::EEG::Derivation`)

Because the acquisition reference cancels whenever the re-reference weights sum
to 1, BNE can be computed directly from the recorded data:

    y_i = x_i − (prop·BN1 + (1−prop)·BN2)      # Avr(C3,C4) cancels

`prop = 0.5` reproduces the "physical balance only" case; other values apply an
additional digital re-balance. `bne()` auto-detects BN1/BN2 and passes DC /
Trigger channels through unchanged.
