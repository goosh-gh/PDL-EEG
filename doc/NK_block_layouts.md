# Nihon Kohden `.EEG` block layouts: `wfmblock` vs `extblock`

Reference for how `PDL::EEG::IO::NihonKohden` locates channels and samples in a
Nihon Kohden `.EEG` file. Every offset below is what the reader actually uses and
has been checked against real recordings and against MNE-Python's
`read_raw_nihon`.

## What the two layouts are

A `.EEG` file stores its waveforms in one of two internal layouts. The names are
the reader's own, chosen after the structure that carries the channel table:

- **`wfmblock`** (waveform block) — the legacy layout. The channel table lives
  *inside* the waveform block, right before the sample stream. Used by the
  EEG-1100 family (and QI-403A, EEG-2100, DAE-2100D).
- **`extblock`** (extended block) — the newer layout. The channel table lives in
  a separate extended-block chain reached through a pointer near the file head;
  the sample stream can also be split into several time segments. Used by the
  EEG-1200 family.

The layout is a property of the **format signature written in the file, not of
the physical recorder model.** An EEG-1290 recorder writes the signature
`EEG-1200A V01.00` (an `extblock` file); an EEG-1214 can write
`EEG-1100C V01.00` (a `wfmblock` file). Dispatch therefore keys on the signature,
never on the model name.

## Layout dispatch (`nk_layout`)

The 16-character signature at file offset `0x0000` decides the layout, in this
order:

1. Exact match in `%FORMAT_LAYOUT` → authoritative (`how = table`).
2. Structural fallback: `ext_address` (`u32LE @ 0x03EE`) is non-zero → `extblock`,
   zero → `wfmblock`. This is a verified discriminant (1100C = 0, 1200A ≠ 0).
3. Name family (`EEG-11*` → `wfmblock`, `EEG-12*` → `extblock`) is used only to
   cross-check the structural result and warn on a mismatch.
4. A signature that does not begin `EEG-`/`QI-`/`DAE-` → `undef` (the reader
   stops rather than risk misreading a non-NK file).

`nk_layout` returns `($signature, $layout, $how)`; `$how` records which rule
fired.

## Common header (both layouts)

- Signature: 16 bytes at `0x0000`.
- `ext_address`: `u32LE @ 0x03EE` (0 ⇒ `wfmblock`, non-zero ⇒ `extblock`).
- Control block: count `u8 @ 0x0091`; address `u32LE @ (0x0092 + i*20)`.
- First data/waveform block address: `u32LE @ (control_addr + 18)`.

## `wfmblock` (EEG-1100 family)

The channel table and samples are inside the waveform block. Relative to the
block address:

| offset        | meaning                                                        |
|---------------|----------------------------------------------------------------|
| `+0x00`       | `0x01` block type                                              |
| `+0x01..+0x10`| ASCII time string (`"TIMEhhmmss..."`, 16 bytes)                |
| `+0x14..+0x19`| BCD timestamp `YY MM DD HH MM SS`                              |
| `+0x1A..+0x1B`| `u16LE`, low 14 bits = sampling rate                           |
| `+0x26`       | `u8` = `n_ch_entries` (valid EEG/physio channels)             |
| `+0x27..`     | channel table: `n_ch_entries` records of 10 bytes each         |
| `+0x27 + n_ch_entries*10` | first data sample (e.g. `+0x171` for 33 channels)  |

The sample stream has `n_ch_entries + 1` channels (a trailing zero-pad channel).

### Channel-table record and the electrode code

Each 10-byte record begins with the **electrode code as its first byte**; the
last two bytes of the record are the fixed footer `0x10 0x05`:

```
+0  ch_code   (1 byte)   <- electrode code
+1 .. +7      (payload)
+8  0x10                 (trailing footer, NOT a header)
+9  0x05
```

So the code for channel `i` is `u8 @ (block + 0x27 + i*10)`, and the sample
stream begins immediately after the last record at `block + 0x27 +
n_ch_entries*10`. The documented data start `+0x171` only reconciles with a
record start of `+0x27` (`0x27 + 33*10 = 0x171`), which pins the record start.

### Fixed off-by-one (labels only)

Earlier code anchored the channel table on the `0x10 0x05` pair as if it were a
*leading* header and read the code from `+0x2F + i*10 + 2`, which equals
`+0x27 + (i+1)*10` — the first byte of the **next** record. Combined with the
`.21E[code-1]` lookup this cancelled for the contiguous 10-20 codes `1..19`, so
those 19 channels looked correct, but every auxiliary channel (where the codes
jump, e.g. `…20, 23, 24, 25, 26, 27, 45…`) landed on the neighbouring
electrode's name: `E` dropped, `X1` shown as `EtCO2`, `DC06` dropped, a spurious
`DC32` appearing, and so on.

Only the *labels* were wrong; the sample columns were always read correctly, so
existing 19-channel 10-20 ERP/LORETA analyses are unaffected. The fix reads the
code from the record start:

```perl
push @ch_indices, _read_u8($fh, $addr + 0x27 + $i * 10) + 1;
```

## `extblock` (EEG-1200 family)

> Deep dive (eb-chain arithmetic, ADC gains, segment headers, worked example):
> see [`NK_extblock.md`](NK_extblock.md).

The channel table lives in the extended-block chain, and the sample stream may
contain several time segments.

Reaching the channel-info block:

- `ext_address = u32LE @ 0x03EE` (non-zero).
- `eb2 = u32LE @ (ext_address + 18)`.
- `eb3 = u32LE @ (eb2 + 20)` — start of the channel-info block.

Inside the channel-info block:

| offset (from `eb3`) | meaning                                              |
|---------------------|------------------------------------------------------|
| `+68`               | `u16LE` = number of valid channels; `n_ch = that + 1` (trailing STIM) |
| `+72..`             | channel table: stride 10 bytes                       |
| `+72 + i*10`        | `u16LE` = electrode code for channel `i`             |

The channel-info block is `hdr_len = 72 + (n_ch-1)*10` bytes long (e.g. 442), and
the first sample follows it at `eb3 + hdr_len`.

### Segments

An `extblock` file is not one contiguous data region. At every recording break
the recorder writes a fresh copy of the channel-info block into the sample
stream (byte-identical to the one at `eb3` except for its timestamp) and then
continues. Each embedded header is sample-aligned and identifiable:

```
+0x00        0x01                      block type
+0x01..+0x04 "TIME"
+0x12,+0x13  0x02 0x02                 version
+0x14..      "YYYYMMDDHHMMSS" (ASCII)  this segment's start time
+0x44        u16 = n_ch - 1            channel count, must agree
+0x48..      the channel table again
```

The reader detects and skips these headers. (Earlier code that assumed one block
to EOF read each embedded header as ~5.8 samples of EEG, slipping the channel
phase by `hdr_len % stride` bytes at every break and over-reporting the sample
count; that is fixed.)

## Shared: electrode code → name (`.21E`)

Both layouts store, per channel, an electrode code that is the **0-based index
into the `.21E` `[ELECTRODE]`/`[REFERENCE]` table**. The reader keeps the code
`+1` (1-based) in `ch_indices` and resolves the name through the shared
`.21E[code-1]` lookup. This is numerically identical to MNE-Python's
`read_raw_nihon`, which uses the raw byte directly against a 0-based table.

`wfmblock` samples are µV throughout. `extblock` uses hardware-fixed ADC gains
(not stored in the file) to reach µV.

## DC trigger numbering by generation

DC channel numbering keys on the **format signature**, not the layout:

- `wfmblock` / EEG-1100 → `DC03`–`DC06`
- `extblock` / EEG-1200 → `DC01`–`DC04`

Both appear in the verification files below and match the vendor exports.

## Reference-channel naming (`$` → `_ref`)

Nihon Kohden's `$`-prefixed reference names (`$A1`, `$A2`) are normalized to a
Perl- and filename-safe `_ref` suffix (`A1_ref`, `A2_ref`). EDFBrowser keeps the
raw `$A1`/`$A2`; the physical channel is the same. This is a naming convention,
not a data difference.

## Verification

- **`wfmblock`** — `YJ0394RY.EEG` (EEG-1100C): after the off-by-one fix the
  channel list matches EDFBrowser's export channel-for-channel
  (`… Pz, E, A1, A2, vEOG, hEOG, X1, DC03–DC06, BN1, BN2 …`).
- **`extblock`** — `JJ0090J6.EEG` (EEG-1200A, EEG-1290 / JE-92NX, 37 ch) matches
  EDFBrowser across all 37 channels
  (`… X1, nose, lm, rm, DC01–DC04, BN1, BN2, $A1/$A2, COM`); `IJ0200S9.EEG`
  (EEG-1200C, multi-segment) reads with the correct channel order.
- Both are consistent with MNE-Python `read_raw_nihon`.

The only remaining difference from EDFBrowser is the reference-name spelling
(`A1_ref`/`A2_ref` vs `$A1`/`$A2`), which is the `_ref` normalization above.
