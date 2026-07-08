# P:EEG08 — NK → EDF, and headbox-independent trigger/label resolution

## What ships here

| File | Package | Role |
|------|---------|------|
| `PDL/EEG/IO/EDF.pm` | `PDL::EEG::IO::EDF` | Write EDF / EDF+ from a `read_nk` record |
| `PDL/EEG/IO/NihonKohden.pm` | `…::NihonKohden` | Reader, **patched**: `label_map` option + `ch_indices` return |
| `PDL/EEG/IO/NihonKohden/PTN.pm` | `…::NihonKohden::PTN` | Parse Neurofax `.PTN` montage files (1100C + 1200A) |
| `PDL/EEG/IO/NihonKohden/Montage.pm` | `…::NihonKohden::Montage` | `.LOG` montage + `.PTN` + signal → `label_map` |
| `PDL/EEG/Signal.pm` | `PDL::EEG::Signal` | **Device-independent** square-pulse/TTL detector |
| `examples/nk_to_edf.pl` | — | NK → EDF CLI |
| `examples/nk_resolve_labels.pl` | — | Resolve trigger labels, optionally write EDF |
| `examples/parse_ptn.pl` | — | Dump one `.PTN` montage |
| `t/04_edf.t` | — | EDF writer round-trip test |

Naming reflects responsibility: the generic TTL detector lives in
`PDL::EEG::Signal` (no vendor knowledge); everything that reads Nihon Kohden
sidecars lives under `PDL::EEG::IO::NihonKohden::*`.

## The core problem (why this is not just "read the file")

Trigger/DC channel names are **not** derivable from the recording format:

- The same trigger line is `DC03–06` on the 1100C headbox (JE-921A) and
  `DC01–04` on the 1200A. Searching for a fixed name is a landmine.
- The bundled `.21e` is a generic template; its DC block does not match every
  headbox's wiring.
- The authoritative display names live in the **montage** (`.PTN`), which labels
  the four TTL lines `TrigBit0/2/4/8`; the electrode table separately calls them
  `DCxx`. Same physical lines, two labels.
- **Which recorded `ch_idx` carries a trigger is only visible in the signal.**
  The `.PTN` gives the trigger *count and names* but stores `G1=0` for them, so
  it does **not** encode their channel index.

No single source is sufficient — the resolver combines all three.

## Pipeline

```
.LOG  ──montage_from_log──▶ "IIA"
                              │  find the .PTN whose NAME == "IIA"
.PTN ──parse_ptn──▶ trigger names [TrigBit0,2,4,8]  (count = 4)
                              │
.EEG data ──detect_square_pulses(n=4)──▶ ch_idx that actually pulse
                              │  (needs the FULL recording; run all_blocks=>1)
                              ▼
        zip names(montage order) ⟷ triggers(ch_idx order)
                              ▼
              label_map { ch_idx => name }   →  read_nk(label_map => …)
```

### Usage

```perl
use PDL::EEG::IO::NihonKohden          qw(read_nk);
use PDL::EEG::IO::NihonKohden::Montage qw(resolve_labels);
use PDL::EEG::IO::EDF                  qw(write_edf);

my $rec = read_nk('YJ0394VB.EEG', all_blocks => 1);
my $r   = resolve_labels($rec, ptn_dir => 'YJ0394VB.PTN');
#  $r->{montage}   => "IIA"
#  $r->{label_map} => { 45=>'TrigBit0', 46=>'TrigBit2', 47=>'TrigBit4', 74=>'TrigBit8' }

my $rec2 = read_nk('YJ0394VB.EEG', all_blocks=>1, label_map => $r->{label_map});
write_edf($rec2, 'YJ0394VB.edf');
```

Or one shot: `perl examples/nk_resolve_labels.pl YJ0394VB.EEG --write YJ0394VB.edf`

To use your own names (physical box labels rather than montage labels):

```perl
resolve_labels($rec, ptn_dir=>'…', names => [qw(DC03 DC04 DC05 DC06)]);
# or pin by hand, highest priority:
read_nk($f, all_blocks=>1, label_map => { 45=>'DC03',46=>'DC04',47=>'DC05',74=>'DC06' });
```

## Honest caveats

- **Full recording required for detection.** Triggers must actually fire to the
  rail to separate cleanly from EEG; run `all_blocks=>1`. Over a short window
  where a bit never toggles, that line won't be detected.
- **Name↔channel order is an assumption.** Montage trigger names (slot order)
  are zipped onto detected triggers sorted by ascending `ch_idx`. That matches
  the acquisition order here, but verify once per headbox; `label_map` overrides.
- **Trigger→ch_idx is not in any file** — it comes from the signal. The `.PTN`
  only supplies names/count.
- **Absolute DC µV is not calibrated here.** Trigger *edges* (hence ERP epoching)
  are gain-independent, so this doesn't matter for event extraction. Matching the
  vendor's 100 mV-scale display would need a one-point calibration.
- **Export-time montage is not recoverable** from the files (`.reg`
  `DISPLAY_USE_LAST_PATTERN=0`). The *recording* montage (from `.LOG`) is, and is
  what this pipeline uses.

## Patch applied to `NihonKohden.pm`

Minimal, additive:
- `label_map => \%h` option (keyed by 1-based `ch_idx`), highest priority above
  `.21e` and `DEFAULT_LABELS`, in both the wfmblock and extblock label loops.
- `ch_indices => \@idx` in the return of both layouts (extblock aliases the
  existing `ch_hw_idx`). Needed to map detector positions ↔ `ch_idx`.
