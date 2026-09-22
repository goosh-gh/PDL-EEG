# PDL::EEG — Nihon Kohden / EDF / BESA EEG toolkit

Read Nihon Kohden Neurofax recordings in PDL, resolve headbox-independent
trigger/channel labels, re-reference (incl. balanced non-cephalic), remove
eye-blink artifacts (GED), and export to EDF/EDF+ or BESA ASCII multiplexed
(`.mul`).

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
| `PDL::EEG::MAP2D` | `plot_topomap` — 2D scalp voltage map for one latency from a voltage vector + an ASA `.elc`. `orientation` selects the round **axial** map (sphere-fit azimuthal-equidistant, nose up, thin-plate-spline clipped to the head disc) or a **sagittal** side view (`sagittal-left` / `sagittal-right`, orthographic projection through the Fz-Cz-Pz plane, thin-plate-spline clipped to a head-profile silhouette). `plot_topomap_panels` draws several views in one figure on a shared colour scale. Renders with `PDL::Graphics::Cairo` (loaded on demand). |
| `PDL::EEG::TFA` | Continuous complex-Morlet time-frequency analysis (CWT), with frequency-domain convolution via `PDL::FFT` (no external wavelet dependency). `tfr_morlet`'s `output` gives total power, inter-trial coherence, or an exact phase-locked (`evoked`) / non-phase-locked (`induced`) power split (`power = evoked + induced`); `tfr_superlet` is the adaptive multiplicative superlet for short HFO bursts; `tfr_stat` is the across-trial reliability `z = mean/SEM`; `apply_baseline` normalises per frequency (`zscore`/`ratio`/`logratio`/`percent`/`mean`). |
| `PDL::EEG::Inverse::MinimumNorm` | L2 minimum-norm distributed source localization on a surface-normal-constrained leadfield (New York Head `V_fem_normal`, or any leadfield of the same shape). `inverse_operator`/`apply_inverse`/`source_estimate`/`source_power` build one data-independent inverse operator and apply it; `method` selects **MNE**, **sLORETA**, or **eLORETA** (same operator, different per-source standardization), `ref` is average (CAR, default — a symmetric-PSD pseudoinverse handles the rank-deficient Gram) or `none`, regularization is `reg_frac`/`alpha`, eLORETA takes `max_iter`/`tol`. `forward_project` is `b = K j`; `avg_reference` re-references an electrode subset for montage studies. `source_power` returns the standardized (dimensionless) source statistic. Pure PDL + `PDL::MatrixOps`; the real leadfield loads via `PDL::IO::NYHead`. |
| `PDL::EEG::GED` | GED (generalized eigenvalue decomposition) spatial filter for eye-blink removal. `ged_cov`/`ged_operator`/`apply_ged` are the generic core: solve `S w = lambda R w` (LAPACK `sygvd` via `PDL::LinearAlgebra`; eigenvectors are R-orthonormal, so patterns `A = R W` are biorthogonal to filters `W` and removing a component subset is the exact oblique projection `X - A_sel (W_sel' X)`), with shrinkage on `R`. `detect_blinks`/`blink_free`/`blink_evoked`/`fit_blink`/`remove_blinks` are the blink workflow: build `S` from a vEOG-detected, peak-normalised blink-averaged evoked and `R` from blink-free background **taken from the recording being cleaned**, pick the blink component by vEOG-evoked correlation, and remove it by back-projection. `remove_blinks` returns `(subtracted, removed, operator)` with `subtracted + removed = input`. Single blink-component removal; data layout `(nch, nt)`. Pure PDL + `PDL::LinearAlgebra`. |
| `PDL::EEG::ICA` | Symmetric FastICA in pure PDL (no LAPACK: whitening and the symmetric orthogonalization use `PDL::MatrixOps` `eigens_sym`). `ica_decompose` returns unmixing/mixing matrices and component activations; `identify_by_reference` scores each component against a reference channel (e.g. vEOG); `apply_ica` reconstructs with a chosen component subset removed. Data layout `(nt, nch)`. Built to compare ICA-based ocular cleaning against `PDL::EEG::GED`. |
| `PDL::EEG::Spans` | Time-interval annotations shared between the viewer and the artifact remover. A span is `{kind,t0,t1}` in seconds; `read_spans`/`write_spans` persist them to a `<file>.spans.tsv` sidecar; `spans_to_mask` rasterises one `kind` to a per-sample byte mask; `kinds_present` lists the kinds present. `kind` is free-form (`blink`, `hsaccade`, `bad`, or any GED class name). |
| `PDL::EEG::Regress` | Temporal nuisance / EOG regression. `regress_out($X, @sources)` removes, from every channel of `$X` (nch,nt), the least-squares fit on a set of source time courses (GED or ICA component activations, or raw EOG channels — with EOG channels it is classic Gratton–Coles regression). Sources are orthonormalised first, so correlated sources are partialled out jointly and the removal depends only on the span of the sources, not on labelling any one vertical/horizontal. Channel DC is preserved. Pure PDL. |
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
| `examples/topomap_demo.pl` | Worked `plot_topomap` example: reads a montage `.elc`, builds a synthetic average, writes a topomap PNG. Needs `PDL::Graphics::Cairo`. |
| `examples/nose_from_nyhead.pl` | Extract the nose-tip vertex from the New York Head skin surface (via `PDL::IO::NYHead`), co-register it to a montage frame over the shared 10-20 electrodes, and print an ASA `.elc` position line for a `nose` electrode. |
| `examples/silhouette_from_nyhead.pl` | Extract the mid-sagittal head-profile silhouette from the New York Head skin surface, co-register it to a montage frame, and write a `Y Z` polyline (`.poly`) for the sagittal views of `plot_topomap`. |
| `examples/sep_hfo_tfa.pl` | SEP high-frequency-oscillation time-frequency map. Reads EEGLAB `.set`+`.fdt` with BIDS sidecars (`_eeg.json`/`_channels.tsv`/`_events.tsv`; a v7.3 HDF5 `.set` via `PDL::IO::HDF5`), concatenates runs, epochs around the stimulus, auto-picks the contralateral central channel, runs the Morlet (or `--superlet`) transform, and renders with `PDL::Graphics::Cairo`. Views: `--itc`, `--stat`, `--decomp` (total/evoked/induced), `--sig` contours; `--excise <ms>` removes the stimulus artifact in-data; `--bl-min/--bl-max`, `--ymin/--ymax`, `--ytick/--xtick`, `--cmap`; `--demo` runs on synthetic data. Needs `PDL::Graphics::Cairo` to render. |
| `examples/nyhead_inverse.pl` | New York Head source-localization demo / montage study: seed a cortical source → forward `b = K j` → sLORETA/eLORETA inverse → peak vertex + distance-to-seed (mm) + Harvard–Oxford area. `--demo` runs on a synthetic leadfield (no `.mat`); `--mat sa_nyhead.mat` uses the real leadfield via `PDL::IO::NYHead`. Flags: `--method`, `--ref`, `--snr <dB>`, `--alpha`/`--reg-frac`, `--seed`/`--seed-area`, `--montage`, `--out powers.dat` (the `Ns`-row source power in cortex75K vertex order, for the GS3D cortical overlay). |
| `examples/sep_n20_inverse.pl` | Real **averaged-SEP N20** source localization on the New York Head. Reads a plain text `ch × time` matrix + one-label-per-line file (dump `evoked.data` / `ch_names` from MNE-Python), maps the recording montage to the 231-electrode leadfield by label (old→new 10-20 aliases; non-scalp channels drop out), auto-picks the N20 latency by GFP (or `--latency`), CAR-re-references the electrode subset, and runs the inverse → `powers.dat` (cortex75K order). Diagnostics (`--diag`): N20 topography, a whole-cortex single-dipole `corr` scan (`best 1-dip`, method-independent), top vertices, best-Postcentral rank. `--method`, `--reg-frac`, `--side`, `--expect`. Needs the NY Head `.mat` via `PDL::IO::NYHead`. |
| `examples/sep_n20_sweep.pl` | Latency sweep. Builds the (data-independent) inverse operator **once** and applies it across a latency window in a single batched solve, writing `powers.<lat>.dat` per latency plus `peak_by_latency.tsv` (latency → peak vertex / MNI / Harvard–Oxford area / Postcentral %). `--lat-min/--lat-max/--lat-step`, `--outdir`, `--summary-only` (tsv only, fastest), `--dump "20.1,20.5"` (write only listed latencies). |
| `examples/sep_gof_sweep.pl` | **Best-single-dipole goodness-of-fit** swept over latency, written as a two-column `latency_ms<TAB>gof` file (feeds the movie's `--gof`). No inverse solve: at each latency it takes the scalp topography and the maximum `|corr|` against every cortical leadfield column (the best fixed-orientation single dipole) — the same quantity `sep_n20_inverse.pl` prints as `best 1-dip`. Electrode-to-leadfield mapping and CAR are as in `sep_n20_sweep.pl`. `--metric corr` (default) or `r2` (variance explained), `--lat-min/--lat-max/--lat-step`, `--avg-ms` (window-average the topography). Reads the same MNE `np.savetxt` `ch × time` dump; needs the NY Head `.mat` via `PDL::IO::NYHead`. |
| `examples/sep_nonHFO_movie.3panel.pl` | **Three-panel latency movie** — waveforms · a 2D scalp topomap · the ECD goodness-of-fit curve, one frame per latency, assembled with `ffmpeg` into an mp4/gif. The three columns are equal width; each waveform trace is coloured by its **polarity at the N20 latency** (cool = negative, warm = positive, discrete palettes) to match the topomap. Reads the same `--text` `ch × time` dump + labels + `--montage` `.elc` (the private `eeg.pm` reader is not used), takes the GoF curve from `sep_gof_sweep.pl` via `--gof`, and draws the topomap with `PDL::EEG::MAP2D`. `--wave-chans`, `--gof-min-ms/--gof-max-ms`, `--anim-min-ms/--anim-max-ms/--step-ms`, `--polarity-ms`, `--neg-up`, `--head-extent`, `--figw/--figh`. Needs `PDL::Graphics::Cairo` **0.2+** (axis-off colour bar + content-aware margins) and `ffmpeg`. |
| `examples/sweep_usda_anim.pl` | Writes a latency window as **one animated USD** (`color3f[] primvars:displayColor.timeSamples`; the mesh, axes and colour-bar are static, so it can be paused and rotated). Computes the source power internally (no intermediate files); time code = latency × 10. `--norm global|frame`, `--cmap`/`--threshold`/`--base-grey`, `--axes`, `--colorbar`; `upAxis="Z"`; threshold recorded in `customLayerData`. cortex10K recommended. |
| `examples/powers_to_usda.pl` | Small wrapper: reads `powers.<lat>.dat` (from `sep_n20_sweep.pl`) and writes **one static USD per latency** — cortex coloured by the source map, `upAxis="Z"`, with optional XYZ axes and an independent colour-bar prim. `--norm global|frame|fixed`, `--lat-min/--lat-max`, `--indir/--outdir`. Uses `PDL::IO::NYHead` for the mesh + `cortex75K → cortexNK` vertex map. |
| `examples/avg_loreta_usda.pl` | **Averaged evoked-response** (ERP or SEP) source localization written straight to **one animated USD**, on the same MinimumNorm + New York Head stack. Reads a Nihon Kohden averaged file through the author's `eeg.pm` reader (like `topomap2d.pl`), or with `--text` a plain `ch × time` matrix + one-label-per-line file (`--names`, `--sfreq`, `--tmin`/`--pretrigger`) — the same MNE `np.savetxt` dump `sep_n20_inverse.pl` reads, so the ERP path can be cross-checked against the SEP N20 result. Maps the montage to the 231-electrode leadfield by label (10-20 aliases), CAR-re-references the subset, builds the inverse operator per method (**sLORETA and eLORETA by default**, `--method sloreta,eloreta`), and sweeps `START..END` ms, colouring the cortex by source power per latency (`color3f[] primvars:displayColor.timeSamples`). One USD per method, `${INFILE}${START}_${END}_<method>.usda`: cortex mesh + optional XYZ axes + colour-bar prim, `upAxis="Z"`, `metersPerUnit`, `defaultPrim`, and `doc`/provenance in `customLayerData`; time code = latency ms with `timeCodesPerSecond=1000` and `framesPerSecond=sfreq`, so usdview steps one sample (1 ms) per frame. `--method` (comma list), `--reg-frac`/`--alpha`, `--res`, `--step`, `--cmap`/`--threshold`/`--grey`, `--norm`, `--fps`, `--axes`/`--colorbar`. Needs `PDL::IO::NYHead` + `sa_nyhead.mat`; the default (non-`--text`) reader is the author's private `eeg.pm`. |
| `examples/includeORexcludeEEC.pl` | **Proxy-electrode A/B source study.** Runs the eLORETA inverse **without** a channel and **with** it, where the extra channel's measured potential is paired with a nearby modelled electrode's leadfield row (a *proxy*), and reports the peak vertex / Harvard–Oxford area / MNI, peak shift (mm), whole-cortex correlation of the two source maps, and normalized max difference; writes `proxy_without.dat` / `proxy_with_<proxy>.dat` (cortex75K order) for the GS3D overlay. `--proxy` (comma list of leadfield electrodes), `--latency` (or GFP auto-pick), `--reg-frac`, `--ear-label`. `--self-test` runs the whole with/without comparison on a synthetic leadfield — no data or `.mat` — to validate the engine path. The real-data reader is the author's private `eeg.pm` (as in `avg_loreta_usda.pl`). |
| `examples/sep_gof_sweep.pl` | **Best-single-dipole goodness-of-fit** swept over latency, written as a two-column `latency_ms<TAB>gof` file (feeds the movie's `--gof`). No inverse solve: at each latency it takes the scalp topography and the maximum `|corr|` against every cortical leadfield column (the best fixed-orientation single dipole) — the same quantity `sep_n20_inverse.pl` prints as `best 1-dip`. Electrode-to-leadfield mapping and CAR are as in `sep_n20_sweep.pl`. `--metric corr` (default) or `r2` (variance explained), `--lat-min/--lat-max/--lat-step`, `--avg-ms` (window-average the topography). Reads the same MNE `np.savetxt` `ch × time` dump; needs the NY Head `.mat` via `PDL::IO::NYHead`. |
| `examples/sep_nonHFO_movie.3panel.pl` | **Three-panel latency movie** — waveforms · a 2D scalp topomap · the ECD goodness-of-fit curve, one frame per latency, assembled with `ffmpeg` into an mp4/gif. The three columns are equal width; each waveform trace is coloured by its **polarity at the N20 latency** (cool = negative, warm = positive, discrete palettes) to match the topomap. Reads the same `--text` `ch × time` dump + labels + `--montage` `.elc` (the private `eeg.pm` reader is not used), takes the GoF curve from `sep_gof_sweep.pl` via `--gof`, and draws the topomap with `PDL::EEG::MAP2D`. `--wave-chans`, `--gof-min-ms/--gof-max-ms`, `--anim-min-ms/--anim-max-ms/--step-ms`, `--polarity-ms`, `--neg-up`, `--head-extent`, `--figw/--figh`. Needs `PDL::Graphics::Cairo` **0.2+** (axis-off colour bar + content-aware margins) and `ffmpeg`. |
| `examples/blink_ged_clean.pl` | Eye-blink removal on a continuous EDF via `PDL::EEG::GED`. Reads `TASK.edf` (plus an optional `VOLUNTARY_BLINKS.edf` whose blinks define the pattern), detects blinks on vEOG, fits the GED blink filter (`S` from the voluntary session's peak-normalised blink evoked, `R` from the task's blink-free background), removes the blink component, and writes the result -- the extension picks the format: `.edf` via `PDL::EEG::IO::EDF::write_edf`, `.mul` via `PDL::EEG::IO::BESA::ASCII::write_mul` (loaded on demand). By default every signal channel is filtered (EOG included) and the **whole recording** is written (filtered channels replaced, DC/Event/Marker and `--exclude` channels passed through, start time preserved), so DC trigger channels survive for epoching/averaging. `--montage-only` writes only the filtered channels; `--out`/`--out-removed` (cleaned / removed-artifact EDF or `.mul`), `--exclude LABEL,...`, `--trim-sec`, `--n`, `--reg`, `--veog`. |
| `examples/raw_plot.pl` | Interactive span annotator over `giza-server` (`show_interactive`). Scrolls an EDF/NK recording, marks time intervals by mouse (left-click start, left-click end to commit; right-click to delete) under a free-form `--kind` (blink, hsaccade, bad, or any GED class name), and appends them to a `<file>.spans.tsv` sidecar read back by `artifact_remove.pl`. `--chans`, `--aux` auxiliary-channel scaling, per-kind colours. Needs `PDL::Graphics::Cairo` + `giza-server`. |
| `examples/artifact_remove.pl` | Remove ocular / stereotyped artifacts by GED spatial filter (`PDL::EEG::GED`) or temporal regression (`PDL::EEG::Regress`) and write the recording back (`--output mul|edf`). Classes are declared with one `--ged NAME[:CHANS[:SOURCE]]` per class, in removal order: `NAME` is a kind in the spans file; `CHANS` is an optional signed channel list used only to give the axis a meaning (correlation report, and the `--regress` source), never to build the filter; `SOURCE` is `span` (covariance of the NAME segments, default) or `blink` (peak-normalised blink-evoked covariance, from `--blink-from` or the target's blink segments). Each class builds R from the segments that are not that artifact and not `bad`, keeping every other declared artifact in R so the classes do not disturb one another. Spatial removal is one `apply_ged` oblique projection per class applied in order; `--regress` uses `regress_out` on the class source time courses. `bad` segments are restored verbatim (excluded later at averaging). Before writing, a review image (original-vs-cleaned waveform overlay + one 2D topography per class, needs `--montage`) is rendered and the write confirmed (`--no-review`/`--yes`). `--rank` is a diagnostic-only mode: it ranks every component of each class by eigenvalue with r(vEOG)/r(hEOG)/kurtosis/low-frequency-ratio/`|pattern|` and auto eye-flags, renders a per-component figure (topography + representative-channel waveforms + component time course), and can `--dump-components` the selected component scores as channels. `--eeg/--extra/--eog/--veog/--heog`, `--blink-thresh/-refr/-half`, `--shrink`, `--gap-min`, `--eye-corr`, `--rank-top`. Needs `PDL::Graphics::Cairo` + `PDL::EEG::MAP2D` for the images. |
| `examples/task_ica_mul.pl` | ICA-based blink + horizontal removal on a continuous recording → full-record `.mul`, via `PDL::EEG::ICA`. |
| `examples/ged_vs_ica_blink.pl` | Quantitative GED-vs-ICA comparison on blink removal (component identification and residual). |
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

## Scalp topography (2D)

`PDL::EEG::MAP2D::plot_topomap` draws a scalp voltage map for one latency.
Electrode positions come from an ASA `.elc` (read by `PDL::EEG::IO::ASA`);
voltages are a per-channel vector, or `avg[chan,time]` plus `time`.

```perl
use PDL::EEG::MAP2D qw(plot_topomap plot_topomap_panels);

plot_topomap(
    avg     => $avg,                 # (n_ch, n_time); or values => <n_ch vector>
    time    => $sample,
    labels  => \@channel_names,      # row order of $avg
    montage => 'standard_1020_eog_nose.elc',
    clim    => 15,                   # ± µV (omit for auto, from scalp channels)
    contours=> 6,
    names   => 1,
    title   => 'wp1 160 ms',
    outfile => 'topo.png',           # or device => 'gs' for a giza-server window
);
```

### Axial (round) map — `orientation => 'axial'` (default)

Scalp sensors are sphere-fitted, projected by an azimuthal-equidistant map (nose
up, right ear right), recentred on the scalp centroid, and interpolated with a
thin-plate spline clipped to the head disc (which extends `overshoot` past the
drawn head circle). Old 10-20 names `T3 T4 T5 T6` map to `T7 T8 P7 P8`;
`LM`/`RM` map to `M1`/`M2` (alias lookup is case-insensitive, so `lm`/`rm`/`t3`
resolve too); fiducials are skipped. Periocular sensors (`/EOG/`, or any sensor
below `periph_deg`, e.g. an aliased `X1`) are kept out of the sphere fit, scale
and colour limits, and by default are drawn but not interpolated
(`eog_interp => 1` includes them). A `nose` electrode is kept out of the head-
circle geometry (so its far-forward position can't distort the circle) but its
value feeds the colour scale and interpolation by default (`nose_interp => 0`
draws it without letting it colour the map). Layout knobs: `margin` (sensor
inset), `overshoot` (colour past the circle), `ear_dy` (ear height),
`head_extent` (half-width of the axes in data units, default 0.6; the head circle
is radius 0.5, so smaller values grow the circle toward the axes edges, e.g.
0.55 leaves little margin for a tightly packed panel). With `PDL::Graphics::Cairo`
0.2+ an `axis('off')` topomap panel can carry its own colour bar (see below).

### Sagittal (side) view — `orientation => 'sagittal-left'` / `'sagittal-right'`

An orthographic projection onto the mid-sagittal (Fz-Cz-Pz) plane. The left view
faces left (Fp1 at screen-left, left hemisphere shown); the right view faces
right (Fp2 at screen-right, right hemisphere shown) — the far hemisphere is
dropped. `nose`, mastoid (`lm`/`rm`) and EOG channels project onto the profile
and take part like any near-side sensor. The head outline is a side-profile
silhouette supplied as a polyline via `silhouette => 'file.poly'` (one `Y Z`
pair per line, montage frame, mm — produced by
`examples/silhouette_from_nyhead.pl`); without one, a convex hull of the
projected sensors is used. Whatever the silhouette, it is grown just enough to
enclose every projected sensor so no electrode falls outside the coloured region
(`fit_silhouette => 0` disables this). `midline_tol` sets how far off the
midline a sensor may sit and still appear in both views.

### Several views in one figure — `plot_topomap_panels`

```perl
plot_topomap_panels(
    avg     => $avg, time => $sample, labels => \@channel_names,
    montage => 'standard_1020_eog_nose.elc',
    silhouette => 'sagittal_silhouette.poly',
    clim    => 15,
    names   => 1,
    title   => 'wp1 160 ms',
    outfile => 'panels.png',         # or device => 'gs'
);
```

Draws the views in `views` (default `['sagittal-left','axial','sagittal-right']`)
side by side on one shared, symmetric colour scale with a single colour bar.
Per-panel labels come from `titles` (default `left` / `axial` / `right`); the
overall `title` is the figure heading. Per-panel fonts are smaller than the
single-map defaults because the maps are packed.

### Text and colour-bar knobs

`names => 1` labels the sensors; `name_size` (default 10) sizes those labels,
`title_size` the titles (12 for panel labels, 13 for the figure heading),
`cbar_size` (default 9) the colour-bar numbers. A label that would overrun the
right edge is placed to the left of its marker instead, so far-lateral names
(`A2`, `rm`, `M2`) stay clear of the colour bar. `cbar_label_x` shifts the
colour-bar numbers, `unit` sets the colour-bar caption (default `uV`).

### Output

`outfile` writes a PNG (the default `device`). `device => 'gs'` (aliases `osx`,
`aqua`, `cocoa`, `giza`) opens a giza-server window; `device => 'gnuplot'`
(aliases `qt`, `x11`) uses the gnuplot backend. With `outfile` set the call
returns the filename; otherwise it returns the figure.

`PDL::Graphics::Cairo` is required only for rendering (loaded on demand);
`project_positions` / `interpolate_topo` and `t/15_map2d.t` need only PDL. Full
options in `perldoc PDL::EEG::MAP2D`; worked example in
`examples/topomap_demo.pl`.

## Time-frequency analysis (HFO)

`PDL::EEG::TFA` is a renderer- and IO-independent time-frequency engine: a
complex Morlet continuous wavelet transform whose convolution is done in the
frequency domain via `PDL::FFT`, so there is no external wavelet dependency. It
was built for somatosensory-evoked high-frequency oscillations, but nothing in
it is HFO-specific.

```perl
use PDL::EEG::TFA qw(tfr_morlet apply_baseline);

# $epochs (nepoch, ntime) real; $sfreq Hz; $times seconds
my $freqs = 70 + 5 * sequence(167);                    # 70..900 Hz
my $power = tfr_morlet($epochs, $sfreq, $freqs, n_cycles => 7);
my $z     = apply_baseline($power, $times, [-0.05, -0.004], mode => 'zscore');
```

`tfr_morlet`'s `output` selects the total power (default, mean `|coef|^2`),
inter-trial coherence (`itc`), the phase-locked evoked power (`evoked`,
`|mean coef|^2`), or the non-phase-locked induced power (`induced`,
`total - evoked`); the split is exact. `n_cycles` is a scalar or a per-frequency
piddle (a ramped-cycle grid suits a wide band). `tfr_superlet` is the adaptive
multiplicative superlet, which concentrates short high-frequency bursts more
sharply in both time and frequency. `tfr_stat` returns the across-trial
reliability `z = mean / SEM`, computed streaming (no per-trial storage) and
distinct from the temporal-baseline z of `apply_baseline` — with many trials it
is approximately Gaussian. `apply_baseline` normalises each frequency against a
baseline window (`zscore`, `ratio`, `logratio`, `percent`, `mean`).

### Worked example — `examples/sep_hfo_tfa.pl`

Reads EEGLAB SEP data in the BIDS layout (a continuous `.set` with a companion
`.fdt`, plus `_eeg.json` / `_channels.tsv` / `_events.tsv`; a v7.3 HDF5 `.set`
is read via `PDL::IO::HDF5` when there is no `.fdt`), concatenates runs, epochs
around the stimulus, auto-picks the contralateral central channel, runs the
transform, z-scores against the pre-stimulus baseline, and renders the heatmap.

```
perl -Ilib examples/sep_hfo_tfa.pl RUN1_eeg.set RUN2_eeg.set -o hfo.png
perl -Ilib examples/sep_hfo_tfa.pl --demo -o hfo_demo.png       # no data needed
```

Views: `--superlet` (superlet power), `--itc` (inter-trial coherence), `--stat`
(across-trial reliability z), `--decomp` (three panels total / evoked / induced
on a diverging, 0-centred, displayed-band scale). `--sig "5,10"` overlays
significance contours at those `|z|` levels (for `--itc`, use 0–1 levels such as
`"0.1,0.2"`). Artifact and baseline control: `--excise <ms>` interpolates the
stimulus artifact in the data before the transform (so the wavelet cannot spread
it into neighbouring latencies), and `--bl-min`/`--bl-max` set the baseline
window. Display: `--ymin`/`--ymax`, `--ytick`/`--xtick`, `--cmap`.
`PDL::Graphics::Cairo` is required only for rendering (loaded on demand); the
engine and `t/16_tfa.t` need only PDL.

## Distributed source localization (inverse)

`PDL::EEG::Inverse::MinimumNorm` solves the EEG inverse problem on a
surface-normal-constrained leadfield (built for the New York Head bundled
leadfield `V_fem_normal`, but takes any leadfield of the same shape). It is the
L2 minimum-norm family: MNE, sLORETA and eLORETA share **one** data-independent
inverse operator and differ only in how each source row is standardized. Pure
PDL — depends on PDL core and `PDL::MatrixOps` only (no `PDL::LinearAlgebra` /
LAPACK). The real New York Head leadfield is loaded through `PDL::IO::NYHead`.

The leadfield is stored `(Ns, Ne)` = (source, electrode); the New York Head
`leadfield()` returns `(Ne, Ns)`, so pass `->transpose`. Every entry point
asserts `Ns > Ne` and rejects a transposed leadfield. The signal model is
`b = K j`.

```perl
use PDL::EEG::Inverse::MinimumNorm
    qw(inverse_operator apply_inverse source_power);

# one-shot: leadfield K (Ns,Ne) + scalp topography b (Ne) -> source power (Ns)
my $op = inverse_operator($K, method => 'eloreta', ref => 'car');
my $pw = source_power($op, $b);

# reuse the operator across many time points (b as (Ne,Nt))
my $J  = apply_inverse($op, $b);        # (Ns,Nt)
```

`inverse_operator` options: `method => 'mne'|'sloreta'|'eloreta'`, `ref =>
'car'` (default) `| 'none'`, `reg_frac => 0.05` (Tikhonov `α =
reg_frac·trace(KKᵀ)/Ne`) or an explicit `alpha`, and for eLORETA `max_iter =>
100`, `tol => 1e-10`. `source_power` returns the standardized source power (the
sLORETA/eLORETA statistic) — a dimensionless localization quantity, not a
physical current density.

An average-referenced leadfield has `1ᵀK = 0`, so the Gram `C = KKᵀ + αH` is
rank `Ne−1`. The operator uses a symmetric-PSD pseudoinverse (eigendecomposition,
`rcond` relative to the largest eigenvalue) that drops the null direction —
correct for both the CAR and the `ref => 'none'` full-rank cases. `V_fem_normal`
is average-referenced over its 231 electrodes; an electrode subset is
re-referenced over the chosen electrodes with `avg_reference`, so a montage is
simulated by taking those leadfield rows.

### Worked example — `examples/nyhead_inverse.pl`

Seeds a known cortical source, forward-projects it to a scalp topography,
inverts, and reports the peak vertex, its distance to the seed (mm), and its
Harvard–Oxford area; `--montage` compares electrode subsets.

```
# synthetic, runs anywhere (no data)
perl -Ilib examples/nyhead_inverse.pl --demo --method eloreta \
     --montage all --montage 19 --montage 0,1,2,3

# real New York Head leadfield (needs PDL::IO::NYHead + sa_nyhead.mat)
perl -Ilib examples/nyhead_inverse.pl --mat sa_nyhead.mat \
     --method eloreta --seed-area "Postcentral Gyrus" --out powers.dat
```

Noiseless, every montage localizes the seed exactly (0.0 mm) — the analytic
exact-localization result; `--snr <dB>` adds noise and the error grows as
electrodes are removed. `--out` writes the `Ns`-row source power in cortex75K
vertex order, which the GS3D New York Head viewer (`PDL::Graphics::Cairo`,
`demo_gs3d_nyhead.pl --source-power`) colours onto the cortical mesh.

The operators are cross-checked element-wise against an independent NumPy
implementation (`examples/verify_inverse_numpy.py`) to machine precision, and
eLORETA reproduces the single-point exact-zero-localization result of
Pascual-Marqui.

## Artifact removal — GED spatial filter, regression, diagnostics

`examples/raw_plot.pl` and `examples/artifact_remove.pl`, with the modules
`PDL::EEG::Spans`, `PDL::EEG::Regress`, `PDL::EEG::ICA` and the generic
`ged_operator` / `apply_ged` core of `PDL::EEG::GED`, remove ocular and other
stereotyped artifacts from a continuous recording and write it back. The
workflow is: mark intervals, estimate one spatial filter per artifact, remove,
keep `bad` for later exclusion.

1. **Mark** with `raw_plot.pl` (interactive, over `giza-server`). Each mouse
   drag adds a time interval under a free-form `--kind`, appended to a
   `<file>.spans.tsv` sidecar. `bad` marks non-stationary artifact (sweat,
   movement) that GED cannot filter.

   ```
   perl -Ilib examples/raw_plot.pl task.edf --kind hsaccade
   perl -Ilib examples/raw_plot.pl task.edf --kind bad
   ```

2. **Remove** with `artifact_remove.pl`. Declare one class per artifact with
   `--ged NAME[:CHANS[:SOURCE]]`, in removal order:

   ```
   perl -Ilib examples/artifact_remove.pl task.edf \
       --ged blink:vEOG:blink --blink-from voluntary_blinks.edf \
       --ged hsaccade:hEOG \
       --spans task.edf.spans.tsv --montage standard_1020_eog_nose.elc \
       --output edf --out task-clean.edf
   ```

   - `NAME` is a kind in the spans file. `SOURCE` is `span` (covariance of the
     NAME segments, default) or `blink` (peak-normalised blink-evoked
     covariance, from `--blink-from` or the target's own blink segments).
   - `CHANS` is an optional signed channel list (e.g. `Fp1,Fp2,-A1,-A2`) used
     **only** to give the estimated axis a meaning — a correlation report, and
     the source for `--regress`. It does **not** build the filter, so a wrong
     or missing `CHANS` never changes what is removed.
   - Each class builds its reference covariance R from the segments that are
     **not** that artifact and **not** `bad`, keeping every **other** declared
     artifact in R. A GED filter is orthogonal (in the R metric) only to what
     lives in R, so keeping the other artifacts in R is what stops one class's
     removal from disturbing another.
   - Spatial removal (default) is one exact `apply_ged` oblique projection per
     class, applied in the given order. `--regress` instead removes the class
     source time courses with `PDL::EEG::Regress::regress_out` (stronger, but
     also removes brain that co-varies in time).
   - `bad` segments are restored to the recorded samples after removal, to be
     excluded later at the averaging stage using the same spans file.
   - Each operator prints its generalized-eigenvalue gap (`lambda_top /
     lambda_next`); a gap near 1 means the marked segments did not isolate one
     axis and the filter is unreliable (`--gap-min` warns).
   - Before writing, a review image (original-vs-cleaned waveform overlay plus
     one 2D scalp topography per class, with `--montage`) is rendered and the
     write confirmed on the terminal; `--no-review` / `--yes` skip it.

3. **Diagnose** with `--rank` (does not remove or write). For each class it
   ranks every component by eigenvalue with r(vEOG), r(hEOG), kurtosis, a
   low-frequency-power ratio, `|pattern|` and automatic eye-flags
   (`--eye-corr`), renders one row per selected component (topography +
   representative-channel waveforms + component time course), and can
   `--dump-components mul|edf` the selected component scores as channels
   (arbitrary units, not µV). Identification is deliberately not by eigenvalue
   rank — the artifact axis is not guaranteed to be the top eigenvalue — but by
   correlation with a reference and by inspecting the topography.

`PDL::EEG::Regress::regress_out` is usable on its own for classic EOG
regression (pass the EOG channels as the sources). `PDL::EEG::ICA` provides the
ICA route for the same ocular cleaning, for comparison against GED.

Layouts: `PDL::EEG::GED`, `Spans`, `Regress` and `artifact_remove.pl` use
`(nch, nt)`; `PDL::EEG::ICA` uses `(nt, nch)`.


## Tests

```
make test        # 20 files (t/06 reserved/skipped), 416 subtests
```

`t/01_nihonkohden` `t/02_edf` `t/03_ptn` `t/04_signal` `t/05_montage`
`t/06_reserved` (skip: reserved) `t/07_blocks` (block_extents + multi-segment
extblock regression) `t/08_epoch` (event placement + wall-clock→sample mapping)
`t/09_i18n` `t/10_besa_ascii` `t/11_edf_to_mul` `t/12_derivation`
`t/13_edf_roundtrip` (µV round-trip incl. DC, per-signal EDF dimension)
`t/14_asa` (ASA `.elc` reader: parse, fiducials, name lookup, Advent shim,
against a 28-point real-coordinate fixture).
`t/15_map2d` (MAP2D: azimuthal projection orientation, scalp recentering,
thin-plate-spline grid; render-free, needs only PDL).
`t/16_tfa` (TFA: Morlet wavelet energy, tone/burst localisation, ITC, baseline
modes, superlet, across-trial reliability, evoked/induced decomposition;
render-free, needs only PDL).
`t/17_minimumnorm` (MinimumNorm inverse: storage convention, forward linearity,
transpose rejection, MNE operator vs the closed form, sLORETA/eLORETA exact
single-source localization, CAR multi-source tracking; render-free, needs only
PDL).
`t/18_regress` (Regress: both artifact sources removed, brain preserved where
the artifact is weak, span-invariance/label-free, channel DC preserved,
cleaned + removed reconstructs the input; render-free, needs only PDL).
`t/19_spans` (Spans: TSV sidecar round-trip, `spans_to_mask` rasterisation,
offset/clip edge cases; needs only PDL).
`t/20_ica` (ICA: whitening, symmetric FastICA source recovery, remove-none
identity, mixing reconstruction; render-free, needs only PDL).

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
- **GED blink removal is a single spatial component.**
  `examples/blink_ged_clean.pl` / `PDL::EEG::GED` remove the one blink component
  (the top vEOG-correlated generalized eigenvector). Frontal channels (Fp1/Fp2)
  clean well; vEOG and especially hEOG reduce less, because their averaged
  deflection also carries non-blink eye activity (vertical drift, horizontal
  saccades) that a vertical-blink filter correctly leaves in place. Removing more
  components (`--n 2`/`3`) reduces the EOG channels further but starts eroding
  real signal, so top-1 is the default. Fit `R` (the background) from the
  recording being cleaned, not a different session, to avoid over-removal.
