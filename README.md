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
| `PDL::EEG::IO::WFDB` | PhysioNet WFDB records. `read_wfdb($stem)` reads the `.hea` and the format-16 signals of a single-segment record from a `.dat` file or a Matlab v4 `.mat` (the PhysioNet / Computing in Cardiology Challenge 2018 layout), keeping the samples as int16; `chan` returns one channel in physical units; a signal file cut short is read up to its last whole sample with a warning. `read_wfdb_annotations` reads MIT-format annotation files (signed SKIP intervals, aux texts); `wfdb_stages_events` turns the aux texts into a sleep stage per sample (W, N1, N2, N3, R, unscored) and a list of `(name ... name)` events. Other signal formats and Matlab v5 / v7.3 files are not read. |
| `PDL::EEG::IO::BrainVision` | `read_brainvision($vhdr, channels => [...], skip_missing => 1)` reads a BrainVision recording (`.vhdr` + binary `.eeg` + `.vmrk`) into the `read_nk` layout: `data [n_ch,n_samp]` in µV, `fs`, `labels`, `units`, `events` (marker onsets in seconds). INT_16, UINT_16, INT_32 and IEEE_FLOAT_32, multiplexed or vectorized; each channel is scaled by its resolution and converted to µV from µV/mV/V/nV (other units are left as they are). With `channels`, a multiplexed file is read in blocks of samples and only those channels are kept. |
| `PDL::EEG::Derivation` | `derive` (general linear derivation `y = M·x`), `bne` (balanced non-cephalic re-reference), `rereference` (single/linked/average). |
| `PDL::EEG::Signal` | Device-independent square-pulse / TTL detector. |
| `PDL::EEG::MAP2D` | `plot_topomap` — 2D scalp voltage map for one latency from a voltage vector + an ASA `.elc`. `orientation` selects the round **axial** map (sphere-fit azimuthal-equidistant, nose up, thin-plate-spline clipped to the head disc) or a **sagittal** side view (`sagittal-left` / `sagittal-right`, orthographic projection through the Fz-Cz-Pz plane, thin-plate-spline clipped to a head-profile silhouette). `plot_topomap_panels` draws several views in one figure on a shared colour scale. Renders with `PDL::Graphics::Cairo` (loaded on demand). |
| `PDL::EEG::TFA` | Continuous complex-Morlet time-frequency analysis (CWT), with frequency-domain convolution through `PDL::EEG::FFT::Dispatcher` (`PDL::FFTDispatcher` when installed (`PDL::PocketFFT`, or the library named by the environment variable `PDL_FFTDISPATCHER_BACKEND`), otherwise `PDL::FFTW3` when installed, otherwise `PDL::FFT` for the lengths it transforms correctly; no external wavelet dependency). `tfr_morlet`'s `output` gives total power, inter-trial coherence, or an exact phase-locked (`evoked`) / non-phase-locked (`induced`) power split (`power = evoked + induced`); `tfr_superlet` is the adaptive multiplicative superlet for short HFO bursts; `tfr_stat` is the across-trial reliability `z = mean/SEM`; `apply_baseline` normalises per frequency (`zscore`/`ratio`/`logratio`/`percent`/`mean`). |
| `PDL::EEG::Inverse::MinimumNorm` | L2 minimum-norm distributed source localization on a surface-normal-constrained leadfield (New York Head `V_fem_normal`, or any leadfield of the same shape). `inverse_operator`/`apply_inverse`/`source_estimate`/`source_power` build one data-independent inverse operator and apply it; `method` selects **MNE**, **sLORETA**, or **eLORETA** (same operator, different per-source standardization), `ref` is average (CAR, default — a symmetric-PSD pseudoinverse handles the rank-deficient Gram) or `none`, regularization is `reg_frac`/`alpha`, eLORETA takes `max_iter`/`tol`. `forward_project` is `b = K j`; `avg_reference` re-references an electrode subset for montage studies. `source_power` returns the standardized (dimensionless) source statistic. Pure PDL; the symmetric eigendecomposition goes through `PDL::EEG::Linalg` (LAPACK `syevd` when `PDL::LinearAlgebra` 0.436 or later is installed, `PDL::MatrixOps` `eigens_sym` otherwise); the real leadfield loads via `PDL::IO::NYHead`. |
| `PDL::EEG::GED` | GED (generalized eigenvalue decomposition) spatial filter for eye-blink removal. `ged_cov`/`ged_operator`/`apply_ged` are the generic core: solve `S w = lambda R w` (LAPACK `sygvd` via `PDL::LinearAlgebra`; eigenvectors are R-orthonormal, so patterns `A = R W` are biorthogonal to filters `W` and removing a component subset is the exact oblique projection `X - A_sel (W_sel' X)`), with shrinkage on `R`. `detect_blinks`/`blink_free`/`blink_evoked`/`fit_blink`/`remove_blinks` are the blink workflow: build `S` from a vEOG-detected, peak-normalised blink-averaged evoked and `R` from blink-free background **taken from the recording being cleaned**, pick the blink component by vEOG-evoked correlation, and remove it by back-projection. `remove_blinks` returns `(subtracted, removed, operator)` with `subtracted + removed = input`. Single blink-component removal; data layout `(nch, nt)`. Pure PDL + `PDL::LinearAlgebra`. |
| `PDL::EEG::ICA` | Symmetric FastICA in pure PDL (whitening and the symmetric orthogonalization use `PDL::EEG::Linalg`: LAPACK `syevd` when `PDL::LinearAlgebra` 0.436 or later is installed, `PDL::MatrixOps` `eigens_sym` otherwise). `ica_decompose` returns unmixing/mixing matrices and component activations; `identify_by_reference` scores each component against a reference channel (e.g. vEOG); `apply_ica` reconstructs with a chosen component subset removed. Data layout `(nt, nch)`. Built to compare ICA-based ocular cleaning against `PDL::EEG::GED`. |
| `PDL::EEG::Linalg` | `eigh_sym` — symmetric eigendecomposition (values ascending, vectors in columns) through LAPACK `syevd` (`PDL::LinearAlgebra` 0.436 or later) when installed, `PDL::MatrixOps` `eigens_sym` otherwise, with one warning per process; `linalg_backend` names the solver in use. Used by `PDL::EEG::ICA` and `PDL::EEG::Inverse::MinimumNorm`. |
| `PDL::EEG::Spans` | Time-interval annotations shared between the viewer and the artifact remover. A span is `{kind,t0,t1}` in seconds; `read_spans`/`write_spans` persist them to a `<file>.spans.tsv` sidecar; `spans_to_mask` rasterises one `kind` to a per-sample byte mask; `kinds_present` lists the kinds present. `kind` is free-form (`blink`, `hsaccade`, `bad`, or any GED class name). |
| `PDL::EEG::SynthNoise` | Noise for synthetic test data. `synth_noise($nt, nch=>, kind=>'white'|'ar1', ar=>, dist=>'gauss'|'uniform', sd=>, seed=>)` returns (nt) or (nt, nch) independent channels, each zero-mean with SD exactly `sd`; AR(1) starts from the stationary distribution. `seed_rng` seeds PDL's generator (`srandom` in PDL 2.090 and later, `srand` before) and Perl's `rand` together; `gaussian`/`uniform` draw `grandom`/`random` with PDL's automatic threading off, so a seed gives the same numbers whatever the number of CPU threads. White uniform noise stays below a `2.5 x 1.4826 x MAD` threshold (as in `detect_blinks`) by itself; Gaussian noise does not. Pure PDL. |
| `PDL::EEG::ElectrodeLocation::Polyline` | Polylines in 2-D or 3-D, as (n, d): `arc_length` (cumulative length at every point), `point_at_length` / `point_at_fraction` (linear interpolation along the line, zero-length segments skipped, clipped to the ends), `project_point` (length along the line and distance from it of the nearest point), `midline_labels('10-20' | '10-10')` (Nasion, Fpz, Fz, Cz, Pz, Oz, Inion at 0, 10, 30, 50, 70, 90, 100 %; 10-10 adds AFz, FCz, CPz, POz). Pure PDL. |
| `PDL::EEG::ElectrodeLocation::TraceLine` | The ordered centre line of a line drawn on an image: `binarize` (threshold on the mean of the colour channels, `invert` for a dark line), `thin_zhang_suen` (Zhang-Suen thinning; only the line pixels are visited, by flat indices and fixed neighbour offsets; equal to a per-pixel implementation in `t/35`), `skeleton_endpoints`, and `centerline` (the longest geodesic path of the largest skeleton component, so side branches drop out; oriented from the chosen end; smoothed by a moving average with the end points fixed). The pixel chain counts slanted stretches up to 8 % long (unit and sqrt(2) steps); the smoothed line is within 0.5 % in `t/35`. Image coordinates: x right, y down, row 0 at the top. |
| `PDL::EEG::ElectrodeLocation::AtlasView` | Side views of an anatomy atlas (screenshots of a 3-D model, the face to the left, e.g. Open Anatomy Atlas / BodyParts3D): `read_png` ((3,nx,ny) 0 .. 1, row 0 at the top), `foreground` (differs from the colour of the top-left pixel), `shift_between` (integer shift between two views from the same camera by phase correlation through `PDL::EEG::FFT::Dispatcher`), `back_outline` (the most posterior body pixel of every row, averaged over 2h+1 rows), `body_outline` (boundary pixels), `refine_tip` (most posterior bone pixel near a guess), `atlas_levels` (inion and spinous-process tips on a bone view, their nearest points on the back outline and lengths from the inion, scaled by inion -> C7 or a given mm per pixel) and `level_of` (the level interval of a length). Tested in `t/36`. |
| `PDL::EEG::Regress` | Temporal nuisance / EOG regression. `regress_out($X, @sources)` removes, from every channel of `$X` (nch,nt), the least-squares fit on a set of source time courses (GED or ICA component activations, or raw EOG channels — with EOG channels it is classic Gratton–Coles regression). Sources are orthonormalised first, so correlated sources are partialled out jointly and the removal depends only on the span of the sources, not on labelling any one vertical/horizontal. Channel DC is preserved. Pure PDL. |
| `PDL::EEG::Heartbeat` | R-peak detection and R-locked averaging. `detect_r($ecg,$fs)` returns R-peak sample indices and R polarity (5–20 Hz band, percentile threshold, 0.30 s minimum distance, ±40 ms refinement on the 0.5 Hz high-passed signal); `beat_epochs` cuts `(nch,nep,nbeat)`; `beat_average` gives the R-locked mean or median, across-beat SD and per-channel beat count, with optional per-beat baseline and per-channel peak-to-peak rejection; `surrogate_r` displaces every R index by a seeded random offset (non-heartbeat time points); `bn_balance($bn1,$bn2)` returns the weight `w` minimising the energy of `w*bn1 + (1-w)*bn2`. Layout `(nch, nt)`. Filtering by `PDL::EEG::Filter`. |
| `PDL::EEG::Heartbeat::Respiration` | Respiratory phase from beat-to-beat series of the heartbeat (ECG-derived respiration). `rr_series($r,$fs)` gives the RR-interval series; `r_amplitude($X,$r,$fs)` the R amplitude of one or several leads at every beat, relative to a window before R; `resp_phase($t,$y)` resamples such a series to 4 Hz, band-passes it (0.20–0.45 Hz by default) and takes its analytic signal, returning also the spectral peak inside the band and the share of power in the band; `phase_at` gives the phase at any time (0 at a maximum of the band-passed series, ±π at a minimum; which of them is inspiration is not determined); `phase_locking($p,$q)` gives the phase-locking value and mean difference of two phase series. One value per heartbeat goes in, so a breathing rate above half the heart rate is not represented. `phase_template($tb,$y,$phi0)` re-estimates the phase at the beats against the shape of one breath (a Fourier template with 2 harmonics fitted on the current phase; per beat, a grid search over the beats within ±2.5 s with a free offset and gain; 3 rounds), using the beat samples themselves so that a harmonic above half the beat rate does not fold into the phase; `phase_template_multi` does the same for several series sharing one phase. `phase_agreement($tb,$p,$q)` gives, at each beat, the circular SD of the difference of two phase estimates over a 120 s window (a constant difference does not count). `fourier_fit`/`fourier_eval`, `wrap_phase`/`unwrap_phase` and `circ_sd` are exported. |
| `PDL::EEG::Heartbeat::BeatFeatures` | `beat_features($x,$fs,r=>$r)` — per-beat shape features of one ECG lead at any sampling rate: each beat (−300..+500 ms) minus the mean of the beats within ±10 s is projected on the average beat, its shift and its stretch in the QRS (−60..+70 ms) and T (+150..+480 ms) windows, giving size (%), shift (ms) and stretch (%) of QRS and T, plus the R amplitude relative to its running mean (%) and the RR intervals ending / starting at the beat (ms). Running means use cumulative sums, so whole-night records need no beats × beats matrix. |
| `PDL::EEG::Filter` | Filters along the time axis, for `(nt)` or `(nch, nt)`. `zp_butter` is a zero-phase Butterworth high-/low-/band-pass applied in the frequency domain (squared-magnitude response of the digital Butterworth filter, odd-reflection padding). `analytic_signal` returns the signal and its Hilbert transform (channel mean removed, zeros added up to a length with prime factors 2, 3, 5; `exact => 1` transforms at the record's own length, which is `scipy.signal.hilbert`); `envelope` and `inst_phase` are its magnitude and angle. The FFTs go through `PDL::EEG::FFT::Dispatcher` (`PDL::FFTDispatcher` when installed (`PDL::PocketFFT`, or the library named by the environment variable `PDL_FFTDISPATCHER_BACKEND`), otherwise `PDL::FFTW3` when installed, otherwise `PDL::FFT` for the lengths it transforms correctly). |
## Command-line tools

| Tool | Role |
|------|------|
| `examples/read_nihonkohden.pl` | Interactive viewer (`--block/--sec/--nch/--chans/--aux`, optional Cairo plot); dispatch is inside `read_nk`, so it needs no format knowledge |
| `examples/nk_to_edf.pl` | NK `.EEG` → EDF/EDF+ (`--subject`, `--equipment`, `--allblocks`, `--gapsamples`; `--rename OLD=NEW,...` or a two-column TSV renames channels before writing, in every mode, and stops on a name that is not in the file or a new name that equals a channel that stays). With `--allblocks`, block starts after the first are written as `Block N` annotations. |
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
| `xt/71_wfdb_real.t` | Real-data check of `PDL::EEG::IO::WFDB` on Challenge 2018 records (pass record stems after `::`): header and signal file agree, the same matrix as `PDL::IO::Matlab` when installed, annotation times in order and within the record |
| `examples/heartbeat_bne_survey.pl` | R-locked survey of a recording with two balanced non-cephalic electrodes (`--bn1/--bn2`, each against the system reference). R peaks are detected on BN1−BN2. Prints (1) P/R/T amplitude ± SE of every channel at the latencies where BN1−BN2 peaks, (2) the weight r minimising `r*BN1 + (1-r)*BN2` in the P, QRS and T windows (with odd/even-beat values and the residual each r leaves in each window), (3) for six references of the scalp channels (system, average, A1+A2, BN r=0.5, BN r=r_QRS, BN r=r_T) the RMS over scalp channels of the QRS peak-to-peak and of the mean amplitude in the T, 200–300, 300–450 and 450–600 ms windows, next to the same measures at `--nsurr` sets of non-heartbeat time points. Writes `PREFIX_avg_<ref>.tsv` and, with `PDL::Graphics::Cairo` + `PDL::EEG::MAP2D`, `PREFIX_bn.png`, `PREFIX_scalp.png` and P/R/T scalp maps for the system and average references. `--skip`, `--band`, `--reject`, `--seed`, `--montage`, `--no-fig`. |
| `examples/heartbeat_hep_scalp.pl` | R-locked scalp maps in the reference the channels were recorded in, tested against non-heartbeat time points. R peaks on BN1−BN2. For the R sample and each `--windows` latency window it prints the map in µV and as z (amplitude / SD over `--nsurr` non-heartbeat averages), Σz² over scalp channels with its rank among the non-heartbeat sets, max \|z\| against the non-heartbeat 95th percentile of max \|z\| and the channels above it, the correlation of the z maps from the first and second half of the beats, the cosine with the map at R, and the share of the map fitted by `a_c*BN1 + b_c*BN2` with coefficients estimated per channel in the QRS window or in `--bnwin` (with the share that fit takes from non-heartbeat averages, and Σz² of what is left). Writes `PREFIX_maps.tsv`, `PREFIX_time.tsv` and, with `PDL::Graphics::Cairo` + `PDL::EEG::MAP2D`, `PREFIX_wave.png` (all scalp channels, RMS against the non-heartbeat 95th percentile, Fz/Cz/Pz with BN1−BN2) and `PREFIX_maps.png` (maps per window: µV, z, µV after each BN fit). |
| `examples/heartbeat_joint.pl` | Draw a tab-separated average (`time_ms` + one column per channel, e.g. the `PREFIX_avg_<ref>.tsv` of `heartbeat_bne_survey.pl`) with `PDL::EEG::MAP2D::plot_joint`: all channels coloured by electrode position with scalp maps above them at the latencies given by `--times` (or by `@TIMES` at the top of the script) and, inside the waveform panel, a small head whose dots give the colour of each electrode's trace. `--xlim`, `--ylim`, `--clim`, `--positive-up`, `--nave`, `--title`, `--montage`, `--out`. |
| `examples/heartbeat_resp_phase.pl` | Per-beat respiratory phase from BN1−BN2 of EDF recordings (`--bn1/--bn2`): the 0.15–0.5 Hz phase of the R size (`resp_phase`), the `phase_template` phase started from it, their `phase_agreement` (q_agree, 120 s window) and a class from q_agree (`--cut 20,25,30` deg); one TSV per recording and a summary line. |
| `examples/heartbeat_resp_vs_belt.pl` | For one Challenge 2018 record (`read_wfdb`): the ECG lead is processed as BN1−BN2 is, the respiratory phase is estimated from the ECG alone (band phase, `phase_template` with several windows and rounds, `phase_template_multi`) and scored against the phase of a breathing belt (`--belt ABD`; inspiration onsets at belt minima): offset, circular SD of the error, share of errors > 45 / 90°, offset SD across 10-min windows, warp over 12 belt-phase bins, belt variance explained. |
| `examples/heartbeat_resp_vs_belt_stage.pl` | The same scores per sleep stage of a Challenge 2018 training record (`.arousal`), in stretches of one stage without a change of the 30 s average QRS (position change flagged from the ECG alone); beats in annotated events and 10 s after them are not scored. |
| `examples/heartbeat_resp_quality_vs_belt.pl` | In 120 s windows of a Challenge 2018 record: the true error of the `phase_template` phase against the belt beside ECG-only quality indices (fold R², fold amplitude / residual, band power share, q_agree, RR variance explained by the phase), with the stage and the share of beats in events per window; TSV per window and a summary. |
| `examples/heartbeat_resp_quality_pool.pl` | Pools the window TSVs of several records (`--maxev` drops windows with a larger share of beats in events): Spearman correlation of each index with the true error over all windows and within records, true error by q_agree bin, and by stage. |
| `examples/heartbeat_sim_bn_resp.pl` | Synthetic BN1/BN2 EDF with a known breathing curve (sine, end-expiratory pause or end-inspiratory hold; RSA, depth modulation, baseline, per-feature lags and lag drift; presets) and a per-beat truth TSV; all random numbers from `PDL::EEG::SynthNoise`. |
| `examples/avg_grid_scrub.pl` | Show a tab-separated average (`time_ms` + one column per channel) as a 5×5 grid of panels laid out like the 10-20 system, all with the same plot area, tick numbers on one panel, and a scalp map (`plot_topomap`, colour bar on its right) at the top centre. In a giza-server window (macOS Cocoa build), a click or drag in a panel chooses the latency: a vertical line in every panel and the map redrawn for it; a double-click opens that channel, or the map, enlarged. The figure title gives the data file and the montage file. Without `--clim` the colour scale of the map follows the map at each latency, rounded up to a multiple of `--cstep` (spacing of the numbers on the colour bar, default 1 µV). `--xlim`, `--ylim`, `--clim`, `--cstep`, `--mark`, `--positive-up`, `--montage`, `--no-map`, `--size`; `--png FILE` writes the grid without a window. Needs a PDL::Graphics::Cairo with `main_scrub`, `cell_margins`, `label_only` and `colorbar(cmap => ..., tick_step => ...)`. Several TSVs are drawn in the same panels with a key (`--names`, `--colors`; the map is that of the first file), and `--cell NAME:ROW:COL` places channels without a 10-20 cell (e.g. BN1, BN2) in free cells. |
| `examples/heartbeat_hep_resp_grid.pl` | R-locked average of the head channels with the cardiac field removed by regression on the averages of cardiac leads (BN1, BN2, and with `--cardiac` more channels or bipolar pairs such as `A1-A2`), with and without a per-beat respiratory variable from the R size. Reference `--ref bne`, `sys`, `a1`, `a2`, `ears`, `chans` or `avg`; `--leads reref` re-references the leads like the head; `--fit A,B` sets the window in which the cardiac field is fitted. Reads `.EEG` (all blocks; beats within 2 s of a recording break are left out), EDF (`Block N` annotations as breaks) and BrainVision (`--ecg`, `--aux`, `--pulse` for a PPG channel). Prints the residual means in 200–300, 300–450 and 450–600 ms with z against averages at uniformly drawn non-heartbeat time points, the slopes on the respiratory variable, and the beat counts at every step (with the channels that throw beats out); writes the averages as TSVs for `avg_grid_scrub.pl` and `plot_joint` PNGs. |
| `examples/hep_reference_t_check.pl` | Runs `heartbeat_hep_resp_grid.pl` on one recording in five conditions (ears with the QRS, P+QRS and whole-beat fit windows, the average reference with P+QRS, and P+QRS with `A1-A2` as a left-right lead) and prints one summary: the residual A2−A1 in each, whether it is the same for the ear and the average reference, and the residual means of chosen channels side by side. Options after `--` go to every run; logs of finished runs are reused. |
| `examples/heartbeat_neck_field_rank.pl` | How many independent time courses the cardiac field has in sets of electrodes (default: chest = ClavL, ClavR, BN2, BN1, Th12; back of the neck = Th12, BN1, Cv45, Cv12; head = the scalp channels + A1, A2), each re-referenced to its own average: for every latency window of the R-locked average, the share of power in the first k principal components and the generalized eigenvalues against averages at non-heartbeat time points. Then, window by window, the share of the head average explained by the first k principal time courses of the chest and neck sets (with the chance level from the non-heartbeat averages and the head's own first k components as the upper bound), the best single body channels and pairs, and the |cos| between the head topographies of different windows. Reads EDF (`Block N` annotations are recording breaks; a block whose `REC START` annotation says `CAL` is left out). Writes the patterns, their time courses and the head's PC patterns per window as TSVs. |
| `examples/nyhead_dipole_gof.pl` | For each scalp topography in a TSV (channel name, then one number per column; `heartbeat_neck_field_rank.pl` writes them as `PREFIX_head_patterns.tsv`), the share of its power explained by a single cortical dipole of the New York Head (`V_fem_normal`, and with `--free` the free-orientation `V_fem`), on the channels that are NYHead electrodes (`--map OLD=NEW`, default T3=T7, T4=T8, T5=P7, T6=P8; others are listed and left out), with leadfield and topography re-referenced to the average of those channels: the best sources with their MNI position and Harvard–Oxford area, the share of the cortex within 90 % of the best fit, and the distribution of the best fit for random spatially white topographies on the same channels. Stops with the column and channel when a value is not a number. `--selftest` checks the fit on a synthetic leadfield with planted dipoles. Needs `PDL::IO::NYHead` + `sa_nyhead.mat`. |
| `examples/nyhead_alpha_sim.pl` | A pure posterior alpha rhythm (one sine per hemisphere, no noise) on the New York Head, written as EDF, one file per reference (`--refs`, default AVG, Cz, Cv12, Cv45, Nk2, BN1, Th12; the reference channel is left out of its file). Sources: equivalent dipoles (free-orientation `V_fem`) at the cortical vertex nearest to the centroid of each `--sites` entry in both hemispheres (default parieto-occipital sulcus, intracalcarine cortex, occipital pole, lateral occipital cortex inferior and superior; the right dipoles mirror the left ones); with one site the orientation is swept, with several the moments are drawn (seeded, `PDL::EEG::SynthNoise`) in batches until `--target-kept` configurations pass the test of a posterior alpha map on the 19 10-20 channels against `--crit-ref` (largest at O1/O2, O1 O2 P3 P4 Pz of one sign, no increase towards the front along O-P-C-F-Fp, Pz-Cz-Fz and P7-T7-F7 on both sides). Own points (`--points`, e.g. the `--points-out` file of PDL-IO-NYHead `nyhead_neck_electrodes.pl`) get their values by `PDL::IO::NYHead::Interp` from the electrode values. Prints, over the kept configurations and for two halves of them, each channel against `--far` as a share of O1 - `--far`, and the differences of candidate references; for the EDF configuration, the amplitude of the main channels in every reference. `--source patch` uses cortical patches with normal current instead. `--falloff PREFIX` writes, for every channel and each `--falloff-refs` reference (default Fpz, BN1, Cv45), the alpha amplitude (`--falloff-y amp`, default) or power relative to O1 in that reference over the kept configurations (median and 10 .. 90 %) against three distances (straight from d, the electrode nearest to the left parieto-occipital dipole; straight from the nearest dipole; along the midline of the head surface from d, continued down the neck by the profile rows of `--points`) to `PREFIX_falloff.tsv`, and draws them in `PREFIX_falloff_d.png`, `_src.png` and `_midline.png` (linear axes; the lower row enlarges the y axis to `--falloff-zoom`; left/right pairs at the same spot share one label). `--selftest` checks the EDF round trip, the patch selection, the dipole maps against the formula on a sphere, and the falloff distances and ratios on a sphere. Needs `PDL::IO::NYHead` 0.04 + `sa_nyhead.mat`. |
| `examples/trace_midline.pl` | Midline electrode positions (10-20, or 10-10 with `--system 10-10`) along a line drawn on a PNG (a white line on black, or `--invert`), e.g. the skin of the midline traced over a sagittal view: the line is thinned, ordered from the nasion end (`--start left|right|top|bottom`) to the inion end, smoothed over `--smooth` points (default 5), and the positions are put at fractions of its length. Prints and writes to `PREFIX_midline.tsv` the image coordinates (sub-pixel), the length from the nasion end in pixels (and mm with `--length-mm`, the measured nasion-inion length), and how far each position moves when measured along the pixel chain instead; draws `PREFIX_midline.png` over the line or over `--background` (an image of the same size), with `--no-labels`, `--no-marks`, `--no-line`. `--selftest` runs on a synthetic semicircle with a side branch. Uses `PDL::IO::PNG` (else `PDL::IO::Pic`) and, for the PNG, `PDL::Graphics::Cairo`. |
| `examples/atlas_neck_levels.pl` | The length along the skin of the back of the neck from the inion to each spinous process (C1 .. T3) on two atlas side views from the same camera (SKIN with the skin, BONE with the bones only, cropped differently): the views are aligned by phase correlation, the back outline is taken from SKIN, the inion (`--inion-row`) and the tips (`--tips` guesses moved to the most posterior bone pixel) from BONE, and the lengths are scaled by `--c7-mm` (inion -> C7 measured on a subject, default 165) or `--mm-per-px`. Prints the tips, lengths and depths under the skin, and with `--mark 'NAME=MM,...'` the level interval of each length; writes `PREFIX_neck_profile.tsv` (the back outline in mm from the inion and the levels, read by PDL-IO-NYHead `nyhead_neck_electrodes.pl --extend-with`) and `PREFIX_neck_levels.png`. `--selftest` on synthetic views. |
| `examples/nyhead_atlas_sagittal.pl` | The New York Head laid over the same atlas side views in the midsagittal plane. The head: Iz and Nz (or Fpz) on the atlas inion and nasion (`--head-fit IN`, default; `--nasion` in SKIN pixels), or the `/sa/head` cut at x = 0 above Iz fitted to the atlas outline by iterated closest points (`all`; `back` for the back half with the scale held); all three are printed with the distances of Iz, Nz and Exxz from the atlas landmarks and back outline. In each body the angle between the inion-nasion line and a back-of-the-neck line (NYHead `--ny-neck`, default from the height of Exxz to Nk2; atlas `--atlas-neck`, default 85 -> 125 mm at the subject scale), and the NYHead neck turned (about the upper end of its chord) to make the two neck lines parallel. With `--neck arc` (default) the NYHead neck electrodes and the `--points` sites are put on the atlas back outline at their length along the NYHead back midline and their vertebral level is read at the subject scale (`--c7-mm`) and at the head scale. Writes `PREFIX_atlas_in_nyhead.tsv` (the atlas levels, tips and back outline in NYHead MNI mm) and `PREFIX_sagittal.png` (`--no-labels` for the marks without text). `--selftest` checks the fits, the inverse of the neck turn and the angle finding on synthetic outlines. Needs `PDL::IO::NYHead` + `sa_nyhead.mat`. |

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

### Waveforms with maps at chosen latencies — `plot_joint`

```perl
use PDL::EEG::MAP2D qw(plot_joint);
plot_joint(
    avg => $avg,                 # (nchan, ntime)
    t   => $t_ms,                # (ntime)
    labels => \@labels, montage => $elc,
    times  => [200, 250, 300, 350, 400],
    nave   => 293, outfile => 'joint.png',
);
```

Draws every channel that has a position in the montage as one trace, coloured by
electrode position (right → red, anterior → blue, superior → green), with a
scalp map for each entry of `times` above the waveforms and a line from each map
to its latency. The maps share one symmetric colour scale (`clim`, default the
largest |amplitude| at the requested latencies); a small head inside the waveform
panel shows which colour is which electrode. `negative_up` (default 1), `xlim`,
`ylim`, `unit`, `title`, `time_format`, `key => 0`, `color` (one colour for all
traces), `width`/`height`. Without `outfile` the figure is returned. Needs
`PDL::Graphics::Cairo` with `Figure::add_connection` and `add_text`.

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
frequency domain via `PDL::EEG::FFT::Dispatcher` (`PDL::FFTDispatcher` when installed (`PDL::PocketFFT`, or the library named by the environment variable `PDL_FFTDISPATCHER_BACKEND`), otherwise `PDL::FFTW3` when installed, otherwise `PDL::FFT` for the lengths it transforms correctly), so there is no external wavelet dependency. It
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
PDL. The symmetric eigendecomposition uses LAPACK `syevd` through
`PDL::LinearAlgebra` 0.436 or later when it is installed, and `PDL::MatrixOps`
`eigens_sym` otherwise (`PDL::EEG::Linalg`; one warning per process in that case). The real New York Head leadfield is loaded through `PDL::IO::NYHead`.

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
make test        # 36 files (t/06 reserved/skipped)
```

`t/00_testdata` (writes the synthetic fixtures with `t/mk_synthetic_nk.pl` when
any is missing) `t/01_nihonkohden` `t/02_edf` `t/03_ptn` `t/04_signal` `t/05_montage`
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
`t/22_heartbeat` (Heartbeat: `detect_r` against the reference R samples in both polarities,
`beat_epochs`/`beat_average` against explicit sums incl. baseline, rejection and
median, `surrogate_r` seeding and bounds, `bn_balance` closed form; deterministic,
needs only PDL).
`t/23_joint` (MAP2D `plot_joint`: renders with and without `outfile`, one
connection per latency, panel count, argument checks; skipped without
`PDL::Graphics::Cairo` providing `Figure::add_connection` and `add_text`).
`t/24_filter` (Filter: `zp_butter` against scipy `filtfilt` values and the
analytic gain; `analytic_signal` with `exact` against scipy `hilbert` values for
an even and an odd length, edges included; `envelope` and `inst_phase` against
the amplitude and phase of a synthetic modulated carrier, with an unmodulated
carrier as control; `(nch,nt)` input, a prime length, a channel offset;
deterministic, needs only PDL).
`t/25_respiration` (Heartbeat::Respiration: a synthetic recording with a known
breathing phase drifting from 0.27 to 0.35 Hz that modulates the RR interval and
the R amplitude of two leads; `rr_series` and `r_amplitude` against the values
built in, the phase from each series against the true phase, the lag built in
between the two, an RR series without breathing as control, `phase_locking`
closed forms, argument checks; deterministic, needs only PDL).
`t/31_resp_template` (Heartbeat::Respiration: beats with a breathing rate drifting
from 0.27 to 0.31 Hz and a non-sinusoidal breath whose second harmonic lies above
half the beat rate; `phase_template` and `phase_template_multi` against the true
phase and below the band-phase error, the template's harmonic ratio, no second
harmonic for a sinusoidal breath, a pure-noise series gated out, `phase_agreement`
closed forms, `fourier_fit`/`circ_sd`/`wrap_phase` checks; seeded through
SynthNoise, needs only PDL).
`t/32_beat_features` (Heartbeat::BeatFeatures: a synthetic lead at 500 and
1000 Hz whose QRS height and width and T height and shift change by known amounts
over a 5-beat cycle; each feature against its built-in change, unchanged
features small, RR; deterministic, needs only PDL).
`t/33_wfdb` (IO::WFDB: records written in the test as Matlab v4 `.mat` and `.dat`,
a file cut short, refused formats; an annotation file with a negative SKIP, aux
texts, a NUM word and a code-0 annotation; stages and events; needs only PDL).
`t/35_electrode_location` (ElectrodeLocation: polyline lengths, interpolation past
a zero-length segment, projection; thinning against a per-pixel Zhang-Suen; the
centre line of a semicircle with a side branch, an ellipse and a curve with a flat
front and a round back: length within 1 %, positions within 0.4 % of the length,
smoothing closer than the pixel chain on the uneven curve).
`t/36_atlas_view` (ElectrodeLocation::AtlasView on synthetic side views: the
body of a rectangle and its boundary, the back outline of a slanted edge against
the mean of the last body column, a tip found from a guess and none in an empty
box, the shift of a textured image and a crop of it, the inion, levels and lengths
along a slanted back within 2 px, `level_of`, the `--c7-mm` scale, and a PNG read
back with `PDL::IO::PNG` when it is installed).

`xt/70_real_data.t` is a real-data regression (not part of `make test`; needs
private recordings). Pass `.EEG` paths and it checks event placement on real
`extblock` and `wfmblock` files — every `REC START` on a block boundary, every
event inside a segment:

```
prove -lv xt/70_real_data.t :: /path/A.EEG /path/B.EEG
```

`xt/71_wfdb_real.t` checks `PDL::EEG::IO::WFDB` on real Challenge 2018 records
(not part of `make test`):

```
prove -lv xt/71_wfdb_real.t :: /data/tr03-0005 /data/te03-0024
```

The synthetic Nihon Kohden fixtures (`t/data/test01.eeg`, `test01.LOG`,
`test02.eeg`, `test03.eeg`) are not in the repository: `t/00_testdata` writes them
with `perl t/mk_synthetic_nk.pl` when any is missing. `t/03_ptn` checks a real
`.PTN` only when `t/data/Pattern_032.PTN` is present (not distributed). `--long[=SEC]` also
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
