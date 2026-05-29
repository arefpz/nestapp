# Third-Party Notices

nestapp ships built-in templates that wrap or vendor code from external research projects. This file enumerates those dependencies, their licenses, and what was copied.

When you publish results produced with one of nestapp's built-in templates, cite the corresponding paper from the table in [`README.md`](README.md#citing-nestapp). The citation is also printed to the batch log on every run (see `templateCitation.m`).

## Pipeline-defining packages

These external packages each define the algorithm and parameters of a built-in nestapp template. They are treated as first-class citations: every batch run that uses one of these templates logs the citation at the start of the run, and every publication that uses one of these templates must cite the corresponding paper.

### TESA (TMS-EEG Signal Analyser)

- **Templates that use it:** `TMS-EEG / TEP (TESA)`, `TMS-EEG / TEP (TESA + Quality Gates)`, *plus most steps of* `TMS-EEG / ARTIST` *and* `TMS-EEG / AARATEP`.
- **Upstream:** https://nigelrogasch.github.io/TESA/ — code at https://github.com/nigelrogasch/TESA
- **License:** GPL-3.0
- **Installed via:** EEGLAB Plugin Manager (bundled with nestapp under `eeglab2026.0.0/plugins/`); not vendored under `third_party/`.
- **What nestapp uses:** `pop_tesa_findpulse` (TMS pulse detection), `pop_tesa_removedata` (artifact removal), `pop_tesa_interpdata` (cubic interpolation), `pop_tesa_filtbutter` (zero-phase Butterworth filtering, also used for the 60 Hz notch in ARTIST and AARATEP templates), `pop_tesa_fastica` (FastICA wrapper), `pop_tesa_compselect` (six-detector TMS-EEG IC classifier — also serves as the FLDA fallback in our ARTIST template), `pop_tesa_sound` (SOUND algorithm wrapper, used by AARATEP template), `pop_tesa_sspsir` (SSP-SIR), `pop_tesa_peakanalysis`/`pop_tesa_peakoutput` (TEP peak extraction).
- **Template design provenance:** the `TMS-EEG / TEP (TESA)` template's step order is annotated against specific steps of the TESA User Manual (Rogasch's published preprocessing recipe) — see `src/buildTemplates.m` lines 19–79 for the manual-step cross-references.
- **Cite as:** Rogasch N.C., Sullivan C., Thomson R.H., Rose N.S., Bailey N.W., Fitzgerald P.B., Farzan F., Hernandez-Pavon J.C. (2017). Analysing concurrent transcranial magnetic stimulation and electroencephalographic data: a review and introduction to the open-source TESA software. *NeuroImage* 147:934-951. doi:[10.1016/j.neuroimage.2017.06.014](https://doi.org/10.1016/j.neuroimage.2017.06.014)
- **Additional citation for SOUND:** Mutanen T.P., Metsomaa J., Liljander S., Ilmoniemi R.J. (2018). Automatic and robust noise suppression in EEG and MEG: The SOUND algorithm. *NeuroImage* 166:135-151. doi:[10.1016/j.neuroimage.2017.10.021](https://doi.org/10.1016/j.neuroimage.2017.10.021) — cite this *in addition to* the TESA paper whenever your pipeline includes `Remove Recording Noise (SOUND)` (this includes the `TMS-EEG / AARATEP` template).

## Vendored code (copies under `third_party/`)

### AARATEPPipeline

- **Path in this repo:** `third_party/aaratep/`
- **Upstream:** https://github.com/chriscline/AARATEPPipeline
- **Pinned commit:** `be75262af689d4e8e5053c05aaa4ed3be258350a` (2025-08-29)
- **License:** MIT (see `third_party/aaratep/LICENSE`)
- **Copyright:** © 2021 Chris Cline
- **What was copied:** the entire repository at the pinned commit, including the `Common/` subtree (~280 files). Nothing was modified in place.
- **What invokes it:** dispatch cases in `src/processOneFile.m` call `c_TMSEEG_fitAndRemoveDecayArtifact` and `c_EEG_ReplaceEpochTimeSegment` directly. `src/ensureAaratepOnPath.m` adds the vendored tree to the MATLAB path on first call.
- **Derivative work:** `src/aaratepMuscleClassifier.m` ports a 12-line block from `c_TMSEEG_Preprocess_AARATEPPipeline.m` (lines 400-412) under the same MIT license. The header comment in `aaratepMuscleClassifier.m` credits the origin.
- **Cite as:** Cline C.C. et al. (2021). Advanced Artifact Removal for Automated TMS-EEG Data Processing. *IEEE NER*. doi:[10.1109/NER49283.2021.9441147](https://doi.org/10.1109/NER49283.2021.9441147)

## Paper-derived implementations (not vendored)

### ARTIST

- **Reference:** Wu W., Keller C.J., Rogasch N.C., Longwell P., Shpigel E., Rolle C.E., Etkin A. (2018). ARTIST: A fully automated artifact rejection algorithm for single-pulse TMS-EEG data. *Hum Brain Mapp* 39(4):1607-1625. doi:[10.1002/hbm.23938](https://doi.org/10.1002/hbm.23938)
- **Public source:** the paper's distribution URL (`etkinlab.stanford.edu/toolboxes/ARTIST/`) is no longer reachable. No GitHub mirror with an open-source license has been located. nestapp does **not** vendor ARTIST code.
- **What nestapp ships:** three pure-MATLAB helpers (`src/artistRejectBadTrials.m`, `src/artistBadChannelsRansac.m`, `src/artistFlagDecayICs.m`) that implement the deterministic stages of the paper (z-score epoch rejection, RANSAC bad-channel rejection, decay-IC magnitude threshold). Numeric defaults are taken from Wu 2018 §2.2.1–§2.2.3 and the source verses are cited in each function header.
- **What nestapp does NOT ship:** the 23-feature Fisher LDA IC classifier from §2.2.1. Reproducing it requires labeled training data not distributed with the paper. The `TMS-EEG / ARTIST` template substitutes TESA's `pop_tesa_compselect` for stage 3 — closer in spirit than vanilla ICLabel because both are hand-crafted TMS-EEG-aware classifiers, but not numerically equivalent to FLDA. This fallback is called out in the template name and in `templateCitation.m`.

#### ARTIST template fidelity gaps

The `TMS-EEG / ARTIST` template follows Wu 2018 §2.2 verbatim where the registry supports it. Known deltas:

| Where paper says | nestapp template does | Why / how to close |
|---|---|---|
| Round-2 ICA: PCA reduction to retain 99.9% variance, discard ICs with < 0.2% variance (§2.2.1) | Runs runica without PCA dimensionality reduction | The `Run ICA` registry step does not yet expose a `pca` parameter. Adding it would mean extending `stepRegistry.m`, the dispatch case in `processOneFile.m`, and the parameter type metadata. Defer until needed. |
| Stage-3 classifier: 23-feature Fisher LDA on Infomax ICs (§2.2.1) | TESA `pop_tesa_compselect` with all detectors on | Original FLDA weights are not publicly distributed (etkinlab.stanford.edu URL is dead). Contact `wwumed@stanford.edu` to obtain, or train from labeled data. |
| Notch: "60 Hz zero-phase FIR notch filter" (§2.2.2) | TESA Butterworth bandstop 58–62 Hz, order 2 (zero-phase via filtfilt) | EEGLAB's pop_cleanline spans trial boundaries on epoched data and prompts the user; TESA bandstop is the same algorithm class (zero-phase, ~symmetric stop band) and safe on epoched data. |

### AARATEP template fidelity gaps

The `TMS-EEG / AARATEP` template's step order matches `c_TMSEEG_Preprocess_AARATEPPipeline.m` line-for-line. Known parameter deltas:

| Where upstream does | nestapp template does | Why / how to close |
|---|---|---|
| Bad channels via `c_TMSEEG_detectBadChannels` with method `{'PREP_deviation','TESA_DDWiener_PerTrial'}` ensemble | nestapp's `Remove Bad Channels` (pop_rejchan kurtosis, threshold 10) | Different algorithm. Same threshold value. To close: add a `Remove Bad Channels (AARATEP)` step dispatching to the vendored `c_TMSEEG_detectBadChannels` helper, similar to how `Remove Decay Artifact` dispatches to `c_TMSEEG_fitAndRemoveDecayArtifact`. |
| Notch via `c_EEG_filter_butterworth` bandstop at `[58, 62] Hz` plus optional harmonics from `lineNoiseNumHarmonics` | TESA Butterworth bandstop 58–62 Hz, order 2, single harmonic | Functionally equivalent for the default `lineNoiseNumHarmonics = 1` case. Multi-harmonic support would need either iterated TESA bandstop calls or a registry extension. |
| HPF via `c_TMSEEG_applyModifiedBandpassFilter` (TMS-aware filter with custom artifact handling) | EEGLAB `pop_eegfiltnew` (firfilt) via `Frequency Filter` | The "modified" prefix in upstream suggests special TMS-edge handling we don't reproduce. Practical impact is mild because AR-Blend interpolation already replaces the artifact samples before filtering. |

**Required MATLAB toolbox for AARATEP:** the `Remove Decay Artifact` step calls MATLAB's `fit()` function (Curve Fitting Toolbox) with constrained nonlinear exponential decay models (see `c_TMSEEG_fitAndRemoveDecayArtifact.m` lines 101 / 107). **Curve Fitting Toolbox must be installed** to run the AARATEP template end-to-end. If it isn't, the pre-flight check (`checkStepDependencies.m`) blocks the run with an install message before the pipeline starts. Workaround: remove `Remove Decay Artifact` from the pipeline (AARATEP will still run but the decay-fit cleanup step is skipped — results will differ from the published pipeline).

## Other bundled EEGLAB stack

These packages are not vendored under `third_party/` but are bundled with the nestapp distribution under `eeglab2026.0.0/`. They are installed and updated through EEGLAB's own plugin manager. Cite each according to its own conventions when used. (TESA is documented above under "Pipeline-defining packages" because it defines two of the built-in templates.)

| Package | Upstream | License | Cite as |
|---|---|---|---|
| EEGLAB | https://eeglab.org | BSD | Delorme & Makeig (2004). *J Neurosci Methods* 134(1):9-21. doi:10.1016/j.jneumeth.2003.10.009 |
| ICLabel | https://sccn.ucsd.edu/wiki/ICLabel | BSD | Pion-Tonachini, Kreutz-Delgado, Makeig (2019). *NeuroImage* 198:181-197. doi:10.1016/j.neuroimage.2019.05.026 |
| CleanLine | EEGLAB plugin manager | BSD | Mullen T. (2012). NITRC: CleanLine. |
| firfilt | EEGLAB plugin manager | GPL-2.0 | Widmann A., Schröger E., Maess B. (2015). *J Neurosci Methods* 250:34-46. doi:10.1016/j.jneumeth.2014.08.002 |
| clean_rawdata | EEGLAB plugin manager | BSD | Mullen et al. (2015). *IEEE Trans Biomed Eng* 62(11):2553-2567. doi:10.1109/TBME.2015.2481482 |
| FastICA | http://research.ics.aalto.fi/ica/fastica/ | GPL-2.0 | Hyvärinen & Oja (2000). *Neural Networks* 13(4-5):411-430. doi:10.1016/S0893-6080(00)00026-5 |
