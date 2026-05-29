# nestapp

![version](https://img.shields.io/badge/version-1.0.0-blue)
![MATLAB](https://img.shields.io/badge/MATLAB-R2023b%2B-orange)
![license](https://img.shields.io/badge/license-GPL--3.0-green)

**nestapp** is a MATLAB GUI for cleaning and analysing TMS-EEG recordings. It wraps [EEGLAB](https://eeglab.org/), [TESA](https://nigelrogasch.gitbook.io/tesa-user-manual/), and [FastICA](http://research.ics.aalto.fi/ica/fastica/) into a point-and-click pipeline builder that lets researchers process data without writing code.

**Contributing** → [CONTRIBUTING.md](CONTRIBUTING.md) · **Changes** → [CHANGELOG.md](CHANGELOG.md) · **Cite** → [CITATION.cff](CITATION.cff) · **Architecture** → [docs/architecture.md](docs/architecture.md)

---

## Requirements

| Dependency | Version |
|---|---|
| MATLAB | R2023b or later (developed on R2025b) |
| EEGLAB | 2025.0.0 |
| TESA | 1.1.1 |
| FastICA | 2.5 |
| Signal Processing Toolbox | any recent |
| Statistics and Machine Learning Toolbox | any recent |
| Curve Fitting Toolbox | required only for the `TMS-EEG / AARATEP` template (`Remove Decay Artifact` step) |

**EEGLAB** — download from [eeglab.org](https://eeglab.org/) and unzip anywhere on your machine. nestapp will add it to the MATLAB path automatically once you set the folder in Preferences (see below).

**TESA** — install from inside EEGLAB via **File → Manage EEGLAB extensions**. Search for TESA and install. TESA lives inside the EEGLAB plugins folder, so no separate path setup is needed after that.

**FastICA** — download from [research.ics.aalto.fi/ica/fastica](http://research.ics.aalto.fi/ica/fastica/) and add the folder to your MATLAB path (Home tab → Set Path → Add Folder). FastICA must be on the path before launching nestapp.

**AARATEP** *(only needed for the `TMS-EEG / AARATEP` template)* — the AARATEP pipeline's helper functions are not bundled with this repository. Clone them into `third_party/aaratep/` from the project root:
```bash
cd third_party
git clone --depth 1 https://github.com/chriscline/AARATEPPipeline.git aaratep
```
nestapp adds this tree to the MATLAB path automatically when an AARATEP step runs. The AARATEP template also requires the **Curve Fitting Toolbox** (for the `Remove Decay Artifact` step). See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for the pinned commit and licensing.

---

## Getting started

**First-time setup (do this once):**

1. Install EEGLAB, TESA, and FastICA as described above.
2. Add FastICA to your MATLAB path via **Home → Set Path → Add Folder**, then save.
3. Clone or download this repository and open MATLAB in the project root.
4. Run the entry point:
   ```matlab
   run_nestapp
   ```
5. Open **Settings → Preferences** and set the EEGLAB installation folder. nestapp will add EEGLAB to the path automatically on every launch after this.

**Typical workflow:**

1. **File → Load Template** to start from a ready-made pipeline, or build one from scratch in the Cleaning tab.
2. Select your data files and press **Run Analysis**.
3. Review the per-file summary in the **Reports** tab.
4. Load the processed files in the **Visualizing** tab to inspect TEPs and topographies.
5. Switch to the **Analysis** tab to view detected TEP components and export peak data to CSV.

---

## Cleaning tab

### Building a pipeline

- The left panel lists all available processing steps. Click a step to see a description and parameter details in the **Info** panel.
- Use **Add** (or double-click a step) to append it to the **Selected Steps** list.
- Use **Remove**, **Move Up**, and **Move Down** to edit the order.
- Click any step in **Selected Steps** to view and edit its parameters in the table below. **Default Value** resets parameters for that step.
- **Re/Start Steps** clears the pipeline after a confirmation prompt.

### Pipeline templates

**File → Load Template** provides ready-to-use starting points:

| Template | Description |
|---|---|
| TMS-EEG / TEP (TESA) | Full TESA two-round-ICA artifact removal for single-pulse TMS (Rogasch 2017) |
| TMS-EEG / TEP (TESA + Quality Gates) | TESA workflow with four Quality Gate checkpoints |
| TMS-EEG / ARTIST | Paper-faithful Wu 2018 pipeline; TESA compselect substitutes for the unreleased FLDA classifier |
| TMS-EEG / AARATEP | Cline 2021 pipeline (SOUND + decay-fit removal + AR-blend interp + TMS-aware muscle classifier); helpers vendored under `third_party/aaratep/` |
| Resting-State EEG | PREP + Delorme 2023 cleaning with ICLabel |
| Minimal (Delorme 2023) | HPF + bad channels + ICA, minimum-handling philosophy |

### Citing nestapp

When you publish results produced with a built-in template, cite the corresponding paper. The template name itself is logged at the start of every batch run alongside the DOI and full reference, so the citation ends up in the run log next to your data.

| Template | Primary citation | Additional |
|---|---|---|
| TMS-EEG / TEP (TESA) | Rogasch et al. (2017). NeuroImage 147:934-951. doi:10.1016/j.neuroimage.2017.06.014 | — |
| TMS-EEG / TEP (TESA + Quality Gates) | Rogasch et al. (2017). doi:10.1016/j.neuroimage.2017.06.014 | — |
| TMS-EEG / ARTIST | Wu et al. (2018). Hum Brain Mapp 39(4):1607-1625. doi:10.1002/hbm.23938 | Rogasch et al. (2017) — TESA compselect substitutes for ARTIST's unpublished FLDA classifier |
| TMS-EEG / AARATEP | Cline et al. (2021). IEEE NER. doi:10.1109/NER49283.2021.9441147 | Mutanen et al. (2018) NeuroImage 166:135-151 doi:10.1016/j.neuroimage.2017.10.021 — SOUND algorithm. Rogasch et al. (2017) — TESA helpers vendored from AARATEPPipeline depend on TESA. |
| Resting-State / Minimal | Delorme (2023). Sci Rep 13:2372. doi:10.1038/s41598-023-27528-0 | — |

See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for the complete attribution of vendored and bundled dependencies, including the EEGLAB stack (Delorme & Makeig 2004), ICLabel (Pion-Tonachini et al. 2019), CleanLine, firfilt, clean_rawdata, and FastICA citations.

### Saving and loading pipelines

- **File → Save Pipeline** saves the current step list and parameters as a `.mat` file.
- **File → Load Pipeline** restores a previously saved pipeline.
- Recent pipelines are listed under **File → Recent Pipelines**.

### Running the pipeline

1. Select the data files to process (`.set`, `.vhdr`, `.cnt`, `.cdt`).
2. Press **Run Analysis**. A progress dialog shows the current file and step; processing can be cancelled between steps.

After the run, the **Reports** tab opens automatically with a per-file summary including channel counts, trial retention, and ICA statistics.

### Pipeline provenance

Every processed file has the full pipeline — steps, parameters, and a timestamp — written into `EEG.history` and preserved inside the saved `.set` file. After a run, the final processed dataset is available in the MATLAB base workspace. To inspect the processing record:

```matlab
EEG.history      % full log as a string
eegh             % browse interactively using EEGLAB's history viewer
```

---

## Visualizing tab

The Visualizing tab works on processed `.set` files — either output from a nestapp pipeline run or any EEGLAB-compatible epoched dataset.

- Select one or more `.set` files.
- Click electrode buttons to define the region of interest (ROI).
- **PLOT TEP** plots the trial-averaged waveform with a shaded SEM band. Multiple files can be overlaid on the same axes.
- **Show Components** detects and overlays the six canonical TEP components (N15, P30, N45, P60, N100, P180) with latency and amplitude labels.
- **TOPOPLOT** plots the scalp topography at a selected time point and window.
- **Export TEP Figure** saves the current plot as PNG, PDF, or `.fig`.
- The **TEP Window** slider sets the time range shown in the plot.

---

## Analysis tab

The Analysis tab works on the files and ROI selected in the Visualizing tab — set those up first, then switch here.

### TEP component table

After clicking **PLOT TEP** in the Visualizing tab, the table populates automatically with the detected latency and amplitude of the six canonical TEP components (N15, P30, N45, P60, N100, P180). Components not found in the current waveform are shown as —.

**Edit Component Windows** opens a dialog to adjust the search window for each component. **Reset Defaults** restores the canonical windows from Beck et al. (2024, *Hum Brain Mapp*, 45:e70048). Changes are applied immediately to the current plot.

### Exporting data

**Export TEP to Workspace** saves the grand-mean ROI waveform as a MATLAB variable. Set the variable name in the field next to the button before clicking.

**Extract Peaks → CSV** runs batch peak extraction across all selected files and saves a long-format CSV suitable for import into R, SPSS, or Excel. Each row is one file × component combination, with columns for latency, amplitude, search window, trial count, and ROI channel list. Any extraction warnings (missing electrodes, failed detections) are listed in an alert after the run.

---

## Reports tab

After each pipeline run, a summary report is added to the **Reports** tab. Each report shows channel counts, trial retention, and ICA component statistics for one file. Reports can be browsed, copied as a ready-to-paste methods paragraph, and exported as CSV for multi-subject batch summaries.

---

## Preferences

Open **Settings → Preferences** to configure:

- EEGLAB installation path
- Default data and pipeline folders
- Whether to show the Reports tab automatically after each run
- Whether to require confirmation before clearing a pipeline

---

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for the development setup, how to run the test suite, the pipeline-step extension recipe, and the pull-request process. New to the codebase? Start with [docs/architecture.md](docs/architecture.md) for a map of what each module does and where to make changes.

## Citing nestapp

If you use nestapp in published work, please cite it using the metadata in [CITATION.cff](CITATION.cff) (GitHub's "Cite this repository" button generates APA/BibTeX for you). Also cite the pipeline papers for any built-in template you run — these are listed per template under **Pipeline templates → Citing nestapp** above and are logged to the run output. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for all bundled and vendored dependencies.

## Authors

**Aref Pariz** — original application (2023), developed at the Royal Institute for Mental Health in Dr. Sara Tremblay's lab ([NESTLAB](https://www.nest-lab.ca/)) and Dr. Jeremie Lefebvre's Lab, University of Ottawa.

**Wesley Dunne** — pipeline engine, progress reporting and provenance, automated reports with ICA tracking, TEP visualisation and batch peak extraction, pipeline templates, quality-control gates, and test suite. Led the 1.0.0 open-source release.
