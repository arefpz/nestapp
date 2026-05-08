# nestapp

**nestapp** is a MATLAB GUI for cleaning and analysing TMS-EEG recordings. It wraps [EEGLAB](https://eeglab.org/), [TESA](https://nigelrogasch.gitbook.io/tesa-user-manual/), and [FastICA](http://research.ics.aalto.fi/ica/fastica/) into a point-and-click pipeline builder that lets researchers process data without writing code.

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

**EEGLAB** — download from [eeglab.org](https://eeglab.org/) and unzip anywhere on your machine. nestapp will add it to the MATLAB path automatically once you set the folder in Preferences (see below).

**TESA** — install from inside EEGLAB via **File → Manage EEGLAB extensions**. Search for TESA and install. TESA lives inside the EEGLAB plugins folder, so no separate path setup is needed after that.

**FastICA** — download from [research.ics.aalto.fi/ica/fastica](http://research.ics.aalto.fi/ica/fastica/) and add the folder to your MATLAB path (Home tab → Set Path → Add Folder). FastICA must be on the path before launching nestapp.

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

**File → Load Template** provides three ready-to-use starting points:

| Template | Description |
|---|---|
| TMS-EEG (TESA) | Full TESA artifact-removal workflow for single-pulse TMS |
| Resting-State EEG | Standard resting-state cleaning with ICA |
| Minimal (HPF + bad channels) | High-pass filter and bad channel removal only |

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

## Contributors

**Aref Pariz** — original application (v1.0, 2023), developed at the Royal Institute for Mental Health in Dr. Sara Tremblay's lab ([NESTLAB](https://www.nest-lab.ca/)) and Dr. Jeremie Lefebvre's Lab, University of Ottawa.

**Wesley Dunne** — v2.0 (2025): pipeline engine, progress reporting and provenance, automated reports with ICA tracking, TEP visualisation and batch peak extraction, pipeline templates, test suite.
