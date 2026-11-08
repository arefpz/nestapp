# nestapp

![version](https://img.shields.io/badge/version-2.2.1-blue)
![MATLAB](https://img.shields.io/badge/MATLAB-R2023b%2B-orange)
![license](https://img.shields.io/badge/license-GPL--3.0-green)

**nestapp** is a MATLAB app for cleaning and analysing TMS-EEG recordings. It wraps [EEGLAB](https://eeglab.org/), [TESA](https://nigelrogasch.gitbook.io/tesa-user-manual/) and the surrounding plugin stack in a point-and-click pipeline builder, so a researcher can preprocess a cohort, check what happened to every file, and get publication-ready TEP figures and peak measurements without writing code.

Three tabs, in the order you use them:

- **Cleaning** builds a preprocessing pipeline, or loads a published one as a template, and runs it over a batch of files.
- **Reports** shows what the run did to each file: channels dropped, trials kept, ICA components removed, and a methods paragraph you can paste into a manuscript.
- **Explore** supports the comparison most TMS-EEG studies are designed around, pre versus post or one cohort against another, with TEP waveforms, topographies, and per-window measures exported to CSV.

Every processed file carries its full pipeline inside `EEG.history` in the saved `.set`: the steps, their parameters, and a timestamp. A result can therefore always be traced back to how it was made.

---

## Installation

About ten minutes, most of it downloads. Everything below is free except
MATLAB itself.

**Before you start:** you need MATLAB R2023b or later with the **Signal
Processing** and **Statistics and Machine Learning** toolboxes. Type `ver` at
the MATLAB prompt to see what you have; add a missing toolbox from
**Home → Add-Ons → Get Add-Ons**.

### 1. Install nestapp

Download the `.mltbx` file from the
[latest release](https://github.com/arefpz/nestapp/releases/latest) and
double-click it (or drag it into the MATLAB Command Window). MATLAB installs
nestapp, keeps it on the path, and lists it under **Add-Ons**, where you can
update or uninstall it later.

If the installer offers to download EEGLAB, say yes and skip step 2.

### 2. Get EEGLAB

Download [`eeglab_current.zip`](https://sccn.ucsd.edu/eeglab/currentversion/eeglab_current.zip)
(about 95 MB; older versions are on
[EEGLAB's download page](https://eeglab.org/others/How_to_download_EEGLAB.html))
and unzip it into your MATLAB folder: `Documents/MATLAB`, or wherever `userpath`
at the prompt says. nestapp finds it there by itself.

Unzipped somewhere else? Set the folder once in nestapp's
**Settings → Preferences**.

### 3. Add the EEGLAB extensions

Open EEGLAB: in MATLAB, go into the EEGLAB folder you just unzipped (for
example `cd(fullfile(userpath, 'eeglab2026.1.0'))`) and type `eeglab`. Then:

- **File → Manage EEGLAB extensions**, and install **TESA**, **PICARD** and
  **CleanLine**.
- **FastICA** is not in the extension manager. Download
  [`FastICA_2.5.zip`](https://research.ics.aalto.fi/ica/fastica/code/FastICA_2.5.zip)
  and unzip it into EEGLAB's `plugins` folder.

ICLabel, clean_rawdata and firfilt come with EEGLAB.

### 4. Launch nestapp

```matlab
nestapp
```

Then run **Help → Check My Install**. It should say there are no problems.
If anything is missing, it names it and how to get it.

### Optional: the AARATEP template

The `TMS-EEG / AARATEP` template needs helper functions from a separate project.
In nestapp, choose **Help → Install AARATEP Helpers...** (about a second), or
type `installAaratep` at the prompt. The helpers are kept outside the nestapp
folder, so updating nestapp does not remove them. The **Curve Fitting Toolbox**
is also needed for that template's `Remove Decay Artifact` step.

nestapp installs a **pinned release** (currently v2.1.1), not whatever is
newest, so that the same template produces the same result for everyone.
If a newer AARATEP exists you are told, but upgrading is a deliberate change.
See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for the pinned commit
and licensing. Already have your own clone? Point nestapp at it instead:

```matlab
setpref('nestapp', 'aaratepPath', '/path/to/AARATEPPipeline')
```

### Optional: TESA 1.2

The extension manager installs TESA 1.1.1, which runs every built-in template.
TESA 1.2 adds six more steps (Robust Detrend, Robust Demean, Modified Bandpass
Filter, Detect Bad Channels, Fit Artifact Model and Interactive Channel Reject,
all TESA versions). 1.2 has not been released yet, so it comes from TESA's
GitHub repository:

1. Download [TESA's current code](https://github.com/nigelrogasch/TESA/archive/refs/heads/master.zip).
2. Unzip it into EEGLAB's `plugins` folder, and rename the folder it creates
   (`TESA-master`) to `TESA1.2.0`. The development code does not state its
   version, so nestapp reads it from this folder name.
3. Delete the old `TESA1.1.1` folder beside it, so only one TESA is found.
4. Restart nestapp. The six steps now appear in the picker.

**Use at your own risk.** TESA 1.2 is unreleased development code: it may
contain bugs, it can change from one download to the next, and we have not
verified its results. For published or production analyses, stay on TESA 1.1.1.

### From source (for development, or to track `develop`)

```matlab
% clone or download the repo, then from the project root:
run_nestapp
```

`run_nestapp` adds `src/` to the path for that session and launches the app.
Steps 2 and 3 above still apply. See
[CONTRIBUTING.md](.github/CONTRIBUTING.md) for the development setup.

---

## Requirements

The installation steps above cover all of this. This is the reference list.

**Required.** nestapp will not run without these:

| | Version | Used for |
|---|---|---|
| MATLAB | R2023b or later | — |
| Signal Processing Toolbox | any recent | filtering and spectra (TESA needs it too) |
| Statistics and Machine Learning Toolbox | any recent | ICA quality metrics |
| EEGLAB | 2025.0 or later (tested with 2026.0 and 2026.1) | data loading and most steps |

**Required for most pipelines.** Every built-in template uses these, and steps
that need one are hidden until it is installed:

| | Used for | Where it comes from |
|---|---|---|
| TESA 1.1.1+ | TMS pulse detection, artifact removal, interpolation, filtering, SSP-SIR, TEP peaks | extension manager (1.2 adds six steps, [see above](#optional-tesa-12)) |
| FastICA | the FastICA decomposition engine | [Aalto University](https://research.ics.aalto.fi/ica/fastica/) |
| PICARD | the Picard decomposition engine | extension manager |
| CleanLine | line-noise removal | extension manager |
| ICLabel | automatic IC classification | comes with EEGLAB |
| clean_rawdata | ASR and bad-channel detection | comes with EEGLAB |
| firfilt | FIR filtering | comes with EEGLAB |

**Optional.** Each unlocks something specific, and nestapp works without them:

| | What it enables | Without it |
|---|---|---|
| Parallel Computing Toolbox | **Parallel Processing** checkbox: processes a batch across workers | The run proceeds serially and says why it skipped |
| Curve Fitting Toolbox | `Remove Decay Artifact` (AARATEP template) | That one step is unavailable |
| [AARATEP helpers](https://github.com/chriscline/AARATEPPipeline) | the `TMS-EEG / AARATEP` template | Installed from inside the app, [see above](#optional-the-aaratep-template) |
| `bva-io` | reading BrainVision `.vhdr` | Those files cannot be loaded |
| `loadcnt` | reading Neuroscan `.cnt` | " |
| `curry` | reading Curry `.cdt` | " |

Some EEGLAB plugins ask for the **Image Processing** or **Wavelet** toolboxes
for their own features. nestapp calls neither directly; **Check My Install**
reports them so you can tell whether a plugin is missing one.

**You do not have to get this list right up front.** nestapp only *offers* a
step whose dependencies are actually installed, and the pre-flight check blocks
a run that references a missing one and names what to install.

### If something looks wrong

**Help → Check My Install** reports what is missing and how to get it. The full
report (MATLAB version, every plugin found and its version) is
**Help → Copy Diagnostics to Clipboard**, or `nestappDoctor` at the prompt.

---

## First run

1. **File → Load Template** and pick a published pipeline, or build one in the **Cleaning** tab.
2. **File → Open Data...** to choose files (`.set`, `.vhdr`, `.cnt`, `.cdt`).
3. **Run Analysis**. A progress dialog names the current file and step, and the run can be cancelled between steps.
4. The **Reports** tab opens with a per-file summary. Copy the methods paragraph, or export the batch as CSV.
5. Load the cleaned files in the **Explore** tab to compare conditions, then export a figure or the per-subject measures.

| Template | Pipeline | Cite |
|---|---|---|
| TMS-EEG / TEP (TESA) | TESA two-round ICA for single-pulse TMS | Rogasch et al. 2017, *NeuroImage* 147:934-951 |
| TMS-EEG / AARATEP | SOUND + decay-fit removal + AR-blend interpolation + TMS-aware muscle classifier | Cline et al. 2021, *IEEE NER* |
| Resting-State EEG | PREP-style robust referencing and bad-channel detection, with ICLabel | Bigdely-Shamlo et al. 2015, *Front Neuroinform* 9:16 (also Delorme 2023) |
| Minimal ERP | High-pass, bad channels, ICA; minimum-handling | Delorme 2023, *Sci Rep* 13:2372 |

Each template logs its own citation and DOI at the start of every batch run, so the reference lands in the run log beside the data.

---

## Documentation

| | |
|---|---|
| [**User guide**](docs/user-guide.md) | Every tab, control and export route in detail |
| [Architecture](docs/architecture.md) | How the code is laid out and where to change things |
| [Contributing](.github/CONTRIBUTING.md) | Development setup, the test suite, adding a pipeline step |
| [Style guide](docs/STYLE.md) | Code and test conventions |
| [Changelog](CHANGELOG.md) | What changed, and when |
| [Third-party notices](THIRD_PARTY_NOTICES.md) | Attribution for every bundled and vendored dependency |

---

## Citing

If you use nestapp in published work, cite it with the metadata in [`CITATION.cff`](CITATION.cff). GitHub's **Cite this repository** button generates APA and BibTeX for you.

**Also cite the pipeline paper for whichever template you ran**, from the table above. nestapp writes these into the run log so the reference is recorded with the data rather than reconstructed later.

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](.github/CONTRIBUTING.md) to get started, and [`docs/architecture.md`](docs/architecture.md) for a map of the codebase.

## Authors

**Aref Pariz** wrote the original application in 2023, at the Royal Institute for Mental Health in Dr. Sara Tremblay's lab ([NESTLAB](https://www.nest-lab.ca/)) and in Dr. Jeremie Lefebvre's lab at the University of Ottawa.

**Wesley Dunne** contributed the pipeline engine, progress reporting and provenance, automated reports with ICA tracking, TEP visualisation and batch peak extraction, pipeline templates, quality-control gates, and the test suite, and led the 1.0.0 open-source release.

## License

GPL-3.0-or-later. See [`LICENSE`](LICENSE).
