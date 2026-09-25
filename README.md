# nestapp

![version](https://img.shields.io/badge/version-2.2.0-blue)
![MATLAB](https://img.shields.io/badge/MATLAB-R2023b%2B-orange)
![license](https://img.shields.io/badge/license-GPL--3.0-green)

**nestapp** is a MATLAB app for cleaning and analysing TMS-EEG recordings. It wraps [EEGLAB](https://eeglab.org/), [TESA](https://nigelrogasch.gitbook.io/tesa-user-manual/) and the surrounding plugin stack in a point-and-click pipeline builder, so a researcher can preprocess a cohort, check what happened to every file, and get publication-ready TEP figures and peak measurements without writing code.

Three tabs, in the order you use them:

- **Cleaning** builds a preprocessing pipeline, or loads a published one as a template, and runs it over a batch of files.
- **Reports** shows what the run did to each file: channels dropped, trials kept, ICA components removed, and a methods paragraph you can paste into a manuscript.
- **Explore** supports the comparison most TMS-EEG studies are designed around, pre versus post or one cohort against another, with TEP waveforms, topographies, and per-window measures exported to CSV.

Every processed file carries its full pipeline inside `EEG.history` in the saved `.set`: the steps, their parameters, and a timestamp. A result can therefore always be traced back to how it was made.

---

## Requirements

**Required.** nestapp will not run without these:

| | Version | How you get it |
|---|---|---|
| MATLAB | R2023b or later | — |
| Signal Processing Toolbox | any recent | MATLAB installer |
| Statistics and Machine Learning Toolbox | any recent | MATLAB installer |
| EEGLAB | 2026.0.0 (2025.x works) | [eeglab.org](https://eeglab.org/) |

**Required for most pipelines.** Every built-in template uses these, and steps
that need one are hidden until it is installed. All come from EEGLAB's own
extension manager:

| | Used for |
|---|---|
| TESA 1.2+ | TMS pulse detection, artifact removal, interpolation, filtering, SSP-SIR, TEP peaks |
| FastICA | the FastICA decomposition engine |
| ICLabel | automatic IC classification |
| PICARD | the Picard decomposition engine |
| CleanLine | line-noise removal |
| clean_rawdata | ASR and bad-channel detection |
| firfilt | FIR filtering |

**Optional.** Each unlocks something specific, and nestapp works without them:

| | What it enables | Without it |
|---|---|---|
| Parallel Computing Toolbox | **Parallel Processing** checkbox — processes a batch across workers | The run proceeds serially and says why it skipped |
| Curve Fitting Toolbox | `Remove Decay Artifact` (AARATEP template) | That one step is unavailable |
| [AARATEP helpers](https://github.com/chriscline/AARATEPPipeline) | the `TMS-EEG / AARATEP` template | Installed from inside the app — see [below](#installing-the-aaratep-helpers) |
| `bva-io` | reading BrainVision `.vhdr` | Those files cannot be loaded |
| `loadcnt` | reading Neuroscan `.cnt` | " |
| `curry` | reading Curry `.cdt` | " |

Some EEGLAB plugins ask for the **Image Processing** or **Wavelet** toolboxes
for their own features. nestapp calls neither directly; `nestappDoctor` reports
them so you can tell whether a plugin is missing one.

**You do not have to get this list right up front.** nestapp only *offers* a
step whose dependencies are actually installed, and the pre-flight check blocks
a run that references a missing one and names what to install.

---

## Installation

### The easy way: install the toolbox

1. **Download the `.mltbx` file** from the
   [latest release](https://github.com/goingloud/nestapp/releases/latest).
2. **Double-click it.** MATLAB installs nestapp, puts it on the path
   permanently, and lists it under **Add-Ons** (where you can later update or
   uninstall it). If you do not have EEGLAB, the installer offers to fetch it.
3. **Install the EEGLAB plugins.** Run `eeglab` at the MATLAB prompt, open
   **File → Manage EEGLAB extensions**, and install TESA, FastICA, ICLabel,
   PICARD, CleanLine, clean_rawdata and firfilt. They go into EEGLAB's own
   `plugins/` folder; nothing needs adding to your path by hand.
4. **Launch nestapp**:
   ```matlab
   nestapp
   ```
5. **Point it at EEGLAB.** Open **Settings → Preferences** and set the EEGLAB
   installation folder. nestapp puts EEGLAB and its plugins on the path on
   every launch from then on.

There is no `cd` to remember and no path to re-add each session.

### From source (for development, or to track `develop`)

```matlab
% clone or download the repo, then from the project root:
run_nestapp
```

`run_nestapp` adds `src/` to the path for that session and launches the app.
Steps 3 and 5 above still apply. See
[CONTRIBUTING.md](.github/CONTRIBUTING.md) for the development setup.

### If something looks wrong

**Help → Check My Install** reports your MATLAB version, which plugins were
found and their versions, and what any unavailable step is waiting for. Most
problems are a missing or out-of-date plugin, and it names which. The same
check is `nestappDoctor` at the prompt.

### Installing the AARATEP helpers

Only needed for the `TMS-EEG / AARATEP` template. Its helper functions are a
separate project. nestapp fetches them on request rather than committing a copy
into this repository, to keep the two codebases from entangling — the licence
(MIT) would permit redistribution:

**Help → Install AARATEP Helpers...** takes about a second and needs no other tools.

Or at the MATLAB prompt:

```matlab
installAaratep
```

nestapp installs a **pinned release** (currently v2.1.1), not whatever is
newest, so that the same template produces the same result for everyone.
If a newer AARATEP exists you are told, but upgrading is a deliberate change.
See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for the pinned commit
and licensing.

Already have your own clone? Point nestapp at it instead:

```matlab
setpref('nestapp', 'aaratepPath', '/path/to/AARATEPPipeline')
```

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
