# nestapp user guide

Every tab, control and export route. For installation and a first run, see the
[README](../README.md); for the code layout, [architecture.md](architecture.md).

- [Cleaning tab](#cleaning-tab)
  - [Parallel processing](#parallel-processing)
- [Reports tab](#reports-tab)
- [Explore tab](#explore-tab)
- [Preferences](#preferences)
- [When something does not work](#when-something-does-not-work)

---

## Cleaning tab

### Building a pipeline

The left panel is a tree of every processing step, grouped by the stage it
belongs to. Click a step to read what it does and what its parameters mean in
the **Info** panel.

- **Add** (or double-click) appends the selected step to **Selected Steps**.
- **Remove**, **Move Up** and **Move Down** edit the order.
- Clicking a step in **Selected Steps** shows its parameters in the table
  below, where they can be edited. **Default Value** resets that step.
- **Re/Start Steps** clears the pipeline, after a confirmation prompt.

**Only installed steps are offered.** The picker asks what each step requires
and withholds any whose plugin is missing or too old. The alternative would be
building a pipeline around a step and finding out at the pre-flight. A saved
pipeline that references an unavailable step still *loads*; the pre-flight then
blocks the run and names both the step and what it needs.

### Templates

**File → Load Template** starts you from a published pipeline rather than a
blank list. See the table in the [README](../README.md#first-run) for what each
one is and the paper to cite. The template's citation and DOI are written to the
run log at the start of every batch, so the reference stays with the data.

### Saving and loading pipelines

- **File → Save Pipeline** writes the step list and every parameter to a `.mat`.
- **File → Load Pipeline** restores one.
- **File → Recent Pipelines** lists the last few.
- **File → Copy Pipeline Description** puts a prose description of the current
  pipeline on the clipboard. This is useful for a methods section before you
  have run anything.

### Running

1. **File → Open Data...** selects the input files (`.set`, `.vhdr`, `.cnt`,
   `.cdt`). Reading the non-EEGLAB formats needs the matching plugin:
   `bva-io` for BrainVision, `loadcnt` for Neuroscan, `curry` for Curry.
2. **Run Analysis** starts the batch. A progress dialog names the current file
   and step, and the run can be cancelled between steps.

A file that fails does not stop the batch: the error is captured with a
diagnostic bundle, the remaining files continue, and the failure is listed in
the report with the step that raised it.

### Parallel processing

**Parallel Processing** on the Cleaning tab spreads a batch across workers. It
needs the Parallel Computing Toolbox, and worker count comes from
**Settings → Preferences** (default 4, never more than the number of files).

It is **requested, not guaranteed.** The run falls back to serial and says so in
the status bar when:

| | Why |
|---|---|
| only one file is selected | there is nothing to parallelise |
| the Parallel Computing Toolbox is not licensed | no workers available |
| the pipeline contains an interactive step | a step that opens a dialog cannot run on a worker |

That last one is the common surprise: adding one interactive step to an
otherwise batch pipeline turns off parallelism for the whole run. The status bar
names the offending steps.

### Provenance

Every processed file gets its full pipeline written into `EEG.history` inside
the saved `.set`: the steps, their parameters, and a timestamp. After a run the final
dataset is also in the MATLAB base workspace, so:

```matlab
EEG.history      % the whole record as text
eegh             % browse it in EEGLAB's history viewer
```

This is deliberate: a `.set` from a nestapp run can always be traced back to how
it was made, even years later and without this repository.

---

## Reports tab

One report per processed file, added as the run proceeds. Each shows:

- channel counts: how many were present, rejected and interpolated, and which;
- trial retention: how many epochs survived, and what removed them;
- ICA statistics: components per round, what each was classified as, and how
  much variance was removed;
- quality-gate verdicts, where thresholds were configured (**Pass**,
  **Marginal** or **Fail**, with the metric and value that decided it).

Three things to do with a report:

- **Copy Methods Text** gives a prose description of what was actually done to
  *this* file, with the numbers filled in, ready to paste into a manuscript.
- **Export Metrics CSV** writes one row per file, for a multi-subject batch
  summary.
- **Export PDF...** writes the report as a document, for one file or all of them.

Reports are also written to disk beside the processed data, so a batch can be
reviewed later without re-running it.

---

## Explore tab

Explore works on processed `.set` files, either the output of a nestapp run or
any EEGLAB-compatible epoched dataset. It supports the comparison most TMS-EEG
studies are designed around: pre versus post, or one cohort against another.

### Defining the comparison

- **Add group...** assigns recordings to a named condition. A group is defined
  over **subjects**, so repeat recordings of one person are averaged before the
  group estimate is formed, and the reported *n* counts participants rather
  than files. Sessions are weighted by their trial counts, because two sessions of
  one person are not equal evidence if one held 10 trials and the other 190.
- **Files, subjects, groups...** shows and corrects which file belongs to whom.
  This is where *n* comes from, so it is worth checking.
- **Design** is set explicitly as paired or unpaired. Paired is offered only
  when every group holds the same subjects, and it drops any subject without a
  complete set, naming them rather than excluding them silently.
- **Region of interest** is chosen from a scalp diagram or a named preset,
  defaulting to `AF3 F1 F3 FC1 FC3`. A file on a non-conforming montage is
  excluded and named.
- **Windows** define the intervals measured. Each carries a peak polarity
  (`auto`, `pos` or `neg`) deciding which way its peak is read. The same table
  switches to **results**, reporting mean amplitude, peak latency and peak
  amplitude for the group selected in the groups list.

Measures are taken **per subject and then averaged**, not read off the group
average. For a peak the two differ systematically, because a peak read off an
average is flattened by latency jitter between subjects.

### Plots

The plot picker offers the ROI waveform, global and local mean field power, the
difference wave, scalp topographies, a TEP-topo grid, and a per-window bar
chart. A plot that cannot be drawn with the current settings is **listed with
the reason** rather than hidden. Plot-specific settings are behind
**Options...**.

Confidence intervals are computed across subjects, with the paired design using
the Cousineau–Morey correction. The interval level actually drawn is stamped
into the figure caption, so a band can never be labelled with a level it was not
computed at.

### Exporting

- **Figure...** composes the plot at a chosen column width and resolution,
  stamping the groups, sample sizes, design and ROI into a caption line so the
  image stays interpretable once separated from the session.
- **Measures → CSV...** writes one row per subject, group and window.
- **Results → MATLAB...** writes the complete result, curves included.
  **File → Load Analysis** reads it back, restoring groups, subject
  assignments, ROI, windows, design and plot selection.

---

## Preferences

**Settings → Preferences**, in the order the dialog presents them.

**Quality Screening**
- auto-generate QC images at each Quality Gate;
- auto-detect the TMS pulse window from EEG events, or set it in ms;
- skip the remaining pipeline steps when a Quality Gate fails;
- auto-save a PDF report per file (text plus checkpoint images);
- keep AARATEP intermediate datasets, which costs roughly **3x the result on
  disk**;
- attribute mode.

**EEGLAB**
- the installation path;
- suppress EEGLAB's processing dialogs (you are still warned about overwrites
  before a run);
- hide the EEGLAB window during processing.

**Default Locations**
- data, pipeline and report folders;
- output root, which can be left blank to write next to the inputs;
- layout, grouping results either by artifact type or per input file.

**Parallel Processing**
- max workers: a cap on simultaneous files when Parallel is on. Greyed out with
  *"Not available - Parallel Computing Toolbox not licensed"* if you do not
  have the toolbox. See [Parallel processing](#parallel-processing) for when it
  is used and when it is skipped.

**Behaviour**
- switch to the Reports tab after each run;
- confirm before clearing a pipeline;
- overwrite existing report files rather than timestamping them.

---

## When something does not work

**Help → Check My Install** first. It says what is missing, how many steps
need it, and where to get it, and lists optional extras you have not installed.
Most problems are a missing or out-of-date plugin, and it names which.

(The full report, with your MATLAB version and every plugin's version, is
`nestappDoctor` at the MATLAB prompt.)

If you need to report a problem, **Help → Collect Support Bundle...** gathers
the versions, the pipeline and the logs into one file to attach. Use
**Help → Copy Diagnostics to Clipboard** instead if you only want to paste them
into an issue.

| Symptom | Usually |
|---|---|
| An AARATEP step is unavailable | The helper functions are not installed — **Help → Install AARATEP Helpers...** |
| A step is greyed out in the picker | Its plugin is missing or older than the step needs — **Check My Install** names it |
| `Undefined function 'pop_...'` | nestapp could not find EEGLAB — set its folder in **Settings → Preferences** |
| A saved pipeline will not run | The pre-flight found an unavailable step; the message names it |
| A file failed mid-batch | The report names the step that raised it, and a diagnostic bundle is saved beside the output |
