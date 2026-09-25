# Changelog

All notable changes to nestapp are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
The version here must match `src/nestappVersion.m` and the release git tag.

## [Unreleased]

## [2.2.0] - 2026-09-24

### Added

- **Subjects from filenames** (Explore > Files, subjects and groups). Replaces
  "Guess from filenames", which ran a fixed per-file heuristic and could not be
  told which part of a name is the person. The new dialog shows one example
  filename cut into clickable parts, plus up to three folders above it, and
  under each part how many distinct values it takes across the cohort. You click
  the part that is the person. The file -> subject table and the change in n
  update on every click, so nothing applies unseen. The rule is positional and
  deterministic: count from the start or the end, optionally split letters from
  digits, and combine several parts. `S1`, `S01` and `s01` count as one subject,
  and a file lacking the chosen part stays its own subject rather than being
  merged. The preselection reads the whole cohort, skips dates and prefers the
  part with the most distinct values; on 40 generated naming schemes it
  recovers the subject in 38. Scripts reproduce the GUI's subjects with
  `exploreDataset(..., struct('subjectRule', rule))`.

### Fixed

- **Paired design was locked off by any missing session.** It required every
  subject in every group, although the analysis already handles the incomplete
  case: it keeps the complete subjects and names the rest. Paired is now offered
  with two or more groups and at least two subjects in every group, and the
  note says how many subjects it will use.
- **Methods paragraph omitted most data-changing steps.** 23 of the 29 steps
  that change the data produced no sentence (RANSAC, Robust Detrend, Remove
  Flagged ICA Components, SSP-SIR, EDM, Manual Command, ...). Every step is now
  either described, declared as never changing the data, or given a plain
  generic sentence, and a test fails for any registry step without a case.
  Clauses now state the thresholds and windows that decide what is removed
  (RANSAC criteria, ICLabel probability ranges, epoch-rejection limits,
  compselect criteria, ASR settings). Where upstream picks the value, the value
  it actually uses is stated. The paragraph ends by pointing to the run's
  `spec.mat`.
- **Repeated steps in the methods paragraph** follow one rule: the same settings
  read "... again" (a second ICA or the final re-reference used to vanish), and
  changed settings read "later, ...". The TMS window change is judged from both
  ends of the window. A notch filter is no longer described as a band-pass.
- **Aggregate counts could visibly fail to add up** ("72 of 76 epochs retained
  (5 rejected)"). Each mean was rounded to an integer on its own; they are now
  given to one decimal. The per-file counts were always consistent.
- **Explore results showed one group's measures without saying whose.** The
  group is now a column, one row per group x window.
- **Scalp maps turned the whole figure pale blue**, on screen and in exported
  figures, because topoplot applies EEGLAB's background colour to the figure.
  The figure colour is now restored.
- **Step Info text credited plugins' functions to EEGLAB** (clean_rawdata,
  ICLabel, firfilt, ...). Every Info string now names the toolbox that supplies
  the implementation.
- **AARATEP install hints said "Bundled with nestapp".** They are fetched, and
  the hint now names the menu item that installs them.

## [2.1.0] - 2026-09-03

### Note on the version number

**This release jumps from 1.0.0 to 2.1.0, and the jump is deliberate.** Two
version lineages had collided. The original application was tagged and released
as `v1.01` through `v1.05` and then `v2.0` between 2024 and 2025; the
open-source re-baseline in May 2026 restarted at `1.0.0` on the reasoning that
1.0.0 is the conventional first public release.

The consequence was that GitHub's **"Latest release" was `v2.0`** - a commit 382
behind, predating the pipeline engine, the Explore tab and everything since - so
a newcomer clicking the download on the repository front page got a 2025
snapshot, and any `1.x` release would have sorted below it permanently. For an
app whose audience is researchers who are not looking closely at tag numbers,
that is worse than publishing no releases at all.

So numbering continues from the line that is already public. `1.0.0` remains in
this file as what it was: the open-source re-baseline of the 2.x application.

### Added

- **`Help > Install AARATEP Helpers...`** - installs the AARATEP helper
  functions from inside the app: 0.6 MB, about a second, no git or shell
  required. It fetches a **pinned release** (v2.1.1) rather than whatever is
  newest, because AARATEP defines the behaviour of a published pipeline and two
  people running the same nestapp template must get the same helpers. A newer
  upstream release is reported but never installed silently. Point nestapp at
  an existing clone with `setpref('nestapp', 'aaratepPath', ...)`.
- **Installable packaging** - `tools/package_toolbox.m` builds
  `nestapp-<version>.mltbx`. Double-clicking it installs nestapp, sets the path
  permanently, and lists it under Add-Ons for update or uninstall, replacing
  "unzip somewhere, cd there, run a script, re-do the cd every session".
  EEGLAB is declared as required additional software, so MATLAB's installer
  offers to fetch it too.

### Changed
- **The test suite has been replaced rather than trimmed.** 916 cases /
  11,201 executable lines across 116 files became 210 / 1,644 across 15 (plus
  19 shared helpers); the test-to-source ratio went 0.72 to 0.106 and a full
  run 174.5 s to 35.3 s. The reason for replacing rather than pruning was that
  nothing consumed the old suite: `run_tests` counted skipped tests and then
  ignored them, so a run in which 100% of tests skipped exited green; 125 of
  916 cases were executed by no automated job at all; and 11 files in the
  documented "no EEGLAB, no GUI" fast suite needed one or the other, so its
  stated contract was false.

  Suites are now **folders**, encoding as a cross-product the only two things
  that gate a test, EEGLAB and a display: `tests/pure/`, `tests/eeglab/`,
  `tests/gui/`, `tests/eeglab_gui/`. A test's requirements are visible in its
  path, so it cannot be misfiled unnoticed. Three rules are enforced: **a skip
  is a failure** (there are now zero `assumeFail` sites, down from 38), an
  empty or missing suite is a failure, and the path is restored on exit.

  The conventions are executable rather than documented. `SuiteHygieneTest`
  holds nine of them, including that every test sits in the folder its
  dependencies require, that no test rolls its own path setup or temp
  directory, that source-scraping is opt-in with a named exception list, and
  that every shared helper has a caller.

  The 10 step characterization goldens were carried across unchanged (now in
  `tests/golden/`) and never re-recorded during the work; re-recording one to
  make a red test green would discard the only protection the step layer has.
  Every regression the old `regression/` cases and `docs/param-audit-findings.md`
  pinned was converted to a behavioural test, and each was verified by
  reintroducing the defect and requiring the suite to go red. See
  `docs/test-rewrite-ledger.md`, which also records two known gaps and the
  mutation-testing results behind them.

- **Documentation rewritten, and several stale claims corrected.** The README
  was doing two jobs - install guide and user manual - so the tab-by-tab detail
  moved to `docs/user-guide.md` and the README is now setup and orientation.
  Verified against the running code rather than assumed, which turned up: the
  documented EEGLAB (2025.0.0) and TESA (1.1.1) versions were *below* what the
  registry requires (2026.0.0, TESA 1.2); FastICA was documented as a separate
  manual download when it is now an EEGLAB plugin; the requirements table
  listed 3 plugins where the registry declares 13; the workflow told users to
  "switch to the Analysis tab", which Stage 7 removed; and three Reports
  buttons were named wrong. `THIRD_PARTY_NOTICES.md` also described the EEGLAB
  stack as "bundled with the nestapp distribution" when a clone contains none
  of it - reworded to say what is true.

### Fixed
- **Adding a group in Explore failed with `Undefined function 'pop_loadset'`.**
  `ensureEeglabReady` judged readiness by the `PLUGINLIST` global, which the
  plugin scan writes and nothing ever clears - so a `restoredefaultpath`, a
  test using `hideFromPath`, or a pathdef reset left it reporting "ready" with
  EEGLAB gone from the path. Two Explore paths also read recordings through
  `pop_loadset` without asking whether EEGLAB was up at all.
  Readiness is now decided by probing the function actually needed, and both
  paths ask. This is the same failure `ensureAaratepOnPath` hit below and
  `tesaVersion` was one path change away from, so the policy moved into a
  shared `pathMemo`: a cached answer about the path is valid only while its
  sentinel resolves to the same file. That also invalidates on a *swapped*
  install - a second TESA, a repointed EEGLAB - which none of the three
  hand-rolled caches noticed.

- `ensureAaratepOnPath` trusted a persistent flag, so anything that removed
  the tree from the path after the first call - a test using `hideFromPath`,
  a `restoredefaultpath` - left it believing the tree was still there, and
  every later AARATEP step in that session failed with a bare "Undefined
  function". It now verifies a sentinel resolves and re-adds if not.
- **An AARATEP run reported almost nothing.** The orchestrator returns only
  the cleaned EEG; its provenance, QC images and intermediates were files in
  its own folder that nestapp never looked at. So the report showed 0 ICA
  components (and `methodsNarrative`, which gates the ICA sentence on
  `nComponents > 0`, left ICA out of the methods paragraph entirely, despite
  two ICA passes having run), 0 bad channels (AARATEP interpolates them back,
  so the count never moves), and no figures at all - the AARATEP template has
  no Quality Gate, so nestapp rendered none either. The new `aaratepHarvest`
  reads the metadata struct saved beside the result and registers the QC
  images as report figures, so they reach the Reports tab and the PDF.
- **Load from Folder was slow.** Three causes. The folder picker opened on
  `lastDataFolder`, which points at raw input data and is typically a network
  share, so the shell had to enumerate a remote path before the dialog could
  paint; it now opens on the folder reports were last loaded from, falling
  back to the output root, and every candidate is checked with `isfolder`
  first so a stale path cannot stall it. Report discovery walked the entire
  tree - a batch root also holds `data/` and `qc/`, so it stat'd hundreds of
  large files to find the few that are reports; it now tries the folder and
  its `reports/` subfolder first and only recurses if neither hits, which
  still lets a parent of several runs work. And each loaded report was passed
  through `exportReport`, which is the text builder *plus* a `save()` of the
  whole struct, so loading N reports wrote N throwaway `.mat` files to
  tempdir; it now calls `buildReportText` directly (4x faster on 34 reports,
  before counting the network).
- **The application window walked up the screen.** `UIFigureSizeChanged`
  enforced its minimum size by assigning `UIFigure.Position(3:4)`, which
  anchors `[left bottom]` - so restoring a height from a raised bottom edge
  lifted the whole window rather than holding the top edge the user set. That
  write re-fired the same callback (and the `drawnow` inside let the re-entry
  run), so over a stream of resize events the window marched off the top of
  the monitor at a constant size. The minimum is now applied by writing the
  full `Position` with the top edge pinned, only when the size is actually
  short, and the handler refuses to re-enter itself.
- `batch/session_summary.csv` now records the retention counts every run -
  `chans_original`, `chans_final`, `chans_interpolated`, `trials_original`,
  `trials_final`, `ica_removed` - alongside the existing operational columns.
  Those numbers previously existed only inside the per-file reports, so a
  methods table or a cross-file quality screen meant pressing a button in the
  Reports tab. Counts a report does not carry are written as NaN, never 0, so
  a missing value can never be averaged in as a real measurement (new
  `reportCounts` helper).
- **Export PDF** now asks whether to export the selected report or every
  listed report, instead of silently doing one file at a time - which made a
  whole batch tedious whenever Auto-export PDF was off, and it is off by
  default. One report failing no longer abandons the rest; failures are
  collected and reported together.
- The session summary is now written for **every** run that produced at least
  one report. It was gated on two or more successful files, so a single-file
  batch got per-file reports and no overall report, in the Reports tab or on
  disk. A cancelled run also writes its batch folder now (spec, CSV, dashboard
  and a `session_summary.txt` headed `*** RUN CANCELLED - PARTIAL RESULTS ***`)
  instead of leaving per-file reports with nothing to read them as a set, and
  the app says where those partial results landed.
- `batch/dashboard.png` was never written for a clean batch:
  `runPipelineCore` passed raw report structs to `anyReportHasGates`, which
  expects Reports-tab entries and so always returned false. The PNG only
  appeared when a file had failed.
- The Plot Type radio group no longer overlaps the topoplot axes. Its base
  geometry ran to x=335 while `UIAxes2` starts at x=340, and since every
  component scales by the same factor that base overlap scaled with the
  window. The group is now 185 px wide (was 195), so the gap stays positive
  at every window size.
- The Cleaning tab's Add / Remove / Move Up / Move Down buttons now span the
  full width of the Selected Steps column instead of huddling on its left.
- Rendering the topoplot no longer disturbs other figures. It called `gcf`
  (which *creates* a figure when none is open) and finished with a bare
  `close`, shutting whatever figure happened to be current - possibly one of
  the user's own plots or an EEGLAB window.
- The topoplot time lookup used exact equality against the time vector, so any
  latency not landing exactly on a sample yielded an empty index and errored.
  It now takes the nearest sample.
- The topoplot averaging-window ("Win") field had no callback and silently did
  nothing until the time spinner or TOPOPLOT was touched.
- **Export TEP Figure** no longer additionally saves the whole application
  window as a `.fig` on every export. Choosing `.fig` now writes a real,
  reopenable plot figure.

### Added
- Step: **Find TMS Pulses (AARATEP)**, wrapping upstream's own
  `c_TMSEEG_findTMSPulses`. It detects from the artifact across channels
  (requiring more than a quarter to cross threshold together) rather than from
  a single reference electrode, so no coil-adjacent channel has to be chosen.
  It sits beside the TESA detector in the picker, and the AARATEP template now
  uses it - keeping detection and cleaning within one toolbox. Its event label
  defaults to `TMS`, not upstream's `Pulse`, so the orchestrator's
  `pulseEvents` matches.
- Setting: **Keep AARATEP intermediate datasets**, on by default. AARATEP
  saves the dataset before SOUND, before decay removal and before ICA
  rejection as well as the result - four full datasets per file, roughly
  325 MB at the template defaults, and those three saves are guarded by
  `if true` upstream rather than by its debug flag. Untick to delete them
  once the step finishes; the vendored code is not patched.
- AARATEP's final `.mat` is now removed when a later Save New Set writes the
  same dataset as a `.set`, so the result exists once, in nestapp's format.
  Its metadata is read into the report first, and it is never deleted when
  that read fails or when no Save New Set follows.
- **Open TEP in Figure** and **Open Topo in Figure** buttons on the Visualizing
  tab. Each copies the plot into a standard MATLAB figure with a classic axes,
  so the plot editor, Property Inspector, Export Setup and `File > Save As` all
  work - axes, limits, colours and labels can be edited by hand before export.
  Backed by the new `popOutAxes` helper.
- The topoplot now carries a **labelled uV colorbar** with symmetric colour
  limits, so the map can actually be read. The limits topoplot computed were
  previously discarded along with its scratch axes, leaving the destination
  autoscaled. The colormap is now a diverging blue-white-red map
  (`divergingColormap`) with white pinned at 0 uV, replacing the cyclic,
  non-perceptual `hsv`. Drawing is shared by the in-app axes and the pop-out
  via the new `drawScalpTopo` helper.
- `nestappDoctor`, environment diagnostics that validate MATLAB/EEGLAB/
  toolbox versions, every plugin the step registry requires, and flag
  shadowed functions. Surfaced in the GUI via **Help → Copy Diagnostics to
  Clipboard** and referenced from the bug-report template.
- `describePipeline`, a readable summary of the current pipeline (steps and
  the parameters that differ from defaults), for methods sections and bug
  reports. Surfaced via **File → Copy Pipeline Description**.
- Debug log (`debugLog` preference), which tees a run's full step-by-step trace
  to a file in the batch output folder (`nestDebugLog` + `nestLog`).
- Metadata-only debug bundles: on a step failure, `saveErrorBundle` writes
  the error, environment, pipeline, and EEG **metadata** (never recordings)
  to `<batch>/debug/`. **Help → Collect Support Bundle** produces the same
  on demand, and **Help → Check My Install** runs the fast test suite.

### Removed
- The Session Quality Dashboard's files x gates **verdict heatmap**. It was
  unreadable past roughly ten files and carried nothing the Failed / Marginal
  table does not already say, one row per flagged file with its actual
  reasons. The table takes the freed space: full width (so the Reasons column
  is readable rather than truncated) and up to the header, which buys visible
  rows.

### Changed
- The step picker's amber dot now marks a step that **waits for you** rather
  than its provenance. Provenance was true of every step once expanded, so a
  colour for it would have marked all 51 rows and said nothing; the tree
  already shows the provider structurally wherever there is a choice between
  variants. Blocking is the exception (6 of 51) and was previously only
  discoverable *after* building a pipeline and pressing Run Analysis, when the
  parallel-processing warning fires. The set is derived from the registry via
  the new `canStepBlock`, so an interactive step marks itself.
- The Steps tree carries a hover tip explaining the dot, shown only after the
  pointer has rested on the tree for three seconds - any movement hides it and
  restarts the clock, so it stays out of the way. `uitreenode` has no Tooltip
  property, so a per-node tip is not possible, and the native tooltip fires on
  a schedule that cannot be delayed; the tip is a small floating panel driven
  by a singleShot timer restarted from `WindowButtonMotionFcn`. Selecting a
  step also shows what the dot means for it, plus which toolbox supplies it,
  in the Info panel above the description.
- **Export Metrics Table** is now **Export Metrics CSV**, and its tooltip
  leads with "Write a CSV file". The previous name did not say what came out,
  and the other exports on that tab produce PDFs or MATLAB objects.
- The Reports tab's **Refresh** button is now **Clear List**: it empties the
  report list instead of re-reading the load folder. The list accumulates
  across every run in a session (each run appends a Session Summary plus one
  entry per file), so what it needed was a way to get back to empty, not a
  reload. It clears session and disk-loaded reports alike, asks first because
  there is no in-session undo, and deletes nothing on disk - Load from Folder
  brings saved reports back. The Quality Dashboard's own Refresh button, which
  re-renders that panel, is unaffected.
- Reports tab buttons renamed and re-tooltipped so it is clear what each
  produces: **Export CSV** is now **Export Metrics Table** (one row per file,
  spanning every listed report including ones loaded from disk, so it can
  cover several batch runs), and **Export PDF** is now **Export PDF...**.
  The **Copy Methods Text** tooltip now says that it copies the *full*
  parameterized paragraph for a single file - longer than the one-sentence
  note shown in the report body - or the cross-file aggregate for the
  Session Summary.
- CI now runs on free GitHub-hosted runners via `matlab-actions/setup-matlab`
  (no self-hosted runner or license needed for the public repo). The fast
  suite, lint, and docs build run in CI; the EEGLAB-dependent integration
  suite is run locally.
- AARATEP path setup (`ensureAaratepOnPath`) no longer lets the bundled
  FastICA shadow the user's own (normally EEGLAB's) FastICA. The bundled
  copy is used only as a fallback when no other FastICA is on the path, and
  a one-time warning (`nestapp:aaratepFastICAMismatch`) is printed when the
  user's FastICA version differs from the one AARATEP was tested with.

## 1.0.0 - 2026-05-29

First public open-source release.

### Added
- **Pipeline builder** for drag-and-drop construction of EEG/TMS-EEG cleaning
  pipelines from a registry of steps (`src/stepRegistry.m`), executed by a
  batch engine (`src/runPipelineCore.m`, `src/processOneFile.m`) with
  serial and parallel (PCT) modes.
- **Built-in pipeline templates** (`src/buildTemplates.m`, `src/templates/`):
  TMS-EEG / TEP (TESA), TESA + Quality Gates, ARTIST (Wu 2018), AARATEP
  (Cline 2021), Resting-State, and Minimal (Delorme 2023).
- **Quality control**: Quality Gate step with absolute and batch (median +
  MAD) thresholds, skip-on-fail, auto-generated QC images, and a Session
  Quality Dashboard on the Reports tab.
- **TEP analysis**: ROI waveform plotting, topographies, TESA-based peak
  detection with a polarity guard, and batch peak extraction to CSV
  (`src/batchTEPExtract.m`, `src/tepPeakFinder.m`).
- **Reports & provenance**: per-file reports, methods-paragraph export, and
  full pipeline provenance written to `EEG.history`.
- **Citations**: built-in templates log their primary reference per run
  (`src/templateCitation.m`); `THIRD_PARTY_NOTICES.md` documents all
  bundled and vendored dependencies.
- **Test suite**: unit, regression, and integration tests with a
  `tests/run_tests.m` harness and CI (`.github/workflows/tests.yml`).

### Notes
- Requires MATLAB R2023b+, EEGLAB, TESA, and FastICA. The AARATEP template
  additionally requires the Curve Fitting Toolbox and the AARATEP helpers
  cloned into `third_party/aaratep/` (see README).

[Unreleased]: https://github.com/goingloud/nestapp/compare/v2.2.0...HEAD
[2.2.0]: https://github.com/goingloud/nestapp/releases/tag/v2.2.0
[2.1.0]: https://github.com/goingloud/nestapp/releases/tag/v2.1.0
