# nestapp architecture

A map of how nestapp is put together and **where to make common changes**.
For the per-function contract of any file, see the generated
[API reference](../site/index.md) (built from the source header comments by
`tools/gen_docs.m`).

## Big picture

nestapp has two halves that share one engine:

```
            ┌─────────────────────────────┐
   GUI ───► │  src/@nestapp/nestapp.m      │   five tabs:
            │  (App Designer class, plain  │   Cleaning · Visualizing ·
            │   text — edit directly)      │   Analysis · Reports · Settings
            └──────────────┬──────────────┘
                           │ builds a pipeline "spec" (ordered steps + params)
                           ▼
            ┌─────────────────────────────┐
   headless │  src/runPipelineCore.m       │   batch engine (serial / parallel)
   callers  │      └─► src/processOneFile.m│   runs one file through the spec
   (tests,  └──────────────┬──────────────┘
   scripts)                │ switch over step name → EEGLAB / TESA / helpers
                           ▼
                 EEGLAB / TESA / FastICA / vendored AARATEP
```

The same `spec` (a struct array of `{name, params}`) drives both the GUI run
button and headless callers like `batchTEPExtract` and the test suite.

## Module map (`src/`)

| Area | Files | Responsibility |
|---|---|---|
| **GUI** | `@nestapp/nestapp.m`, `@nestapp/createComponents.m`, `@nestapp/rescaleComponents.m` | The App Designer class. `nestapp.m` holds all state, callbacks, and tab logic. `createComponents`/`rescaleComponents` build and lay out the UI. **Edit `nestapp.m` directly — never the `.mlapp`.** |
| **Step registry** | `stepRegistry.m`, `makePipelineStep.m`, `checkStepDependencies.m` | The catalogue of pipeline steps: each step's name, default params, UI metadata, and dependency requirements. The single source of truth for "what steps exist." |
| **Execution** | `runPipelineCore.m`, `processOneFile.m`, `paramsToVarin.m`, `varinToStruct.m`, `stripVarinKeys.m`, `stripEmptyVarin.m`, `nestLog.m` | The batch engine and the per-file dispatch `switch`. Each step name maps to an EEGLAB/TESA call or a nestapp helper here. |
| **Templates** | `buildTemplates.m`, `templates/*.mat`, `templateCitation.m`, `specFromSaved.m` | Built-in pipelines. `buildTemplates.m` is the source; the `.mat` files are **generated artifacts** (see gotchas). |
| **Step helpers** | `aaratepMuscleClassifier.m`, `artist*.m`, `ensureAaratepOnPath.m`, `computeICAActivation.m`, `tepPeakFinder.m`, `batchTEPExtract.m`, `defaultTEPComponentDefs.m`, `tepPeakFinder.m` | Algorithm implementations behind specific steps and analyses. |
| **Quality control** | `qa/*.m` | Quality Gate scoring, batch verdicts, QC images, dashboard, attribute matrices. |
| **Reporting / IO** | `buildReportText.m`, `initPipelineReport.m`, `exportReport.m`, `summarizeReports.m`, `buildHistoryEntry.m`, `io/*.m` | Per-file reports, methods paragraphs, provenance, and output-path layout. |
| **Version** | `nestappVersion.m` | Single source of truth for the app version (SemVer). |

## "If you want to change X, edit Y"

| Goal | Where | Notes |
|---|---|---|
| Add a processing step | `stepRegistry.m` (register) + `processOneFile.m` (dispatch `case`) | Recipe in [CONTRIBUTING.md](../CONTRIBUTING.md#adding-a-pipeline-step). |
| Change/add a built-in template | `buildTemplates.m`, then run `buildTemplates()` and commit the regenerated `templates/*.mat` | The `.mat` is generated — never edit it directly. |
| Add a citation for a template | `templateCitation.m` | Logged per run by `runPipelineCore.m`. |
| Change a Quality Gate metric | `qa/qualityGate.m` (+ `qa/finalizeBatchVerdicts.m` for batch mode) | Step params live in `stepRegistry.m`. |
| Change TEP peak detection | `tepPeakFinder.m` (interactive) / `batchTEPExtract.m` (CSV) | Both feed the smoothed waveform to `tepPeakFinder`. |
| Change a tab's UI/behaviour | `@nestapp/nestapp.m` (callbacks) + `@nestapp/createComponents.m` (layout) | Plain-text class; diffable. |
| Change report contents | `buildReportText.m`, `initPipelineReport.m` | |
| Bump the version | `nestappVersion.m` + `CHANGELOG.md` (+ git tag) | A CI check keeps the three in sync. |

## Data flow of a run

1. The GUI (or a headless caller) assembles a `spec`: an ordered struct array
   of `{name, params}` (`makePipelineStep` builds each from the registry).
2. `runPipelineCore(spec, files, opts)` sets up output paths, logging,
   citations, and serial/parallel execution.
3. For each file, `processOneFile` walks the spec; its `switch` maps every
   step `name` to its implementation (mostly `pop_*` EEGLAB/TESA calls, plus
   nestapp helpers). EEGLAB state lives in globals reset per worker.
4. Quality Gates score the data at checkpoints; reports and QC images are
   written; provenance is appended to `EEG.history`.

## Gotchas (also in CONTRIBUTING)

- **Generated templates.** `templates/*.mat` come from `buildTemplates()`.
- **Allowlist `.gitignore`.** New files in new dirs are invisible to git until
  allowlisted.
- **`nestapp.m`, not the `.mlapp`.** The `.mlapp` would overwrite hand edits.
- **EEGLAB and `third_party/` are not committed** — external dependencies.
- **`processOneFile` uses EEGLAB globals** (`EEG`, `ALLEEG`, …); headless
  callers should expect shared state to be reset per worker.

## Tests

`tests/run_tests.m` is the harness: `run_tests` (fast: unit + regression),
`run_tests('all')` (adds integration, needs EEGLAB/TESA). Unit tests avoid
EEGLAB; integration tests `assumeFail` (skip) when it's absent. See
`tests/unit/test_newStepDispatch.m` for the conventions.
