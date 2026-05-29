# Contributing to nestapp

Thanks for your interest in improving nestapp. This guide covers how to set up
a development environment, the conventions we follow, and how to get a change
merged. New to the codebase? Read [docs/architecture.md](docs/architecture.md)
first — it maps what each module does and where to make common changes.

## Development setup

1. **MATLAB R2023b or later** (developed on R2025b+). Required toolboxes:
   Signal Processing, Statistics and Machine Learning. The `TMS-EEG / AARATEP`
   template additionally needs the Curve Fitting Toolbox.
2. **EEGLAB**, **TESA**, **FastICA** — install per the [README](README.md#requirements).
   These are *not* bundled with the repo.
3. **AARATEP helpers** (only if working on that template) — clone into
   `third_party/aaratep/`:
   ```bash
   cd third_party && git clone --depth 1 https://github.com/chriscline/AARATEPPipeline.git aaratep
   ```
4. **Clone and launch**:
   ```bash
   git clone https://github.com/goingloud/nestapp.git
   ```
   ```matlab
   run_nestapp            % launches the app
   ```
5. **Run the tests** (from the repo root):
   ```matlab
   cd tests
   run_tests              % fast: unit + regression, no EEGLAB needed
   run_tests('all')       % adds integration tests (EEGLAB + TESA required)
   ```
   `run_tests` prints a summary and errors on any failure. Integration tests
   `assumeFail` (skip) gracefully when EEGLAB/TESA are absent.

## Branching & commits

- Branch from **`develop`**. `main` is the released line; `develop` is the
  integration branch. Open PRs against `develop`.
- Use **Conventional Commits**: `feat:`, `fix:`, `docs:`, `test:`, `chore:`,
  `refactor:`, with an optional scope, e.g. `feat(templates): add ...`. The
  `CHANGELOG.md` is curated from these.
- Keep commits focused; one logical change per commit.

## Pull requests

- Ensure `run_tests` passes. CI runs the fast suite (unit + regression) on a
  GitHub-hosted runner — no setup required. The integration suite needs
  EEGLAB/TESA, which CI does not provide, so **run `run_tests('all')`
  locally** when your change touches pipeline execution or anything
  EEGLAB-dependent.
- Add or update tests for the behaviour you change. When fixing a bug, add a
  regression test that fails before your fix (the "if you liked it, put a test
  on it" rule).
- Update `CHANGELOG.md` under `## [Unreleased]`.
- A maintainer (see `.github/CODEOWNERS`) reviews and merges.

## Code style

See [STYLE.md](STYLE.md) for the full conventions. In short: camelCase
functions and variables, 4-space indent, a header comment block on every
function (see the docstring contract below), and named intermediate variables
over dense one-liners. CI runs `tools/run_lint.m` (a `checkcode` wrapper) — fix
reported errors before merging.

### Function docstring contract

Every function file in `src/` starts with a header block in this shape (the
generated API docs are built from it, so keep it accurate):

```matlab
function out = myFunction(in)
% MYFUNCTION  One-line summary (the H1 line).
%
%   Syntax
%     out = myFunction(in)
%
%   Inputs
%     in   <type> - what it is.
%
%   Outputs
%     out  <type> - what it is.
%
%   Side effects   (only if any: globals, EEG mutation, files written)
%   Requires       (only if any: EEGLAB/TESA/toolbox functions)
%   See also: relatedFunction
```

## Adding a pipeline step

Steps are data-driven. To add one (worked detail in
[docs/architecture.md](docs/architecture.md)):

1. **Register it** in `src/stepRegistry.m` — append a block with `.name`,
   `.defaults`, `.info`, `.params` (per-field metadata), and optional
   `.requires` (plugin/toolbox dependencies surfaced by the pre-flight check).
2. **Dispatch it** in `src/processOneFile.m` — add a `case 'Your Step Name'`
   to the switch that maps the step to its implementation.
3. **(Optional) Use it in a template** — edit `src/buildTemplates.m`, then
   regenerate: run `buildTemplates()` and commit the updated `src/templates/*.mat`.
4. Add a test under `tests/unit/` (see `test_newStepDispatch.m` for the pattern).

## Diagnosing environment problems

If steps fail to run, the app won't launch, or a dependency seems missing, run:

```matlab
nestappDoctor          % prints a validated environment report
```

It checks MATLAB/EEGLAB/toolbox versions, every plugin the step registry
requires, and flags shadowed functions (a common cause of wrong-version
calls). The same report is available in the GUI via **Help → Copy Diagnostics
to Clipboard** and should be attached to bug reports.

For a deeper trace, enable the debug log before a run:

```matlab
setpref('nestapp', 'debugLog', true)   % writes a full run trace to the batch folder
```

If a run fails mid-pipeline, nestapp also writes a **metadata-only** error
bundle (no recordings) under `<batch>/debug/`. **Help → Collect Support
Bundle** produces the same environment + pipeline summary on demand. These
artifacts are what to attach when reporting a hard-to-reproduce problem.

## Gotchas (read before your first change)

- **Allowlist `.gitignore`.** The repo ignores everything by default and
  allowlists specific paths. **A new file in a new directory is invisible to
  git** until you add a matching `!path/` rule to `.gitignore`. Check
  `git status` shows your new file before committing.
- **Templates are generated binaries.** `src/templates/*.mat` are produced by
  `buildTemplates()`. Editing a template means editing `buildTemplates.m`,
  re-running it, and committing the regenerated `.mat` — not editing the `.mat`.
- **Edit `nestapp.m` directly, never the `.mlapp`.** The GUI is a plain-text
  class (`src/@nestapp/nestapp.m`). Opening and saving the App Designer
  `.mlapp` would overwrite hand-edited methods.
- **EEGLAB and `third_party/` are not committed.** They are external
  dependencies, like in many MATLAB projects. Don't add them to git.

## Reporting bugs & requesting features

Use the GitHub issue templates. For security issues, do **not** open a public
issue — see [SECURITY.md](SECURITY.md).

## License

By contributing, you agree your contributions are licensed under the project's
GPL-3.0 license.
