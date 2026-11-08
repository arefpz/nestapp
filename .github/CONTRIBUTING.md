# Contributing to nestapp

Thanks for your interest in improving nestapp. This guide covers how to set up
a development environment, the conventions we follow, and how to get a change
merged. New to the codebase? Read [docs/architecture.md](../docs/architecture.md)
first. It maps what each module does and where to make common changes.

## Development setup

1. **MATLAB R2023b or later** (developed on R2026a). Required: Signal
   Processing and Statistics and Machine Learning toolboxes. Optional: Parallel
   Computing (for batch parallelism) and Curve Fitting (for the AARATEP
   template's decay step). See the
   [README](../README.md#requirements) for the full required/optional split.
2. **EEGLAB plus its plugins** (TESA, FastICA, ICLabel, PICARD, CleanLine,
   clean_rawdata, firfilt). Install these per the
   [README](../README.md#installation). For development, install TESA 1.2
   ([README, Optional: TESA 1.2](../README.md#optional-tesa-12)) so the six
   1.2-only steps can be exercised. **None of these are redistributed with
   this repo**; `eeglab2026.0.0/` in a development checkout is untracked.
3. **AARATEP helpers** (only if working on that template). Run
   `installAaratep` at the MATLAB prompt, which fetches the pinned release into
   `prefdir/nestapp/aaratep`. A `git clone` of the same tag into
   `third_party/aaratep/` works too and is found first; the installer stamps
   the version it wrote, so a hand-cloned tree is reported as "version unknown"
   rather than assumed to be the pin.
4. **Clone and launch**:
   ```bash
   git clone https://github.com/arefpz/nestapp.git
   ```
   ```matlab
   run_nestapp            % launches the app
   ```
5. **Run the tests** (from the repo root):
   ```matlab
   cd tests
   run_tests                 % pure: no EEGLAB, no display, ~4 s
   run_tests('eeglab')       % needs EEGLAB (includes the step goldens)
   run_tests('gui')          % needs a display; opens windows
   run_tests('eeglab_gui')   % needs both; launches the app
   run_tests('all')          % everything, ~40 s
   ```
   Folders are the suites, and they encode the only two things that gate a
   test: **does it need EEGLAB, does it need a display.** Both are binary, so
   the cross-product is four folders and nothing else is needed to select
   tests. A test's requirements are visible in its path, so it cannot be
   misfiled unnoticed. `tests/pure/SuiteHygieneTest.m` checks that convention,
   and eight others, executably.

   `run_tests` enforces three rules:

   - **A skip is a failure.** The folder already declares what a test needs, so
     no test has any business deciding it cannot run. There are zero
     `assumeFail` sites; if you find yourself wanting one, the test is in the
     wrong folder.
   - **An empty or missing suite is a failure**, so a mistyped or renamed
     folder cannot pass by running nothing.
   - **The path is restored** on the way out, so runs don't leak into each
     other.

## Branching & commits

- Branch from **`develop`**. `main` is the released line; `develop` is the
  integration branch. Open PRs against `develop`.
- Use **Conventional Commits**: `feat:`, `fix:`, `docs:`, `test:`, `chore:`,
  `refactor:`, with an optional scope, e.g. `feat(templates): add ...`. The
  `CHANGELOG.md` is curated from these.
- Keep commits focused; one logical change per commit.

## Pull requests

- Ensure `run_tests('all')` passes locally. CI runs the `pure` suite and the
  linter on every push to `main` and `develop` and on every pull request, but
  it cannot run the `eeglab`, `gui` or `eeglab_gui` suites: EEGLAB is not
  tracked in this repository and a hosted runner has no display. So CI proves
  the pure suite and static analysis; the other three are on you.
- When you add a test, put it in the folder its dependencies require and let
  the suite run with no skips. When you add a *helper*, make sure something
  calls it. An unused helper is a second thing to keep in step, and the
  hygiene test will fail on it.
- Add or update tests for the behaviour you change. When fixing a bug, add a
  regression test that fails before your fix (the "if you liked it, put a test
  on it" rule).
- Update `CHANGELOG.md` under `## [Unreleased]`.
- A maintainer (see `.github/CODEOWNERS`) reviews and merges.

## Code style

See [STYLE.md](../docs/STYLE.md) for the full conventions. In short: camelCase
functions and variables, 4-space indent, a header comment block on every
function (see the docstring contract below), and named intermediate variables
over dense one-liners. CI runs `tools/run_lint.m` (a `checkcode` wrapper); fix
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
[docs/architecture.md](../docs/architecture.md)):

1. **Register it** in `src/registry/stepRegistry.m`, appending a block with `.name`,
   `.defaults`, `.info`, `.params` (per-field metadata), and optional
   `.requires` (plugin/toolbox dependencies surfaced by the pre-flight check).
2. **Dispatch it** in `src/processOneFile.m`, adding a `case 'Your Step Name'`
   to the switch that maps the step to its implementation.
3. **(Optional) Use it in a template.** Edit `src/util/buildTemplates.m`, then
   regenerate: run `buildTemplates()` and commit the updated `src/templates/*.mat`.
4. Add a test. A step's dispatch is covered by `tests/eeglab/DispatchContractTest.m`
   automatically, because it sweeps the whole registry, so what a new step
   needs is a *behavioural* test of what it does. If its output is worth pinning
   exactly, add a case to `tests/helpers/characterizationCases.m` and record a
   golden with `recordGoldens('Your Step Name')`.

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
  re-running it, and committing the regenerated `.mat`, never editing the `.mat` itself.
- **Edit `nestapp.m` directly, never the `.mlapp`.** The GUI is a plain-text
  class (`src/@nestapp/nestapp.m`). Opening and saving the App Designer
  `.mlapp` would overwrite hand-edited methods.
- **EEGLAB and `third_party/` are not committed.** They are external
  dependencies, like in many MATLAB projects. Don't add them to git.

## Using AI coding assistants

LLM-assisted contributions are welcome, but the quality of the tool matters.
This is a research codebase with non-obvious cross-file invariants (the step
registry ↔ dispatch ↔ template contract, EEGLAB/TESA quirks, the allowlist
`.gitignore`), and a weak assistant will confidently get these wrong.

- **Use coding-grade agents on advanced models at high effort**, for example Claude
  **Opus 4.8** or OpenAI **Codex** (or the latest equivalents) in a high-effort /
  extended-reasoning mode, running an agent designed for software engineering.
- **Do not** use general-purpose chat assistants or small, underpowered models
  with limited context windows for non-trivial changes. They miss the
  cross-file contracts above and tend to produce plausible-but-wrong edits.
- You own what you submit: read the diff, make sure it follows the conventions
  here, and ensure `run_tests` passes. "The AI wrote it" is not a review.

## Reporting bugs & requesting features

Use the GitHub issue templates. For security issues, do **not** open a public
issue; see [SECURITY.md](SECURITY.md).

## License

By contributing, you agree your contributions are licensed under the project's
GPL-3.0 license.
