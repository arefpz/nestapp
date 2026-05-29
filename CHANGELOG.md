# Changelog

All notable changes to nestapp are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
The version here must match `src/nestappVersion.m` and the release git tag.

## [Unreleased]

## [1.0.0] - 2026-05-29

First public open-source release.

### Added
- **Pipeline builder** — drag-and-drop construction of EEG/TMS-EEG cleaning
  pipelines from a registry of steps (`src/stepRegistry.m`), executed by a
  batch engine (`src/runPipelineCore.m`, `src/processOneFile.m`) with
  serial and parallel (PCT) modes.
- **Built-in pipeline templates** (`src/buildTemplates.m`, `src/templates/`):
  TMS-EEG / TEP (TESA), TESA + Quality Gates, ARTIST (Wu 2018), AARATEP
  (Cline 2021), Resting-State, and Minimal (Delorme 2023).
- **Quality control** — Quality Gate step with absolute and batch (median +
  MAD) thresholds, skip-on-fail, auto-generated QC images, and a Session
  Quality Dashboard on the Reports tab.
- **TEP analysis** — ROI waveform plotting, topographies, TESA-based peak
  detection with a polarity guard, and batch peak extraction to CSV
  (`src/batchTEPExtract.m`, `src/tepPeakFinder.m`).
- **Reports & provenance** — per-file reports, methods-paragraph export, and
  full pipeline provenance written to `EEG.history`.
- **Citations** — built-in templates log their primary reference per run
  (`src/templateCitation.m`); `THIRD_PARTY_NOTICES.md` documents all
  bundled and vendored dependencies.
- **Test suite** — unit, regression, and integration tests with a
  `tests/run_tests.m` harness and CI (`.github/workflows/tests.yml`).

### Notes
- Requires MATLAB R2023b+, EEGLAB, TESA, and FastICA. The AARATEP template
  additionally requires the Curve Fitting Toolbox and the AARATEP helpers
  cloned into `third_party/aaratep/` (see README).

[Unreleased]: https://github.com/goingloud/nestapp/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/goingloud/nestapp/releases/tag/v1.0.0
