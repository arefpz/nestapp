# Security Policy

## Supported versions

nestapp follows semantic versioning. Fixes are applied to the latest released
minor version. There is no long-term support for older lines during the 1.x
series.

| Version | Supported |
|---------|-----------|
| 1.0.x   | yes       |
| < 1.0   | no        |

## Scope

nestapp is a desktop research tool that runs locally in MATLAB and processes
EEG/TMS-EEG data files the user already trusts. It is not a network service and
does not handle credentials. The most relevant concerns are therefore:

- Crafted data files (e.g. `.set`, `.vhdr`) that could trigger unsafe behavior
  when loaded.
- The `Manual Command` pipeline step, which executes user-supplied MATLAB code
  by design — only run pipelines you trust.
- Third-party dependencies (EEGLAB, TESA, FastICA, AARATEP); vulnerabilities in
  those should be reported upstream, but let us know if a nestapp default
  exposes one.

## Reporting a vulnerability

Please do **not** open a public GitHub issue for security problems.

Email **dunne.wesley@gmail.com** with:

- a description of the issue and its impact,
- steps to reproduce (a minimal example if possible),
- the nestapp version (`nestappVersion()` or Help → About) and MATLAB version.

We will acknowledge receipt within a reasonable timeframe, work with you on a
fix, and credit you in the release notes unless you prefer otherwise.
