# Security Policy

## Supported versions

nestapp follows semantic versioning. Fixes are applied to the latest released
minor version only; there is no long-term support for older lines.

| Version | Supported |
|---------|-----------|
| The latest release (see [CHANGELOG](../CHANGELOG.md)) | yes |
| Anything older | no |

`1.0.0` is not a lower version than `2.0` despite the number: it was a separate
release line, superseded by the 2.x releases. Always take the highest version.

## Scope

nestapp is a desktop research tool that runs locally in MATLAB and processes
EEG/TMS-EEG data files the user already trusts. It is not a network service and
does not handle credentials. The most relevant concerns are therefore:

- Crafted data files (e.g. `.set`, `.vhdr`) that could trigger unsafe behavior
  when loaded.
- The `Manual Command` pipeline step, which executes user-supplied MATLAB code
  by design, so only run pipelines you trust.
- Third-party dependencies: EEGLAB and its plugins (TESA, FastICA, ICLabel,
  PICARD, CleanLine, clean_rawdata, firfilt) and the AARATEP helpers.
  Vulnerabilities in those belong upstream, but tell us if a nestapp default
  exposes one. Note that nestapp **redistributes none of them**: EEGLAB and its
  plugins are installed by the user through EEGLAB's extension manager, and
  `Help > Install AARATEP Helpers...` downloads a pinned upstream release over
  HTTPS from GitHub.
- The toolbox package (`.mltbx`) attached to a release, which MATLAB installs
  into the user's add-ons folder. It contains only files tracked in this
  repository, because the packaging script takes its file list from `git ls-files`
  rather than globbing the working tree, so an untracked local dependency
  cannot be shipped inside it.

## Reporting a vulnerability

Please do **not** open a public GitHub issue for security problems.

Email **dunne.wesley@gmail.com** with:

- a description of the issue and its impact,
- steps to reproduce (a minimal example if possible),
- environment details: run `nestappDoctor` (or Help → Copy Diagnostics)
  and include the output.

We will acknowledge receipt within a reasonable timeframe, work with you on a
fix, and credit you in the release notes unless you prefer otherwise.
