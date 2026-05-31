# CI/CD setup — recreating the GitHub Actions

This is a complete, reproducible runbook for standing up nestapp's continuous
integration and docs deployment in a fresh GitHub repository (e.g. when these
changes are merged into a repo that has never run the workflows). It documents
every workflow under `.github/workflows/`, the repository settings each one
needs, and the one-time commands that are easy to miss.

Everything here works on a **public** repository with **no paid runners and no
secrets** — the MATLAB Actions install MATLAB on the free GitHub-hosted runner
with no license, and Pages/CI use the automatic `GITHUB_TOKEN`.

There are three workflows:

| Workflow file | Name | Purpose | Needs repo setup? |
|---|---|---|---|
| `tests.yml` | **Tests** | Fast unit + regression suite on every push/PR | Just enable Actions |
| `lint.yml` | **Lint** | `checkcode` static analysis on every push/PR | Just enable Actions |
| `jekyll-gh-pages.yml` | **Docs** | Build + deploy the API-reference site to GitHub Pages | **Pages + environment policy (below)** |

---

## 0. Enable Actions

**Settings → Actions → General**

- *Actions permissions*: **Allow all actions and reusable workflows** (the
  workflows use `actions/*`, `matlab-actions/*`, and `actions/jekyll-build-pages`).
- *Workflow permissions*: **Read repository contents** is enough for Tests/Lint;
  the Docs workflow requests the elevated scopes it needs in-file (see §3), so
  you do **not** need to widen the default here.

No secrets are required for any workflow. Do not add any.

---

## 1. Tests (`tests.yml`)

Runs the **fast** suite (`run_tests('fast')` = unit + regression) on
`ubuntu-latest` for every push to `main`/`develop` and every pull request.

```yaml
env:
  FORCE_JAVASCRIPT_ACTIONS_TO_NODE24: true   # silences the Node 20 deprecation warning
...
- uses: matlab-actions/setup-matlab@v2
  with:
    products: >
      Signal_Processing_Toolbox
      Statistics_and_Machine_Learning_Toolbox
      Parallel_Computing_Toolbox
- uses: matlab-actions/run-command@v2
  with:
    command: |
      cd(fullfile(getenv('GITHUB_WORKSPACE'), 'tests'));
      results = run_tests('fast');
      assert(~any([results.Failed]), 'run_tests: %d test(s) failed', sum([results.Failed]));
```

Notes for whoever inherits this:

- **No license, no self-hosted runner.** `matlab-actions/setup-matlab@v2`
  installs MATLAB on the free GitHub-hosted runner; this is free for **public**
  repos. (On a private repo you would need a MATLAB batch-licensing token.)
- **Toolboxes** are limited to the three nestapp actually uses (Signal
  Processing, Statistics, Parallel Computing). Image Processing / Wavelet /
  Optimization are intentionally omitted.
- **Integration tests are NOT run in CI.** They need EEGLAB + TESA on disk, which
  the runner does not have. Run them locally before merging anything that touches
  pipeline execution: `run_tests('all')` (see `.github/CONTRIBUTING.md`). There is
  deliberately no `[run-integration]` trigger.
- `FORCE_JAVASCRIPT_ACTIONS_TO_NODE24: true` is set because some actions still
  ship a Node 20 entrypoint that GitHub now warns about; this opts them into the
  installed Node 24.

**Setup required: none beyond §0.** It runs as soon as Actions is enabled.

---

## 2. Lint (`lint.yml`)

`checkcode` static analysis via `tools/run_lint.m`, on every push to
`main`/`develop` and every PR. Base MATLAB only (no toolboxes).

```yaml
- uses: matlab-actions/setup-matlab@v2
- uses: matlab-actions/run-command@v2
  with:
    command: |
      addpath(fullfile(getenv('GITHUB_WORKSPACE'), 'tools'));
      run_lint('error');
```

`run_lint('error')` fails the job on any `checkcode` message at severity *error*
(MDOTM/SYNER/UNDEF/DEFNU and friends); warnings are reported but non-fatal.

**Setup required: none beyond §0.**

---

## 3. Docs → GitHub Pages (`jekyll-gh-pages.yml`)

This is the workflow that needs real configuration. It runs `tools/gen_docs.m`
to generate the API reference from source headers, builds it with Jekyll, and
deploys to GitHub Pages. It triggers on **push to `main` or `develop`**, on
**tags `v*`**, and on manual dispatch.

It is **optional** — the docs site is a convenience, not part of validating code.
If you do not want a Pages site, simply delete this workflow. If you do, follow
all of the steps below or deploys will fail.

### 3a. Branch trigger

The workflow triggers on `push` to **`main` and `develop`** (and `v*` tags), so it
runs out of the box on a `main`-default repo — no branch renaming needed. If you
only want to publish from one branch, trim the `branches:` list (and then you only
need that branch in the environment policy below, §3d).

### 3b. Enable Pages with the Actions source

**Settings → Pages → Build and deployment → Source: GitHub Actions.**
(Not "Deploy from a branch" — this workflow publishes an artifact via the
`actions/deploy-pages` path.)

### 3c. Workflow permissions (already in the file)

The workflow declares the scopes Pages deployment requires; no UI change needed,
but verify they survive any edits:

```yaml
permissions:
  contents: read
  pages: write
  id-token: write
concurrency:
  group: "pages"
  cancel-in-progress: false
```

### 3d. The `github-pages` environment deployment-branch policy (the easy-to-miss one)

GitHub auto-creates a protected `github-pages` **environment** the first time you
enable Pages. By default it only allows deployments from the **default branch**.
Because this workflow deploys from `develop` (a non-default branch), you must add
`develop` to the environment's allowed branches, or every deploy fails with:

> Branch "develop" is not allowed to deploy to github-pages due to environment
> protection rules.

One-time fix (replace `OWNER/REPO`, e.g. `arefpz/nestapp`):

```bash
# allow the default branch (usually already allowed) AND develop:
gh api -X POST repos/OWNER/REPO/environments/github-pages/deployment-branch-policies -f name=main
gh api -X POST repos/OWNER/REPO/environments/github-pages/deployment-branch-policies -f name=develop
```

(Equivalently: **Settings → Environments → github-pages → Deployment branches →
add `develop`.**) If you changed the trigger to `main` in §3a, you only need the
`main` policy.

### 3e. What the build actually does (for reference)

1. `tools/gen_docs.m` walks `src/`, extracts each file's H1 + help block, and
   writes Markdown into `site/` (one page per folder + an index + the
   step-extension reference from `stepRegistry.m`).
2. A `site/_config.yml` is written (theme `jekyll-theme-cayman`, plugin
   `jekyll-relative-links`), `docs/architecture.md` is copied in, and an empty
   front-matter block is prepended to any generated page that lacks one (Jekyll
   only renders Markdown that starts with `---`).
3. `actions/jekyll-build-pages` builds `./site` → `./_site`, which
   `actions/upload-pages-artifact` + `actions/deploy-pages` publish.

The generated HTML is **never committed** — it is a derived artifact rebuilt on
every deploy, the same way the template `.mat` files are regenerated from
`buildTemplates.m`.

---

## Quick checklist

- [ ] Settings → Actions → General → allow all actions.
- [ ] (Tests, Lint) nothing else — they run on enable.
- [ ] (Docs, optional) Settings → Pages → Source: **GitHub Actions**.
- [ ] (Docs) it triggers on `main` and `develop` already — no branch change needed.
- [ ] (Docs) add the deploying branch(es) to the `github-pages` environment policy
      (`gh api ... deployment-branch-policies -f name=<branch>`).
- [ ] No secrets. Public repo. Free runners.
