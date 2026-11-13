
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [report, info] = nestappDoctor(varargin)
% NESTAPPDOCTOR  Collect and validate the nestapp runtime environment.
%
%   Syntax
%     nestappDoctor                      % print a diagnostics report
%     report = nestappDoctor(...)        % also return the Markdown report
%     [report, info] = nestappDoctor(...)% also return the raw data struct
%     nestappDoctor('Copy', true)        % copy the report to the clipboard
%     nestappDoctor('Quiet', true)       % do not print to the command window
%
%   Inputs (name-value)
%     Copy   logical - copy the Markdown report to the system clipboard
%                      (default false).
%     Quiet  logical - suppress printing to the command window (default false).
%
%   Outputs
%     report  char   - a Markdown diagnostics report, ready to paste into a
%                      GitHub issue (see .github/ISSUE_TEMPLATE/bug_report.yml).
%     info    struct - the collected data: .nestapp, .matlab, .toolboxes,
%                      .eeglab, .dependencies, .parallel, .prefs, .shadows,
%                      .problems (cellstr: what stops a default install
%                      working) and .notes (cellstr: optional extras absent).
%
%   nestappDoctor diagnoses, it does not only dump: the report ends with a
%   "Problems detected" section listing actionable issues (missing EEGLAB,
%   absent required toolboxes, shadowed functions, unsupported MATLAB
%   version, ...). The dependency list is derived from stepRegistry, so it
%   stays in sync as steps are added. Runs with or without EEGLAB on the path.
%
%   Side effects   Optionally writes to the system clipboard; prints to the
%                  command window unless Quiet is set.
%   See also: nestappVersion, checkStepDependencies, stepRegistry

p = inputParser;
p.addParameter('Copy',  false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Quiet', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
doCopy  = logical(p.Results.Copy);
beQuiet = logical(p.Results.Quiet);

info = struct();
info.problems = {};
info.notes    = {};

% Do what the app does on launch first - EEGLAB from Preferences (or found in
% its usual places), the AARATEP helpers from where installAaratep puts them -
% so the report describes what the app will see. Without this, a check run at
% the prompt reported EEGLAB and AARATEP missing when both were installed.
ensureEeglabReady();
try
    ensureAaratepOnPath();
catch
    % Not installed: reported below as a note.
end

info.nestapp      = collectNestapp();
info.matlab       = collectMatlab();
info.toolboxes    = collectToolboxes();
info.eeglab       = collectEeglab();
info.dependencies = collectDependencies();
info.parallel     = collectParallel();
info.prefs        = collectPrefs();
info.shadows      = collectShadows();

[info.problems, info.notes] = diagnose(info);

report = buildReport(info);

if ~beQuiet
    fprintf('%s\n', report);
end
if doCopy
    try
        clipboard('copy', report);
        if ~beQuiet
            fprintf('(diagnostics copied to clipboard)\n');
        end
    catch ME
        warning('nestappDoctor:clipboard', ...
            'Could not copy to clipboard: %s', ME.message);
    end
end
if nargout == 0
    clear report   % typed without a semicolon, it would print a second time as ans
end
end

% ── collectors ──────────────────────────────────────────────────────────────

function s = collectNestapp()
s.version = nestappVersion();
s.gitCommit = '';
try
    [st, out] = system(sprintf('git -C "%s" rev-parse --short HEAD', nestappRoot()));
    if st == 0
        s.gitCommit = strtrim(out);
    end
catch
end
end

function s = collectMatlab()
s.version = version;
try
    s.release = version('-release');
catch
    s.release = '';
end
s.computer = computer;
end

function s = collectToolboxes()
% Names nestapp steps care about; presence is checked, version recorded.
s.required = {'Signal Processing Toolbox', 'Statistics and Machine Learning Toolbox'};
s.optional = {'Curve Fitting Toolbox', 'Parallel Computing Toolbox', ...
              'Image Processing Toolbox', 'Wavelet Toolbox'};
s.installed = struct('name', {}, 'version', {});
v = ver;
for i = 1:numel(v)
    s.installed(end+1) = struct('name', v(i).Name, ...
        'version', strtrim([v(i).Version ' ' v(i).Release]));
end
end

function tf = hasToolbox(info, name)
tf = any(strcmp({info.toolboxes.installed.name}, name));
end

function s = collectEeglab()
s.found   = ~isempty(which('eeglab'));
s.path    = '';
s.version = '';
if s.found
    s.path = fileparts(which('eeglab'));
end
if ~isempty(which('eeg_getversion'))
    try
        s.version = eeg_getversion();
    catch
    end
end
end

function deps = collectDependencies()
% Derive the dependency set from the step registry so it never goes stale.
% Some requirements are optional - their absence is reported but is not a
% problem, because a default install is not expected to have them:
%   - format loaders (tied to a file extension, e.g. bva-io for .vhdr);
%   - newer-version extras (a minVersion, e.g. the TESA 1.2 steps, which the
%     extension manager cannot install yet);
%   - the AARATEP helpers, which only the AARATEP template uses;
%   - MATLAB toolboxes (a license feature), which diagnose judges from the
%     toolbox list instead, so they are not reported twice.
deps = struct('plugin', {}, 'fn', {}, 'found', {}, 'path', {}, ...
              'optional', {}, 'minVersion', {}, 'installNote', {}, 'nSteps', {});
try
    reg = stepRegistry();
catch
    return
end
for s = 1:numel(reg)
    r = reg(s).requires;
    for i = 1:numel(r)
        k = find(strcmp({deps.fn}, r(i).fn), 1);
        if ~isempty(k)
            deps(k).nSteps = deps(k).nSteps + 1;
            continue
        end
        wh = which(r(i).fn);
        isOptional = ~isempty(r(i).fileExt) || ~isempty(r(i).minVersion) || ...
                     ~isempty(r(i).feature) || startsWith(r(i).plugin, 'AARATEP');
        deps(end+1) = struct('plugin', r(i).plugin, 'fn', r(i).fn, ...
            'found', ~isempty(wh), 'path', wh, 'optional', isOptional, ...
            'minVersion', r(i).minVersion, 'installNote', r(i).installNote, ...
            'nSteps', 1); %#ok<AGROW>
    end
end
end

function s = collectParallel()
s.licensed = license('test', 'Distrib_Computing_Toolbox');
s.poolSize = 0;
if s.licensed
    try
        pool = gcp('nocreate');
        if ~isempty(pool)
            s.poolSize = pool.NumWorkers;
        end
    catch
    end
end
end

function s = collectPrefs()
s = struct();
if ispref('nestapp')
    try
        s = getpref('nestapp');
    catch
    end
end
end

function shadows = collectShadows()
% Functions known to collide between MATLAB, EEGLAB, and the vendored
% AARATEP forks. A genuine shadow is two or more plain function files on the
% path with the same name. Class/package methods (paths under @ or +) and
% the canonical built-in are not shadows, so they are filtered out.
names = {'version', 'topoplot', 'epoch', 'pop_loadbv', 'pop_resample', ...
         'fastica', 'runica', 'eeglab', 'finputcheck'};
shadows = struct('name', {}, 'paths', {});
for i = 1:numel(names)
    hits = which(names{i}, '-all');
    if ~iscell(hits); continue; end
    isReal = cellfun(@(p) isempty(regexp(p, '[\\/][@+]', 'once')) && ...
                          ~startsWith(p, 'built-in'), hits);
    real = hits(isReal);
    if numel(real) > 1
        shadows(end+1) = struct('name', names{i}, 'paths', {real}); %#ok<AGROW>
    end
end
end

% ── diagnosis ─────────────────────────────────────────────────────────────────

function [problems, notes] = diagnose(info)
% problems: things that stop a default install working. notes: optional
% extras that are absent, each with what it would add. A fresh install with
% the documented steps must come back with no problems.
problems = {};
notes    = {};

% MATLAB version floor (R2023b).
rel = info.matlab.release;   % e.g. '2026a'
yr  = str2double(regexp(rel, '^\d+', 'match', 'once'));
if ~isnan(yr) && (yr < 2023 || (yr == 2023 && contains(rel, 'a')))
    problems{end+1} = sprintf(['MATLAB %s is below the supported floor ' ...
        '(R2023b). Some App Designer features may not work.'], rel);
end

% Required toolboxes.
for t = info.toolboxes.required
    if ~hasToolbox(info, t{1})
        problems{end+1} = sprintf('Required toolbox missing: %s.', t{1}); %#ok<AGROW>
    end
end
if ~hasToolbox(info, 'Curve Fitting Toolbox')
    notes{end+1} = ['Curve Fitting Toolbox not installed - only the AARATEP ' ...
        'template''s "Remove Decay Artifact" step needs it.'];
end

deps = info.dependencies;
if ~info.eeglab.found
    % Every plugin lives inside EEGLAB, so listing each missing function
    % here would be one cause reported dozens of times.
    problems{end+1} = ['EEGLAB not found. Download it from ' ...
        'https://eeglab.org/others/How_to_download_EEGLAB.html, unzip it, and ' ...
        'set its folder in Settings > Preferences. EEGLAB plugins are checked ' ...
        'once it is found.'];
    deps = deps(startsWith({deps.plugin}, 'AARATEP'));   % the only non-EEGLAB plugin
end

% One line per plugin, however many of its functions are missing. A plugin
% that is missing outright is a problem; only when it is present do its
% optional extras (TESA's 1.2 steps) get a note of their own.
missing = deps(~[deps.found]);
plugins = unique({missing.plugin}, 'stable');
for k = 1:numel(plugins)   % not "for p = plugins": an empty list would iterate once
    plugin = plugins(k);
    m   = missing(strcmp({missing.plugin}, plugin{1}));
    req = m(~[m.optional]);
    opt = m([m.optional]);
    if ~isempty(req)
        problems{end+1} = sprintf('%s not found - %d step(s) need it. %s', ...
            plugin{1}, sum([req.nSteps]), req(1).installNote); %#ok<AGROW>
    elseif ~isempty(opt(1).minVersion)
        notes{end+1} = sprintf('%d step(s) need %s %s or later and are hidden. %s', ...
            sum([opt.nSteps]), plugin{1}, opt(1).minVersion, opt(1).installNote); %#ok<AGROW>
    elseif startsWith(plugin{1}, 'AARATEP')
        notes{end+1} = sprintf(['AARATEP helpers not installed - only the ' ...
            'TMS-EEG / AARATEP template needs them. %s'], opt(1).installNote); %#ok<AGROW>
    end
    % Format loaders (bva-io, ...) are shown in the table only: they matter
    % just for that file type, and loading such a file names them.
end

% Shadowed functions.
for i = 1:numel(info.shadows)
    problems{end+1} = sprintf(['Function "%s" is shadowed (%d copies on the ' ...
        'path) - this can cause wrong-version calls. See the Shadows section.'], ...
        info.shadows(i).name, numel(info.shadows(i).paths)); %#ok<AGROW>
end
end

% ── report ─────────────────────────────────────────────────────────────────────

function md = buildReport(info)
L = {};
L{end+1} = '### nestapp diagnostics';
L{end+1} = '';
L{end+1} = '```text';
L{end+1} = sprintf('nestapp:  %s%s', info.nestapp.version, ...
    ternary(isempty(info.nestapp.gitCommit), '', [' (' info.nestapp.gitCommit ')']));
L{end+1} = sprintf('MATLAB:   %s', info.matlab.version);
L{end+1} = sprintf('Computer: %s', info.matlab.computer);
L{end+1} = sprintf('EEGLAB:   %s', describeEeglab(info.eeglab));
L{end+1} = sprintf('Parallel: %s', describeParallel(info.parallel));
L{end+1} = '';
L{end+1} = 'Dependencies:';
for i = 1:numel(info.dependencies)
    d = info.dependencies(i);
    suffix = '';
    if d.optional; suffix = '  (optional)'; end
    L{end+1} = sprintf('  [%s] %-22s (%s)%s', tick(d.found), d.plugin, d.fn, suffix); %#ok<AGROW>
end
L{end+1} = '';
L{end+1} = 'Toolboxes:';
for t = info.toolboxes.required
    L{end+1} = sprintf('  [%s] %s', tick(hasToolbox(info, t{1})), t{1}); %#ok<AGROW>
end
for t = info.toolboxes.optional
    L{end+1} = sprintf('  [%s] %s  (optional)', tick(hasToolbox(info, t{1})), t{1}); %#ok<AGROW>
end
L{end+1} = '```';
L{end+1} = '';

if ~isempty(info.shadows)
    L{end+1} = '<details><summary>Shadowed functions</summary>';
    L{end+1} = '';
    L{end+1} = '```text';
    for i = 1:numel(info.shadows)
        L{end+1} = sprintf('%s:', info.shadows(i).name); %#ok<AGROW>
        for k = 1:numel(info.shadows(i).paths)
            L{end+1} = sprintf('  %s', info.shadows(i).paths{k}); %#ok<AGROW>
        end
    end
    L{end+1} = '```';
    L{end+1} = '</details>';
    L{end+1} = '';
end

L{end+1} = '#### Problems detected';
if isempty(info.problems)
    L{end+1} = '';
    L{end+1} = 'None - environment looks healthy.';
else
    L{end+1} = '';
    for i = 1:numel(info.problems)
        L{end+1} = sprintf('- %s', info.problems{i}); %#ok<AGROW>
    end
end
if ~isempty(info.notes)
    L{end+1} = '';
    L{end+1} = '#### Optional extras not installed';
    L{end+1} = '';
    for i = 1:numel(info.notes)
        L{end+1} = sprintf('- %s', info.notes{i}); %#ok<AGROW>
    end
end

md = strjoin(L, newline);
end

function s = describeEeglab(e)
if ~e.found
    s = 'NOT FOUND';
else
    s = sprintf('%s  %s', firstNonEmpty(e.version, '(version unknown)'), e.path);
end
end

function s = describeParallel(p)
if ~p.licensed
    s = 'PCT not licensed';
elseif p.poolSize > 0
    s = sprintf('PCT licensed, pool active (%d workers)', p.poolSize);
else
    s = 'PCT licensed, no active pool';
end
end

function m = tick(tf)
if tf; m = 'x'; else; m = ' '; end
end

function out = ternary(cond, a, b)
if cond; out = a; else; out = b; end
end

function out = firstNonEmpty(a, b)
if ~isempty(a); out = a; else; out = b; end
end
