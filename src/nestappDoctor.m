
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
%                      and .problems (cellstr of detected issues).
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

info.nestapp      = collectNestapp();
info.matlab       = collectMatlab();
info.toolboxes    = collectToolboxes();
info.eeglab       = collectEeglab();
info.dependencies = collectDependencies();
info.parallel     = collectParallel();
info.prefs        = collectPrefs();
info.shadows      = collectShadows();

info.problems = diagnose(info);

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
end

% ── collectors ──────────────────────────────────────────────────────────────

function s = collectNestapp()
s.version = nestappVersion();
s.gitCommit = '';
try
    root = fileparts(fileparts(mfilename('fullpath')));   % src/ -> repo root
    [st, out] = system(sprintf('git -C "%s" rev-parse --short HEAD', root));
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
% Format-specific loaders (those tied to a file extension, e.g. bva-io for
% .vhdr) are marked optional - their absence is reported but not a problem.
deps = struct('plugin', {}, 'fn', {}, 'found', {}, 'path', {}, 'optional', {});
try
    reg  = stepRegistry();
    reqs = [reg.requires];
catch
    reqs = struct('fn', {}, 'plugin', {}, 'fileExt', {});
end
seen = {};
for i = 1:numel(reqs)
    key = reqs(i).fn;
    if any(strcmp(seen, key)); continue; end
    seen{end+1} = key; %#ok<AGROW>
    isOptional = isfield(reqs, 'fileExt') && ~isempty(reqs(i).fileExt);
    deps(end+1) = struct('plugin', reqs(i).plugin, 'fn', key, ...
        'found', ~isempty(which(key)), 'path', which(key), ...
        'optional', isOptional); %#ok<AGROW>
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

function problems = diagnose(info)
problems = {};

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
    problems{end+1} = ['Curve Fitting Toolbox not found - the TMS-EEG / ' ...
        'AARATEP template''s "Remove Decay Artifact" step requires it.'];
end

% EEGLAB.
if ~info.eeglab.found
    problems{end+1} = ['EEGLAB not found on the MATLAB path. Set the EEGLAB ' ...
        'folder in Settings > Preferences (most steps need it).'];
end

% Key plugins (optional format-loaders are reported in the table, not here).
for i = 1:numel(info.dependencies)
    d = info.dependencies(i);
    if ~d.found && ~d.optional
        problems{end+1} = sprintf('%s not found (function "%s" missing).', ...
            d.plugin, d.fn); %#ok<AGROW>
    end
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
L{end+1} = 'Required toolboxes:';
allReq = [info.toolboxes.required, {'Curve Fitting Toolbox'}];
for t = allReq
    L{end+1} = sprintf('  [%s] %s', tick(hasToolbox(info, t{1})), t{1}); %#ok<AGROW>
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
