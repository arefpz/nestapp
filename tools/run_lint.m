
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function results = run_lint(failOn)
% RUN_LINT  Static-analysis lint of all nestapp source via checkcode.
%
%   Syntax
%     run_lint                 % report; error if any 'error'-severity message
%     run_lint('error')        % same as default
%     run_lint('warning')      % stricter: error on warnings too
%     results = run_lint(...)  % return findings, do not throw
%
%   Inputs
%     failOn  char - severity that makes lint fail: 'error' (default) or
%                    'warning'. Ignored when an output is captured.
%
%   Outputs
%     results  struct array of findings: .file, .line, .id, .severity, .message.
%              When called with no output, run_lint prints a summary and
%              errors if any finding at or above failOn is present (for CI).
%
%   Scans every .m file under src/ with MATLAB's built-in checkcode. Mirrors
%   the static analysis used during development.
%
%   See also: checkcode, run_tests

if nargin < 1 || isempty(failOn)
    failOn = 'error';
end

toolsRoot = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(toolsRoot);
srcRoot   = fullfile(repoRoot, 'src');

files = dir(fullfile(srcRoot, '**', '*.m'));
results = struct('file', {}, 'line', {}, 'id', {}, 'severity', {}, 'message', {});

for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    msgs  = checkcode(fpath, '-id', '-struct');
    for k = 1:numel(msgs)
        sev = severityOf(msgs(k).id);
        results(end+1) = struct( ...
            'file',     relpath(repoRoot, fpath), ...
            'line',     msgs(k).line(1), ...
            'id',       msgs(k).id, ...
            'severity', sev, ...
            'message',  msgs(k).message); %#ok<AGROW>
    end
end

nErr  = sum(strcmp({results.severity}, 'error'));
nWarn = sum(strcmp({results.severity}, 'warning'));

fprintf('\nrun_lint: %d file(s) scanned, %d error(s), %d warning(s)\n', ...
    numel(files), nErr, nWarn);
for i = 1:numel(results)
    fprintf('  [%s] %s:%d  %s (%s)\n', upper(results(i).severity), ...
        results(i).file, results(i).line, results(i).message, results(i).id);
end

if nargout == 0
    if strcmpi(failOn, 'warning')
        nBad = nErr + nWarn;
    else
        nBad = nErr;
    end
    if nBad > 0
        error('run_lint: %d lint finding(s) at or above "%s" severity.', nBad, failOn);
    end
end
end

% ── helpers ───────────────────────────────────────────────────────────────

function sev = severityOf(id)
% checkcode does not expose severity directly; treat a small set of
% correctness-relevant message ids as errors and the rest as style warnings.
errorIds = { ...
    'MDOTM', ...   % file/function name mismatch
    'SYNER', ...   % syntax error
    'UNDEF', ...   % undefined function or variable (where detectable)
    'DEFNU'};      % function defined but never used (likely a mistake)
if any(strcmp(id, errorIds))
    sev = 'error';
else
    sev = 'warning';
end
end

function r = relpath(root, fpath)
r = strrep(fpath, [root filesep], '');
end
