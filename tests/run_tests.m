%% run_tests  Run the full nestapp test suite (or a subset).
%
%   USAGE
%     run_tests            % unit + regression (no EEGLAB required)
%     run_tests('all')     % unit + regression + integration (EEGLAB required)
%     run_tests('unit')    % unit tests only
%     run_tests('regression') % regression tests only
%     run_tests('integration') % integration tests only (EEGLAB required)
%
%   OUTPUT
%     results  (optional) matlab.unittest.TestResult array. When called with
%              no output (e.g. interactively), run_tests errors on any
%              failure so the failure is impossible to miss. When the caller
%              captures results (e.g. CI), it inspects [results.Failed]
%              itself and run_tests does not throw.
%     Also prints a pass/fail summary to the MATLAB command window.
%
%   See also: runtests

% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.

function results = run_tests(suite)
if nargin < 1
    suite = 'fast';   % unit + regression, no EEGLAB
end

testRoot = fileparts(mfilename('fullpath'));
repoRoot = fileparts(testRoot);
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'src')));
addpath(fullfile(testRoot, 'helpers'));

switch lower(suite)
    case 'fast'
        suites = {fullfile(testRoot, 'unit'), fullfile(testRoot, 'regression')};
    case 'unit'
        suites = {fullfile(testRoot, 'unit')};
    case 'regression'
        suites = {fullfile(testRoot, 'regression')};
    case 'integration'
        suites = {fullfile(testRoot, 'integration')};
    case 'all'
        suites = {fullfile(testRoot, 'unit'), ...
                  fullfile(testRoot, 'regression'), ...
                  fullfile(testRoot, 'integration')};
    otherwise
        error('run_tests: unknown suite "%s". Valid: fast, unit, regression, integration, all', suite);
end

% Per-test progress: log every test's start/end to the console and a file, so
% a hanging or slow test is easy to pinpoint (the last START with no matching
% END is the culprit, and the file survives a killed session). Built with an
% explicit TestRunner so the logger plugin can attach to the standard text
% output.
logFile = fullfile(tempdir, 'nestapp_test_progress.log');
if exist(logFile, 'file'); delete(logFile); end
fprintf('Per-test progress log: %s\n\n', logFile);

% Sweep stale scratch folders from prior runs before starting. Integration
% tests clean up their own temp dirs on teardown, but a run that is killed
% mid-way (e.g. a frozen session closed by hand) leaves them orphaned in
% tempdir. Removing them here keeps test artifacts from piling up over time.
sweepStaleTestArtifacts();

runner = matlab.unittest.TestRunner.withTextOutput;
runner.addPlugin(TestProgressLogger(logFile));

results = matlab.unittest.TestResult.empty;
for i = 1:numel(suites)
    if ~exist(suites{i}, 'dir')
        warning('run_tests: suite directory not found: %s', suites{i});
        continue
    end
    s = matlab.unittest.TestSuite.fromFolder(suites{i});
    r = runner.run(s);
    results = [results, r]; %#ok<AGROW>
end

%% Summary
nPass = sum([results.Passed]);
nFail = sum([results.Failed]);
nInc  = sum([results.Incomplete]);
nTot  = numel(results);

fprintf('\n');
fprintf('══════════════════════════════════════════════\n');
fprintf('  nestapp test suite — %s\n', suite);
fprintf('══════════════════════════════════════════════\n');
fprintf('  Total:      %3d\n', nTot);
fprintf('  Passed:     %3d\n', nPass);
if nFail > 0
    fprintf('  FAILED:     %3d  ← fix before committing\n', nFail);
else
    fprintf('  Failed:     %3d\n', nFail);
end
if nInc > 0
    fprintf('  Incomplete: %3d\n', nInc);
end
fprintf('══════════════════════════════════════════════\n\n');

%% Print failed test names for quick diagnosis
if nFail > 0
    fprintf('Failed tests:\n');
    for i = 1:numel(results)
        if results(i).Failed
            fprintf('  ✗ %s\n', results(i).Name);
        end
    end
    fprintf('\n');
end

if nFail > 0 && nargout == 0
    % In CI contexts: exit with non-zero status
    % (MATLAB does not have a direct exit code, but the caller can check)
    error('run_tests: %d test(s) failed.', nFail);
end
end

function sweepStaleTestArtifacts()
% Remove orphaned integration-test scratch folders left in tempdir by runs
% that did not tear down cleanly. Only nestapp's own test-folder prefixes are
% touched; each prefix matches the unique-UUID dirs the integration tests
% create (e.g. tempdir/nestapp_qc_<uuid>). Best-effort: a folder still locked
% by another process is skipped rather than erroring the run.
prefixes = {'nestapp_qc_', 'qg_skip_', 'qg_advisory_', 'qg_autoname_', 'qg_batch_'};
for pi = 1:numel(prefixes)
    stale = dir(fullfile(tempdir, [prefixes{pi} '*']));
    for k = 1:numel(stale)
        if ~stale(k).isdir; continue; end
        try
            rmdir(fullfile(stale(k).folder, stale(k).name), 's');
        catch
            % Locked or already gone - leave it for the next sweep.
        end
    end
end
end
