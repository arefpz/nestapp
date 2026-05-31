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
%     Prints pass/fail summary to the MATLAB command window.
%     Non-zero exit status on any failure (useful for CI).
%
%   See also: runtests

function run_tests(suite)
if nargin < 1
    suite = 'fast';   % unit + regression, no EEGLAB
end

testRoot = fileparts(mfilename('fullpath'));
repoRoot = fileparts(testRoot);
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
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

results = matlab.unittest.TestResult.empty;
for i = 1:numel(suites)
    if ~exist(suites{i}, 'dir')
        warning('run_tests: suite directory not found: %s', suites{i});
        continue
    end
    r = runtests(suites{i});
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
