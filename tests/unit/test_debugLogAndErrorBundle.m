
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_debugLogAndErrorBundle
% TEST_DEBUGLOGANDERRORBUNDLE  Unit tests for the debug log tee
%   (nestDebugLog + nestLog) and the metadata-only error bundle
%   (saveErrorBundle). No EEGLAB required.
%
%   Run: runtests('tests/unit/test_debugLogAndErrorBundle')
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));
end

function teardown(testCase) %#ok<INUSD>
% Always close any debug log a test left open.
try; nestDebugLog('stop'); catch; end
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

% ── debug log tee ─────────────────────────────────────────────────────────────

function test_debugLogCapturesNestLogLines(testCase)
folder = tempname; mkdir(folder);
restore = onCleanup(@() rmdirSafe(folder));
logPath = nestDebugLog('start', folder);
nestLog('TEST', 'marker %d', 7);
nestDebugLog('stop');
testCase.verifyTrue(isfile(logPath), 'debug log file must be created');
txt = fileread(logPath);
testCase.verifyTrue(contains(txt, 'marker 7'), ...
    'nestLog lines must be teed into the debug log file');
end

function test_debugLogStopClearsGlobal(testCase)
folder = tempname; mkdir(folder);
restore = onCleanup(@() rmdirSafe(folder));
nestDebugLog('start', folder);
nestDebugLog('stop');
% After stop, nestLog must not error and nothing more is written.
testCase.verifyWarningFree(@() nestLog('TEST', 'after stop'));
end

function test_nestLogWorksWithNoDebugFile(testCase)
% With no debug log open, nestLog must behave normally (no error).
nestDebugLog('stop');   % ensure closed
testCase.verifyWarningFree(@() nestLog('TEST', 'plain %s', 'line'));
end

% ── error bundle (metadata only) ──────────────────────────────────────────────

function bd = makeBundle()
try
    error('nestapp:test', 'simulated failure');
catch e
end
EEG = struct('nbchan', 61, 'trials', 79, 'pnts', 2000, 'srate', 1000, ...
    'xmin', -1, 'xmax', 1, 'data', zeros(61, 50, 3), ...
    'history', 'pop_a; pop_b;', 'chanlocs', struct('labels', {'Fz', 'Cz'}));
spec = struct('name', {'Load Data', 'Re-Sample'}, ...
    'params', {struct(), struct('freq', 500)});
bd = saveErrorBundle(tempname, struct('err', e, 'EEG', EEG, 'spec', spec, ...
    'stepName', 'Re-Sample', 'stepIndex', 2, 'fileName', 'sub01.set', ...
    'pipelineName', 'Test'));
end

function test_bundleWritesExpectedFiles(testCase)
bd = makeBundle();
restore = onCleanup(@() rmdirSafe(bd));
testCase.verifyTrue(isfolder(bd), 'bundle folder must be created');
for f = {'error.txt', 'eeg_metadata.txt', 'environment.md', 'pipeline.md'}
    testCase.verifyTrue(isfile(fullfile(bd, f{1})), ...
        sprintf('bundle must contain %s', f{1}));
end
end

function test_bundleRecordsErrorAndStep(testCase)
bd = makeBundle();
restore = onCleanup(@() rmdirSafe(bd));
txt = fileread(fullfile(bd, 'error.txt'));
testCase.verifyTrue(contains(txt, 'simulated failure'), 'error message must be recorded');
testCase.verifyTrue(contains(txt, 'Re-Sample'), 'failing step must be recorded');
end

function test_bundleMetadataExcludesRawData(testCase)
% The metadata file must contain the data SHAPE but not the raw values.
bd = makeBundle();
restore = onCleanup(@() rmdirSafe(bd));
meta = fileread(fullfile(bd, 'eeg_metadata.txt'));
testCase.verifyTrue(contains(meta, 'data size'), 'metadata must record the data shape');
testCase.verifyTrue(contains(meta, 'channels: Fz Cz'), 'metadata must list channel labels');
% Heuristic: the raw zeros matrix (61x50x3 = 9150 values) is not dumped;
% the file stays small.
info = dir(fullfile(bd, 'eeg_metadata.txt'));
testCase.verifyLessThan(info.bytes, 4096, ...
    'metadata file is suspiciously large - raw data may have leaked in');
end

function test_bundleNoErrorReturnsEmpty(testCase)
bd = saveErrorBundle(tempname, struct());   % no .err
testCase.verifyEmpty(bd, 'with no error in ctx, no bundle should be written');
end

% ── on-demand support bundle ──────────────────────────────────────────────────

function test_supportBundleWritesEnvironmentAndPipeline(testCase)
spec = struct('name', {'Load Data', 'Re-Sample'}, ...
    'params', {struct(), struct('freq', 500)});
bd = collectSupportBundle(tempname, spec);
restore = onCleanup(@() rmdirSafe(bd));
testCase.verifyTrue(isfolder(bd), 'support bundle folder must be created');
testCase.verifyTrue(isfile(fullfile(bd, 'environment.md')), 'must include environment.md');
testCase.verifyTrue(isfile(fullfile(bd, 'pipeline.md')),    'must include pipeline.md');
testCase.verifyTrue(isfile(fullfile(bd, 'README.md')),      'must include README.md');
end

function test_supportBundleWithoutPipeline(testCase)
bd = collectSupportBundle(tempname);   % no spec
restore = onCleanup(@() rmdirSafe(bd));
testCase.verifyTrue(isfile(fullfile(bd, 'environment.md')), ...
    'environment must be captured even with no pipeline');
testCase.verifyFalse(isfile(fullfile(bd, 'pipeline.md')), ...
    'no pipeline.md when no spec is given');
end

% ── helpers ───────────────────────────────────────────────────────────────────

function rmdirSafe(d)
try
    if isfolder(d); rmdir(d, 's'); end
catch
end
end
