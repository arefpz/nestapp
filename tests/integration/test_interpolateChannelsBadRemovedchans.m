function tests = test_interpolateChannelsBadRemovedchans
% TEST_INTERPOLATECHANNELSBADREMOVEDCHANS  Regression test for chanlocs loss.
%
%   Reproduces the failure seen on S08_1_pre_singlepulse_CHAN.set: the input
%   carried stale entries in EEG.chaninfo.removedchans with numeric-string
%   labels and EMPTY coordinates. The old "Interpolate Channels" step fed
%   ALL of removedchans to pop_interp, which appended the bogus channels to
%   chanlocs (without matching data rows). eeg_checkset then discarded the
%   whole chanlocs struct on the size mismatch (chanlocs -> []), and the
%   next "Re-Reference" step crashed dereferencing EEG.chanlocs with:
%       "Dot indexing is not supported for variables of type double."
%
%   After the fix, "Interpolate Channels" ignores removed-channel records
%   that are not genuinely absent or that lack coordinates, so chanlocs
%   survives and Re-Reference succeeds.
%
%   Run: runtests('tests/integration/test_interpolateChannelsBadRemovedchans')
tests = functiontests(localfunctions);
end

% ── setup ─────────────────────────────────────────────────────────────────

function setupOnce(testCase)
if ~exist('eeglab', 'file')
    testCase.assumeFail('EEGLAB not on path — skipping integration test');
end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..', 'src')));
global EEG ALLEEG CURRENTSET ALLCOM  %#ok<GVMIS>
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');
sf = findSampleData();
if isempty(sf)
    testCase.assumeFail('EEGLAB sample data not found — skipping integration test');
end
testCase.TestData.sampleFile = sf;
end

function filePath = findSampleData()
eeglabRoot = fileparts(which('eeglab'));
candidates = {
    fullfile(eeglabRoot, 'sample_data', 'eeglab_data.set'), ...
    fullfile(eeglabRoot, 'sample_data', 'EEG.set')
};
filePath = '';
for i = 1:numel(candidates)
    if exist(candidates{i}, 'file')
        filePath = candidates{i};
        return
    end
end
end

% ── test ──────────────────────────────────────────────────────────────────

function test_staleRemovedchansDoesNotDropChanlocs(testCase)
global EEG  %#ok<GVMIS>

% Load sample data and bake in a bogus "removed channel": numeric-string
% label, no coordinates - exactly the shape that destroyed chanlocs.
[p, n, e] = fileparts(testCase.TestData.sampleFile);
EEG = pop_loadset([n e], p);

bad = EEG.chanlocs(1);
for f = fieldnames(bad)'
    bad.(f{1}) = [];
end
bad.labels = '99';                       % not a real electrode name
EEG.chaninfo.removedchans = bad;         % single coordinate-less record

tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() rmdir(tmpDir, 's'));
fn = 'badremoved.set';
EEG = pop_saveset(EEG, 'filename', fn, 'filepath', tmpDir);

% Confirm the bogus record survived the save/load round trip, otherwise
% the test would silently pass for the wrong reason.
EEGcheck = pop_loadset(fn, tmpDir);
testCase.assertEqual(numel(EEGcheck.chaninfo.removedchans), 1, ...
    'Setup precondition: bogus removedchans record must persist in the .set file');

% Minimal pipeline that exercises the real step code paths.
reg  = stepRegistry();
spec = [ makePipelineStep('Load Data', reg), ...
         makePipelineStep('Interpolate Channels', reg), ...
         makePipelineStep('Re-Reference', reg) ];
spec(3).params.ref = '[]';   % average reference -> no interactive prompt

opts = struct('pipelineName', 'reproTest', 'fileIndex', 1);
[~, stepLog] = processOneFile(spec, fullfile(tmpDir, fn), opts);

% chanlocs must survive interpolation as a populated struct.
testCase.verifyTrue(isstruct(EEG.chanlocs), ...
    'chanlocs must remain a struct after Interpolate Channels');
testCase.verifyNotEmpty(EEG.chanlocs, ...
    'chanlocs must not be discarded by the chanlocs/data size mismatch');

% Re-Reference must have run without the dot-indexing crash.
rerefRec = stepLog(strcmp({stepLog.step}, 'Re-Reference'));
testCase.verifyNotEmpty(rerefRec, 'Re-Reference step must have executed');
testCase.verifyEqual(rerefRec(end).error, '', ...
    'Re-Reference must complete without error');
end
