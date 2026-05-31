function tests = test_pipeline_integration
% TEST_PIPELINE_INTEGRATION  Smoke tests using EEGLAB built-in sample data.
%
%   These tests verify that core EEGLAB signal processing operations
%   (load, filter, epoch, re-reference, save/reload) behave as expected
%   with the version of EEGLAB bundled in the project directory.
%
%   REQUIREMENTS:
%     - EEGLAB 2025.0.0 must be on the MATLAB path (eeglab('nogui') succeeds)
%     - The standard EEGLAB sample data must be present in the EEGLAB
%       'sample_data' subdirectory (bundled with EEGLAB)
%
%   Run: runtests('tests/integration/test_pipeline_integration')
tests = functiontests(localfunctions);
end

% ── setup ─────────────────────────────────────────────────────────────────

function setupOnce(testCase)
% Initialise EEGLAB without GUI; locate sample data.
if ~exist('eeglab', 'file')
    testCase.assumeFail('EEGLAB not on path — skipping integration tests');
end
global EEG ALLEEG CURRENTSET ALLCOM  %#ok<GVMIS>
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');
testCase.TestData.sampleFile = findSampleData();
if isempty(testCase.TestData.sampleFile)
    testCase.assumeFail('EEGLAB sample data not found — skipping integration tests');
end
end

function filePath = findSampleData()
% Locate the eeglab_data.set sample file in any EEGLAB installation on path.
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

% ── tests ─────────────────────────────────────────────────────────────────

function test_eeglabInitialises(testCase)
global EEG ALLEEG CURRENTSET  %#ok<GVMIS>
testCase.verifyTrue(isstruct(ALLEEG), 'ALLEEG must be a struct after eeglab(''nogui'')');
testCase.verifyTrue(isnumeric(CURRENTSET), 'CURRENTSET must be numeric');
end

function test_loadSampleData(testCase)
global EEG  %#ok<GVMIS>
[fPath, fName, fExt] = fileparts(testCase.TestData.sampleFile);
EEG = pop_loadset([fName, fExt], fPath);
testCase.verifyGreaterThan(EEG.nbchan, 0, 'Loaded EEG must have channels');
testCase.verifyGreaterThan(EEG.srate, 0,  'Loaded EEG must have positive srate');
testCase.verifyFalse(isempty(EEG.data),   'Loaded EEG must have data');
end

function test_resampleStep(testCase)
global EEG  %#ok<GVMIS>
[fPath, fName, fExt] = fileparts(testCase.TestData.sampleFile);
EEG          = pop_loadset([fName, fExt], fPath);
origSrate    = EEG.srate;
targetSrate  = origSrate / 2;
EEGr = pop_resample(EEG, targetSrate);
testCase.verifyEqual(EEGr.srate, targetSrate, 'Resampled srate must match target');
testCase.verifyEqual(EEGr.nbchan, EEG.nbchan, 'Resample must not change channel count');
end

function test_bandpassFilterStep(testCase)
global EEG  %#ok<GVMIS>
[fPath, fName, fExt] = fileparts(testCase.TestData.sampleFile);
EEG  = pop_loadset([fName, fExt], fPath);
EEGf = pop_eegfiltnew(EEG, 1, 40);   % 1–40 Hz bandpass
testCase.verifyEqual(EEGf.nbchan, EEG.nbchan, 'Filter must not change channel count');
testCase.verifyEqual(EEGf.srate,  EEG.srate,  'Filter must not change srate');
end

function test_rereferenceStep(testCase)
global EEG  %#ok<GVMIS>
[fPath, fName, fExt] = fileparts(testCase.TestData.sampleFile);
EEG  = pop_loadset([fName, fExt], fPath);
EEGr = pop_reref(EEG, []);   % average reference
testCase.verifyEqual(size(EEGr.data, 2), size(EEG.data, 2), ...
    'Re-reference must not change time dimension');
end

function test_epochingStep(testCase)
global EEG  %#ok<GVMIS>
[fPath, fName, fExt] = fileparts(testCase.TestData.sampleFile);
EEG = pop_loadset([fName, fExt], fPath);
% Only epoch if there are events
if isempty(EEG.event); return; end
eventTypes = unique({EEG.event.type});
EEGe = pop_epoch(EEG, eventTypes, [-0.5 0.5]);
testCase.verifyGreaterThan(EEGe.trials, 1, 'Epoching must produce multiple trials');
end

function test_saveAndReloadSet(testCase)
global EEG  %#ok<GVMIS>
[fPath, fName, fExt] = fileparts(testCase.TestData.sampleFile);
EEG  = pop_loadset([fName, fExt], fPath);
tmpPath = [tempname, '.set'];
[tmpDir, tmpName, ~] = fileparts(tmpPath);
pop_saveset(EEG, 'filename', [tmpName, '.set'], 'filepath', tmpDir);
EEGr = pop_loadset([tmpName, '.set'], tmpDir);
testCase.verifyEqual(EEGr.nbchan, EEG.nbchan, 'Round-trip save/load preserves channel count');
testCase.verifyEqual(EEGr.srate,  EEG.srate,  'Round-trip save/load preserves srate');
% Clean up
delete(tmpPath);
fdt = strrep(tmpPath, '.set', '.fdt');
if exist(fdt, 'file'); delete(fdt); end
end
