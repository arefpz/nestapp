function tests = test_qualityGateBatchMode
% TEST_QUALITYGATEBATCHMODE  End-to-end test of batch-relative verdicts.
%   Five synthetic files run through a Quality Gate in 'batch' mode.
%   Four have low trial amplitude; one is 20x larger -> batch MAD
%   detects the outlier and flags it Fail. skipOnQualityFail is on
%   but ignored (batch verdicts arrive after the run; logged warning).
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
r = repoRoot();
addpath(r);
addpath(genpath(fullfile(r, 'src')));

if isempty(which('eeglab'))
    testCase.assumeFail('EEGLAB is not on the MATLAB path.');
end
evalc('eeglab(''nogui'')');

if ispref('nestapp', 'skipOnQualityFail')
    testCase.TestData.prevSkip = getpref('nestapp', 'skipOnQualityFail');
    testCase.TestData.prevSkipSet = true;
else
    testCase.TestData.prevSkipSet = false;
end
end

function teardownOnce(testCase)
if testCase.TestData.prevSkipSet
    setpref('nestapp', 'skipOnQualityFail', testCase.TestData.prevSkip);
else
    if ispref('nestapp', 'skipOnQualityFail')
        rmpref('nestapp', 'skipOnQualityFail');
    end
end
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

% -- tests ----------------------------------------------------------------

function test_outlier_file_flagged_after_batch_resolution(testCase)
tmpDir = fullfile(tempdir, ['qg_batch_', char(matlab.lang.internal.uuid())]);
mkdir(tmpDir);
testCase.addTeardown(@() rmdir(tmpDir, 's'));

% Five files; first four small-amplitude (no saturated channels),
% file 5 huge-amplitude (channel 1 saturates). Seed RNG so the
% random draws don't tip a "low" file over the saturation threshold.
rng(42);
amps = [0.1 0.1 0.1 0.1 10];
setPaths = cell(1, 5);
for k = 1:5
    baseName = sprintf('sub_%d', k);
    setPath  = fullfile(tmpDir, [baseName, '.set']);
    EEG      = makeSyntheticEEG(amps(k));
    evalc(sprintf('pop_saveset(EEG, ''filename'', ''%s.set'', ''filepath'', tmpDir);', baseName));
    setPaths{k} = setPath;
end

setpref('nestapp', 'skipOnQualityFail', true);  % deliberately on - should be ignored

spec(1).name   = 'Load Data';
spec(1).params = struct();
spec(2).name   = 'Quality Gate';
spec(2).params = struct( ...
    'gateLabel',      'batch-amplitude', ...
    'thresholdMode',  'batch', ...
    'marginalSlack',  0.8, ...
    'maxSatChans',    1, ...   % nonzero -> nSatChans is collected
    'outlierSigmas',  3);

opts = struct( ...
    'uiFigure',     [], ...
    'pipelineName', 'qg-batch-test', ...
    'statusBar',    [], ...
    'parallel',     false, ...
    'chanLocFile',  '');

[allReports, ~] = runPipelineCore(spec, setPaths, opts);

testCase.verifyLength(allReports, 5, ...
    'All 5 files should complete in batch mode (no inline skip).');

verdicts = cellfun(@(r) r.quality.worstVerdict, allReports, 'UniformOutput', false);
% File 5 (the outlier) should be Fail; the others Pass or Marginal.
testCase.verifyEqual(verdicts{5}, 'Fail');
for k = 1:4
    testCase.verifyTrue(any(strcmp(verdicts{k}, {'Pass', 'Marginal'})), ...
        sprintf('Non-outlier file %d should not be Fail; got %s', k, verdicts{k}));
end
end

% -- helpers --------------------------------------------------------------

function EEG = makeSyntheticEEG(amplitude)
nChan   = 16;
nPnts   = 1000;
srate   = 500;
EEG.setname  = 'synth';
EEG.filename = '';
EEG.filepath = '';
EEG.nbchan   = nChan;
EEG.pnts     = nPnts;
EEG.srate    = srate;
EEG.trials   = 1;
EEG.xmin     = 0;
EEG.xmax     = (nPnts - 1) / srate;
EEG.times    = 1000 * (0:nPnts-1) / srate;
EEG.data     = single(amplitude * randn(nChan, nPnts));
% Force one channel to be very noisy so maxBadChanPct picks it up.
EEG.data(1, :) = single(amplitude * 50 * randn(1, nPnts));
EEG.event    = [];
EEG.epoch    = [];
EEG.icawinv   = [];
EEG.icaweights = [];
EEG.icasphere  = [];
EEG.icaact     = [];
EEG.icachansind = [];
EEG.ref       = 'common';
EEG.history   = '';
EEG.urevent   = [];
EEG.chanlocs  = struct('labels', {}, 'X', {}, 'Y', {}, 'Z', {});
for c = 1:nChan
    EEG.chanlocs(c).labels = sprintf('Ch%d', c);
    EEG.chanlocs(c).X = cos(2*pi*(c-1)/nChan);
    EEG.chanlocs(c).Y = sin(2*pi*(c-1)/nChan);
    EEG.chanlocs(c).Z = 0;
end
EEG = eeg_checkset(EEG);
end
