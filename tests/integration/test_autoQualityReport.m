
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_autoQualityReport
% TEST_AUTOQUALITYREPORT  Integration test for the autoQualityReport flow.
%
%   End-to-end test of the visual quality screening (Phase 1):
%   - Sets the autoQualityReport preference on
%   - Saves a tiny synthetic EEG to disk as a .set file
%   - Runs a one-step pipeline (Load Data, which is a checkpoint)
%   - Asserts the <batchRoot>/qc/<filebase>/NN_StepName.png hierarchy
%     was created and that report.quality.figures records it.
%     (The batch root itself is named <yyyyMMdd_HHmmss>_<pipelineSlug>
%     by buildBatchContext - see io/outputPaths.m for the layout map.)
%
%   Requires EEGLAB on the MATLAB path. assumeFail if absent.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
r = repoRoot();
addpath(r);
addpath(genpath(fullfile(r, 'src')));

if isempty(which('eeglab'))
    testCase.assumeFail( ...
        'EEGLAB is not on the MATLAB path. Add EEGLAB before running this test.');
end

% Initialize EEGLAB silently so pop_saveset / pop_loadset are usable.
evalc('eeglab(''nogui'')');

% Snapshot the autoQualityReport pref so we can restore it.
if ispref('nestapp', 'autoQualityReport')
    testCase.TestData.prevPref     = getpref('nestapp', 'autoQualityReport');
    testCase.TestData.prevPrefSet  = true;
else
    testCase.TestData.prevPrefSet  = false;
end
end

function teardownOnce(testCase)
if ~isfield(testCase.TestData, 'prevPrefSet')
    return   % setupOnce assumeFail'd (no EEGLAB) before snapshotting prefs
end
if testCase.TestData.prevPrefSet
    setpref('nestapp', 'autoQualityReport', testCase.TestData.prevPref);
elseif ispref('nestapp', 'autoQualityReport')
    rmpref('nestapp', 'autoQualityReport');
end
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

% -- tests ----------------------------------------------------------------

function test_qc_folder_and_png_are_created(testCase)
% Tiny synthetic .set file in a temp folder.
tmpDir = fullfile(tempdir, ['nestapp_qc_', char(matlab.lang.internal.uuid())]);
mkdir(tmpDir);
testCase.addTeardown(@() rmdir(tmpDir, 's'));

baseName = 'tiny_synth';
setPath  = fullfile(tmpDir, [baseName, '.set']);
EEG      = makeSyntheticEEG();
evalc('pop_saveset(EEG, ''filename'', [baseName, ''.set''], ''filepath'', tmpDir);');

% Enable auto quality report. QC images are now captured after every
% Quality Gate step (no hardcoded checkpoint list anymore), so the
% recipe ends with a permissive Quality Gate to trigger the render.
setpref('nestapp', 'autoQualityReport', true);

spec(1).name   = 'Load Data';
spec(1).params = struct();
spec(2).name   = 'Quality Gate';
spec(2).params = struct('gateLabel', 'auto-qc-smoke');

opts = struct( ...
    'uiFigure',     [], ...
    'pipelineName', 'qc-integration-test', ...
    'statusBar',    [], ...
    'parallel',     false, ...
    'chanLocFile',  '', ...
    'outputRoot',   tmpDir);   % isolate from any user outputRoot pref

% Pop a UIFigure-free dlg: runPipelineCore tolerates an empty uiFigure by
% creating a standalone progress figure. We close it via the dlg cleanup.
[allReports, ~] = runPipelineCore(spec, {setPath}, opts);

% --- Assertions ---
testCase.verifyNotEmpty(allReports, 'Pipeline produced no reports.');
report = allReports{1};

testCase.verifyTrue(isfield(report, 'quality'), 'report missing quality field');
testCase.verifyTrue(isfield(report.quality, 'figures'), ...
    'report.quality missing figures field');
testCase.verifyNotEmpty(report.quality.figures, ...
    'report.quality.figures should list at least one PNG.');

pngPath = report.quality.figures{1};
testCase.verifyTrue(exist(pngPath, 'file') == 2, ...
    sprintf('Recorded QC PNG was not written to disk: %s', pngPath));

% Layout under the batch root (typeBased default):
%   <batchRoot>/qc/<stem>/NN_StepName.png
% so fileparts(pngPath)         -> <batchRoot>/qc/<stem>
%    fileparts(fileparts(...))  -> <batchRoot>/qc
%    fileparts(...3 levels up)  -> <batchRoot>
perFileDir = fileparts(pngPath);
qcDir      = fileparts(perFileDir);
batchRoot  = fileparts(qcDir);

[~, qcFolder] = fileparts(qcDir);
testCase.verifyEqual(qcFolder, 'qc', ...
    sprintf('Top-level QC folder should be "qc": got %s', qcFolder));

[~, perFileFolder] = fileparts(perFileDir);
testCase.verifyEqual(perFileFolder, baseName, ...
    'Per-file subfolder should match input basename.');

% PNG filename should be NN_<sanitized step>.png.
[~, pngBase] = fileparts(pngPath);
testCase.verifyMatches(pngBase, '^\d{2}_[A-Za-z0-9_]+$', ...
    sprintf('PNG name does not match NN_StepName pattern: %s', pngBase));

% Batch root should be <yyyyMMdd_HHmmss>_<pipelineSlug>.
[~, batchFolder] = fileparts(batchRoot);
testCase.verifyMatches(batchFolder, '^\d{8}_\d{6}_[a-z0-9_]+$', ...
    sprintf('Batch folder name does not match expected pattern: %s', batchFolder));

% Batch-level artifacts (spec.mat is the cheapest one to assert).
specPath = fullfile(batchRoot, 'batch', 'spec.mat');
testCase.verifyTrue(exist(specPath, 'file') == 2, ...
    'Expected <batchRoot>/batch/spec.mat to be written.');
end

% -- helpers --------------------------------------------------------------

function EEG = makeSyntheticEEG()
nChan   = 8;
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
EEG.data     = randn(nChan, nPnts, 'single');
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
