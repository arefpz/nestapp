function tests = test_qualityGateSkip
% TEST_QUALITYGATESKIP  End-to-end test of the skip-on-fail short circuit.
%   Recipe: Load Data -> Quality Gate(expectedChans=999) -> Save New Set
%   With skipOnQualityFail=true and an 8-channel synthetic file, the
%   gate fails and the rest of the pipeline is skipped. The failure
%   is tagged kind='skipped' by parseFailure and surfaces in the
%   post-run summary as a Quality Gate skip (not an error).
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
evalc('eeglab(''nogui'')');

if ispref('nestapp', 'skipOnQualityFail')
    testCase.TestData.prevPref = getpref('nestapp', 'skipOnQualityFail');
    testCase.TestData.prevPrefSet = true;
else
    testCase.TestData.prevPrefSet = false;
end
end

function teardownOnce(testCase)
if testCase.TestData.prevPrefSet
    setpref('nestapp', 'skipOnQualityFail', testCase.TestData.prevPref);
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

function test_skipOnFail_short_circuits_pipeline(testCase)
% Two files. The gate insists on 999 channels; both fail (each has 8).
% Save New Set in step 3 should NOT run for any file. Both files end
% up in the failure list with kind='skipped'.
tmpDir = fullfile(tempdir, ['qg_skip_', char(matlab.lang.internal.uuid())]);
mkdir(tmpDir);
testCase.addTeardown(@() rmdir(tmpDir, 's'));

setPaths = cell(1, 2);
for k = 1:2
    baseName = sprintf('tiny_%d', k);
    setPath  = fullfile(tmpDir, [baseName, '.set']);
    EEG      = makeSyntheticEEG();
    evalc(sprintf('pop_saveset(EEG, ''filename'', ''%s.set'', ''filepath'', tmpDir);', baseName));
    setPaths{k} = setPath;
end

setpref('nestapp', 'skipOnQualityFail', true);
setpref('nestapp', 'autoQualityReport', false);

spec(1).name   = 'Load Data';
spec(1).params = struct();
spec(2).name   = 'Quality Gate';
spec(2).params = struct( ...
    'gateLabel',     'must-have-999-chans', ...
    'thresholdMode', 'absolute', ...
    'marginalSlack', 0.8, ...
    'expectedChans', 999, ...
    'outlierSigmas', 3);
spec(3).name   = 'Save New Set';
spec(3).params = struct('savenew', 'should_never_run', 'includeFileName', 'yes');

opts = struct( ...
    'uiFigure',     [], ...
    'pipelineName', 'qg-skip-test', ...
    'statusBar',    [], ...
    'parallel',     false, ...
    'chanLocFile',  '', ...
    'outputRoot',   tmpDir);

% Both files fail -> runPipelineCore throws nestapp:allFilesFailed.
% That's the expected end-state; the value of this test is the
% side-effect that Save New Set never ran.
testCase.verifyError(@() runPipelineCore(spec, setPaths, opts), ...
    'nestapp:allFilesFailed');

for k = 1:2
    [~, baseName] = fileparts(setPaths{k});
    saveTarget = fullfile(tmpDir, [baseName, '_should_never_run.set']);
    testCase.verifyFalse(exist(saveTarget, 'file') == 2, ...
        sprintf('Save New Set should not have run for %s.', baseName));
end
end

function test_skipOff_lets_pipeline_continue_after_fail(testCase)
% Same recipe, but skipOnQualityFail=false (advisory): the gate
% records Fail in the report; the rest of the pipeline still runs.
tmpDir = fullfile(tempdir, ['qg_advisory_', char(matlab.lang.internal.uuid())]);
mkdir(tmpDir);
testCase.addTeardown(@() rmdir(tmpDir, 's'));

baseName = 'tiny_synth';
setPath  = fullfile(tmpDir, [baseName, '.set']);
EEG      = makeSyntheticEEG();
evalc('pop_saveset(EEG, ''filename'', [baseName, ''.set''], ''filepath'', tmpDir);');

setpref('nestapp', 'skipOnQualityFail', false);
setpref('nestapp', 'autoQualityReport', false);

spec(1).name   = 'Load Data';
spec(1).params = struct();
spec(2).name   = 'Quality Gate';
spec(2).params = struct( ...
    'gateLabel',     'advisory', ...
    'thresholdMode', 'absolute', ...
    'marginalSlack', 0.8, ...
    'expectedChans', 999, ...
    'outlierSigmas', 3);

opts = struct( ...
    'uiFigure',     [], ...
    'pipelineName', 'qg-advisory-test', ...
    'statusBar',    [], ...
    'parallel',     false, ...
    'chanLocFile',  '', ...
    'outputRoot',   tmpDir);

allReports = runPipelineCore(spec, {setPath}, opts);

testCase.verifyNotEmpty(allReports);
report = allReports{1};
testCase.verifyTrue(isfield(report, 'quality'));
testCase.verifyEqual(report.quality.worstVerdict, 'Fail');
testCase.verifyEqual(report.quality.gates{1}.label, 'advisory');
testCase.verifyEqual(report.quality.gates{1}.verdict, 'Fail');
testCase.verifyNotEmpty(report.quality.gates{1}.reasons);
end

function test_unset_gate_labels_get_step_indexed_names(testCase)
% Two Quality Gates in the same pipeline, neither with an explicit
% gateLabel. processOneFile must rewrite them to gate-NN using the
% step index so they don't collide in the dashboard / batch finalizer.
tmpDir = fullfile(tempdir, ['qg_autoname_', char(matlab.lang.internal.uuid())]);
mkdir(tmpDir);
testCase.addTeardown(@() rmdir(tmpDir, 's'));

baseName = 'tiny_synth';
setPath  = fullfile(tmpDir, [baseName, '.set']);
EEG      = makeSyntheticEEG();
evalc('pop_saveset(EEG, ''filename'', [baseName, ''.set''], ''filepath'', tmpDir);');

setpref('nestapp', 'skipOnQualityFail', false);
setpref('nestapp', 'autoQualityReport', false);

% Three-step recipe: Load Data, two unnamed Quality Gates. Both gates
% are permissive (empty params -> all checks disabled -> Pass) so the
% test runs to completion and we can inspect the labels.
spec(1).name   = 'Load Data';
spec(1).params = struct();
spec(2).name   = 'Quality Gate';
spec(2).params = struct();          % no gateLabel set
spec(3).name   = 'Quality Gate';
spec(3).params = struct('gateLabel', 'gate');   % literal default value

opts = struct( ...
    'uiFigure',     [], ...
    'pipelineName', 'qg-autoname-test', ...
    'statusBar',    [], ...
    'parallel',     false, ...
    'chanLocFile',  '', ...
    'outputRoot',   tmpDir);

allReports = runPipelineCore(spec, {setPath}, opts);
report = allReports{1};
gates = report.quality.gates;
testCase.verifyLength(gates, 2);

% Step 2 and step 3 -> gate-02 and gate-03. Distinct labels so the
% dashboard renders one row per gate, not one collapsed row.
testCase.verifyEqual(gates{1}.label, 'gate-02');
testCase.verifyEqual(gates{2}.label, 'gate-03');
testCase.verifyNotEqual(gates{1}.label, gates{2}.label);
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
