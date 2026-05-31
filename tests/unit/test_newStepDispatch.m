
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_newStepDispatch
% TEST_NEWSTEPDISPATCH  Regression tests for the ARTIST + AARATEP step
%   dispatch wiring. These would have caught:
%     - Calling helpers with a struct positional arg when the helper
%       declares opts.foo (e.g. "Invalid argument at position 2.
%       Function requires exactly 1 positional input(s).")
%     - addpath of the vendored AARATEP tree shadowing EEGLAB pop_*
%       functions (Common/ThirdParty/FromEEGLab/).
%
%   Run: runtests('tests/unit/test_newStepDispatch')
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));
addpath(fullfile(r, 'tests', 'helpers'));   % hideFromPath
end

function r = repoRoot()
r = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

% ── helper invocation: each new helper must accept the name-value pairs
%    that paramsToVarin produces from its registry defaults ──────────────────

function test_artistRejectBadTrials_acceptsDispatchedArgs(testCase)
EEG = makeFakeEpochedEEG();
vars = defaultsFor('Reject Bad Trials (ARTIST)');
fn = @() artistRejectBadTrials(EEG, vars{:});
testCase.verifyWarningFree(fn, ...
    'Helper must accept the name-value pairs produced by paramsToVarin.');
end

function test_artistFlagDecayICs_acceptsDispatchedArgs(testCase)
EEG = makeFakeEpochedEEGWithICA();
vars = defaultsFor('Flag ICA Components (ARTIST Decay)');
fn = @() artistFlagDecayICs(EEG, vars{:});
testCase.verifyWarningFree(fn, ...
    'Helper must accept the name-value pairs produced by paramsToVarin.');
end

function test_aaratepMuscleClassifier_acceptsDispatchedArgs(testCase)
EEG = makeFakeEpochedEEGWithICA();
vars = defaultsFor('Flag ICA Components (AARATEP Muscle)');
fn = @() aaratepMuscleClassifier(EEG, vars{:});
testCase.verifyWarningFree(fn, ...
    'Helper must accept the name-value pairs produced by paramsToVarin.');
end

function test_artistBadChannelsRansac_acceptsDispatchedArgs(testCase)
EEG = makeFakeEpochedEEG();
% pop_select needs to exist to run the helper to completion, but the
% arg-shape check happens before that. We assert the helper at least
% accepts the args without an arg-shape error.
vars = defaultsFor('Remove Bad Channels (ARTIST)');
fn = @() artistBadChannelsRansac(EEG, vars{:});
err = tryAndCatch(fn);
testCase.verifyFalse(isInvalidArgError(err), ...
    sprintf(['Helper rejected dispatched args with a signature error:' ...
             newline '%s'], err));
end

% ── dispatch case wiring: simulate what processOneFile actually does for
%    each new step name ──────────────────────────────────────────────────────

function test_dispatchPattern_artistRejectBadTrials(testCase)
verifyDispatchPattern(testCase, 'Reject Bad Trials (ARTIST)', ...
    'artistRejectBadTrials', makeFakeEpochedEEG());
end

function test_dispatchPattern_artistFlagDecayICs(testCase)
verifyDispatchPattern(testCase, 'Flag ICA Components (ARTIST Decay)', ...
    'artistFlagDecayICs', makeFakeEpochedEEGWithICA());
end

function test_dispatchPattern_aaratepMuscleClassifier(testCase)
verifyDispatchPattern(testCase, 'Flag ICA Components (AARATEP Muscle)', ...
    'aaratepMuscleClassifier', makeFakeEpochedEEGWithICA());
end

% ── path hygiene: ensureAaratepOnPath must not shadow EEGLAB pop_* funcs ──

function test_ensureAaratepOnPath_doesNotShadowEEGLabForks(testCase)
% Verify the FromEEGLab fork directory is NOT added to the MATLAB path
% when ensureAaratepOnPath runs. Start from a clean state so we measure
% only what ensureAaratepOnPath itself adds (the session may already have
% the dir on path from prior runs).
shadowDir = fullfile(repoRoot(), 'third_party', 'aaratep', ...
                     'Common', 'ThirdParty', 'FromEEGLab');
clear ensureAaratepOnPath;
prePath = path;
restorePath = onCleanup(@() path(prePath));

% Scrub any pre-existing FromEEGLab entries inherited from the session.
pathParts = strsplit(path, pathsep);
toDrop = pathParts(startsWith(pathParts, shadowDir));
for k = 1:numel(toDrop)
    rmpath(toDrop{k});
end

ensureAaratepOnPath();

pathParts = strsplit(path, pathsep);
shadowed  = any(startsWith(pathParts, shadowDir));

testCase.verifyFalse(shadowed, sprintf( ...
    ['ensureAaratepOnPath added the forked-EEGLAB subtree to the path. ' ...
     'This shadows EEGLAB''s topoplot, epoch, pop_loadbv, pop_resample. ' ...
     'Shadow root: %s'], shadowDir));
end

function test_ensureAaratepOnPath_addsCoreHelpers(testCase)
clear ensureAaratepOnPath;
prePath = path;
restorePath = onCleanup(@() path(prePath));

ensureAaratepOnPath();

testCase.verifyNotEmpty(which('c_EEG_ReplaceEpochTimeSegment'), ...
    'AR-Blend interp helper not on path after ensureAaratepOnPath.');
testCase.verifyNotEmpty(which('c_TMSEEG_fitAndRemoveDecayArtifact'), ...
    'Decay-fit helper not on path after ensureAaratepOnPath.');
end

function test_ensureAaratepOnPath_doesNotShadowBundledFastICA(testCase)
% The bundled FastICA must not be added to the path when another fastica
% (e.g. EEGLAB's) is already resolvable - users get the one they expect.
testCase.assumeNotEmpty(which('fastica'), ...
    'No fastica on the path - cannot verify the bundled one is not preferred.');
fasticaDir = fullfile(repoRoot(), 'third_party', 'aaratep', ...
                      'Common', 'ThirdParty', 'FastICA');
clear ensureAaratepOnPath;
prePath = path;
restorePath = onCleanup(@() restoreAaratepState(prePath));

ensureAaratepOnPath();

onPath = any(strcmp(strsplit(path, pathsep), fasticaDir));
testCase.verifyFalse(onPath, ...
    'ensureAaratepOnPath added the bundled FastICA dir, shadowing the user''s.');
resolved = which('fastica');
testCase.verifyFalse(strcmpi(resolved, fullfile(fasticaDir, 'fastica.m')), ...
    'fastica resolved to the bundled copy instead of the user''s install.');
end

function test_ensureAaratepOnPath_warnsOnFastICAMismatch(testCase)
% A user fastica whose version differs from the bundled one must trigger a
% one-time mismatch warning when the AARATEP tree is activated.
bundled = fullfile(repoRoot(), 'third_party', 'aaratep', ...
                  'Common', 'ThirdParty', 'FastICA', 'fastica.m');
testCase.assumeTrue(isfile(bundled), 'Bundled FastICA absent - cannot test mismatch.');

fakeDir = tempname; mkdir(fakeDir);
fid = fopen(fullfile(fakeDir, 'fastica.m'), 'w');
fprintf(fid, 'function varargout = fastica(varargin)\n%% FastICA version 99.9\nend\n');
fclose(fid);

clear ensureAaratepOnPath;
prePath = path;
restorePath = onCleanup(@() restoreAaratepState(prePath, fakeDir));
addpath(fakeDir);   % addpath prepends, so the fake resolves first

lastwarn('');
ws = warning('off', 'all');
ensureAaratepOnPath();
[~, wid] = lastwarn();
warning(ws);

testCase.verifyEqual(wid, 'nestapp:aaratepFastICAMismatch', ...
    'A FastICA version mismatch must raise nestapp:aaratepFastICAMismatch.');
end

function test_ensureAaratepOnPath_warnsOnFastICAFallback(testCase)
% With NO fastica on the path, ensureAaratepOnPath must add the bundled copy
% and warn nestapp:aaratepFastICAFallback. Hide any installed fastica so this
% branch runs deterministically even on a machine that has FastICA.
bundled = fullfile(repoRoot(), 'third_party', 'aaratep', ...
                  'Common', 'ThirdParty', 'FastICA', 'fastica.m');
testCase.assumeTrue(isfile(bundled), 'Bundled FastICA absent - cannot test fallback.');

cleanupHide = hideFromPath('fastica'); %#ok<NASGU>
testCase.assertEmpty(which('fastica'), 'Could not hide fastica from the path.');

clear ensureAaratepOnPath;
prePath = path;
restorePath = onCleanup(@() restoreAaratepState(prePath)); %#ok<NASGU>

lastwarn('');
ws = warning('off', 'all');
ensureAaratepOnPath();
[~, wid] = lastwarn();
warning(ws);

testCase.verifyEqual(wid, 'nestapp:aaratepFastICAFallback', ...
    'With no fastica installed, the bundled FastICA fallback warning must fire.');
end

% ── fixtures and helpers ─────────────────────────────────────────────────────

function restoreAaratepState(prePath, varargin)
% Restore the path and reset ensureAaratepOnPath's persistent so a
% path-mutating test cannot leave later tests with a stale path state.
path(prePath);
clear ensureAaratepOnPath;
for i = 1:numel(varargin)
    if isfolder(varargin{i})
        try; rmdir(varargin{i}, 's'); catch; end
    end
end
end

function EEG = makeFakeEpochedEEG()
% Minimal epoched EEG struct: 8 channels x 200 samples x 20 trials at 1 kHz,
% times in ms covering [-100, 100). Channel locations populated so RANSAC
% rejector doesn't bail out. Includes the housekeeping fields that
% eeg_checkset / eeg_getica expect (setname, comments, ref, urevent, ...).
nChan = 8; nTime = 200; nTrial = 20;
EEG = struct();
EEG.setname = 'fixture';
EEG.filename = '';
EEG.filepath = '';
EEG.subject = '';
EEG.group = '';
EEG.condition = '';
EEG.session = [];
EEG.comments = '';
EEG.ref = 'common';
EEG.data = randn(nChan, nTime, nTrial) * 5;
EEG.srate = 1000;
EEG.nbchan = nChan;
EEG.pnts = nTime;
EEG.trials = nTrial;
EEG.xmin = -0.1;
EEG.xmax = 0.099;
EEG.times = linspace(-100, 99.5, nTime);
EEG.chanlocs = struct('labels', {}, 'theta', {}, 'radius', {});
for c = 1:nChan
    EEG.chanlocs(c).labels = sprintf('Ch%d', c);
    EEG.chanlocs(c).theta  = 360 * c / nChan;
    EEG.chanlocs(c).radius = 0.5;
end
EEG.chaninfo = struct();
EEG.urchanlocs = EEG.chanlocs;
EEG.event = struct([]);
EEG.urevent = struct([]);
EEG.epoch = struct([]);
EEG.epochdescription = {};
EEG.reject = struct();
EEG.stats = struct();
EEG.specdata = [];
EEG.specicaact = [];
EEG.splinefile = '';
EEG.icasplinefile = '';
EEG.dipfit = [];
EEG.history = '';
EEG.saved = 'no';
EEG.etc = struct();
EEG.icaact = [];
EEG.icaweights = [];
EEG.icasphere = [];
EEG.icawinv = [];
EEG.icachansind = [];
end

function EEG = makeFakeEpochedEEGWithICA()
% As above plus an identity ICA decomposition so IC-flagging steps run.
EEG = makeFakeEpochedEEG();
nChan = EEG.nbchan;
EEG.icaweights = eye(nChan);
EEG.icasphere  = eye(nChan);
EEG.icawinv    = eye(nChan);
EEG.icachansind = 1:nChan;
EEG.reject.gcompreject = false(1, nChan);
end

function vars = defaultsFor(stepName)
% Build a name-value cell {key, val, ...} from the registry defaults for
% one step - exactly what paramsToVarin produces and what dispatch
% forwards via vars{:}.
reg  = stepRegistry();
step = makePipelineStep(stepName, reg);
vars = paramsToVarin(step.params);
end

function verifyDispatchPattern(testCase, stepName, helperName, EEG)
% Smoke-test that the helper-call form `helper(EEG, vars{:})` survives the
% arg-count check. We do not assert the algorithm completes; just that
% the function signature accepts the dispatched arguments.
vars = defaultsFor(stepName);
helperFn = str2func(helperName);
err = tryAndCatch(@() helperFn(EEG, vars{:}));
testCase.verifyFalse(isInvalidArgError(err), sprintf(...
    ['Dispatch wiring for step "%s" -> %s broken: %s' newline ...
     'This is the regression that caused' newline ...
     '  "Invalid argument at position 2. Function requires exactly 1 ' ...
     'positional input(s)."'], ...
    stepName, helperName, err));
end

function msg = tryAndCatch(fn)
try
    fn();
    msg = '';
catch ME
    msg = ME.message;
end
end

function tf = isInvalidArgError(msg)
% Match MATLAB's arguments-block mismatch messages.
tf = ~isempty(msg) && ( ...
    contains(msg, 'Invalid argument at position', 'IgnoreCase', true) || ...
    contains(msg, 'requires exactly',             'IgnoreCase', true) || ...
    contains(msg, 'requires at least',            'IgnoreCase', true) || ...
    contains(msg, 'requires at most',             'IgnoreCase', true) || ...
    contains(msg, 'Unexpected input argument',    'IgnoreCase', true) || ...
    contains(msg, 'Unknown name-value argument',  'IgnoreCase', true));
end
