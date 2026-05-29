
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_pipelineTemplates
% TEST_PIPELINETEMPLATES  Unit tests for the built-in pipeline template .mat files.
%
%   Verifies that each template has the correct data shape, that all step
%   names exist in stepRegistry, and that key ordering constraints hold.
%   Templates are stored as .mat files in src/templates/ — this test suite
%   loads them directly without EEGLAB.
%
%   Run: runtests('tests/unit/test_pipelineTemplates')
tests = functiontests(localfunctions);
end

% ── setup ─────────────────────────────────────────────────────────────────

function setupOnce(testCase) %#ok<INUSD>
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

function templates = loadTemplates(testCase)
% Load all .mat files from src/templates/ and return a struct array.
% v3 format: pipelineName (string), spec (struct array with .name/.params), version.
r        = repoRoot();
matFiles = dir(fullfile(r, 'src', 'templates', '*.mat'));
testCase.verifyFalse(isempty(matFiles), 'No template .mat files found in src/templates/');
templates = struct('name', {}, 'steps', {}, 'spec', {});
for i = 1:numel(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    templates(i).name  = data.pipelineName;
    templates(i).steps = {data.spec.name};
    templates(i).spec  = data.spec;
end
end

function allNames = validStepNames()
steps    = stepRegistry();
allNames = {steps.name};
end

% ── .mat file shape ────────────────────────────────────────────────────────

function test_templatesCanBeLoaded(testCase)
templates = loadTemplates(testCase);
testCase.verifyGreaterThan(numel(templates), 0, 'Must have at least one template');
end

function test_allTemplatesHaveRequiredFields(testCase)
templates = loadTemplates(testCase);
for i = 1:numel(templates)
    testCase.verifyFalse(isempty(templates(i).name), ...
        sprintf('Template %d: pipelineName must not be empty', i));
    testCase.verifyFalse(isempty(templates(i).steps), ...
        sprintf('Template %d: spec must not be empty', i));
end
end

function test_sixTemplatesExist(testCase)
templates = loadTemplates(testCase);
testCase.verifyEqual(numel(templates), 6, 'Expected exactly 6 templates');
end

% ── template names ────────────────────────────────────────────────────────

function test_tmsEEGTemplateExists(testCase)
templates = loadTemplates(testCase);
testCase.verifyTrue(any(strcmp({templates.name}, 'TMS-EEG / TEP (TESA)')), ...
    'Must have a TMS-EEG template');
end

function test_restingStateTemplateExists(testCase)
templates = loadTemplates(testCase);
testCase.verifyTrue(any(contains({templates.name}, 'Resting')), ...
    'Must have a Resting-State template');
end

function test_minimalTemplateExists(testCase)
templates = loadTemplates(testCase);
testCase.verifyTrue(any(contains({templates.name}, 'Minimal')), ...
    'Must have a Minimal template');
end

% ── step validity ─────────────────────────────────────────────────────────

function test_allTmsEEGStepsInRegistry(testCase)
templates = loadTemplates(testCase);
t = templates(strcmp({templates.name}, 'TMS-EEG / TEP (TESA)'));
allNames = validStepNames();
for i = 1:numel(t.steps)
    testCase.verifyTrue(ismember(t.steps{i}, allNames), ...
        sprintf('TMS-EEG step "%s" not found in stepRegistry', t.steps{i}));
end
end

function test_allRestingStateStepsInRegistry(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'Resting'));
allNames = validStepNames();
for i = 1:numel(t.steps)
    testCase.verifyTrue(ismember(t.steps{i}, allNames), ...
        sprintf('Resting-State step "%s" not found in stepRegistry', t.steps{i}));
end
end

function test_allMinimalStepsInRegistry(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'Minimal'));
allNames = validStepNames();
for i = 1:numel(t.steps)
    testCase.verifyTrue(ismember(t.steps{i}, allNames), ...
        sprintf('Minimal step "%s" not found in stepRegistry', t.steps{i}));
end
end

% ── step counts ──────────────────────────────────────────────────────────

function test_tmsEEGHasAtLeast5Steps(testCase)
templates = loadTemplates(testCase);
t = templates(strcmp({templates.name}, 'TMS-EEG / TEP (TESA)'));
testCase.verifyGreaterThanOrEqual(numel(t.steps), 5, ...
    'TMS-EEG template must have at least 5 steps');
end

function test_restingStateHasAtLeast5Steps(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'Resting'));
testCase.verifyGreaterThanOrEqual(numel(t.steps), 5, ...
    'Resting-State template must have at least 5 steps');
end

function test_minimalHasAtLeast3Steps(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'Minimal'));
testCase.verifyGreaterThanOrEqual(numel(t.steps), 3, ...
    'Minimal template must have at least 3 steps');
end

% ── key ordering constraints ─────────────────────────────────────────────

function test_tmsEEGLoadDataIsFirst(testCase)
templates = loadTemplates(testCase);
t = templates(strcmp({templates.name}, 'TMS-EEG / TEP (TESA)'));
testCase.verifyEqual(t.steps{1}, 'Load Data', ...
    'TMS-EEG template must start with Load Data');
end

function test_tmsEEGSaveNewSetIsLast(testCase)
templates = loadTemplates(testCase);
t = templates(strcmp({templates.name}, 'TMS-EEG / TEP (TESA)'));
testCase.verifyEqual(t.steps{end}, 'Save New Set', ...
    'TMS-EEG template must end with Save New Set');
end

function test_tmsEEGFindPulsesBeforeRemove(testCase)
templates = loadTemplates(testCase);
t = templates(strcmp({templates.name}, 'TMS-EEG / TEP (TESA)'));
findIdx   = find(strcmp(t.steps, 'Find TMS Pulses (TESA)'),      1);
removeIdx = find(strcmp(t.steps, 'Remove TMS Artifacts (TESA)'), 1);
testCase.verifyTrue(~isempty(findIdx) && ~isempty(removeIdx), ...
    'Must have both Find TMS Pulses and Remove TMS Artifacts steps');
testCase.verifyLessThan(findIdx, removeIdx, ...
    'Find TMS Pulses must come before Remove TMS Artifacts');
end

function test_restingStateHasLoadDataFirst(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'Resting'));
testCase.verifyEqual(t.steps{1}, 'Load Data', ...
    'Resting-State template must start with Load Data');
end

% ── parameter values ───────────────────────────────────────────────────────

function test_restingStateHasFrequencyFilter(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'Resting'));
filterIdx = find(strcmp(t.steps, 'Frequency Filter'), 1);
testCase.verifyFalse(isempty(filterIdx), 'Resting-State must have Frequency Filter step');
locVal = t.spec(filterIdx).params.locutoff;
testCase.verifyGreaterThan(locVal, 0, 'Resting-State HPF locutoff must be > 0 Hz');
end

function test_minimalHasFrequencyFilter(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'Minimal'));
filterIdx = find(strcmp(t.steps, 'Frequency Filter'), 1);
testCase.verifyFalse(isempty(filterIdx), 'Minimal must have Frequency Filter step');
locVal = t.spec(filterIdx).params.locutoff;
testCase.verifyGreaterThan(locVal, 0, 'Minimal HPF locutoff must be > 0 Hz');
end

% ── ARTIST template ───────────────────────────────────────────────────────

function test_artistTemplateExists(testCase)
templates = loadTemplates(testCase);
testCase.verifyTrue(any(contains({templates.name}, 'ARTIST')), ...
    'Must have an ARTIST template');
end

function test_allArtistStepsInRegistry(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'ARTIST'));
allNames = validStepNames();
for i = 1:numel(t.steps)
    testCase.verifyTrue(ismember(t.steps{i}, allNames), ...
        sprintf('ARTIST step "%s" not found in stepRegistry', t.steps{i}));
end
end

function test_artistLoadDataIsFirst(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'ARTIST'));
testCase.verifyEqual(t.steps{1}, 'Load Data', ...
    'ARTIST template must start with Load Data');
end

function test_artistSaveNewSetIsLast(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'ARTIST'));
testCase.verifyEqual(t.steps{end}, 'Save New Set', ...
    'ARTIST template must end with Save New Set');
end

function test_artistEpochWindowMatchesPaper(testCase)
% Wu 2018 §2.2.2: "-500 to +1500 ms by default" -> [-0.5, 1.5] s.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'ARTIST'));
epochIdx = find(strcmp(t.steps, 'Epoching'), 1);
testCase.verifyEqual(t.spec(epochIdx).params.timelim, [-0.5, 1.5], ...
    'ARTIST epoch window must be [-0.5, 1.5] s per Wu 2018.');
end

function test_artistDownsampleMatchesPaper(testCase)
% Wu 2018 §2.2.1: "downsampled to 1 kHz".
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'ARTIST'));
rsIdx = find(strcmp(t.steps, 'Re-Sample'), 1);
testCase.verifyEqual(t.spec(rsIdx).params.freq, 1000, ...
    'ARTIST must downsample to 1000 Hz per Wu 2018.');
end

function test_artistBaselineMatchesPaper(testCase)
% Wu 2018 §2.2.3: "-300 to -100 ms baseline by default". This is the
% final pre-save TEP baseline. The first (early) baseline is a
% full-epoch demean for ICA stability and is verified separately by
% test_artistFirstBaselineIsFullEpoch.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'ARTIST'));
blIdx = find(strcmp(t.steps, 'Remove Baseline'), 1, 'last');
testCase.verifyEqual(t.spec(blIdx).params.timerange, [-300, -100], ...
    'ARTIST final baseline window must be [-300, -100] ms per Wu 2018 §2.2.3.');
end

% ── AARATEP template ──────────────────────────────────────────────────────

function test_aaratepTemplateExists(testCase)
templates = loadTemplates(testCase);
testCase.verifyTrue(any(contains({templates.name}, 'AARATEP')), ...
    'Must have an AARATEP template');
end

function test_allAaratepStepsInRegistry(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'AARATEP'));
allNames = validStepNames();
for i = 1:numel(t.steps)
    testCase.verifyTrue(ismember(t.steps{i}, allNames), ...
        sprintf('AARATEP step "%s" not found in stepRegistry', t.steps{i}));
end
end

function test_aaratepLoadDataIsFirst(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'AARATEP'));
testCase.verifyEqual(t.steps{1}, 'Load Data', ...
    'AARATEP template must start with Load Data');
end

function test_aaratepSaveNewSetIsLast(testCase)
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'AARATEP'));
testCase.verifyEqual(t.steps{end}, 'Save New Set', ...
    'AARATEP template must end with Save New Set');
end

function test_aaratepSoundLambdaMatchesSource(testCase)
% AARATEP source: SOUNDlambda = 10^-1.5.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'AARATEP'));
soundIdx = find(strcmp(t.steps, 'Remove Recording Noise (SOUND)'), 1);
testCase.verifyEqual(t.spec(soundIdx).params.lambdaValue, 10^-1.5, ...
    'AbsTol', 1e-6, ...
    'AARATEP SOUND lambda must be 10^-1.5 per c_TMSEEG_Preprocess_AARATEPPipeline.m');
end

function test_aaratepDownsampleMatchesSource(testCase)
% AARATEP source: downsampleTo = 1000.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'AARATEP'));
rsIdx = find(strcmp(t.steps, 'Re-Sample'), 1);
testCase.verifyEqual(t.spec(rsIdx).params.freq, 1000, ...
    'AARATEP must downsample to 1000 Hz per upstream defaults.');
end

function test_aaratepArtifactWindowMatchesSource(testCase)
% AARATEP source: artifactTimespan = [-0.002, 0.012] s -> [-2, 12] ms.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'AARATEP'));
arIdx = find(strcmp(t.steps, 'Interpolate Missing Data (AR-Blend)'), 1);
testCase.verifyEqual(t.spec(arIdx).params.artifactStartMs, -2, ...
    'AARATEP AR-Blend start must be -2 ms.');
testCase.verifyEqual(t.spec(arIdx).params.artifactEndMs, 12, ...
    'AARATEP AR-Blend end must be 12 ms.');
end

% ── CleanLine ordering: never put CleanLine after Epoching ───────────────
% CleanLine's default 4-s window spans trial boundaries on epoched data
% and prompts the user with 'y/n' or creates discontinuity artifacts. Every
% built-in template that uses CleanLine must run it on continuous data.

function test_cleanlineNeverAfterEpoching(testCase)
templates = loadTemplates(testCase);
for i = 1:numel(templates)
    t = templates(i);
    epochIdx = find(strcmp(t.steps, 'Epoching'), 1);
    cleanIdx = find(strcmp(t.steps, 'Frequency Filter (CleanLine)'), 1);
    if isempty(epochIdx) || isempty(cleanIdx)
        continue
    end
    testCase.verifyLessThan(cleanIdx, epochIdx, sprintf( ...
        ['Template "%s" runs CleanLine after Epoching. CleanLine ' ...
         'uses sliding-window spectral fits and prompts when its ' ...
         'window spans trial boundaries. Use Frequency Filter (TESA) ' ...
         'bandstop 58-62 Hz on epoched data instead.'], t.name));
end
end

function test_artistDoesNotUseCleanlineAfterEpoch(testCase)
% Explicit ARTIST guard - documents that we use TESA bandstop 58-62 Hz
% as the paper-faithful substitute for the FIR notch on epoched data.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'ARTIST'));
testCase.verifyTrue(any(strcmp(t.steps, 'Frequency Filter (TESA)')), ...
    'ARTIST template must include Frequency Filter (TESA) as the 60 Hz notch.');
testCase.verifyFalse(any(strcmp(t.steps, 'Frequency Filter (CleanLine)')), ...
    ['ARTIST template must NOT use CleanLine - it would span trial ' ...
     'boundaries since ARTIST epochs before filtering.']);
end

% ── Paper-fidelity audit guards ──────────────────────────────────────────

function test_artistFiltersBeforeStage2(testCase)
% Wu 2018 §2.2.2 lists filtering as part of "Preprocessing" alongside
% epoch/baseline; §2.2.1 places stage 2 (reject) AFTER preprocessing.
% Filtering must therefore come before bad-trial / bad-channel rejection.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'ARTIST'));
filtIdx   = find(strcmp(t.steps, 'Frequency Filter'),                1);
notchIdx  = find(strcmp(t.steps, 'Frequency Filter (TESA)'),         1);
rejTrIdx  = find(strcmp(t.steps, 'Reject Bad Trials (ARTIST)'),      1);
rejChIdx  = find(strcmp(t.steps, 'Remove Bad Channels (ARTIST)'),    1);
testCase.verifyLessThan(filtIdx,  rejTrIdx, ...
    'ARTIST: bandpass must run before stage-2 trial rejection.');
testCase.verifyLessThan(notchIdx, rejTrIdx, ...
    'ARTIST: 60 Hz notch must run before stage-2 trial rejection.');
testCase.verifyLessThan(filtIdx,  rejChIdx, ...
    'ARTIST: bandpass must run before stage-2 channel rejection.');
testCase.verifyLessThan(notchIdx, rejChIdx, ...
    'ARTIST: 60 Hz notch must run before stage-2 channel rejection.');
end

function test_artistICAUsesRunica(testCase)
% Wu 2018 §2.2.1: "Infomax algorithm" -> runica in EEGLAB. Both rounds.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'ARTIST'));
icaIdx = find(strcmp(t.steps, 'Run ICA'));
testCase.verifyNumElements(icaIdx, 2, ...
    'ARTIST must have two Run ICA steps (round 1 decay, round 2 classify).');
for k = 1:numel(icaIdx)
    testCase.verifyEqual(t.spec(icaIdx(k)).params.icatype, 'runica', ...
        sprintf('ARTIST Run ICA #%d must use runica (Infomax) per Wu 2018.', k));
end
end

function test_artistFirstBaselineIsFullEpoch(testCase)
% Pre-ICA baseline should be a full-epoch demean for ICA stability
% (matches template 1's first baseline). The final pre-save baseline
% [-300, -100] ms is the TEP analysis window from Wu 2018 §2.2.3.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'ARTIST'));
blIdx = find(strcmp(t.steps, 'Remove Baseline'));
testCase.verifyEqual(numel(blIdx), 2, ...
    'ARTIST must have two Remove Baseline occurrences (early + final).');
testCase.verifyEqual(t.spec(blIdx(1)).params.timerange, [-500, 1500], ...
    'First baseline must be full-epoch demean.');
testCase.verifyEqual(t.spec(blIdx(2)).params.timerange, [-300, -100], ...
    'Final baseline must be [-300, -100] ms (Wu 2018 §2.2.3).');
end

function test_aaratepRerefBeforeEarlyEyeICA(testCase)
% Upstream c_TMSEEG_Preprocess_AARATEPPipeline.m line 244:
% pop_reref(EEG, []) is called BEFORE the early eye-IC ICA so the early
% decomposition is on average-referenced data.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'AARATEP'));
rerefIdx = find(strcmp(t.steps, 'Re-Reference'));
icaIdx   = find(strcmp(t.steps, 'Run ICA'));
testCase.verifyNotEmpty(rerefIdx, 'AARATEP must include Re-Reference.');
testCase.verifyNotEmpty(icaIdx,   'AARATEP must include Run ICA.');
testCase.verifyLessThan(rerefIdx(1), icaIdx(1), ...
    'AARATEP early Re-Reference must precede the first ICA round.');
end

function test_aaratepEarlyEyeICAFlagsEyeOnly(testCase)
% Upstream line 254: muscleComponentThreshold = NaN,
% brainComponentThreshold = NaN, otherComponentThreshold = NaN. Only Eye
% is flagged in the early pass. The registry default for Muscle is
% [0.9, 1] so the template MUST override it to [NaN, NaN].
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'AARATEP'));
flagIdx = find(strcmp(t.steps, 'Flag ICA Components for Rejection'));
p = t.spec(flagIdx(1)).params;
testCase.verifyTrue(all(isnan(p.Muscle)), ...
    ['AARATEP early Flag step must have Muscle = [NaN, NaN] (upstream ' ...
     'muscleComponentThreshold = NaN).']);
testCase.verifyTrue(all(isnan(p.Brain)), ...
    'AARATEP early Flag step must have Brain = [NaN, NaN].');
testCase.verifyTrue(all(isnan(p.Heart)), ...
    'AARATEP early Flag step must have Heart = [NaN, NaN].');
testCase.verifyEqual(p.Eye, [0.9, 1], ...
    'AARATEP early Flag step must reject ICs with Eye prob >= 0.9.');
end

function test_aaratepMuscleFlagAfterICLabelFlag(testCase)
% pop_icflag in "Flag ICA Components for Rejection" ASSIGNS to
% EEG.reject.gcompreject (replacing prior flags). The AARATEP muscle
% classifier ORs into gcompreject. Therefore Flag (ICLabel) must run
% BEFORE Flag (AARATEP Muscle), or pop_icflag clobbers the muscle flags.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'AARATEP'));
labelIdx = find(strcmp(t.steps, 'Label ICA Components'), 1, 'last');
flagIdx  = find(strcmp(t.steps, 'Flag ICA Components for Rejection'), 1, 'last');
muscleIdx = find(strcmp(t.steps, 'Flag ICA Components (AARATEP Muscle)'), 1);
testCase.verifyNotEmpty(muscleIdx, 'AARATEP must include the muscle classifier step.');
testCase.verifyLessThan(labelIdx, flagIdx, ...
    'Label ICA Components must precede Flag ICA Components for Rejection.');
testCase.verifyLessThan(flagIdx, muscleIdx, ...
    ['AARATEP muscle flag must come AFTER the ICLabel flag step - ' ...
     'pop_icflag replaces gcompreject and would clobber muscle flags.']);
end

function test_aaratepDoesNotUseCleanline(testCase)
% Upstream uses c_EEG_filter_butterworth bandstop, not CleanLine. Our
% substitute is Frequency Filter (TESA) bandstop 58-62 Hz.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'AARATEP'));
testCase.verifyFalse(any(strcmp(t.steps, 'Frequency Filter (CleanLine)')), ...
    'AARATEP must NOT use CleanLine; use Frequency Filter (TESA) bandstop instead.');
notchIdx = find(strcmp(t.steps, 'Frequency Filter (TESA)'), 1);
testCase.verifyNotEmpty(notchIdx, 'AARATEP must include a Frequency Filter (TESA) step.');
end

function test_aaratepFinalRereferenceIsAverage(testCase)
% Upstream line 498: pop_reref(EEG, []) (average reference) before save.
templates = loadTemplates(testCase);
t = templates(contains({templates.name}, 'AARATEP'));
rerefIdx = find(strcmp(t.steps, 'Re-Reference'));
testCase.verifyEqual(t.spec(rerefIdx(end)).params.ref, '[]', ...
    'AARATEP final Re-Reference must be average ([]).');
end
