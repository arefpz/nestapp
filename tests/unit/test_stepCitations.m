
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_stepCitations
% TEST_STEPCITATIONS  Citations are derived from the pipeline steps used.
%
%   stepCitations(stepNames) returns one reference per method whose trigger
%   steps appear in stepNames - so any constellation of steps cites exactly
%   the methods it used, with no spurious "if you also used X" entries.
%
%   Run: runtests('tests/unit/test_stepCitations')
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(r);
addpath(fullfile(r, 'src'));
end

function refs = refsFor(stepNames)
c = stepCitations(stepNames);
refs = {c.reference};
end

function tf = anyContains(refs, needle)
tf = any(cellfun(@(r) contains(r, needle), refs));
end

% -- output shape ---------------------------------------------------------

function test_emptyInputReturnsEmpty(testCase)
c = stepCitations({});
testCase.verifyEmpty(c);
testCase.verifyTrue(isfield(c, 'reference') && isfield(c, 'doi'), ...
    'Empty result must still have the reference/doi fields');
end

function test_fieldsAreReferenceAndDoiOnly(testCase)
c = stepCitations({'Load Data'});
testCase.verifyEqual(sort(fieldnames(c)), sort({'reference'; 'doi'}), ...
    'Output must expose only reference and doi (the match predicate is internal)');
end

% -- per-method triggers --------------------------------------------------

function test_loadDataCitesEEGLAB(testCase)
testCase.verifyTrue(anyContains(refsFor({'Load Data'}), 'Delorme'));
end

function test_tesaStepCitesRogasch(testCase)
testCase.verifyTrue(anyContains(refsFor({'Remove ICA Components (TESA)'}), 'Rogasch'));
end

function test_soundStepCitesMutanen(testCase)
testCase.verifyTrue(anyContains(refsFor({'Remove Recording Noise (SOUND)'}), 'Mutanen'));
end

function test_artistStepCitesWu(testCase)
testCase.verifyTrue(anyContains(refsFor({'Reject Bad Trials (ARTIST)'}), 'Wu W.'));
end

function test_aaratepStepCitesCline(testCase)
testCase.verifyTrue(anyContains(refsFor({'Remove Decay Artifact'}), 'Cline'));
end

function test_iclabelStepCitesPionTonachini(testCase)
testCase.verifyTrue(anyContains(refsFor({'Label ICA Components'}), 'Pion-Tonachini'));
end

% -- the point: union by constellation, no spurious entries ---------------

function test_combinedStepsCiteBothMethods(testCase)
% The user's example: TESA compselect + ARTIST trial rejection -> both refs.
refs = refsFor({'Remove ICA Components (TESA)', 'Reject Bad Trials (ARTIST)'});
testCase.verifyTrue(anyContains(refs, 'Rogasch'), 'TESA must be cited');
testCase.verifyTrue(anyContains(refs, 'Wu W.'),   'ARTIST must be cited');
end

function test_noSpuriousSoundWhenSoundNotUsed(testCase)
% A TESA pipeline without a SOUND step must NOT cite Mutanen/SOUND.
refs = refsFor({'Load Data', 'Remove ICA Components (TESA)', 'Save New Set'});
testCase.verifyFalse(anyContains(refs, 'Mutanen'), ...
    'SOUND/Mutanen must not be cited without a SOUND step');
end

function test_duplicateTriggersCiteMethodOnce(testCase)
% Several (TESA) steps must yield a single TESA citation, not one per step.
refs = refsFor({'Frequency Filter (TESA)', 'Run TESA ICA', ...
                'Remove ICA Components (TESA)'});
nRogasch = sum(cellfun(@(r) contains(r, 'Rogasch'), refs));
testCase.verifyEqual(nRogasch, 1, 'TESA must be cited exactly once across multiple TESA steps');
end

function test_uncitedStepsYieldNothing(testCase)
testCase.verifyEmpty(stepCitations({'Remove Baseline', 'Re-Reference', 'Epoching'}));
end

% -- every built-in template cites at least one method --------------------

function test_everyBuiltInTemplateGetsACitation(testCase)
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
matFiles = dir(fullfile(r, 'src', 'templates', '*.mat'));
testCase.verifyFalse(isempty(matFiles), 'No template .mat files found.');
for i = 1:numel(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    c = stepCitations({data.spec.name});
    testCase.verifyNotEmpty(c, sprintf( ...
        'Template "%s" produced no citations from its steps.', data.pipelineName));
end
end
