
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_describePipeline
% TEST_DESCRIBEPIPELINE  Unit tests for the pipeline description exporter.
%
%   No EEGLAB required - describePipeline uses only stepRegistry defaults.
%
%   Run: runtests('tests/unit/test_describePipeline')
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

function spec = makeSpec()
% Two steps: one left at defaults, one with an override.
reg = stepRegistry();
s1 = makePipelineStep('Load Data', reg);
s2 = makePipelineStep('Re-Sample', reg);
s2.params.freq = 500;   % override the default (1000)
spec = [s1, s2];
end

% ── output contract ─────────────────────────────────────────────────────────

function test_returnsCharAndStruct(testCase)
[text, info] = describePipeline(makeSpec(), 'Quiet', true);
testCase.verifyClass(text, 'char');
testCase.verifyNotEmpty(text);
testCase.verifyEqual(info.nSteps, 2);
end

function test_listsEveryStepName(testCase)
spec = makeSpec();
text = describePipeline(spec, 'Quiet', true);
for i = 1:numel(spec)
    testCase.verifyTrue(contains(text, spec(i).name), ...
        sprintf('description must mention step "%s"', spec(i).name));
end
end

function test_showsOverriddenParam(testCase)
text = describePipeline(makeSpec(), 'Quiet', true);
testCase.verifyTrue(contains(text, 'freq = 500'), ...
    'an overridden parameter must be shown with its value');
end

function test_omitsDefaultParams(testCase)
% Re-Sample's fc/df are left at defaults and must NOT appear as overrides.
[~, info] = describePipeline(makeSpec(), 'Quiet', true);
resampleStep = info.steps(strcmp({info.steps.name}, 'Re-Sample'));
ovKeys = fieldnames(resampleStep.overrides);
testCase.verifyTrue(ismember('freq', ovKeys), 'freq override must be recorded');
testCase.verifyFalse(ismember('fc', ovKeys), 'default fc must not be recorded as an override');
end

function test_emptySpecIsHandled(testCase)
[text, info] = describePipeline(struct('name', {}, 'params', {}), 'Quiet', true);
testCase.verifyEqual(info.nSteps, 0);
testCase.verifyTrue(contains(text, 'empty pipeline'), ...
    'empty spec must produce a clear "empty pipeline" note');
end

function test_versionInHeader(testCase)
text = describePipeline(makeSpec(), 'Quiet', true);
testCase.verifyTrue(contains(text, nestappVersion()), ...
    'description header must include the nestapp version');
end
