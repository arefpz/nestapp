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

function test_threeTemplatesExist(testCase)
templates = loadTemplates(testCase);
testCase.verifyEqual(numel(templates), 3, 'Expected exactly 3 templates');
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
