function tests = test_templateCitation
% TEST_TEMPLATECITATION  Citation parity: every built-in template must
%   have a non-empty citation registered in templateCitation.m, so the
%   batch log includes attribution for each pipeline.
%
%   Run: runtests('tests/unit/test_templateCitation')
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));
end

function r = repoRoot()
r = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

function templates = loadTemplates(testCase)
r = repoRoot();
matFiles = dir(fullfile(r, 'src', 'templates', '*.mat'));
testCase.verifyFalse(isempty(matFiles), 'No template .mat files found.');
templates = struct('name', {}, 'steps', {});
for i = 1:numel(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    templates(i).name  = data.pipelineName;
    templates(i).steps = {data.spec.name};
end
end

% ── every built-in template gets a citation ───────────────────────────────

function test_everyBuiltInTemplateHasCitation(testCase)
templates = loadTemplates(testCase);
for i = 1:numel(templates)
    c = templateCitation(templates(i).name);
    testCase.verifyNotEmpty(c.reference, sprintf( ...
        ['Template "%s" has no citation registered in ' ...
         'templateCitation.m. Add an entry so the batch log ' ...
         'attributes the pipeline.'], templates(i).name));
    testCase.verifyNotEmpty(c.doi, sprintf( ...
        ['Template "%s" has no DOI registered in templateCitation.m.'], ...
        templates(i).name));
end
end

% ── specific citation content sanity checks ───────────────────────────────

function test_tesaCitationPresent(testCase)
c = templateCitation('TMS-EEG / TEP (TESA)');
testCase.verifyEqual(c.doi, '10.1016/j.neuroimage.2017.06.014');
testCase.verifyTrue(contains(c.reference, 'Rogasch'));
end

function test_tesaQualityGatesCitationPresent(testCase)
c = templateCitation('TMS-EEG / TEP (TESA + Quality Gates)');
testCase.verifyEqual(c.doi, '10.1016/j.neuroimage.2017.06.014');
testCase.verifyTrue(contains(c.reference, 'Rogasch'));
end

function test_artistCitationPresent(testCase)
c = templateCitation('TMS-EEG / ARTIST');
testCase.verifyEqual(c.doi, '10.1002/hbm.23938');
testCase.verifyTrue(contains(c.reference, 'Wu'));
end

function test_artistAcknowledgesTESACompselect(testCase)
% ARTIST template uses TESA pop_tesa_compselect as FLDA fallback -
% the citation note must steer users to also cite Rogasch 2017.
c = templateCitation('TMS-EEG / ARTIST');
testCase.verifyTrue(contains(c.notes, 'Rogasch'), ...
    ['ARTIST citation note must direct users to also cite Rogasch 2017 ' ...
     'because pop_tesa_compselect substitutes for the unreleased FLDA classifier.']);
end

function test_aaratepCitationPresent(testCase)
c = templateCitation('TMS-EEG / AARATEP');
testCase.verifyEqual(c.doi, '10.1109/NER49283.2021.9441147');
testCase.verifyTrue(contains(c.reference, 'Cline'));
end

function test_aaratepAcknowledgesSOUNDAndTESA(testCase)
% AARATEP template runs SOUND (Mutanen 2018) via TESA's pop_tesa_sound
% (Rogasch 2017). Both must be co-cited.
c = templateCitation('TMS-EEG / AARATEP');
testCase.verifyTrue(contains(c.notes, 'Mutanen'), ...
    ['AARATEP citation note must mention Mutanen 2018 (SOUND) - the ' ...
     'AARATEP template runs Remove Recording Noise (SOUND).']);
testCase.verifyTrue(contains(c.notes, 'Rogasch'), ...
    ['AARATEP citation note must mention Rogasch 2017 (TESA) - the ' ...
     'vendored AARATEP helpers depend on TESA''s pop_tesa_sound.']);
end

function test_tesaNotesMentionSOUNDForOptionalUse(testCase)
% TESA citation notes should remind users about SOUND co-citation when
% Remove Recording Noise (SOUND) is in the pipeline. The TESA TEP
% template itself does not use SOUND, but the citation function should
% mention the contingent attribution rule.
c = templateCitation('TMS-EEG / TEP (TESA)');
testCase.verifyTrue(contains(c.notes, 'Mutanen') || contains(c.notes, 'SOUND'), ...
    ['TESA citation notes should remind users about Mutanen 2018 (SOUND) ' ...
     'co-citation when Remove Recording Noise (SOUND) is in the pipeline.']);
end
