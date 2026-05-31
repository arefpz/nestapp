function tests = test_pipelineReport
% TEST_PIPELINEREPORT  Unit tests for initPipelineReport and exportReport.
%
%   Covers: struct shape, field defaults, exportReport formatting,
%   methods-summary branch coverage (the main regression target for the
%   channel-count bug fixed in Apr 2026).
%
%   Run: runtests('tests/unit/test_pipelineReport')
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

% ── initPipelineReport — struct shape ────────────────────────────────────

function test_initReportHasRequiredFields(testCase)
report = initPipelineReport('subject01.set');
testCase.verifyTrue(isfield(report, 'inputFile'),   'Missing inputFile');
testCase.verifyTrue(isfield(report, 'processedAt'), 'Missing processedAt');
testCase.verifyTrue(isfield(report, 'steps'),       'Missing steps');
testCase.verifyTrue(isfield(report, 'channels'),    'Missing channels');
testCase.verifyTrue(isfield(report, 'trials'),      'Missing trials');
testCase.verifyTrue(isfield(report, 'ica'),         'Missing ica');
end

function test_initReportChannelFields(testCase)
report = initPipelineReport('subject01.set');
ch = report.channels;
testCase.verifyTrue(isfield(ch, 'original'),       'Missing channels.original');
testCase.verifyTrue(isfield(ch, 'nRejected'),      'Missing channels.nRejected');
testCase.verifyTrue(isfield(ch, 'nInterpolated'),  'Missing channels.nInterpolated');
testCase.verifyTrue(isfield(ch, 'final'),          'Missing channels.final');
testCase.verifyEqual(ch.original,      0);
testCase.verifyEqual(ch.nRejected,     0);
testCase.verifyEqual(ch.nInterpolated, 0);
testCase.verifyEqual(ch.final,         0);
end

function test_initReportTrialFields(testCase)
report = initPipelineReport('subject01.set');
tr = report.trials;
testCase.verifyTrue(isfield(tr, 'original'), 'Missing trials.original');
testCase.verifyTrue(isfield(tr, 'rejected'), 'Missing trials.rejected');
testCase.verifyTrue(isfield(tr, 'final'),    'Missing trials.final');
testCase.verifyEqual(tr.original, 0);
testCase.verifyEqual(tr.rejected, 0);
testCase.verifyEqual(tr.final,    0);
end

function test_initReportIcaFields(testCase)
report = initPipelineReport('subject01.set');
testCase.verifyTrue(isfield(report.ica, 'nComponents'), 'Missing ica.nComponents');
testCase.verifyTrue(isfield(report.ica, 'nRejected'),   'Missing ica.nRejected');
testCase.verifyTrue(isfield(report.ica, 'varRemoved'),  'Missing ica.varRemoved');
testCase.verifyEqual(report.ica.nComponents, 0);
testCase.verifyEqual(report.ica.nRejected,   0);
testCase.verifyTrue(isnan(report.ica.varRemoved), 'varRemoved must init as NaN');
end

function test_initReportStepsIsEmpty(testCase)
report = initPipelineReport('subject01.set');
testCase.verifyEmpty(report.steps, 'steps should be empty at initialisation');
end

function test_initReportStoresFilename(testCase)
report = initPipelineReport('my_eeg_data.set');
testCase.verifyEqual(report.inputFile, 'my_eeg_data.set');
end

function test_initReportIcaCategoriesPresent(testCase)
report = initPipelineReport('test.set');
cats = report.ica.categories;
testCase.verifyTrue(isfield(cats, 'names'),    'Missing ica.categories.names');
testCase.verifyTrue(isfield(cats, 'nRemoved'), 'Missing ica.categories.nRemoved');
testCase.verifyEqual(numel(cats.names), 7, 'Should have 7 ICLabel categories');
testCase.verifyEqual(cats.nRemoved, zeros(1,7), 'nRemoved should start at zero');
end

% ── exportReport — basic formatting ──────────────────────────────────────

function test_exportReportReturnsText(testCase)
report = initPipelineReport('subject01.set');
report.channels.original = 64;
report.channels.final    = 64;
summaryText = exportReport(report, '');
testCase.verifyFalse(isempty(summaryText), 'summaryText must not be empty');
testCase.verifyTrue(ischar(summaryText) || isstring(summaryText), 'must be text');
end

function test_exportReportContainsFilename(testCase)
report = initPipelineReport('subject01.set');
summaryText = exportReport(report, '');
testCase.verifyTrue(contains(summaryText, 'subject01.set'), ...
    'summaryText must contain the input filename');
end

function test_exportReportChannelSummary(testCase)
report = initPipelineReport('test.set');
report.channels.original  = 64;
report.channels.nRejected = 3;
report.channels.final     = 61;
summaryText = exportReport(report, '');
testCase.verifyTrue(contains(summaryText, '64'), 'Must mention original channel count');
testCase.verifyTrue(contains(summaryText, '61'), 'Must mention final channel count');
end

function test_exportReportTrialSummary(testCase)
report = initPipelineReport('test.set');
report.trials.original = 80;
report.trials.rejected = 5;
report.trials.final    = 75;
summaryText = exportReport(report, '');
testCase.verifyTrue(contains(summaryText, '80'), 'Must mention original trial count');
testCase.verifyTrue(contains(summaryText, '75'), 'Must mention final trial count');
end

function test_exportReportSkipsSaveWhenEmptyDir(testCase)
report = initPipelineReport('test.set');
[summaryText, matPath] = exportReport(report, '');
testCase.verifyFalse(isempty(summaryText), 'summaryText must be non-empty');
% matPath may be empty if no valid output dir — that is acceptable
testCase.verifyTrue(ischar(matPath) || isstring(matPath), 'matPath must be text type');
end

% ── exportReport — METHODS SUMMARY branch coverage ───────────────────────
% These tests pin the Apr-2026 fix: methods text now uses nRejected and
% nInterpolated directly rather than computing finCh/origCh fractions.

function test_methodsNone(testCase)
% No rejection, no interpolation → "none rejected"
report = initPipelineReport('test.set');
report.channels.original      = 32;
report.channels.nRejected     = 0;
report.channels.nInterpolated = 0;
report.channels.final         = 32;
txt = exportReport(report, '');
testCase.verifyTrue(contains(txt, 'none rejected'), ...
    'Methods should say "none rejected" when no channels were rejected');
end

function test_methodsRejectionOnly(testCase)
% 3 channels rejected, none interpolated
report = initPipelineReport('test.set');
report.channels.original      = 64;
report.channels.nRejected     = 3;
report.channels.nInterpolated = 0;
report.channels.final         = 61;
txt = exportReport(report, '');
testCase.verifyTrue(contains(txt, '3'), 'Should mention 3 rejected channels');
testCase.verifyTrue(contains(txt, 'identified as bad'), ...
    'Should use "identified as bad" phrasing');
% Must NOT claim all channels retained
testCase.verifyFalse(contains(txt, '100%'), ...
    'Must not claim 100% retained when channels were rejected');
end

function test_methodsRejectionAndInterpolation(testCase)
% 5 rejected, 3 interpolated → "5 removed; 3 interpolated"
report = initPipelineReport('test.set');
report.channels.original      = 64;
report.channels.nRejected     = 5;
report.channels.nInterpolated = 3;
report.channels.final         = 62;
txt = exportReport(report, '');
testCase.verifyTrue(contains(txt, '5'), 'Should mention 5 rejected');
testCase.verifyTrue(contains(txt, '3'), 'Should mention 3 interpolated');
testCase.verifyTrue(contains(txt, 'interpolated'), 'Should use "interpolated"');
end

function test_methodsInterpolationOnly(testCase)
% 0 rejected but 4 interpolated (unusual but valid)
report = initPipelineReport('test.set');
report.channels.original      = 64;
report.channels.nRejected     = 0;
report.channels.nInterpolated = 4;
report.channels.final         = 64;
txt = exportReport(report, '');
testCase.verifyTrue(contains(txt, 'interpolated'), 'Should mention interpolation');
testCase.verifyTrue(contains(txt, '4'), 'Should mention 4 interpolated channels');
end

function test_methodsTrialRetentionPercent(testCase)
% Trial sentence should include a percentage
report = initPipelineReport('test.set');
report.trials.original = 80;
report.trials.rejected = 8;
report.trials.final    = 72;
txt = exportReport(report, '');
testCase.verifyTrue(contains(txt, '90') || contains(txt, '72') || contains(txt, '80'), ...
    'Methods trial sentence should mention counts or percentage');
end

function test_methodsICANoRejection(testCase)
% ICA ran but nothing removed → "none rejected"
report = initPipelineReport('test.set');
report.ica.nComponents = 30;
report.ica.nRejected   = 0;
report.ica.nKept       = 30;
txt = exportReport(report, '');
testCase.verifyTrue(contains(txt, '30'), 'Should mention 30 components');
testCase.verifyTrue(contains(txt, 'none rejected'), ...
    'Should say "none rejected" when no ICA components removed');
end

function test_methodsICAWithVariance(testCase)
% Single round, variance available
report = initPipelineReport('test.set');
report.ica.nComponents = 30;
report.ica.nRejected   = 4;
report.ica.nKept       = 26;
report.ica.varRemoved  = 18.5;
report.ica.varMin      = 2.3;
report.ica.varMax      = 7.1;
txt = exportReport(report, '');
testCase.verifyTrue(contains(txt, '18'), 'Should mention variance removed (~18%)');
end

function test_methodsICAWithCategories(testCase)
% When ICLabel category data is present the report must show the breakdown.
report = initPipelineReport('test.set');
report.ica.nComponents = 25;
report.ica.nRejected   = 5;
report.ica.nKept       = 20;
% Set two non-zero category counts (Muscle and Eye)
report.ica.categories.nRemoved(2) = 3;   % Muscle
report.ica.categories.nRemoved(3) = 2;   % Eye
report.ica.categories.varShare(2) = 12.0;
report.ica.categories.varShare(3) = 6.5;
txt = exportReport(report, '');
testCase.verifyTrue(contains(txt, 'Muscle') || contains(lower(txt), 'muscle'), ...
    'Report must list Muscle category when components were rejected in that category');
testCase.verifyTrue(contains(txt, 'Eye') || contains(lower(txt), 'eye'), ...
    'Report must list Eye category when components were rejected in that category');
end

function test_methodsMultiRoundNoVariance(testCase)
% Multi-round TESA ICA: top-level variance fields are NaN (variance is not
% additive across ICA bases); only per-round detail should carry variance.
report = initPipelineReport('test.set');
report.ica.nComponents = 30;
report.ica.nRejected   = 8;
report.ica.nKept       = 22;
% Simulate two rounds — varRemoved stays NaN at the top level
rnd1.nComponents = 30; rnd1.nRejected = 4; rnd1.varRemoved = 15.0;
rnd1.varMin = 2.0; rnd1.varMax = 6.0;
rnd1.categories = report.ica.categories;
rnd2.nComponents = 26; rnd2.nRejected = 4; rnd2.varRemoved = NaN;
rnd2.varMin = NaN; rnd2.varMax = NaN;
rnd2.categories = report.ica.categories;
report.ica.rounds = {rnd1, rnd2};
txt = exportReport(report, '');
% Must mention multiple rounds
testCase.verifyTrue(contains(txt, '2') && contains(lower(txt), 'round'), ...
    'Multi-round report must mention round count');
% Top-level summary must NOT assert a single combined variance figure
% (because variance across different ICA bases is not additive)
testCase.verifyFalse(contains(txt, sprintf('%.1f%% ICA variance', ...
    rnd1.varRemoved + rnd2.varRemoved)), ...
    'Must not sum variance across rounds in top-level summary');
end
