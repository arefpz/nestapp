
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_icaReporting
% TEST_ICAREPORTING  Unit tests for the unified ICA reporting helpers.
%
%   Covers openICARound / addICARemoval (a round == an ICA decomposition;
%   removals fold into the current round) and markICClass (per-IC category
%   tagging by custom classifiers), which together give all ICA paths
%   TESA-style parity: a correct round count, detected count, and a
%   per-category breakdown - even for a round that removes nothing.
%
%   Run: runtests('tests/unit/test_icaReporting')
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(r);
addpath(genpath(fullfile(r, 'src')));
end

% -- helpers --------------------------------------------------------------

function rmv = removal(nComp, names, counts, varargin)
% Build a removal struct for addICARemoval.
% varargin{1}: [varRemoved varMin varMax]; default NaNs.
v = [NaN NaN NaN];
if nargin >= 4 && ~isempty(varargin{1}); v = varargin{1}; end
rmv = struct('nComponents', nComp, 'nRejected', sum(counts), ...
    'varRemoved', v(1), 'varMin', v(2), 'varMax', v(3));
rmv.categories = struct('names', {names}, 'nRemoved', counts, ...
    'varShare', zeros(1, numel(names)));
end

function n = catCount(report, name)
j = find(strcmp(report.ica.categories.names, name), 1);
if isempty(j); n = 0; else; n = report.ica.categories.nRemoved(j); end
end

% -- openICARound / addICARemoval: a round == an ICA decomposition --------

function test_singleRoundAccumulates(testCase)
report = initPipelineReport('x.set');
report = openICARound(report, 20);
report = addICARemoval(report, removal(20, {'Muscle','Eye'}, [3 2]));
testCase.verifyEqual(numel(report.ica.rounds), 1);
testCase.verifyEqual(report.ica.nComponents, 20);
testCase.verifyEqual(report.ica.nRejected, 5);
testCase.verifyEqual(report.ica.nKept, 15);
testCase.verifyEqual(catCount(report, 'Muscle'), 3);
testCase.verifyEqual(catCount(report, 'Eye'), 2);
end

function test_eachDecompositionIsItsOwnRound(testCase)
report = initPipelineReport('x.set');
report = openICARound(report, 20);
report = addICARemoval(report, removal(20, {'Muscle'}, 2));
report = openICARound(report, 18);
report = addICARemoval(report, removal(18, {'Eye'}, 3));
testCase.verifyEqual(numel(report.ica.rounds), 2);
testCase.verifyEqual(report.ica.nRejected, 5);
testCase.verifyEqual(catCount(report, 'Muscle'), 2);
testCase.verifyEqual(catCount(report, 'Eye'), 3);
end

function test_zeroRemovalRoundStillCounts(testCase)
% Regression (AARATEP): a 2nd ICA decomposition that removes nothing must
% still be a recorded round, and nKept must stay correct (was negative when
% nComponents got overwritten by the last Run ICA).
report = initPipelineReport('x.set');
report = openICARound(report, 15);
report = addICARemoval(report, removal(15, {'Muscle'}, 12));
report = openICARound(report, 3);          % 2nd decomposition, no removal
testCase.verifyEqual(numel(report.ica.rounds), 2, 'A zero-removal round must still count');
testCase.verifyEqual(report.ica.nComponents, 15, 'nComponents stays the first decomposition');
testCase.verifyEqual(report.ica.nRejected, 12);
testCase.verifyEqual(report.ica.nKept, 3, 'nKept must be the final survivors, never negative');
testCase.verifyEqual(report.ica.rounds{2}.nRejected, 0);
end

function test_categoriesMergeAcrossRounds(testCase)
report = initPipelineReport('x.set');
report = openICARound(report, 30);
report = addICARemoval(report, removal(30, {'Muscle'}, 2));
report = openICARound(report, 28);
report = addICARemoval(report, removal(28, {'Muscle','Eye'}, [1 4]));
testCase.verifyEqual(catCount(report, 'Muscle'), 3);
testCase.verifyEqual(catCount(report, 'Eye'), 4);
end

function test_newCategoryNamesAreAppended(testCase)
% A non-ICLabel scheme (TESA / classifier names) appends rather than dropping,
% so schemes coexist and an unused default category never shows up.
report = initPipelineReport('x.set');
report = openICARound(report, 30);
report = addICARemoval(report, removal(30, {'TMS Muscle','Decay'}, [2 1]));
testCase.verifyEqual(catCount(report, 'TMS Muscle'), 2);
testCase.verifyEqual(catCount(report, 'Decay'), 1);
testCase.verifyEqual(catCount(report, 'Brain'), 0);
end

function test_topVarianceFromFirstRoundOnly(testCase)
% Variance is not additive across ICA bases: only the first round sets the
% top-level variance.
report = initPipelineReport('x.set');
report = openICARound(report, 20);
report = addICARemoval(report, removal(20, {'Muscle'}, 2, [12.5 1 6]));
report = openICARound(report, 18);
report = addICARemoval(report, removal(18, {'Eye'}, 1, [99.9 1 99]));
testCase.verifyEqual(report.ica.varRemoved, 12.5);
end

function test_removalWithoutPriorRoundOpensOne(testCase)
% Defensive: a removal with no preceding Run ICA opens a round so totals stay
% sane rather than erroring.
report = initPipelineReport('x.set');
report = addICARemoval(report, removal(10, {'Manual'}, 4));
testCase.verifyEqual(numel(report.ica.rounds), 1);
testCase.verifyEqual(report.ica.nRejected, 4);
testCase.verifyEqual(report.ica.nKept, 6);
end

% -- markICClass ----------------------------------------------------------

function test_markICClassTagsFlaggedComps(testCase)
EEG = struct('etc', struct());
EEG = markICClass(EEG, [false true false true], 'Muscle');
testCase.verifyEqual(EEG.etc.nestappICClass, {'', 'Muscle', '', 'Muscle'});
end

function test_markICClassInitialisesWhenMissing(testCase)
EEG = struct();   % no etc field at all
EEG = markICClass(EEG, logical([1 0 0]), 'Decay');
testCase.verifyEqual(EEG.etc.nestappICClass, {'Decay', '', ''});
end

function test_markICClassLaterLabelWins(testCase)
EEG = struct('etc', struct());
EEG = markICClass(EEG, logical([1 1 0]), 'Muscle');
EEG = markICClass(EEG, logical([0 1 0]), 'Decay');   % IC 2 reclassified
testCase.verifyEqual(EEG.etc.nestappICClass, {'Muscle', 'Decay', ''});
end

% -- icaCategoriesFromFlags (classifier > ICLabel > Manual priority) ------

function test_categoriesFromClassifierLabels(testCase)
% Custom-classifier labels (from markICClass) drive the categories.
p = struct('classLabels', {{'Muscle', '', 'Decay', 'Muscle'}});
cats = icaCategoriesFromFlags(logical([1 0 1 1]), p);
testCase.verifyEqual(catN(cats, 'Muscle'), 2);
testCase.verifyEqual(catN(cats, 'Decay'), 1);
end

function test_categoriesFallBackToICLabelArgmax(testCase)
% With no classifier label, the ICLabel argmax category is used.
probs = [0.9 0.1 0 0 0 0 0;    % comp1 -> Brain (not removed)
         0.1 0.8 0.1 0 0 0 0;  % comp2 -> Muscle
         0   0.1 0.9 0 0 0 0]; % comp3 -> Eye
p = struct('iclabelProbs', probs);
cats = icaCategoriesFromFlags(logical([0 1 1]), p);
testCase.verifyEqual(catN(cats, 'Muscle'), 1);
testCase.verifyEqual(catN(cats, 'Eye'), 1);
testCase.verifyEqual(catN(cats, 'Brain'), 0);   % comp1 not removed
end

function test_categoriesClassifierBeatsICLabel(testCase)
% Classifier label wins over ICLabel for the same component.
probs = [0 0 1 0 0 0 0];   % ICLabel would say Eye
p = struct('classLabels', {{'Muscle'}}, 'iclabelProbs', probs);
cats = icaCategoriesFromFlags(true, p);
testCase.verifyEqual(catN(cats, 'Muscle'), 1);
testCase.verifyEqual(catN(cats, 'Eye'), 0);
end

function test_categoriesDefaultToManual(testCase)
% No classifier and no ICLabel -> 'Manual'.
cats = icaCategoriesFromFlags(logical([1 0 1]), struct());
testCase.verifyEqual(catN(cats, 'Manual'), 2);
end

function test_categoriesVarShareAccumulates(testCase)
% varShare sums the per-component variance within each category.
p = struct('classLabels', {{'Muscle', 'Muscle'}}, 'compVarPct', [4 6]);
cats = icaCategoriesFromFlags(logical([1 1]), p);
testCase.verifyEqual(cats.varShare(strcmp(cats.names, 'Muscle')), 10);
end

function n = catN(cats, name)
j = find(strcmp(cats.names, name), 1);
if isempty(j); n = 0; else; n = cats.nRemoved(j); end
end
