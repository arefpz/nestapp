
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_icaReporting
% TEST_ICAREPORTING  Unit tests for the unified ICA reporting helpers.
%
%   Covers recordICARound (the shared per-round accumulator every ICA removal
%   path feeds) and markICClass (per-IC category tagging by custom
%   classifiers), which together give all ICA paths TESA-style parity:
%   detected count, N rounds, and a per-category breakdown.
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

function rnd = makeRound(roundNum, nComp, names, counts, varargin)
% varargin{1}: varShare vector (default zeros); varargin{2}: varRemoved scalar.
varShare = zeros(1, numel(names));
if nargin >= 5 && ~isempty(varargin{1}); varShare = varargin{1}; end
varRemoved = NaN;
if nargin >= 6; varRemoved = varargin{2}; end
rnd = struct('roundNum', roundNum, 'nComponents', nComp, ...
    'nRejected', sum(counts), 'varRemoved', varRemoved, ...
    'varMin', NaN, 'varMax', NaN);
rnd.categories = struct('names', {names}, 'nRemoved', counts, 'varShare', varShare);
end

function n = catCount(report, name)
j = find(strcmp(report.ica.categories.names, name), 1);
if isempty(j); n = 0; else; n = report.ica.categories.nRemoved(j); end
end

% -- recordICARound -------------------------------------------------------

function test_singleRoundAccumulates(testCase)
report = initPipelineReport('x.set');
report.ica.nComponents = 20;
report = recordICARound(report, makeRound(1, 20, {'Muscle','Eye'}, [3 2]));
testCase.verifyEqual(numel(report.ica.rounds), 1);
testCase.verifyEqual(report.ica.nRejected, 5);
testCase.verifyEqual(report.ica.nKept, 15);
testCase.verifyEqual(catCount(report, 'Muscle'), 3);
testCase.verifyEqual(catCount(report, 'Eye'), 2);
end

function test_twoRoundsAreMultiRoundAndSum(testCase)
report = initPipelineReport('x.set');
report.ica.nComponents = 20;
report = recordICARound(report, makeRound(1, 20, {'Muscle'}, 2));
report = recordICARound(report, makeRound(2, 18, {'Eye'}, 3));
testCase.verifyEqual(numel(report.ica.rounds), 2);
testCase.verifyEqual(report.ica.nRejected, 5);
testCase.verifyEqual(catCount(report, 'Muscle'), 2);
testCase.verifyEqual(catCount(report, 'Eye'), 3);
end

function test_categoriesMergeByName(testCase)
% The same category appearing in two rounds accumulates into one entry.
report = initPipelineReport('x.set');
report.ica.nComponents = 30;
report = recordICARound(report, makeRound(1, 30, {'Muscle'}, 2));
report = recordICARound(report, makeRound(2, 28, {'Muscle','Eye'}, [1 4]));
testCase.verifyEqual(catCount(report, 'Muscle'), 3);
testCase.verifyEqual(catCount(report, 'Eye'), 4);
end

function test_newCategoryNamesAreAppended(testCase)
% A round using a non-ICLabel scheme (e.g. TESA / classifier names) appends
% its categories rather than being dropped, so all schemes coexist.
report = initPipelineReport('x.set');
report.ica.nComponents = 30;
report = recordICARound(report, makeRound(1, 30, {'TMS Muscle','Decay'}, [2 1]));
testCase.verifyEqual(catCount(report, 'TMS Muscle'), 2);
testCase.verifyEqual(catCount(report, 'Decay'), 1);
% Default ICLabel categories that flagged nothing stay zero, not negative.
testCase.verifyEqual(catCount(report, 'Brain'), 0);
end

function test_topVarianceFromFirstRoundOnly(testCase)
% Variance across different ICA bases is not additive: only round 1 sets the
% top-level variance fields.
report = initPipelineReport('x.set');
report.ica.nComponents = 20;
report = recordICARound(report, makeRound(1, 20, {'Muscle'}, 2, 5.0, 12.5));
report = recordICARound(report, makeRound(2, 18, {'Eye'}, 1, 9.0, 99.9));
testCase.verifyEqual(report.ica.varRemoved, 12.5);
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
