
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_icaDispatch
% TEST_ICADISPATCH  Regression tests for the Run ICA and Label ICA
%   Components dispatch cases. These would have caught:
%     - "Output argument 'sphere' (and possibly others) not assigned"
%       from runica when FastICA-specific params (approach, g,
%       stabilization) are forwarded to runica's parser
%     - "Unable to perform assignment with 0 elements" from
%       Label ICA Components when the dispatch looked up the wrong key
%       (the registry uses 'version', not 'iclabelVersion').
%
%   Run: runtests('tests/unit/test_icaDispatch')
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

% ── stripVarinKeys ────────────────────────────────────────────────────────

function test_stripVarinKeys_dropsMatchingPairs(testCase)
vars = {'icatype','runica','approach','symm','g','tanh','stabilization','on'};
out  = stripVarinKeys(vars, {'approach','g','stabilization'});
testCase.verifyEqual(out, {'icatype','runica'});
end

function test_stripVarinKeys_caseInsensitive(testCase)
vars = {'Approach','symm','G','tanh'};
out  = stripVarinKeys(vars, {'approach','g'});
testCase.verifyTrue(isempty(out), ...
    'Case-insensitive strip should drop all matched name-value pairs.');
end

function test_stripVarinKeys_emptyInputsAreSafe(testCase)
testCase.verifyEqual(stripVarinKeys({},          {'foo'}), {});
testCase.verifyEqual(stripVarinKeys({'a',1},     {}),      {'a',1});
testCase.verifyEqual(stripVarinKeys({'a',1},     {'b'}),   {'a',1});
end

function test_stripVarinKeys_preservesOrder(testCase)
vars = {'a',1,'b',2,'c',3,'d',4};
out  = stripVarinKeys(vars, {'b','d'});
testCase.verifyEqual(out, {'a',1,'c',3});
end

% ── Run ICA dispatch shape ────────────────────────────────────────────────

function test_runICA_runicaStripsFastICAParams(testCase)
% Simulates exactly what the Run ICA dispatch does when the registry's
% varin includes FastICA fields but icatype was overridden to runica.
vars = {'icatype','runica','approach','symm','g','tanh','stabilization','on'};
idx = find(strcmpi(vars, 'icatype'), 1);
testCase.verifyEqual(vars{idx+1}, 'runica');
filtered = stripVarinKeys(vars, {'approach','g','stabilization'});
testCase.verifyEqual(filtered, {'icatype','runica'}, ...
    ['Run ICA dispatch must strip FastICA-only params (approach, g, ' ...
     'stabilization) before forwarding to pop_runica when icatype = ' ...
     'runica - else runica.m crashes with "Output argument sphere not assigned".']);
end

function test_runICA_fastIcaKeepsAllParams(testCase)
% When icatype is fastica (the default), all params must pass through.
vars = {'icatype','fastica','approach','symm','g','tanh','stabilization','on'};
idx = find(strcmpi(vars, 'icatype'), 1);
isRunica = strcmpi(vars{idx+1}, 'runica');
testCase.verifyFalse(isRunica);
testCase.verifyEqual(vars, vars);
end

function test_runICA_artistTemplateProducesRunicaArgs(testCase)
% End-to-end check from the ARTIST template's stored params through
% paramsToVarin and the dispatch-side filter. Must produce a varin
% that only contains keys runica understands.
reg = stepRegistry();
artistStep = makePipelineStep('Run ICA', reg);
artistStep.params.icatype = 'runica';   % match the ARTIST template override

vars = paramsToVarin(artistStep.params);
vars = convertContainedStringsToChars(vars);
idx = find(strcmpi(vars, 'icatype'), 1);
if ~isempty(idx) && strcmpi(vars{idx+1}, 'runica')
    vars = stripVarinKeys(vars, {'approach','g','stabilization'});
end

remainingKeys = vars(1:2:end);
testCase.verifyFalse(any(strcmpi(remainingKeys, 'approach')), ...
    'ARTIST Run ICA dispatch must not pass "approach" to runica.');
testCase.verifyFalse(any(strcmpi(remainingKeys, 'g')), ...
    'ARTIST Run ICA dispatch must not pass "g" to runica.');
testCase.verifyFalse(any(strcmpi(remainingKeys, 'stabilization')), ...
    'ARTIST Run ICA dispatch must not pass "stabilization" to runica.');
testCase.verifyTrue(any(strcmpi(remainingKeys, 'icatype')), ...
    'icatype must survive the filter.');
end

% ── Label ICA Components dispatch shape ───────────────────────────────────

function test_labelICA_findsRegistryVersionKey(testCase)
% Reproduces the bug: the registry stores 'version', the dispatch must
% look for 'version' (not 'iclabelVersion'). Empty find -> 0-element
% index -> "Unable to perform assignment with 0 elements" error.
reg = stepRegistry();
labelStep = makePipelineStep('Label ICA Components', reg);
vars = paramsToVarin(labelStep.params);
vars = convertContainedStringsToChars(vars);

ind = find(strcmpi(vars, 'version'), 1);
testCase.verifyNotEmpty(ind, ...
    ['Label ICA Components dispatch must locate the registry "version" ' ...
     'key in varin. Old dispatch looked for "iclabelVersion" and crashed.']);

iclabelVersion = vars{ind+1};
testCase.verifyEqual(iclabelVersion, 'default', ...
    'Default ICLabel classifier version should be "default".');
end

function test_labelICA_missingKeyFallsBackToDefault(testCase)
% Even if a future refactor drops 'version' from a saved template, the
% dispatch should fall back to 'default' rather than crash with the
% "0 elements" error.
vars = {};
ind = find(strcmpi(vars, 'version'), 1);
if isempty(ind)
    iclabelVersion = 'default';
else
    iclabelVersion = vars{ind+1};
end
testCase.verifyEqual(iclabelVersion, 'default');
end
