
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_nestappDoctor
% TEST_NESTAPPDOCTOR  Unit tests for the environment diagnostics command.
%
%   nestappDoctor must run with or without EEGLAB on the path (it reports
%   absence rather than erroring), so these are plain unit tests.
%
%   Run: runtests('tests/unit/test_nestappDoctor')
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

% ── output contract ─────────────────────────────────────────────────────────

function test_returnsCharReportAndStruct(testCase)
[report, info] = nestappDoctor('Quiet', true);
testCase.verifyClass(report, 'char');
testCase.verifyNotEmpty(report);
testCase.verifyClass(info, 'struct');
end

function test_infoHasExpectedFields(testCase)
[~, info] = nestappDoctor('Quiet', true);
for f = {'nestapp','matlab','toolboxes','eeglab','dependencies', ...
         'parallel','prefs','shadows','problems'}
    testCase.verifyTrue(isfield(info, f{1}), ...
        sprintf('info struct missing field "%s"', f{1}));
end
end

function test_problemsIsCellstr(testCase)
[~, info] = nestappDoctor('Quiet', true);
testCase.verifyTrue(iscell(info.problems), 'problems must be a cell array');
testCase.verifyTrue(all(cellfun(@ischar, info.problems)), ...
    'every problem must be a char string');
end

function test_versionAppearsInReport(testCase)
report = nestappDoctor('Quiet', true);
testCase.verifyTrue(contains(report, nestappVersion()), ...
    'report must include the nestapp version');
end

function test_reportHasProblemsSection(testCase)
report = nestappDoctor('Quiet', true);
testCase.verifyTrue(contains(report, 'Problems detected'), ...
    'report must contain a "Problems detected" section');
end

% ── dependency list is derived from the registry (stays in sync) ──────────────

function test_dependencyListDerivedFromRegistry(testCase)
[~, info] = nestappDoctor('Quiet', true);
testCase.verifyNotEmpty(info.dependencies, ...
    'dependency list must be populated from stepRegistry requires');
% Every step requirement plugin should be represented.
plugins = {info.dependencies.plugin};
testCase.verifyTrue(any(contains(plugins, 'TESA')), ...
    'TESA (a registry requirement) must appear in the dependency list');
end

function test_optionalLoadersFlaggedOptional(testCase)
% Format-specific loaders (e.g. bva-io for .vhdr) carry a fileExt in the
% registry and must be marked optional so their absence is not a problem.
[~, info] = nestappDoctor('Quiet', true);
fns = {info.dependencies.fn};
idx = find(strcmp(fns, 'pop_loadbv'), 1);
testCase.assumeNotEmpty(idx, 'pop_loadbv requirement not present in registry');
testCase.verifyTrue(info.dependencies(idx).optional, ...
    'pop_loadbv (a .vhdr loader) must be marked optional');
end

% ── robustness ────────────────────────────────────────────────────────────────

function test_runsQuietlyWithoutError(testCase)
% Quiet mode must not print and must not throw.
fn = @() nestappDoctor('Quiet', true);
testCase.verifyWarningFree(fn);
end
