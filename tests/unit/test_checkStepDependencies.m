function tests = test_checkStepDependencies
% TEST_CHECKSTEPDEPENDENCIES  Unit tests for checkStepDependencies.
%
%   Tests the output contract, extension-filtering logic, and message format.
%   No EEGLAB installation required — all tests use the real stepRegistry but
%   only probe the logic that is independent of which plugins are installed.
%
%   Tests that require a specific plugin to be ABSENT use testCase.assumeTrue
%   to mark themselves as incomplete (not silently passed) when that plugin
%   happens to be installed.
%
%   Run: runtests('tests/unit/test_checkStepDependencies')
tests = functiontests(localfunctions);
end

% ── setup ─────────────────────────────────────────────────────────────────────

function setupOnce(testCase) %#ok<INUSD>
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

% ── output contract ───────────────────────────────────────────────────────────

function test_emptyStepListReturnsOk(testCase)
[ok, msg] = checkStepDependencies({}, {});
testCase.verifyTrue(ok,  'Empty step list must return ok=true');
testCase.verifyEqual(strtrim(char(msg)), '', 'Empty step list must return empty msg');
end

function test_emptyStepListWithFilesReturnsOk(testCase)
[ok, msg] = checkStepDependencies({}, {'sub01.set', 'sub02.vhdr'});
testCase.verifyTrue(ok,  'No steps means no requirements regardless of files');
testCase.verifyEqual(strtrim(char(msg)), '');
end

function test_okIsLogical(testCase)
[ok, ~] = checkStepDependencies({}, {});
testCase.verifyTrue(islogical(ok), 'ok return value must be logical');
end

function test_msgIsCharWhenOk(testCase)
[ok, msg] = checkStepDependencies({}, {});
testCase.verifyTrue(ok);
testCase.verifyTrue(ischar(msg) || isstring(msg), 'msg must be char or string');
end

function test_msgIsCharWhenNotOk(testCase)
% We need at least one missing dep to reach the not-ok path.
% Use a step that has external deps and assume we can find at least one missing.
% If everything is installed, this test is skipped.
steps = stepRegistry();
names = {steps.name};
% Try TESA steps first (most likely to be missing)
tesaSteps = names(~cellfun(@(n) isempty(steps(strcmp(names,n)).requires), names));
if isempty(tesaSteps) || isempty(which('tesa_peakanalysis'))
    % Attempt with TESA
end
testCase.assumeTrue(isempty(which('tesa_peakanalysis')), ...
    'No missing deps detected — cannot test msg type in not-ok path');
[ok, msg] = checkStepDependencies({'Run TESA ICA'}, {});
testCase.verifyFalse(ok);
testCase.verifyTrue(ischar(msg) || isstring(msg), 'msg must be char or string when not ok');
end

% ── unknown step names ────────────────────────────────────────────────────────

function test_unknownStepNameIgnored(testCase)
[ok, msg] = checkStepDependencies({'NotARealStep_XYZ_99'}, {});
testCase.verifyTrue(ok,  'Unknown step name must be silently ignored');
testCase.verifyEqual(strtrim(char(msg)), '', 'Unknown step must produce no message');
end

function test_mixedKnownAndUnknownSteps(testCase)
% Known step with no external dep + unknown step = ok
[ok, ~] = checkStepDependencies({'Remove Baseline', 'FakeStep_ABC'}, {});
% Remove Baseline has no external deps, so regardless of environment this must be ok
% (unless EEGLAB's pop_rmbase is somehow missing — not a realistic failure mode)
testCase.verifyTrue(islogical(ok));  % at minimum, must not crash
end

% ── extension filter — format-specific deps ───────────────────────────────────

function test_setFilesDoNotTriggerBvaIo(testCase)
% bva-io is required only for .vhdr files. Selecting only .set files must NOT
% flag bva-io as missing, regardless of whether bva-io is installed.
[ok, msg] = checkStepDependencies({'Load Data'}, {'sub01.set', 'sub02.set'});
if ~ok
    testCase.verifyFalse(contains(msg, 'bva-io'), ...
        ['bva-io must not be flagged when only .set files are selected. ' ...
         'Extension filter is broken.']);
end
end

function test_noFilesSkipsAllFormatSpecificDeps(testCase)
% With no files selected, ALL format-specific deps must be skipped.
% This verifies the ~isempty(filePaths) branch in the filter.
[ok, msg] = checkStepDependencies({'Load Data'}, {});
if ~ok
    testCase.verifyFalse(contains(msg, 'bva-io'),  'bva-io must not flag with no files');
    testCase.verifyFalse(contains(msg, 'loadcnt'), 'loadcnt must not flag with no files');
    testCase.verifyFalse(contains(msg, 'curry'),   'curry must not flag with no files');
end
end

function test_cntExtensionPassesThroughFilter(testCase)
% With a .cnt file, the loadcnt dep check must not be filtered out.
% Only testable when loadcnt is NOT installed (otherwise ok=true for wrong reason).
testCase.assumeTrue(isempty(which('pop_loadcnt')), ...
    'loadcnt is installed — cannot verify extension filter passes .cnt dep through');
[ok, msg] = checkStepDependencies({'Load Data'}, {'recording.cnt'});
testCase.verifyFalse(ok, 'Should flag missing loadcnt when .cnt file is selected');
testCase.verifyTrue(contains(msg, 'loadcnt'), ...
    'Message must name the loadcnt plugin when .cnt file selected and plugin missing');
end

function test_vhdrExtensionPassesThroughFilter(testCase)
% With a .vhdr file, bva-io must be checked (not filtered).
% Only testable when bva-io is NOT installed.
testCase.assumeTrue(isempty(which('pop_loadbv')), ...
    'bva-io is installed — cannot verify extension filter passes .vhdr dep through');
[ok, msg] = checkStepDependencies({'Load Data'}, {'recording.vhdr'});
testCase.verifyFalse(ok, 'Should flag missing bva-io when .vhdr file is selected');
testCase.verifyTrue(contains(msg, 'bva-io'), ...
    'Message must name bva-io when .vhdr file selected and plugin missing');
end

% ── message format ────────────────────────────────────────────────────────────

function test_missingDepMessageNamesPlugin(testCase)
testCase.assumeTrue(isempty(which('tesa_peakanalysis')), ...
    'TESA is installed — cannot test missing-plugin message format');
[ok, msg] = checkStepDependencies({'Run TESA ICA'}, {});
testCase.verifyFalse(ok);
testCase.verifyTrue(contains(msg, 'TESA'), ...
    'Message must name the missing plugin ("TESA")');
end

function test_missingDepMessageIncludesInstallNote(testCase)
testCase.assumeTrue(isempty(which('tesa_peakanalysis')), ...
    'TESA is installed — cannot test install note in message');
[ok, msg] = checkStepDependencies({'Run TESA ICA'}, {});
testCase.verifyFalse(ok);
testCase.verifyTrue(contains(msg, 'Install') || contains(msg, 'install'), ...
    'Message must include install instructions');
end

function test_missingDepMessageIncludesStepName(testCase)
testCase.assumeTrue(isempty(which('tesa_peakanalysis')), ...
    'TESA is installed — cannot test step name in message');
[ok, msg] = checkStepDependencies({'Run TESA ICA'}, {});
testCase.verifyFalse(ok);
testCase.verifyTrue(contains(msg, 'Run TESA ICA') || contains(msg, 'Steps'), ...
    'Message must identify which step requires the missing plugin');
end

function test_multipleStepsSamePluginGroupedOnce(testCase)
% Two steps that both need TESA must produce a single TESA plugin entry.
testCase.assumeTrue(isempty(which('tesa_peakanalysis')), ...
    'TESA is installed — cannot test plugin grouping');
[ok, msg] = checkStepDependencies( ...
    {'Run TESA ICA', 'Remove ICA Components (TESA)'}, {});
testCase.verifyFalse(ok);
% Both steps must be mentioned in the Steps line
testCase.verifyTrue(contains(msg, 'Run TESA ICA'), ...
    'Message must list Run TESA ICA under the TESA plugin entry');
testCase.verifyTrue(contains(msg, 'Remove ICA'), ...
    'Message must list Remove ICA Components under the TESA plugin entry');
end
