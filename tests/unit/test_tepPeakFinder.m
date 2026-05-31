
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_tepPeakFinder
% TEST_TEPPEAKFINDER  Unit tests for tepPeakFinder that do NOT require TESA.
%
%   The only meaningful unit-testable behaviour is that the function throws
%   the correct error (identifier + message) when TESA is missing.
%   All structural and behavioural tests (output fields, component counts,
%   polarity, latency windows) are in the integration suite:
%
%     tests/integration/test_tepPeakFinder_tesa.m   (requires TESA on path)
%
%   Run: runtests('tests/unit/test_tepPeakFinder')
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));
addpath(fullfile(r, 'tests', 'helpers'));
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

% ── error contract when TESA is absent ───────────────────────────────────────

function test_throwsCorrectErrorWhenTESAAbsent(testCase)
% Exercise the no-TESA error path even when TESA IS installed, by hiding
% tesa_peakanalysis from the path for the duration of the test (restored on
% cleanup). Asserts the hide worked so the test fails loudly rather than
% silently skipping.
cleanup = hideFromPath('tesa_peakanalysis'); %#ok<NASGU>
testCase.assertEmpty(which('tesa_peakanalysis'), ...
    'Could not hide tesa_peakanalysis from the path; cannot test the no-TESA branch.');

times    = -50:2:300;
waveform = zeros(1, numel(times));

testCase.verifyError(@() tepPeakFinder(waveform, times), 'tepPeakFinder:noTESA', ...
    'tepPeakFinder must throw tepPeakFinder:noTESA when TESA is not on the path');
end
