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
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

% ── error contract when TESA is absent ───────────────────────────────────────

function test_throwsCorrectErrorWhenTESAAbsent(testCase)
% This test only runs when TESA is NOT on the path.
% If TESA is present, the test is skipped (incomplete) — not silently passed.
testCase.assumeTrue(isempty(which('tesa_peakanalysis')), ...
    'TESA is installed — this test only exercises the no-TESA error path');

times    = -50:2:300;
waveform = zeros(1, numel(times));

testCase.verifyError(@() tepPeakFinder(waveform, times), 'tepPeakFinder:noTESA', ...
    'tepPeakFinder must throw tepPeakFinder:noTESA when TESA is not on the path');
end
