
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_tepPeakFinder_tesa
% TEST_TEPPEAKFINDER_TESA  Integration tests for tepPeakFinder requiring TESA.
%
%   Tests output structure, component count, polarity, field types, and the
%   custom-compDefs interface.  All tests call tesa_peakanalysis via
%   tepPeakFinder — TESA must be on the MATLAB path.
%
%   Run: runtests('tests/integration/test_tepPeakFinder_tesa')
tests = functiontests(localfunctions);
end

% ── setup ─────────────────────────────────────────────────────────────────────

function setupOnce(testCase)
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));

% Hard requirement — if TESA is absent these tests must fail visibly, not pass.
if isempty(which('tesa_peakanalysis'))
    testCase.assumeFail( ...
        ['TESA (tesa_peakanalysis) is not on the MATLAB path. ' ...
         'Install TESA via the EEGLAB Plugin Manager before running these tests.']);
end
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

% ── helpers ───────────────────────────────────────────────────────────────────

function [waveform, times] = flatWaveform()
times    = -50:2:300;
waveform = zeros(1, numel(times));
end

% ── output struct contract ────────────────────────────────────────────────────

function test_outputIsStructArray(testCase)
[w, t] = flatWaveform();
peaks = tepPeakFinder(w, t);
testCase.verifyTrue(isstruct(peaks),         'Output must be a struct array');
testCase.verifyGreaterThan(numel(peaks), 0,  'Must return at least one entry');
end

function test_outputHasRequiredFields(testCase)
[w, t] = flatWaveform();
peaks = tepPeakFinder(w, t);
for fi = 1:numel({'name','polarity','latencyMs','amplitudeUV','found'})
    field = {'name','polarity','latencyMs','amplitudeUV','found'};
    testCase.verifyTrue(isfield(peaks, field{fi}), ...
        sprintf('Output struct missing field: %s', field{fi}));
end
end

function test_defaultWindowsReturn6Components(testCase)
[w, t] = flatWaveform();
peaks = tepPeakFinder(w, t);
testCase.verifyEqual(numel(peaks), 6, ...
    'Default call must return exactly 6 component entries (N15 P30 N45 P60 N100 P180)');
end

function test_defaultComponentNamesInOrder(testCase)
[w, t] = flatWaveform();
peaks = tepPeakFinder(w, t);
expected = {'N15','P30','N45','P60','N100','P180'};
for i = 1:numel(expected)
    testCase.verifyEqual(peaks(i).name, expected{i}, ...
        sprintf('Component %d should be %s, got %s', i, expected{i}, peaks(i).name));
end
end

function test_polarityStringsAreNegOrPos(testCase)
[w, t] = flatWaveform();
peaks = tepPeakFinder(w, t);
for i = 1:numel(peaks)
    testCase.verifyTrue(ismember(peaks(i).polarity, {'neg','pos'}), ...
        sprintf('Component %d polarity must be "neg" or "pos", got "%s"', ...
            i, peaks(i).polarity));
end
end

function test_latencyAndAmplitudeAreNumeric(testCase)
[w, t] = flatWaveform();
peaks = tepPeakFinder(w, t);
for i = 1:numel(peaks)
    testCase.verifyTrue(isnumeric(peaks(i).latencyMs), ...
        sprintf('%s: latencyMs must be numeric', peaks(i).name));
    testCase.verifyTrue(isnumeric(peaks(i).amplitudeUV), ...
        sprintf('%s: amplitudeUV must be numeric', peaks(i).name));
end
end

function test_foundIsLogicalOrNumeric(testCase)
[w, t] = flatWaveform();
peaks = tepPeakFinder(w, t);
for i = 1:numel(peaks)
    testCase.verifyTrue(islogical(peaks(i).found) || isnumeric(peaks(i).found), ...
        sprintf('%s: found must be logical or numeric', peaks(i).name));
end
end

% ── custom compDefs interface ─────────────────────────────────────────────────

function test_customCompDefsReducesOutputCount(testCase)
[w, t] = flatWaveform();
compDefs = struct( ...
    'name',       {'N100', 'P180'}, ...
    'polarity',   {'neg',  'pos'}, ...
    'nomLatency', {100,    180}, ...
    'winStart',   {80,     140}, ...
    'winEnd',     {140,    260});
peaks = tepPeakFinder(w, t, compDefs);
testCase.verifyEqual(numel(peaks), 2, ...
    '2-component compDefs must yield exactly 2 output structs');
testCase.verifyEqual(peaks(1).name, 'N100');
testCase.verifyEqual(peaks(2).name, 'P180');
end

function test_customCompDefsPreservesNames(testCase)
[w, t] = flatWaveform();
compDefs = struct( ...
    'name',       {'MyPeak'}, ...
    'polarity',   {'pos'}, ...
    'nomLatency', {60}, ...
    'winStart',   {45}, ...
    'winEnd',     {80});
peaks = tepPeakFinder(w, t, compDefs);
testCase.verifyEqual(peaks(1).name, 'MyPeak', ...
    'Custom component name must appear in output exactly');
end

function test_emptyCompDefsUsesDefaults(testCase)
[w, t] = flatWaveform();
peaks = tepPeakFinder(w, t, []);
testCase.verifyEqual(numel(peaks), 6, ...
    'Passing [] as compDefs must fall back to the 6 canonical components');
end

% ── polarity (sign) guard ──────────────────────────────────────────────────────
% ERP/TEP nomenclature: an N component is a negative baseline-relative
% deflection, a P component is positive. tesa_peakanalysis only finds a
% turning point, so the valley between two positive peaks is reported as a
% "negative" peak with positive amplitude. tepPeakFinder must reject that.

function [w, t] = twoPositivePeaksWaveform()
% P30 and P60 positivities with a positive-valued valley (~45 ms) between
% them, and a genuine N100 negativity. Mirrors the reported N45 false
% positive (valley between P30/P60) and a real negative component.
t = -50:1:300;
w = 2.0*exp(-((t-30)/8).^2) ...   % P30
  + 2.0*exp(-((t-60)/9).^2) ...   % P60
  - 3.0*exp(-((t-100)/18).^2);    % N100 (dips well below zero)
end

function test_negativeComponentWithPositiveAmplitudeRejected(testCase)
[w, t] = twoPositivePeaksWaveform();
compDefs = struct('name',{'N45'},'polarity',{'neg'}, ...
    'nomLatency',{45},'winStart',{40},'winEnd',{55});
peaks = tepPeakFinder(w, t, compDefs);
testCase.verifyFalse(logical(peaks(1).found), ...
    ['A negative component whose amplitude is positive (the valley ' ...
     'between P30 and P60) must NOT be reported as found.']);
end

function test_genuineNegativeComponentStillFound(testCase)
[w, t] = twoPositivePeaksWaveform();
compDefs = struct('name',{'N100'},'polarity',{'neg'}, ...
    'nomLatency',{100},'winStart',{80},'winEnd',{140});
peaks = tepPeakFinder(w, t, compDefs);
testCase.verifyTrue(logical(peaks(1).found), ...
    'A genuine negative deflection (amplitude < 0) must still be found.');
testCase.verifyLessThan(peaks(1).amplitudeUV, 0, ...
    'A found N component must have negative amplitude.');
end

function test_positiveComponentWithNegativeAmplitudeRejected(testCase)
% Inverse case: a "positive" component sitting in a negative trough.
t = -50:1:300;
w = -2.0*exp(-((t-30)/8).^2) - 2.0*exp(-((t-60)/9).^2);  % two negativities
compDefs = struct('name',{'P45'},'polarity',{'pos'}, ...
    'nomLatency',{45},'winStart',{40},'winEnd',{55});
peaks = tepPeakFinder(w, t, compDefs);
testCase.verifyFalse(logical(peaks(1).found), ...
    'A positive component whose amplitude is negative must NOT be found.');
end
