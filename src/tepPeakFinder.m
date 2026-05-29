
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function peaks = tepPeakFinder(waveform, times, compDefs)
% TEPPEAKFINDER  Detect TMS-EEG components using TESA's peak analysis.
%
%   peaks = tepPeakFinder(waveform, times)
%   peaks = tepPeakFinder(waveform, times, compDefs)
%
%   Inputs
%     waveform  1xT numeric - ROI-averaged grand-mean TEP (uV)
%     times     1xT numeric - time points (ms)
%     compDefs  (optional) struct array with fields:
%                 .name, .polarity, .nomLatency, .winStart, .winEnd
%
%   Output
%     peaks  1xN struct array:
%              .name, .polarity, .latencyMs, .amplitudeUV, .found
%
%   A peak is reported found only if it is (a) a local extremum of the
%   correct polarity within its window (tesa_peakanalysis) AND (b) its
%   baseline-relative amplitude carries the component's defining sign
%   (negative for N components, positive for P components). The sign guard
%   rejects local minima that are not true negativities - notably the
%   valley between two positive peaks - per ERP/TEP nomenclature
%   (Luck 2014; Rogasch et al. 2017). It assumes the input waveform is
%   baseline-corrected so that zero is the pre-stimulus baseline.
%
%   Callers should pass the same smoothed, trial-averaged ROI waveform
%   they display, so detection runs on the curve the user sees rather than
%   on un-smoothed noise (see plotTEP / batchTEPExtract).
%
%   Requires: TESA toolbox (tesa_peakanalysis) on the MATLAB path.

if isempty(which('tesa_peakanalysis'))
    error('tepPeakFinder:noTESA', ...
        ['TESA toolbox not found on the MATLAB path. ' ...
         'Add TESA to the path via the EEGLAB plugin manager or Preferences.']);
end

if nargin < 3 || isempty(compDefs)
    % Default windows follow Beck et al. 2024 (Hum Brain Mapp, 45:e70048).
    compDefs = struct( ...
        'name',       {'N15',  'P30',  'N45',  'P60',  'N100', 'P180'}, ...
        'polarity',   {'neg',  'pos',  'neg',  'pos',  'neg',  'pos'}, ...
        'nomLatency', {15,     30,     45,     60,     100,    180}, ...
        'winStart',   {10,     20,     40,     50,     70,     150}, ...
        'winEnd',     {20,     40,     55,     70,     150,    240});
end

waveform = waveform(:)';
times    = times(:)';

% Pre-allocate output
blank = struct('name','','polarity','','latencyMs',NaN,'amplitudeUV',NaN,'found',false);
peaks = repmat(blank, 1, numel(compDefs));
for k = 1:numel(compDefs)
    peaks(k).name     = compDefs(k).name;
    peaks(k).polarity = compDefs(k).polarity;
end

% Build minimal EEG stub - only the fields tesa_peakanalysis reads
EEGstub.times          = times;
EEGstub.ROI.R1.tseries = waveform;

isNeg = strcmp({compDefs.polarity}, 'neg');
isPos = strcmp({compDefs.polarity}, 'pos');

negDefs = compDefs(isNeg);
posDefs = compDefs(isPos);

% Call tesa_peakanalysis; evalc suppresses its per-peak fprintf output.
if any(isNeg)
    NEG_PEAKS = [negDefs.nomLatency];
    NEG_WINS  = buildWinMatrix(negDefs);
    EEGstub = tesa_peakanalysis(EEGstub, 'ROI', 'negative', NEG_PEAKS, NEG_WINS);
end
if any(isPos)
    POS_PEAKS = [posDefs.nomLatency];
    POS_WINS  = buildWinMatrix(posDefs);
    EEGstub = tesa_peakanalysis(EEGstub, 'ROI', 'positive', POS_PEAKS, POS_WINS);
end

% Parse results
for k = 1:numel(compDefs)
    compName = compDefs(k).name;
    if ~isfield(EEGstub.ROI.R1, compName)
        continue
    end
    c = EEGstub.ROI.R1.(compName);
    found = strcmp(c.found, 'yes');

    % Polarity guard. ERP/TEP nomenclature defines an N component as a
    % NEGATIVE baseline-relative deflection and a P component as a POSITIVE
    % one (Luck 2014; Rogasch et al. 2017). tesa_peakanalysis only tests for
    % a local extremum (turning point), so it reports the valley between two
    % positive peaks (e.g. between P30 and P60) as a "negative" peak whose
    % amplitude is still positive - a relative minimum, not a negativity.
    % A genuine component must carry its defining sign relative to the
    % (baseline-corrected) pre-stimulus zero.
    if found
        isNeg = strcmpi(compDefs(k).polarity, 'neg');
        if (isNeg && c.amp >= 0) || (~isNeg && c.amp <= 0)
            found = false;
        end
    end

    peaks(k).found       = found;
    peaks(k).latencyMs   = c.lat;
    peaks(k).amplitudeUV = c.amp;
end
end

function W = buildWinMatrix(defs)
W = [[defs.winStart]', [defs.winEnd]'];
end
