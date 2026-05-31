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
    peaks(k).found       = strcmp(c.found, 'yes');
    peaks(k).latencyMs   = c.lat;
    peaks(k).amplitudeUV = c.amp;
end
end

function W = buildWinMatrix(defs)
W = [[defs.winStart]', [defs.winEnd]'];
end
