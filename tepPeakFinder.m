function peaks = tepPeakFinder(waveform, times)
% TEPPEAKFINDER  Detect canonical TMS-EEG components in a TEP waveform.
%
%   peaks = tepPeakFinder(waveform, times) searches for six canonical
%   TMS-EEG components in the 1×T grand-mean waveform (µV) by delegating
%   to TESA's tesa_peakanalysis, which checks that each candidate sample
%   is larger/smaller than its ±5 neighbouring samples (local-peak criterion).
%
%   Inputs
%     waveform  1×T numeric vector — the ROI-averaged grand-mean TEP (µV)
%     times     1×T numeric vector — corresponding time points (ms)
%
%   Output
%     peaks     1×6 struct array with fields:
%                 name        — component label ('N15','P30',…,'P180')
%                 polarity    — 'neg' or 'pos'
%                 latencyMs   — detected peak latency in ms (NaN if not found)
%                 amplitudeUV — detected peak amplitude in µV (NaN if not found)
%                 found       — logical scalar
%
%   Component search windows follow standard TMS-EEG conventions:
%     Rogasch et al. (2013) Clin Neurophysiol 124:2059-72
%     Farzan et al. (2016) Ann NY Acad Sci 1265:1-28
%
%   Implementation: builds a minimal EEG stub (ROI.R1.tseries = waveform)
%   and calls tesa_peakanalysis, so peak detection is identical to what
%   TESA would produce in a standard pipeline.
%
%   Requires: TESA toolbox (tesa_peakanalysis).
%
%   See also: computeTEPTStat, computeSplitHalf, computeCompositeQuality

% Component definitions — grouped by polarity for the two tesa_peakanalysis calls.
% %#ok<NASGU> below: variables used inside evalc strings; analyzer cannot see them.
NEG_PEAKS = [15, 45, 100];              %#ok<NASGU> % nominal latencies (ms)
NEG_WINS  = [5, 28; 38, 65; 80, 140];  %#ok<NASGU> % search windows (ms)
POS_PEAKS = [30, 60, 180];              %#ok<NASGU>
POS_WINS  = [22, 50; 50, 90; 140, 260]; %#ok<NASGU>

% Canonical order for the output struct
COMP_ORDER = {'N15','P30','N45','P60','N100','P180'};
COMP_POL   = {'neg','pos','neg','pos','neg','pos'};

% Pre-allocate output
blank = struct('name','','polarity','','latencyMs',NaN,'amplitudeUV',NaN,'found',false);
peaks = repmat(blank, 1, numel(COMP_ORDER));
for k = 1:numel(COMP_ORDER)
    peaks(k).name     = COMP_ORDER{k};
    peaks(k).polarity = COMP_POL{k};
end

% Minimal EEG stub — only the fields tesa_peakanalysis reads
EEGstub.times        = times(:)';         % 1×T
EEGstub.ROI.R1.tseries = waveform(:)';   % 1×T

% Run TESA peak analysis; evalc suppresses the per-peak fprintf messages.
try
    [~] = evalc('EEGstub = tesa_peakanalysis(EEGstub, ''ROI'', ''negative'', NEG_PEAKS, NEG_WINS)');
    [~] = evalc('EEGstub = tesa_peakanalysis(EEGstub, ''ROI'', ''positive'', POS_PEAKS, POS_WINS)');
catch ME
    warning('tepPeakFinder:tesaFailed', ...
        'tesa_peakanalysis failed (%s) — returning empty peaks.', ME.message);
    return
end

% Parse results from EEG stub into output struct
for k = 1:numel(COMP_ORDER)
    compName = COMP_ORDER{k};
    if ~isfield(EEGstub.ROI.R1, compName)
        continue
    end
    c = EEGstub.ROI.R1.(compName);
    peaks(k).found       = strcmp(c.found, 'yes');
    peaks(k).latencyMs   = c.lat;
    peaks(k).amplitudeUV = c.amp;
end
end
