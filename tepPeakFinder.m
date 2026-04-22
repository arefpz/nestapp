function peaks = tepPeakFinder(waveform, times, compDefs)
% TEPPEAKFINDER  Detect TMS-EEG components in a grand-mean TEP waveform.
%
%   peaks = tepPeakFinder(waveform, times) uses the default canonical
%   component windows (Rogasch et al. 2013; Farzan et al. 2016).
%
%   peaks = tepPeakFinder(waveform, times, compDefs) uses the search
%   windows defined in compDefs, a struct array with fields:
%     .name        — component label string (e.g. 'N15')
%     .polarity    — 'neg' or 'pos'
%     .nomLatency  — nominal peak latency in ms (seed for tesa_peakanalysis)
%     .winStart    — search window start in ms
%     .winEnd      — search window end in ms
%
%   Inputs
%     waveform  1×T numeric vector — the ROI-averaged grand-mean TEP (µV)
%     times     1×T numeric vector — corresponding time points (ms)
%     compDefs  (optional) struct array — custom component definitions
%
%   Output
%     peaks     1×N struct array with fields:
%                 name        — component label
%                 polarity    — 'neg' or 'pos'
%                 latencyMs   — detected peak latency in ms (NaN if not found)
%                 amplitudeUV — detected peak amplitude in µV (NaN if not found)
%                 found       — logical scalar
%
%   Implementation: builds a minimal EEG stub (ROI.R1.tseries = waveform)
%   and calls tesa_peakanalysis, so peak detection is identical to what
%   TESA would produce in a standard pipeline.
%
%   Requires: TESA toolbox (tesa_peakanalysis).
%
%   See also: nestapp, overlayTEPComponents, populateTEPComponentTable

if nargin < 3 || isempty(compDefs)
    % Default canonical windows (Rogasch 2013 Clin Neurophysiol; Farzan 2016 Ann NY Acad Sci)
    compDefs = struct( ...
        'name',       {'N15',  'P30',  'N45',  'P60',  'N100', 'P180'}, ...
        'polarity',   {'neg',  'pos',  'neg',  'pos',  'neg',  'pos'}, ...
        'nomLatency', {15,     30,     45,     60,     100,    180}, ...
        'winStart',   {5,      22,     38,     50,     80,     140}, ...
        'winEnd',     {28,     50,     65,     90,     140,    260});
end

COMP_ORDER = {compDefs.name};
COMP_POL   = {compDefs.polarity};

% Pre-allocate output
blank = struct('name','','polarity','','latencyMs',NaN,'amplitudeUV',NaN,'found',false);
peaks = repmat(blank, 1, numel(COMP_ORDER));
for k = 1:numel(COMP_ORDER)
    peaks(k).name     = COMP_ORDER{k};
    peaks(k).polarity = COMP_POL{k};
end

% Minimal EEG stub — only the fields tesa_peakanalysis reads
EEGstub.times          = times(:)';
EEGstub.ROI.R1.tseries = waveform(:)';

% Separate neg/pos components and build array arguments for tesa_peakanalysis.
% NEG_WINS/POS_WINS are referenced inside evalc strings; %#ok<NASGU> suppresses
% the "variable might be unused" warning that the static analyzer would otherwise raise.
isNeg = strcmp(COMP_POL, 'neg');
isPos = strcmp(COMP_POL, 'pos');

negDefs = compDefs(isNeg);
posDefs = compDefs(isPos);

NEG_PEAKS = [negDefs.nomLatency];
POS_PEAKS = [posDefs.nomLatency];

NEG_WINS = buildWinMatrix(negDefs); %#ok<NASGU>
POS_WINS = buildWinMatrix(posDefs); %#ok<NASGU>

% Run TESA peak analysis; evalc suppresses the per-peak fprintf messages.
try
    if ~isempty(NEG_PEAKS)
        [~] = evalc('EEGstub = tesa_peakanalysis(EEGstub, ''ROI'', ''negative'', NEG_PEAKS, NEG_WINS)');
    end
    if ~isempty(POS_PEAKS)
        [~] = evalc('EEGstub = tesa_peakanalysis(EEGstub, ''ROI'', ''positive'', POS_PEAKS, POS_WINS)');
    end
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

%% ---- helpers ---------------------------------------------------------------

function W = buildWinMatrix(defs)
% Build an N×2 matrix of [winStart, winEnd] rows from a struct array.
W = zeros(numel(defs), 2);
for k = 1:numel(defs)
    W(k, :) = [defs(k).winStart, defs(k).winEnd];
end
end
