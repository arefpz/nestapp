function peaks = tepPeakFinder(waveform, times)
% TEPPEAKFINDER  Detect canonical TMS-EEG components in a TEP waveform.
%
%   peaks = tepPeakFinder(waveform, times) searches for six canonical
%   TMS-EEG components in the 1×T grand-mean waveform (µV) using
%   MATLAB's findpeaks (Signal Processing Toolbox).
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
%   Prominence threshold: max(0.05, 0.10 × waveform range). When multiple
%   qualifying peaks exist in a window the one with highest prominence wins.
%
%   Requires: Signal Processing Toolbox (findpeaks).
%
%   See also: computeTEPTStat, computeSplitHalf, computeCompositeQuality

% Component definitions: name | polarity | [tMin tMax] ms
COMPONENTS = {
    'N15',  'neg', [5   28  ];
    'P30',  'pos', [22  50  ];
    'N45',  'neg', [38  65  ];
    'P60',  'pos', [50  90  ];
    'N100', 'neg', [80  140 ];
    'P180', 'pos', [140 260 ];
};
nComp = size(COMPONENTS, 1);

% Prominence threshold: 10% of waveform range, floored at 0.05 µV
minProm = max(0.05, 0.10 * (max(waveform) - min(waveform)));

% Pre-allocate output struct array
blank = struct('name','','polarity','','latencyMs',NaN,'amplitudeUV',NaN,'found',false);
peaks = repmat(blank, 1, nComp);

for k = 1:nComp
    compName = COMPONENTS{k,1};
    polarity = COMPONENTS{k,2};
    win      = COMPONENTS{k,3};

    peaks(k).name     = compName;
    peaks(k).polarity = polarity;

    inWin = times >= win(1) & times <= win(2);
    if ~any(inWin)
        continue
    end
    segment  = waveform(inWin);
    winTimes = times(inWin);

    % Find peaks with minimum prominence; negate segment for negative components
    if strcmp(polarity, 'neg')
        [pv, locs, ~, proms] = findpeaks(-segment, 'MinProminence', minProm);
        pv = -pv;   % restore original sign (negative amplitudes)
    else
        [pv, locs, ~, proms] = findpeaks(segment, 'MinProminence', minProm);
    end

    if isempty(locs)
        continue
    end

    % Select the peak with highest prominence
    [~, best] = max(proms);
    peaks(k).latencyMs   = winTimes(locs(best));
    peaks(k).amplitudeUV = pv(best);
    peaks(k).found       = true;
end
end
