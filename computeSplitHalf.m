function r = computeSplitHalf(EEG, varargin)
% COMPUTESPLITHS  Compute split-half reproducibility of the TEP waveform.
%
%   r = COMPUTESPLITHS(EEG) computes the Pearson correlation between the
%   GFP waveforms derived from odd and even trial averages in the default
%   post-stimulus window [10 250] ms.  A value near +1 indicates a highly
%   reproducible TEP; values near 0 indicate noise dominates.
%
%   r = COMPUTESPLITHS(EEG, 'signalWin', [10 250]) uses a custom window.
%
%   Returns NaN if fewer than 4 epochs are present or the time window
%   falls outside the data.
%
%   Reference: Groppe et al. (2009) Reliable ICs via split-half.
%              NeuroImage doi:10.1016/j.neuroimage.2008.12.038
%
%   See also: computeTEPTStat, computeCompositeQuality, tepPeakFinder

p = inputParser;
addParameter(p, 'signalWin', [10, 250]);
parse(p, varargin{:});
sigWin = p.Results.signalWin;

r = NaN;

if isempty(EEG) || ~isstruct(EEG)
    return
end
if ~isfield(EEG, 'times') || isempty(EEG.times)
    return
end
if ~isfield(EEG, 'trials') || EEG.trials < 4
    return
end

t       = EEG.times;
sigMask = t >= sigWin(1) & t <= sigWin(2);
if ~any(sigMask)
    return
end

oddIdx  = 1:2:EEG.trials;
evenIdx = 2:2:EEG.trials;

% Average each half then compute GFP (std across channels)
oddAvg  = mean(EEG.data(:, sigMask, oddIdx),  3);  % ch × t
evenAvg = mean(EEG.data(:, sigMask, evenIdx), 3);  % ch × t

gfpOdd  = std(oddAvg,  0, 1);  % 1 × t
gfpEven = std(evenAvg, 0, 1);

% Pearson r between the two GFP waveforms.
% Guard: corrcoef throws a division-by-zero warning when either GFP is
% constant (e.g. a trial half that is entirely zeros after artifact removal).
if std(gfpOdd) == 0 || std(gfpEven) == 0
    return  % r stays NaN
end
corrMat = corrcoef(gfpOdd', gfpEven');
r = corrMat(1, 2);
end
