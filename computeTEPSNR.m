function snr = computeTEPSNR(EEG, varargin)
% COMPUTETEPSNR  Grand-mean peak-GFP to baseline-RMS signal-to-noise ratio.
%
%   snr = COMPUTETEPSNR(EEG) computes the SNR of the TMS-evoked potential
%   from the trial-averaged (grand-mean) waveform:
%
%       snr = max(GFP in signalWin) / RMS(GFP in baselineWin)
%
%   where GFP (global field power) is the standard deviation across channels
%   at each time point of the grand-mean EEG.
%
%   Unlike the t-statistic (which operates per-trial and penalises low trial
%   counts), SNR operates on the grand average and provides a complementary
%   view of amplitude cleaning quality — how large the peak response is
%   relative to the pre-stimulus noise floor.
%
%   snr = COMPUTETEPSNR(EEG, 'signalWin', [10 250], 'baselineWin', [-200 -10])
%   uses custom windows (milliseconds).
%
%   Returns NaN if the data lacks a pre-stimulus baseline, has fewer than
%   4 epochs, has fewer than 2 channels, or either window falls outside
%   the data range.
%
%   See also: computeTEPTStat, computeSplitHalf, computeCompositeQuality

p = inputParser;
addParameter(p, 'signalWin',   [10,  250]);
addParameter(p, 'baselineWin', [-200, -10]);
parse(p, varargin{:});
sigWin = p.Results.signalWin;
basWin = p.Results.baselineWin;

snr = NaN;

if isempty(EEG) || ~isstruct(EEG),               return; end
if ~isfield(EEG, 'times') || isempty(EEG.times), return; end
if ~isfield(EEG, 'data')  || isempty(EEG.data),  return; end
if size(EEG.data, 1) < 2,                        return; end
if ~isfield(EEG, 'trials') || EEG.trials < 4,    return; end

sigMask = EEG.times >= sigWin(1) & EEG.times <= sigWin(2);
basMask = EEG.times >= basWin(1) & EEG.times <= basWin(2);
if ~any(sigMask) || ~any(basMask)
    return
end

% Grand-mean across trials then compute GFP (std across channels at each time point).
grandMean = mean(EEG.data, 3);  % ch × t
gfp = std(grandMean, 0, 1);     % 1 × t

peakSignal = max(gfp(sigMask));
rmsNoise   = rms(gfp(basMask));

if rmsNoise == 0 || isnan(rmsNoise)
    return
end

snr = peakSignal / rmsNoise;
end
