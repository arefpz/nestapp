function tstat = computeTEPTStat(EEG, varargin)
% COMPUTETEPTTSTAT  T-statistic of the TEP peak amplitude across trials.
%
%   tstat = COMPUTETEPTTSTAT(EEG) returns the one-sample t-statistic testing
%   whether the per-trial GFP peak in the post-stimulus window differs from
%   zero:
%
%       tstat = mean(peak_per_trial) / (std(peak_per_trial) / sqrt(n_trials))
%
%   Unlike peak_GFP / baseline_RMS (SNR), this metric rewards both signal
%   cleanliness AND trial retention simultaneously.  A pipeline keeping 200
%   moderate-quality trials can outscore one keeping 20 pristine trials when
%   the extra trials contribute more precision than they cost in noise — which
%   is precisely the regime where statistical power is maximised.
%
%   This makes tstat a direct proxy for the probability of detecting a
%   significant TEP in a one-sample t-test, which is the hypothesis most
%   TMS-EEG studies test.
%
%   tstat = COMPUTETEPTTSTAT(EEG, 'signalWin', [10 250]) uses a custom
%   post-stimulus window (milliseconds).
%
%   Returns NaN if EEG is not epoched, has fewer than 4 epochs, or the
%   signal window falls outside the data range.
%
%   See also: computeSplitHalf, computeCompositeQuality, tepPeakFinder

p = inputParser;
addParameter(p, 'signalWin', [10, 250]);
parse(p, varargin{:});
sigWin = p.Results.signalWin;

tstat = NaN;

if isempty(EEG) || ~isstruct(EEG),                   return; end
if ~isfield(EEG, 'times') || isempty(EEG.times),     return; end
if ~isfield(EEG, 'trials') || EEG.trials < 4,        return; end
if ~isfield(EEG, 'data')  || isempty(EEG.data),      return; end

sigMask = EEG.times >= sigWin(1) & EEG.times <= sigWin(2);
if ~any(sigMask), return; end

% Per-trial GFP in the signal window.
% GFP = std across channels at each time point — a single channel-independent
% amplitude measure that reflects the global signal strength.
% Result: nTime × nTrials
gfp = squeeze(std(EEG.data(:, sigMask, :), 0, 1));

% Handle edge case: single channel (std returns zeros)
if size(EEG.data, 1) < 2, return; end

% Peak GFP per trial — the amplitude the t-test is run on.
% max over the time axis → one value per trial.
if isvector(gfp)
    % Only one time point in window (shouldn't happen in practice)
    peakPerTrial = gfp(:);
else
    peakPerTrial = max(gfp, [], 1)';   % nTrials × 1
end

n  = numel(peakPerTrial);
sd = std(peakPerTrial);
if sd == 0 || isnan(sd), return; end

tstat = mean(peakPerTrial) / (sd / sqrt(n));
end
