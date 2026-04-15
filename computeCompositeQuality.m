function score = computeCompositeQuality(EEG, report, varargin)
% COMPUTECOMPOSITEQUALITY  Composite pipeline quality score (0-1).
%
%   score = COMPUTECOMPOSITEQUALITY(EEG, report) returns a weighted
%   combination of two quality sub-metrics:
%
%     1. Adjusted SNR   (weight 0.65) — grand-mean peak GFP / baseline RMS,
%                                       scaled by sqrt(n_final / n_original)
%     2. Split-half r   (weight 0.35) — Spearman-Brown corrected odd/even
%                                       trial reproducibility
%
%   ADJUSTED SNR combines signal quality with trial retention in one term.
%   Grand-mean noise scales as 1/sqrt(n), so SNR already grows with trial
%   count; multiplying by sqrt(retention) applies an additional explicit
%   penalty for overcleaning: discarding half the trials costs ~30% of the
%   score, matching the penalty structure of a t-statistic while keeping SNR
%   as the base quality measure (which, unlike GFP peak consistency, cannot
%   be inflated by systematic muscle artifacts).
%
%   SPLIT-HALF reproducibility guards against pipelines that achieve high SNR
%   by retaining a consistent artifact rather than genuine neural signal: true
%   TEP waveforms replicate across random trial halves; systematic artifacts
%   may too, but noise-dominated averages do not. The Spearman-Brown
%   correction r_sb = 2r/(1+r) adjusts for the reduced trial count in each
%   half, giving an estimate of reliability at the full trial count.
%
%   score = COMPUTECOMPOSITEQUALITY(EEG, []) omits the retention adjustment
%   (uses raw SNR) when trial count information is unavailable.
%
%   score = COMPUTECOMPOSITEQUALITY(EEG, report, 'snrWeight', 0.7, ...)
%   overrides individual weights. Weights need not sum to 1 — normalised
%   internally.
%
%   See also: computeTEPSNR, computeSplitHalf, tepPeakFinder

% Default weights.
% Adjusted SNR (0.65) is primary: captures both signal quality and the
% statistical cost of trial loss in a single interpretable number.
% Split-half (0.35) provides an independent reproducibility check that
% penalises consistent artifacts the SNR term cannot detect.
W_SNR = 0.65;
W_SH  = 0.35;
SNR_MAX = 8;  % adjusted SNR above this → full score (typical clean TMS-EEG ceiling)

p = inputParser;
addParameter(p, 'snrWeight', W_SNR);
addParameter(p, 'shWeight',  W_SH);
addParameter(p, 'snrMax',    SNR_MAX);
parse(p, varargin{:});

wSNR   = p.Results.snrWeight;
wSH    = p.Results.shWeight;
snrMax = p.Results.snrMax;

%% Adjusted SNR — SNR × sqrt(trial retention), normalised to [0, 1]
snrRaw = computeTEPSNR(EEG);

if isnan(snrRaw)
    % No pre-stimulus baseline — drop SNR weight rather than penalising.
    adjSNRnorm = 0;
    wSNR = 0;
else
    % Apply retention adjustment when trial counts are available.
    if isstruct(report) && isfield(report, 'trials') && ...
            isfield(report.trials, 'original') && report.trials.original > 0
        retention = report.trials.final / report.trials.original;
        adjSNR = snrRaw * sqrt(retention);
    else
        adjSNR = snrRaw;
    end
    adjSNRnorm = min(max(adjSNR / snrMax, 0), 1);
end

%% Split-half — Spearman-Brown corrected, mapped from [-1, 1] to [0, 1]
shRaw = computeSplitHalf(EEG);

if isnan(shRaw)
    shNorm = 0;
    wSH = 0;
else
    % Spearman-Brown correction: estimates reliability at full trial count.
    % Defined for all r; for r < -0.5 the denominator approaches zero so clamp.
    denomSB = 1 + max(shRaw, -0.99);
    shCorrected = (2 * shRaw) / denomSB;
    shNorm = (shCorrected + 1) / 2;  % map [-1, 1] → [0, 1]
end

%% Guard: if no metric is computable the score is meaningless
if wSNR == 0 && wSH == 0
    score = NaN;
    return
end

%% Weighted combination (normalise so weights always sum to 1)
totalW = wSNR + wSH;
score  = (wSNR * adjSNRnorm + wSH * shNorm) / totalW;
score  = max(0, min(score, 1));
end
