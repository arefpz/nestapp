function score = computeCompositeQuality(EEG, report, varargin)
% COMPUTECOMPOSITEQUALITY  Composite pipeline quality score (0-1).
%
%   score = COMPUTECOMPOSITEQUALITY(EEG, report) returns a weighted
%   combination of three quality sub-metrics:
%
%     1. TEP t-stat        (weight 0.55) — mean(peak_per_trial) / SEM across trials
%     2. Split-half r      (weight 0.35) — odd vs even trial reproducibility
%     3. Retention penalty (weight 0.10) — soft penalty when <60% of trials survive
%
%   The t-statistic component replaces the old peak-GFP/baseline-RMS SNR.
%   Unlike SNR, the t-stat rewards both signal cleanliness AND trial count
%   simultaneously: t = mean_peak / (SD_peak / sqrt(n)).  A pipeline that
%   achieves high SNR by aggressively discarding trials will be penalised
%   by the falling sqrt(n) term.  This makes the composite a direct proxy
%   for statistical power — the probability of detecting the TEP in a
%   one-sample t-test.
%
%   The split-half component guards against optimising for muscle-artifact
%   reduction rather than genuine neural signal: true TEP waveforms replicate
%   across random trial halves; burst-like artifacts do not.
%
%   score = COMPUTECOMPOSITEQUALITY(EEG, []) uses only t-stat + split-half
%   (trial retention weight drops to 0 when report is unavailable).
%
%   score = COMPUTECOMPOSITEQUALITY(EEG, report, 'tstatWeight', 0.6, ...)
%   overrides individual weights.  Weights need not sum to 1 — normalised
%   internally.
%
%   See also: computeTEPTStat, computeSplitHalf, tepPeakFinder

% Default weights.
% Rationale: t-stat is the primary metric (0.55) because it captures both
% signal quality and trial quantity in one number.  Split-half (0.35) provides
% an independent reproducibility check.  The retention penalty (0.10) only
% fires below 60% trial survival to catch catastrophic data loss without
% discouraging reasonable cleaning.
W_TSTAT   = 0.55;
W_SH      = 0.35;
W_PENALTY = 0.10;
TSTAT_MAX = 15;          % t-stat above this → full score (prevents outlier dominance)
RETENTION_FLOOR = 0.60;  % penalty kicks in below this trial-retention fraction

p = inputParser;
addParameter(p, 'tstatWeight',    W_TSTAT);
addParameter(p, 'shWeight',       W_SH);
addParameter(p, 'penaltyWeight',  W_PENALTY);
addParameter(p, 'tstatMax',       TSTAT_MAX);
addParameter(p, 'retentionFloor', RETENTION_FLOOR);
parse(p, varargin{:});

wTstat   = p.Results.tstatWeight;
wSH      = p.Results.shWeight;
wPenalty = p.Results.penaltyWeight;
tstatMax = p.Results.tstatMax;
retFloor = p.Results.retentionFloor;

%% TEP t-statistic — normalised to [0, 1]
tstatRaw = computeTEPTStat(EEG);
if isnan(tstatRaw)
    tstatNorm = 0;
else
    tstatNorm = min(max(tstatRaw / tstatMax, 0), 1);
end

%% Split-half correlation — mapped from [-1, 1] to [0, 1]
shRaw = computeSplitHalf(EEG);
if isnan(shRaw)
    shNorm = 0;
else
    shNorm = (shRaw + 1) / 2;
end

% If both primary metrics are uncomputable (continuous data, too few trials, etc.)
% the composite is meaningless — return NaN rather than a misleading score.
if isnan(tstatRaw) && isnan(shRaw)
    score = NaN;
    return
end

%% Trial retention penalty
% Fires only when trial survival drops below retFloor — catches catastrophic
% data loss without double-penalising the reasonable cleaning that the
% t-stat denominator already handles.
retentionPenalty = 0;
if nargin >= 2 && isstruct(report) && isfield(report, 'trials') && ...
        isfield(report.trials, 'original') && report.trials.original > 0
    retention = report.trials.final / report.trials.original;
    retentionPenalty = max(0, retFloor - retention);  % 0 unless below floor
else
    % Report not available — skip penalty rather than penalising unfairly
    wPenalty = 0;
end

%% Weighted combination (normalise weights so they always sum to 1)
totalW = wTstat + wSH + wPenalty;
if totalW == 0
    score = NaN;
    return
end
score = (wTstat * tstatNorm + wSH * shNorm - wPenalty * retentionPenalty) / totalW;

% Clamp to [0, 1] — the penalty subtraction can push below 0 in extreme cases
score = max(0, min(score, 1));
end
