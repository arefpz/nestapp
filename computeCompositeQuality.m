function score = computeCompositeQuality(EEG, report, varargin)
% COMPUTECOMPOSITEQUALITY  Composite pipeline quality score (0-1).
%
%   score = COMPUTECOMPOSITEQUALITY(EEG, report) returns a weighted
%   combination of four quality sub-metrics:
%
%     1. TEP t-stat        (weight 0.45) — mean(peak_per_trial) / SEM across trials
%     2. Split-half r      (weight 0.30) — odd vs even trial reproducibility
%     3. SNR               (weight 0.15) — grand-mean peak GFP / baseline RMS
%     4. Retention penalty (weight 0.10) — soft penalty when <60% of trials survive
%
%   The t-statistic is the primary metric because it rewards both signal
%   cleanliness AND trial retention simultaneously: t = mean_peak /
%   (SD_peak / sqrt(n)).  A pipeline that achieves high SNR by aggressively
%   discarding trials will be penalised by the falling sqrt(n) term.
%
%   SNR (grand-mean peak GFP / baseline RMS) provides a complementary view
%   of amplitude cleaning quality that is independent of trial count.
%   Together, t-stat and SNR capture both statistical power and signal
%   amplitude relative to the noise floor.
%
%   The split-half component guards against optimising for muscle-artifact
%   reduction rather than genuine neural signal: true TEP waveforms replicate
%   across random trial halves; burst-like artifacts do not.
%
%   score = COMPUTECOMPOSITEQUALITY(EEG, []) uses only t-stat + split-half + SNR
%   (trial retention weight drops to 0 when report is unavailable).
%
%   score = COMPUTECOMPOSITEQUALITY(EEG, report, 'tstatWeight', 0.5, ...)
%   overrides individual weights.  Weights need not sum to 1 — normalised
%   internally.
%
%   See also: computeTEPTStat, computeTEPSNR, computeSplitHalf, tepPeakFinder

% Default weights.
% Rationale: t-stat is primary (0.45) — captures statistical power directly.
% Split-half (0.30) provides independent reproducibility check.
% SNR (0.15) adds a grand-mean amplitude quality view complementary to t-stat.
% Retention penalty (0.10) fires only below 60% trial survival to catch
% catastrophic data loss without discouraging reasonable cleaning.
W_TSTAT   = 0.45;
W_SH      = 0.30;
W_SNR     = 0.15;
W_PENALTY = 0.10;
TSTAT_MAX = 15;          % t-stat above this → full score (prevents outlier dominance)
SNR_MAX   = 8;           % SNR above this → full score (typical clean TMS-EEG ceiling)
RETENTION_FLOOR = 0.60;  % penalty kicks in below this trial-retention fraction

p = inputParser;
addParameter(p, 'tstatWeight',    W_TSTAT);
addParameter(p, 'shWeight',       W_SH);
addParameter(p, 'snrWeight',      W_SNR);
addParameter(p, 'penaltyWeight',  W_PENALTY);
addParameter(p, 'tstatMax',       TSTAT_MAX);
addParameter(p, 'snrMax',         SNR_MAX);
addParameter(p, 'retentionFloor', RETENTION_FLOOR);
parse(p, varargin{:});

wTstat   = p.Results.tstatWeight;
wSH      = p.Results.shWeight;
wSNR     = p.Results.snrWeight;
wPenalty = p.Results.penaltyWeight;
tstatMax = p.Results.tstatMax;
snrMax   = p.Results.snrMax;
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

%% SNR — normalised to [0, 1]
snrRaw = computeTEPSNR(EEG);
if isnan(snrRaw)
    % SNR requires a pre-stimulus baseline; drop its weight rather than
    % penalising data that was epoched without one.
    snrNorm = 0;
    wSNR    = 0;
else
    snrNorm = min(max(snrRaw / snrMax, 0), 1);
end

% If all primary metrics are uncomputable (continuous data, too few trials,
% etc.) the composite is meaningless — return NaN rather than a misleading score.
if isnan(tstatRaw) && isnan(shRaw) && isnan(snrRaw)
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
totalW = wTstat + wSH + wSNR + wPenalty;
if totalW == 0
    score = NaN;
    return
end
score = (wTstat * tstatNorm + wSH * shNorm + wSNR * snrNorm - wPenalty * retentionPenalty) / totalW;

% Clamp to [0, 1] — the penalty subtraction can push below 0 in extreme cases
score = max(0, min(score, 1));
end
