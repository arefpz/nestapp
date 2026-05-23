function gate = qualityGate(EEG, params)
% QUALITYGATE  Apply numeric quality thresholds and emit a verdict.
%   gate = QUALITYGATE(EEG, params) measures a battery of quality
%   metrics on EEG and compares each to a threshold from the step's
%   params struct. Returns a struct describing what was measured and
%   how the file scored.
%
%   Disabled checks: any threshold equal to 0 is skipped (the metric
%   is still recorded for batch-mode aggregation).
%
%   params fields (all optional except gateLabel / thresholdMode /
%   marginalSlack / outlierSigmas which always carry defaults):
%     gateLabel          string shown in CSV / report
%     thresholdMode      'absolute' | 'batch'
%     marginalSlack      scalar in (0, 1] - within slack*threshold of
%                        the cut is Marginal; above is Fail
%     expectedChans      exact-match check on EEG.nbchan
%     expectedSrate      exact-match check on EEG.srate
%     minTriggers        EEG.event count must be >= threshold
%     maxFlatChans       count of var ~ 0 channels must be <= threshold
%     maxSatChans        count of |data| > 250 uV channels must be <= threshold
%     minRankRatio       rank(EEG.data) / nbchan must be >= threshold
%     maxBadTrialPct     % outlier trials (SM-based) must be <= threshold
%     maxBadChanPct      % outlier channels (SM-based) must be <= threshold
%     minTrials          EEG.trials must be >= threshold
%     maxEMGFraction     ICA classifier emg fraction must be <= threshold
%     maxElectrodeCount  ICA electrode-artifact count must be <= threshold
%     outlierSigmas      scalar - N for the median + N * 1.4826 * MAD rule
%                        used both for SM outlier detection and (in batch
%                        mode) for cross-file outlier detection
%
%   Output gate:
%     .label, .mode
%     .verdict      'Pass' | 'Marginal' | 'Fail' | 'Pending'
%                     ('Pending' in batch mode; resolved by
%                      finalizeBatchVerdicts after the run completes)
%     .reasons      cellstr - one per check that flagged (empty for Pass)
%     .metrics      struct of raw metric values (one field per enabled
%                     check, used by batch-mode finalization)
%     .thresholds   struct mirroring the absolute thresholds that were
%                     used (for log readability and batch fallback)
%
%   Reuses Phase 1 helpers:
%     computeAttributeMatrix    SM matrix + flat/sat masks + per-axis medians
%     computeICAQualityMetrics  source-aware ICA classification
%
%   See also: finalizeBatchVerdicts, computeAttributeMatrix,
%             computeICAQualityMetrics

if nargin < 2 || ~isstruct(params), params = struct(); end
params = applyDefaults(params);

gate = struct( ...
    'label',      params.gateLabel, ...
    'mode',       params.thresholdMode, ...
    'verdict',    'Pass', ...
    'reasons',    {{}}, ...
    'metrics',    struct(), ...
    'thresholds', struct());

gate.metrics    = collectMetrics(EEG, params);
gate.thresholds = enabledThresholds(params);

if strcmpi(params.thresholdMode, 'batch')
    gate.verdict = 'Pending';
    return
end

[gate.verdict, gate.reasons] = evaluateAbsolute(gate.metrics, params);
end

% -- metric collection ----------------------------------------------------

function m = collectMetrics(EEG, params)
% Compute every metric used by the gate. Cheap operations always run;
% SM-derived metrics only when at least one SM-using check is enabled.

m.nbchan    = getField(EEG, 'nbchan', size(EEG.data, 1));
m.srate     = getField(EEG, 'srate', NaN);
m.nTriggers = numEvents(EEG);
m.nTrials   = getField(EEG, 'trials', max(size(EEG.data, 3), 1));
m.rankRatio = computeRankRatio(EEG);

needsSM = anyEnabled(params, {'maxFlatChans','maxSatChans', ...
    'maxBadTrialPct','maxBadChanPct'});
if needsSM
    [SM, sm] = computeAttributeMatrix(EEG, ...
        struct('attribute', 'minmax_no_tms'));
    m.nFlatChans   = sum(sm.flatChanMask);
    m.nSatChans    = sum(sm.satChanMask);
    m.pctBadTrials = pctOutliers(sm.perTrialMedian, params.outlierSigmas);
    skipMask = sm.flatChanMask | sm.satChanMask;
    chanScores = sm.perChanMedian;
    chanScores(skipMask) = NaN;     % exclude already-flagged channels
    m.pctBadChans  = pctOutliers(chanScores, params.outlierSigmas);
    m.smShape      = size(SM);      % stored for batch-mode diagnostics
else
    m.nFlatChans   = NaN;
    m.nSatChans    = NaN;
    m.pctBadTrials = NaN;
    m.pctBadChans  = NaN;
end

needsICA = anyEnabled(params, {'maxEMGFraction','maxElectrodeCount'});
if needsICA
    icaM = computeICAQualityMetrics(EEG);
    if isempty(icaM)
        m.emgFraction   = NaN;
        m.electrodeCount = NaN;
    else
        labels = {icaM.classification};
        m.emgFraction    = ratioOf(labels, {'EMG','Muscle','TMS Muscle'});
        m.electrodeCount = countOf(labels, ...
            {'Electrode','Elec Noise','Channel Noise'});
    end
else
    m.emgFraction    = NaN;
    m.electrodeCount = NaN;
end
end

% -- absolute-mode evaluation ---------------------------------------------

function [verdict, reasons] = evaluateAbsolute(m, p)
verdict = 'Pass';
reasons = {};

% Exact match checks (no Marginal tier).
if p.expectedChans > 0 && m.nbchan ~= p.expectedChans
    [verdict, reasons] = bump(verdict, reasons, 'Fail', ...
        sprintf('nbchan %d != expected %d', m.nbchan, p.expectedChans));
end
if p.expectedSrate > 0 && abs(m.srate - p.expectedSrate) > eps
    [verdict, reasons] = bump(verdict, reasons, 'Fail', ...
        sprintf('srate %g != expected %g Hz', m.srate, p.expectedSrate));
end

% Min checks: fail if metric < threshold; marginal if threshold <=
% metric < warn cutoff. Warn cutoff = warnAt if set, else threshold /
% slack (symmetric with the max-check default).
[verdict, reasons] = checkMin(verdict, reasons, m.nTriggers,  p.minTriggers, ...
    p.marginalSlack, p.minTriggersWarnAt,  'triggers');
[verdict, reasons] = checkMin(verdict, reasons, m.rankRatio,  p.minRankRatio, ...
    p.marginalSlack, p.minRankRatioWarnAt, 'rank/nbchan');
[verdict, reasons] = checkMin(verdict, reasons, m.nTrials,    p.minTrials, ...
    p.marginalSlack, p.minTrialsWarnAt,    'trials');

% Max checks: fail if metric > threshold; marginal if metric above the
% warn cutoff. Warn cutoff = warnAt if set, else slack * threshold.
[verdict, reasons] = checkMax(verdict, reasons, m.nFlatChans,     p.maxFlatChans, ...
    p.marginalSlack, p.maxFlatChansWarnAt,    'flat channels');
[verdict, reasons] = checkMax(verdict, reasons, m.nSatChans,      p.maxSatChans, ...
    p.marginalSlack, p.maxSatChansWarnAt,     'saturated channels');
[verdict, reasons] = checkMax(verdict, reasons, m.pctBadTrials,   p.maxBadTrialPct, ...
    p.marginalSlack, p.maxBadTrialPctWarnAt,  '% bad trials');
[verdict, reasons] = checkMax(verdict, reasons, m.pctBadChans,    p.maxBadChanPct, ...
    p.marginalSlack, p.maxBadChanPctWarnAt,   '% bad channels');
[verdict, reasons] = checkMax(verdict, reasons, m.emgFraction,    p.maxEMGFraction, ...
    p.marginalSlack, p.maxEMGFractionWarnAt,  'EMG fraction');
[verdict, reasons] = checkMax(verdict, reasons, m.electrodeCount, p.maxElectrodeCount, ...
    p.marginalSlack, p.maxElectrodeCountWarnAt, 'electrode-artifact comps');
end

function [verdict, reasons] = checkMin(verdict, reasons, value, threshold, slack, warnAt, name)
% slack is unused for min checks (see comment below) but kept in the
% signature for parity with checkMax.
%#ok<*INUSD>
if threshold <= 0 || isnan(value), return, end
% warnAt sits above threshold for a min check: the marginal band is
% [threshold, warnAt). With warnAt <= threshold (including the default
% of 0) there is no marginal band - any value below threshold fails
% outright. Min metrics like rankRatio are capped at 1.0, so a sensible
% default cannot be derived from slack alone; users opt into the
% marginal band by setting WarnAt above threshold.
warnCutoff = max(warnAt, threshold);
if value < threshold
    [verdict, reasons] = bump(verdict, reasons, 'Fail', ...
        sprintf('%s %g < %g', name, value, threshold));
elseif value < warnCutoff
    [verdict, reasons] = bump(verdict, reasons, 'Marginal', ...
        sprintf('%s %g near min %g', name, value, warnCutoff));
end
end

function [verdict, reasons] = checkMax(verdict, reasons, value, threshold, slack, warnAt, name)
if threshold <= 0 || isnan(value), return, end
warnCutoff = pickCutoff(warnAt, slack * threshold);
if value > threshold
    [verdict, reasons] = bump(verdict, reasons, 'Fail', ...
        sprintf('%s %g > %g', name, value, threshold));
elseif value > warnCutoff
    [verdict, reasons] = bump(verdict, reasons, 'Marginal', ...
        sprintf('%s %g near max %g', name, value, threshold));
end
end

function c = pickCutoff(warnAt, slackCutoff)
% warnAt > 0 overrides the slack-derived boundary; otherwise fall back.
if warnAt > 0
    c = warnAt;
else
    c = slackCutoff;
end
end

function [verdict, reasons] = bump(verdict, reasons, newVerdict, reason)
verdict = worstVerdict(verdict, newVerdict);
reasons{end+1} = reason;
end

function v = worstVerdict(a, b)
order = {'Pass', 'Marginal', 'Fail'};
ia = find(strcmp(a, order));
ib = find(strcmp(b, order));
if isempty(ia), ia = 0; end
if isempty(ib), ib = 0; end
v = order{max(ia, ib)};
end

% -- small math helpers ---------------------------------------------------

function pct = pctOutliers(values, nSigmas)
% % of finite entries that exceed median + nSigmas * 1.4826 * MAD.
v = values(:);
v = v(~isnan(v));
if isempty(v)
    pct = NaN;
    return
end
med   = median(v);
madV  = median(abs(v - med));
cutoff = med + nSigmas * 1.4826 * madV;
pct    = 100 * sum(v > cutoff) / numel(v);
end

function r = computeRankRatio(EEG)
if ~isfield(EEG, 'data') || isempty(EEG.data) || EEG.nbchan == 0
    r = NaN;
    return
end
data2D = reshape(EEG.data, size(EEG.data, 1), []);
r = rank(double(data2D)) / EEG.nbchan;
end

function n = numEvents(EEG)
if isfield(EEG, 'event') && ~isempty(EEG.event)
    n = numel(EEG.event);
else
    n = 0;
end
end

function frac = ratioOf(labels, targets)
hits = false(size(labels));
for k = 1:numel(targets)
    hits = hits | strcmp(labels, targets{k});
end
frac = sum(hits) / numel(labels);
end

function c = countOf(labels, targets)
hits = false(size(labels));
for k = 1:numel(targets)
    hits = hits | strcmp(labels, targets{k});
end
c = sum(hits);
end

% -- param plumbing -------------------------------------------------------

function p = applyDefaults(p)
defs = struct( ...
    'gateLabel',         'gate', ...
    'thresholdMode',     'absolute', ...
    'marginalSlack',     0.8, ...
    'expectedChans',     0, ...
    'expectedSrate',     0, ...
    'minTriggers',       0, ...
    'maxFlatChans',      0, ...
    'maxSatChans',       0, ...
    'minRankRatio',      0, ...
    'maxBadTrialPct',    0, ...
    'maxBadChanPct',     0, ...
    'minTrials',         0, ...
    'maxEMGFraction',    0, ...
    'maxElectrodeCount', 0, ...
    'outlierSigmas',     3, ...
    'minTriggersWarnAt',       0, ...
    'maxFlatChansWarnAt',      0, ...
    'maxSatChansWarnAt',       0, ...
    'minRankRatioWarnAt',      0, ...
    'maxBadTrialPctWarnAt',    0, ...
    'maxBadChanPctWarnAt',     0, ...
    'minTrialsWarnAt',         0, ...
    'maxEMGFractionWarnAt',    0, ...
    'maxElectrodeCountWarnAt', 0);
fns = fieldnames(defs);
for k = 1:numel(fns)
    if ~isfield(p, fns{k}) || isempty(p.(fns{k}))
        p.(fns{k}) = defs.(fns{k});
    end
end
if ischar(p.gateLabel) || isstring(p.gateLabel)
    p.gateLabel = char(p.gateLabel);
end
end

function t = enabledThresholds(p)
% Mirror of params, but only the threshold-bearing fields, kept for log.
fields = {'expectedChans','expectedSrate','minTriggers','maxFlatChans', ...
    'maxSatChans','minRankRatio','maxBadTrialPct','maxBadChanPct', ...
    'minTrials','maxEMGFraction','maxElectrodeCount'};
t = struct();
for k = 1:numel(fields)
    t.(fields{k}) = p.(fields{k});
end
t.marginalSlack = p.marginalSlack;
t.outlierSigmas = p.outlierSigmas;
end

function tf = anyEnabled(p, fields)
tf = false;
for k = 1:numel(fields)
    if isfield(p, fields{k}) && p.(fields{k}) > 0
        tf = true; return
    end
end
end

function v = getField(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
