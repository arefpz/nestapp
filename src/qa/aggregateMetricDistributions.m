function dist = aggregateMetricDistributions(reports)
% AGGREGATEMETRICDISTRIBUTIONS  Per-metric value arrays across the batch.
%   dist = AGGREGATEMETRICDISTRIBUTIONS(reports)
%
%   For every Quality Gate metric that was *enabled* on at least one
%   gate anywhere in the batch, collect the values that were measured
%   across the batch (one entry per gate-instance that enabled it).
%   Disabled metrics are silently skipped so the dashboard does not
%   plot meaningless histograms.
%
%   Returns a struct array (one entry per enabled metric):
%     .name           field name on gate.metrics (e.g. 'pctBadTrials')
%     .displayName    pretty label for axis / panel title
%     .values         numeric vector (NaN-stripped)
%     .absThresholds  numeric vector of absolute thresholds when used
%                       (empty when no gate used absolute mode)
%     .batchCutoffs   numeric vector of batch cutoffs when used
%                       (empty when no gate used batch mode)
%     .mode           'absolute' | 'batch' | 'mixed'
%
%   The metric-to-threshold mapping mirrors finalizeBatchVerdicts: a
%   metric counts as enabled when the corresponding 'max*' or 'min*'
%   param on the gate's thresholds struct is > 0.

dist = struct('name', {}, 'displayName', {}, 'values', {}, ...
    'absThresholds', {}, 'batchCutoffs', {}, 'mode', {});

if isempty(reports), return, end

% Field <-> threshold-param mapping. Higher-is-worse (max) and
% lower-is-worse (min) families are handled together so the dashboard
% gets a single panel per metric.
fieldMap = { ...
    'nFlatChans',     'maxFlatChans',      'Flat channels',          'max'; ...
    'nSatChans',      'maxSatChans',       'Saturated channels',     'max'; ...
    'pctBadTrials',   'maxBadTrialPct',    '% bad trials',           'max'; ...
    'pctBadChans',    'maxBadChanPct',     '% bad channels',         'max'; ...
    'emgFraction',    'maxEMGFraction',    'EMG / muscle fraction',  'max'; ...
    'electrodeCount', 'maxElectrodeCount', 'Electrode-artifact comps','max'; ...
    'nTriggers',      'minTriggers',       'Trigger count',          'min'; ...
    'nTrials',        'minTrials',         'Trial count',            'min'; ...
    'rankRatio',      'minRankRatio',      'Rank / nbchan',          'min'};

acc = containers.Map();   % field name -> struct accumulator

for ri = 1:numel(reports)
    r = reports{ri};
    if ~isstruct(r) || ~isfield(r, 'quality') ...
            || ~isfield(r.quality, 'gates')
        continue
    end
    for gi = 1:numel(r.quality.gates)
        g = r.quality.gates{gi};
        if ~isfield(g, 'thresholds') || ~isfield(g, 'metrics'), continue, end
        for fi = 1:size(fieldMap, 1)
            metricName = fieldMap{fi, 1};
            paramName  = fieldMap{fi, 2};
            display    = fieldMap{fi, 3};
            if ~isfield(g.thresholds, paramName) ...
                    || g.thresholds.(paramName) <= 0
                continue
            end
            if ~isfield(g.metrics, metricName), continue, end
            v = g.metrics.(metricName);
            if isnan(v), continue, end

            if ~isKey(acc, metricName)
                acc(metricName) = struct( ...
                    'displayName',  display, ...
                    'values',       [], ...
                    'absThresholds',[], ...
                    'batchCutoffs', [], ...
                    'modes',        {{}});
            end
            entry = acc(metricName);
            entry.values(end+1) = v;
            if strcmpi(getOr(g, 'mode', 'absolute'), 'batch')
                entry.modes{end+1} = 'batch';
                if isfield(g, 'batchCutoffs') ...
                        && isfield(g.batchCutoffs, metricName) ...
                        && ~isnan(g.batchCutoffs.(metricName))
                    entry.batchCutoffs(end+1) = g.batchCutoffs.(metricName);
                end
            else
                entry.modes{end+1} = 'absolute';
                entry.absThresholds(end+1) = g.thresholds.(paramName);
            end
            acc(metricName) = entry;
        end
    end
end

% Build the output struct array, preserving the metric order from
% fieldMap (deterministic for downstream layout).
out = {};
for fi = 1:size(fieldMap, 1)
    metricName = fieldMap{fi, 1};
    if ~isKey(acc, metricName), continue, end
    entry = acc(metricName);
    out{end+1} = struct( ...
        'name',          metricName, ...
        'displayName',   entry.displayName, ...
        'values',        entry.values, ...
        'absThresholds', entry.absThresholds, ...
        'batchCutoffs',  entry.batchCutoffs, ...
        'mode',          summariseMode(entry.modes)); %#ok<AGROW>
end
if ~isempty(out)
    dist = [out{:}];
end
end

% -- small helpers --------------------------------------------------------

function m = summariseMode(modes)
hasAbs   = any(strcmp(modes, 'absolute'));
hasBatch = any(strcmp(modes, 'batch'));
if hasAbs && hasBatch
    m = 'mixed';
elseif hasBatch
    m = 'batch';
else
    m = 'absolute';
end
end

function v = getOr(s, field, default)
if isfield(s, field) && ~isempty(s.(field))
    v = s.(field);
else
    v = default;
end
end
