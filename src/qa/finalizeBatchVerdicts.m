function reports = finalizeBatchVerdicts(reports)
% FINALIZEBATCHVERDICTS  Resolve Pending gates with batch-relative thresholds.
%   reports = FINALIZEBATCHVERDICTS(reports)
%
%   reports : cell array of pipeline report structs returned by
%             runPipelineCore. Each report may carry
%             report.quality.gates{*}, some of which may be in
%             'batch' mode with verdict = 'Pending'.
%
%   Gates are grouped across files by (label, stepIndex). For each
%   group, every higher-is-worse metric is compared to:
%       cutoff = median + outlierSigmas * 1.4826 * MAD
%   A file whose metric exceeds the cutoff is flagged 'Fail'. A file
%   within marginalSlack of the cutoff is 'Marginal'. Otherwise 'Pass'.
%
%   Verdicts for non-batch gates and gates that were already resolved
%   (not 'Pending') are left untouched. report.quality.worstVerdict
%   is recomputed on every report that had at least one gate updated.
%
%   See also: qualityGate, runPipelineCore

if isempty(reports), return, end

% Higher-is-worse metric fields evaluated by batch mode, paired with
% the threshold-param name the user uses to opt in. Min-style metrics
% (triggers / rank / trials) are intentionally absent because a
% uniformly low batch would silently pass; those belong to absolute
% mode where the experimenter sets the absolute floor explicitly.
%
% Per-field opt-in: a metric is compared in batch mode only when the
% corresponding 'max*' threshold on the gate's thresholds struct is
% > 0. That keeps incidental metrics (e.g., pctBadTrials gets computed
% whenever maxSatChans is enabled because both share the SM helper)
% from silently driving the verdict.
FIELD_TO_PARAM = { ...
    'nFlatChans',     'maxFlatChans'; ...
    'nSatChans',      'maxSatChans'; ...
    'pctBadTrials',   'maxBadTrialPct'; ...
    'pctBadChans',    'maxBadChanPct'; ...
    'emgFraction',    'maxEMGFraction'; ...
    'electrodeCount', 'maxElectrodeCount'};

% Build an index of pending (reportIdx, gateIdx) pairs grouped by
% (gateLabel, stepIndex).
groups = struct();   % field name -> struct(members, label, stepIndex, sigmas, slack)
touchedReports = false(1, numel(reports));

for ri = 1:numel(reports)
    r = reports{ri};
    if ~isfield(r, 'quality') || ~isfield(r.quality, 'gates'), continue, end
    for gi = 1:numel(r.quality.gates)
        g = r.quality.gates{gi};
        if ~strcmpi(g.mode, 'batch') || ~strcmpi(g.verdict, 'Pending')
            continue
        end
        key = makeGroupKey(g);
        if ~isfield(groups, key)
            groups.(key) = struct( ...
                'label',     g.label, ...
                'stepIndex', getOr(g, 'stepIndex', 0), ...
                'sigmas',    getOr(g.thresholds, 'outlierSigmas', 3), ...
                'slack',     getOr(g.thresholds, 'marginalSlack', 0.8), ...
                'members',   {{}});
        end
        groups.(key).members{end+1} = [ri gi];
    end
end

groupNames = fieldnames(groups);
for k = 1:numel(groupNames)
    grp = groups.(groupNames{k});

    % Determine which fields the user opted into on this gate (any
    % member's thresholds work; they all came from the same step).
    sampleIdx = grp.members{1};
    sampleThresholds = reports{sampleIdx(1)}.quality.gates{sampleIdx(2)}.thresholds;
    enabledFields = {};
    enabledParams = {};
    for fi = 1:size(FIELD_TO_PARAM, 1)
        paramName = FIELD_TO_PARAM{fi, 2};
        if isfield(sampleThresholds, paramName) ...
                && sampleThresholds.(paramName) > 0
            enabledFields{end+1} = FIELD_TO_PARAM{fi, 1}; %#ok<AGROW>
            enabledParams{end+1} = paramName; %#ok<AGROW>
        end
    end

    % Collect each member's metric values for each enabled field.
    cutoffs = struct();
    for fi = 1:numel(enabledFields)
        f = enabledFields{fi};
        vals = zeros(1, numel(grp.members));
        for mi = 1:numel(grp.members)
            idx = grp.members{mi};
            m = reports{idx(1)}.quality.gates{idx(2)}.metrics;
            vals(mi) = safeGet(m, f);
        end
        cutoffs.(f) = madCutoff(vals, grp.sigmas);
    end

    % Apply cutoffs and write back into each member.
    for mi = 1:numel(grp.members)
        idx = grp.members{mi};
        ri = idx(1); gi = idx(2);
        m  = reports{ri}.quality.gates{gi}.metrics;

        verdict = 'Pass';
        reasons = {};
        for fi = 1:numel(enabledFields)
            f         = enabledFields{fi};
            paramName = enabledParams{fi};
            cutoff    = cutoffs.(f);
            value     = safeGet(m, f);
            if isnan(cutoff) || isnan(value), continue, end
            % Per-metric warn override (Phase 4): if *WarnAt > 0 on the
            % gate's thresholds, use it as the Marginal cutoff instead
            % of slack * batchCutoff. The Fail boundary (batch cutoff)
            % is unchanged.
            warnAt     = getOr(sampleThresholds, [paramName 'WarnAt'], 0);
            warnCutoff = pickCutoff(warnAt, grp.slack * cutoff);
            if value > cutoff
                [verdict, reasons] = bump(verdict, reasons, 'Fail', ...
                    sprintf('%s %g > batch cutoff %g', f, value, cutoff));
            elseif value > warnCutoff
                [verdict, reasons] = bump(verdict, reasons, 'Marginal', ...
                    sprintf('%s %g near batch cutoff %g', f, value, cutoff));
            end
        end

        reports{ri}.quality.gates{gi}.verdict      = verdict;
        reports{ri}.quality.gates{gi}.reasons      = reasons;
        reports{ri}.quality.gates{gi}.batchCutoffs = cutoffs;
        touchedReports(ri) = true;
    end
end

% Recompute worstVerdict on every touched report.
for ri = find(touchedReports)
    reports{ri} = refreshWorstVerdict(reports{ri});
end
end

% -- small helpers --------------------------------------------------------

function key = makeGroupKey(g)
% Sanitise label + stepIndex into a valid struct field name.
label  = regexprep(getOr(g, 'label', 'gate'),    '[^A-Za-z0-9_]', '_');
stepIx = getOr(g, 'stepIndex', 0);
key    = sprintf('g_%s__%d', label, stepIx);
end

function r = refreshWorstVerdict(r)
worst = 'Pass';
for gi = 1:numel(r.quality.gates)
    worst = worstVerdict(worst, r.quality.gates{gi}.verdict);
end
r.quality.worstVerdict = worst;
end

function c = madCutoff(vals, nSigmas)
v = vals(:);
v = v(~isnan(v));
if isempty(v)
    c = NaN; return
end
med  = median(v);
madV = median(abs(v - med));
c    = med + nSigmas * 1.4826 * madV;
end

function v = safeGet(s, name)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = NaN;
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

function v = getOr(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function [verdict, reasons] = bump(verdict, reasons, newVerdict, reason)
verdict = worstVerdict(verdict, newVerdict);
reasons{end+1} = reason;
end

function v = worstVerdict(a, b)
order = {'Pass', 'Marginal', 'Fail', 'Pending'};
ia = find(strcmp(a, order));
ib = find(strcmp(b, order));
if isempty(ia), ia = 0; end
if isempty(ib), ib = 0; end
v = order{max(ia, ib)};
end
