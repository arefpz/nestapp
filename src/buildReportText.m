
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function summaryText = buildReportText(report)
% BUILDREPORTTEXT  Pure text formatter for a PipelineReport struct.
%   summaryText = BUILDREPORTTEXT(report)
%   Used by exportReport (for the .mat side + UI display) and by
%   exportFileReportPDF (for the cover page) so both share one
%   canonical formatter.

lines = {};
lines{end+1} = '=== Pipeline Report ===';
lines{end+1} = sprintf('File:      %s', report.inputFile);
lines{end+1} = sprintf('Processed: %s', string(report.processedAt, 'yyyy-MM-dd HH:mm:ss'));

% Quality reports - single folder path (all PNGs share a parent).
if isfield(report, 'quality') && isfield(report.quality, 'figures') ...
        && ~isempty(report.quality.figures)
    qcFolder = fileparts(report.quality.figures{1});
    lines{end+1} = sprintf('QC images: %s  (%d files)', ...
        qcFolder, numel(report.quality.figures));
end

% Quality Gate verdict - one summary line for the whole file.
if isfield(report, 'quality') && isfield(report.quality, 'worstVerdict') ...
        && ~strcmp(report.quality.worstVerdict, 'NotChecked')
    lines{end+1} = sprintf('Quality:   %s', report.quality.worstVerdict);
end

lines{end+1} = '';

% Quality Gate details - per-gate verdict and reasons (if any gates ran).
if isfield(report, 'quality') && isfield(report.quality, 'gates') ...
        && ~isempty(report.quality.gates)
    lines{end+1} = 'QUALITY GATES';
    for gi = 1:numel(report.quality.gates)
        g = report.quality.gates{gi};
        stepIx = '';
        if isfield(g, 'stepIndex') && ~isempty(g.stepIndex)
            stepIx = sprintf(' (step %d)', g.stepIndex);
        end
        lines{end+1} = sprintf('  %s%s: %s', g.label, stepIx, g.verdict); %#ok<AGROW>
        if isfield(g, 'reasons') && ~isempty(g.reasons)
            for ri = 1:numel(g.reasons)
                lines{end+1} = sprintf('    - %s', g.reasons{ri}); %#ok<AGROW>
            end
        end
    end
    lines{end+1} = '';
end

% Channels
lines{end+1} = 'CHANNELS';
origCh = report.channels.original;
rejCh  = report.channels.nRejected;
intpCh = report.channels.nInterpolated;
finCh  = report.channels.final;
lines{end+1} = sprintf('  Original:     %d', origCh);
if rejCh > 0
    lines{end+1} = sprintf('  Rejected:     %d', rejCh);
end
if intpCh > 0
    lines{end+1} = sprintf('  Interpolated: %d', intpCh);
end
lines{end+1} = sprintf('  Final:        %d', finCh);
lines{end+1} = '';

% Trials
lines{end+1} = 'TRIALS';
if report.trials.original > 0
    origTr = report.trials.original;
    rejTr  = report.trials.rejected;
    finTr  = report.trials.final;
    lines{end+1} = sprintf('  Original:  %d', origTr);
    if rejTr > 0
        lines{end+1} = sprintf('  Rejected:  %d', rejTr);
    end
    lines{end+1} = sprintf('  Final:     %d', finTr);
else
    lines{end+1} = '  Not epoched (continuous data)';
end
lines{end+1} = '';

% ICA
lines{end+1} = 'ICA';
if report.ica.nComponents > 0
    nComp = report.ica.nComponents;
    nRej  = report.ica.nRejected;
    % nKept added in M3; fall back for reports saved before that field existed.
    if isfield(report.ica, 'nKept')
        nKept = report.ica.nKept;
    else
        nKept = nComp - nRej;
    end
    lines{end+1} = sprintf('  Identified: %d components', nComp);
    if nRej > 0
        hasVar   = ~isnan(report.ica.varRemoved);
        hasRange = ~isnan(report.ica.varMin);
        multiRound = isfield(report.ica, 'rounds') && numel(report.ica.rounds) > 1;

        if multiRound
            lines{end+1} = sprintf('  Removed:    %d total (%d rounds)', ...
                nRej, numel(report.ica.rounds));
        elseif hasVar && hasRange
            lines{end+1} = sprintf( ...
                '  Removed:    %d  (%.1f%% ICA variance, %.1f-%.1f%% per component)', ...
                nRej, report.ica.varRemoved, report.ica.varMin, report.ica.varMax);
        elseif hasVar
            lines{end+1} = sprintf('  Removed:    %d  (%.1f%% ICA variance)', ...
                nRej, report.ica.varRemoved);
        else
            lines{end+1} = sprintf('  Removed:    %d', nRej);
        end
        lines{end+1} = sprintf('  Kept:       %d', nKept);

        % Per-category summary (totals across all rounds)
        if isfield(report.ica, 'categories') && any(report.ica.categories.nRemoved > 0)
            lines{end+1} = '  By category (all rounds):';
            lines = appendCategoryLines(lines, report.ica.categories, hasVar && ~multiRound);
        end

        % Per-round detail for multi-round TESA
        if multiRound
            for ri = 1:numel(report.ica.rounds)
                rnd = report.ica.rounds{ri};
                rndHasVar = ~isnan(rnd.varRemoved);
                if rndHasVar
                    lines{end+1} = sprintf('  Round %d: %d components, %d removed (%.1f%% ICA var, %.1f-%.1f%% per comp)', ...
                        ri, rnd.nComponents, rnd.nRejected, rnd.varRemoved, rnd.varMin, rnd.varMax);
                else
                    lines{end+1} = sprintf('  Round %d: %d components, %d removed', ...
                        ri, rnd.nComponents, rnd.nRejected);
                end
                lines = appendCategoryLines(lines, rnd.categories, rndHasVar);
            end
        end
    else
        lines{end+1} = sprintf('  Removed:    0  (kept all %d)', nKept);
    end
else
    lines{end+1} = '  ICA not run';
end
lines{end+1} = '';

% Steps run
lines{end+1} = 'STEPS RUN';
for k = 1:numel(report.steps)
    rec = report.steps{k};
    chanNote = '';
    if rec.chansAfter ~= rec.chansBefore
        chanNote = sprintf('  [%d -> %d ch]', rec.chansBefore, rec.chansAfter);
    end
    trialNote = '';
    if rec.trialsAfter ~= rec.trialsBefore && rec.trialsBefore > 1
        trialNote = sprintf('  [%d -> %d trials]', rec.trialsBefore, rec.trialsAfter);
    elseif rec.trialsAfter > 1 && rec.trialsBefore <= 1
        trialNote = sprintf('  [-> %d trials]', rec.trialsAfter);
    end
    lines{end+1} = sprintf('  %2d. %-35s %.1fs%s%s', ...
        k, rec.name, rec.duration, chanNote, trialNote);
end
lines{end+1} = '';

% Methods note - one concise sentence per file. The full cross-file methods
% prose (mean +/- SD across files) lives in the session summary; see
% summarizeReports / methodsParagraphAggregate.
lines{end+1} = 'METHODS';
lines{end+1} = ['  ', methodsParagraph(report)];

% Citation - references for the methods this file's pipeline actually used,
% derived from the steps that ran. Same block runPipelineCore prints to the
% batch log and the session summary shows; rendered here so each per-file
% report (and its PDF) carries it.
citeLines = citationLines(reportStepNames(report));
if ~isempty(citeLines)
    lines{end+1} = '';
    lines = [lines, citeLines];
end

summaryText = strjoin(lines, newline);
end
function lines = appendCategoryLines(lines, cats, showVar)
for ci = 1:numel(cats.names)
    if cats.nRemoved(ci) > 0
        if showVar
            lines{end+1} = sprintf('    %-12s %d  (%.1f%% ICA var)', ...
                [cats.names{ci} ':'], cats.nRemoved(ci), cats.varShare(ci));
        else
            lines{end+1} = sprintf('    %-12s %d', ...
                [cats.names{ci} ':'], cats.nRemoved(ci));
        end
    end
end
end