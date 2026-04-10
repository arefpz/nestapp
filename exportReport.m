function [summaryText, matPath] = exportReport(report, outputDir)
% EXPORTREPORT  Export a PipelineReport to disk and return a formatted summary.
%
%   [summaryText, matPath] = EXPORTREPORT(report, outputDir)
%
%   Saves a .mat file containing the full report struct into outputDir and
%   returns a multi-line text summary suitable for display in the GUI.
%
%   summaryText - formatted char suitable for uitextarea or uialert
%   matPath     - full path to saved .mat file, or '' if save failed
%
%   See also: initPipelineReport, runPipeline

if nargin < 2 || isempty(outputDir)
    outputDir = fileparts(report.inputFile);
    if isempty(outputDir)
        outputDir = pwd;
    end
end

%% Build summary text
lines = {};
lines{end+1} = '=== Pipeline Report ===';
lines{end+1} = sprintf('File:      %s', report.inputFile);
lines{end+1} = sprintf('Processed: %s', datestr(report.processedAt, 'yyyy-mm-dd HH:MM:SS'));
lines{end+1} = '';

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
    nKept = report.ica.nKept;
    lines{end+1} = sprintf('  Identified: %d components', nComp);
    if nRej > 0
        hasVar   = ~isnan(report.ica.varRemoved);
        hasRange = ~isnan(report.ica.varMin);
        if hasVar && hasRange
            lines{end+1} = sprintf( ...
                '  Removed:    %d  (%.1f%% total variance, %.1f-%.1f%% per component)', ...
                nRej, report.ica.varRemoved, report.ica.varMin, report.ica.varMax);
        elseif hasVar
            lines{end+1} = sprintf('  Removed:    %d  (%.1f%% of signal variance)', ...
                nRej, report.ica.varRemoved);
        else
            lines{end+1} = sprintf('  Removed:    %d', nRej);
        end
        lines{end+1} = sprintf('  Kept:       %d', nKept);

        % Per-category breakdown (only show categories with >= 1 removed)
        if isfield(report.ica, 'categories') && any(report.ica.categories.nRemoved > 0)
            lines{end+1} = '  By category (removed):';
            cats = report.ica.categories;
            for ci = 1:numel(cats.names)
                if cats.nRemoved(ci) > 0
                    if hasVar
                        lines{end+1} = sprintf('    %-12s %d  (%.1f%% var)', ...
                            [cats.names{ci} ':'], cats.nRemoved(ci), cats.varShare(ci));
                    else
                        lines{end+1} = sprintf('    %-12s %d', ...
                            [cats.names{ci} ':'], cats.nRemoved(ci));
                    end
                end
            end
        end
    else
        lines{end+1} = sprintf('  Removed:    0  (kept all %d)', nKept);
    end
else
    lines{end+1} = '  ICA not run';
end
lines{end+1} = '';

% TEPs (reserved for M5)
lines{end+1} = 'TEP COMPONENTS';
if isfield(report.teps, 'components') && ~isempty(report.teps.components)
    tepData = report.teps.components;
    for k = 1:size(tepData, 1)
        lines{end+1} = sprintf('  %-6s  %6.1f ms  %8.4f uV  ROI: %s', ...
            tepData{k,1}, tepData{k,2}, tepData{k,3}, tepData{k,4});
    end
else
    lines{end+1} = '  TEP analysis not run';
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

% Methods summary
chStr = '';
if origCh > 0
    chStr = sprintf('%d/%d channels retained', finCh, origCh);
elseif finCh > 0
    chStr = sprintf('%d channels', finCh);
end
trStr = '';
if report.trials.original > 0
    trStr = sprintf('%d/%d trials retained', report.trials.final, report.trials.original);
end
icaStr = '';
if report.ica.nRejected > 0
    icaStr = sprintf('%d ICA component(s) rejected', report.ica.nRejected);
end
parts = {chStr, trStr, icaStr};
parts = parts(~cellfun(@isempty, parts));
lines{end+1} = 'METHODS SUMMARY';
if ~isempty(parts)
    lines{end+1} = sprintf('  %s.', strjoin(parts, ', '));
else
    lines{end+1} = '  No metrics available.';
end

summaryText = strjoin(lines, newline);

%% Save .mat file
[~, baseName] = fileparts(report.inputFile);
if isempty(baseName)
    baseName = 'pipeline';
end
timestamp   = datestr(report.processedAt, 'yyyymmdd_HHMMSS');
matFileName = sprintf('%s_report_%s.mat', baseName, timestamp);
matPath     = fullfile(outputDir, matFileName);

try
    pipelineReport = report;
    save(matPath, 'pipelineReport');
catch
    matPath = '';
end
end
