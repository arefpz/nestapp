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
    % Use the user-specified report folder if set, otherwise data folder
    prefFolder = getpref('nestapp', 'reportFolder', '');
    if ~isempty(prefFolder) && isfolder(prefFolder)
        outputDir = prefFolder;
    else
        outputDir = fileparts(report.inputFile);
        if isempty(outputDir)
            outputDir = pwd;
        end
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
            cats = report.ica.categories;
            for ci = 1:numel(cats.names)
                if cats.nRemoved(ci) > 0
                    if hasVar && ~multiRound
                        lines{end+1} = sprintf('    %-12s %d  (%.1f%% ICA var)', ...
                            [cats.names{ci} ':'], cats.nRemoved(ci), cats.varShare(ci));
                    else
                        lines{end+1} = sprintf('    %-12s %d', ...
                            [cats.names{ci} ':'], cats.nRemoved(ci));
                    end
                end
            end
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
                cats = rnd.categories;
                for ci = 1:numel(cats.names)
                    if cats.nRemoved(ci) > 0
                        if rndHasVar
                            lines{end+1} = sprintf('    %-12s %d  (%.1f%% ICA var)', ...
                                [cats.names{ci} ':'], cats.nRemoved(ci), cats.varShare(ci));
                        else
                            lines{end+1} = sprintf('    %-12s %d', ...
                                [cats.names{ci} ':'], cats.nRemoved(ci));
                        end
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

% TEP quality vector
if isstruct(report.teps) && isfield(report.teps, 'version')
    q = report.teps;
    lines{end+1} = sprintf('TEP QUALITY  [metric %s]', q.version);
    if ~isempty(q.baselineWarning)
        lines{end+1} = sprintf('  WARNING: %s', q.baselineWarning);
    end
    lines = appendAxis(lines, q.retention,         'Retention',          '');
    lines = appendAxis(lines, q.artifactReduction,  'Artifact reduction', '');
    lines = appendAxis(lines, q.bgRestoration,      'Background',         '');
    lines = appendAxis(lines, q.reproducibility,    'Reproducibility',    '');
    lines = appendAxis(lines, q.aepLikeness,        'AEP-likeness',       '[EXPERIMENTAL]');
    lines{end+1} = '';
end

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

% Methods summary — prose suitable for copy-paste into a methods section
lines{end+1} = 'METHODS SUMMARY';
methLines = {};

% Channels sentence
nBad  = report.channels.nRejected;
nIntp = report.channels.nInterpolated;
if origCh > 0
    if nBad == 0 && nIntp == 0
        methLines{end+1} = sprintf('  %d channels recorded; none rejected.', origCh);
    elseif nBad > 0 && nIntp == 0
        methLines{end+1} = sprintf( ...
            '  Of %d channels, %d were identified as bad and removed (%.0f%% retained).', ...
            origCh, nBad, 100*(origCh-nBad)/origCh);
    elseif nBad > 0 && nIntp > 0
        methLines{end+1} = sprintf( ...
            ['  Of %d channels, %d were identified as bad and removed; ' ...
             '%d were subsequently interpolated (%d channels for analysis).'], ...
            origCh, nBad, nIntp, finCh);
    else
        methLines{end+1} = sprintf( ...
            '  %d channels recorded; %d interpolated (%d for analysis).', ...
            origCh, nIntp, finCh);
    end
elseif finCh > 0
    methLines{end+1} = sprintf('  %d channels.', finCh);
end

% Trials sentence
if report.trials.original > 0
    origTr_ = report.trials.original;
    rejTr_  = report.trials.rejected;
    finTr_  = report.trials.final;
    pctRet  = 100 * finTr_ / origTr_;
    if rejTr_ == 0
        methLines{end+1} = sprintf('  %d trials; none rejected.', origTr_);
    else
        methLines{end+1} = sprintf( ...
            '  Of %d trials, %d were rejected (%d retained, %.0f%%).', ...
            origTr_, rejTr_, finTr_, pctRet);
    end
end

% ICA sentence
if report.ica.nComponents > 0
    nComp      = report.ica.nComponents;
    nRej       = report.ica.nRejected;
    hasVar     = ~isnan(report.ica.varRemoved);
    multiRound = isfield(report.ica, 'rounds') && numel(report.ica.rounds) > 1;
    hasCats    = isfield(report.ica, 'categories') && any(report.ica.categories.nRemoved > 0);

    if nRej == 0
        methLines{end+1} = sprintf( ...
            '  ICA decomposition identified %d components; none rejected.', nComp);
    else
        if multiRound
            icaBase = sprintf( ...
                '  ICA decomposition (%d rounds) identified %d components; %d were removed', ...
                numel(report.ica.rounds), nComp, nRej);
        else
            icaBase = sprintf( ...
                '  ICA decomposition identified %d components; %d were removed', nComp, nRej);
        end

        catStr = '';
        if hasCats
            cats = report.ica.categories;
            catParts = {};
            for ci = 1:numel(cats.names)
                if cats.nRemoved(ci) > 0
                    catParts{end+1} = sprintf('%s: %d', cats.names{ci}, cats.nRemoved(ci)); %#ok<AGROW>
                end
            end
            if ~isempty(catParts)
                catStr = sprintf(' (%s)', strjoin(catParts, ', '));
            end
        end

        if hasVar && ~multiRound
            methLines{end+1} = sprintf('%s%s, accounting for %.1f%% of ICA variance.', ...
                icaBase, catStr, report.ica.varRemoved);
        else
            methLines{end+1} = sprintf('%s%s.', icaBase, catStr);
        end
    end
end

if isempty(methLines)
    lines{end+1} = '  No metrics available.';
else
    lines = [lines, methLines];
end

summaryText = strjoin(lines, newline);

%% Save .mat file
[~, baseName] = fileparts(report.inputFile);
if isempty(baseName)
    baseName = 'pipeline';
end
if getpref('nestapp', 'overwriteReports', false)
    matFileName = sprintf('%s_report.mat', baseName);
else
    timestamp   = datestr(report.processedAt, 'yyyymmdd_HHMMSS');
    matFileName = sprintf('%s_report_%s.mat', baseName, timestamp);
end
matPath = fullfile(outputDir, matFileName);

try
    pipelineReport = report;
    save(matPath, 'pipelineReport');
catch
    matPath = '';
end
end

%% ---- helpers ---------------------------------------------------------------

function lines = appendAxis(lines, ax, label, tag)
% Format one quality axis as a single report line.
if isempty(tag)
    prefix = sprintf('  %-22s', label);
else
    prefix = sprintf('  %-22s', [label ' ' tag]);
end
if isnan(ax.value)
    lines{end+1} = sprintf('%s  —  (%s)', prefix, ax.status);
    if ~isempty(ax.status_detail)
        lines{end+1} = sprintf('    %s', ax.status_detail);
    end
else
    lines{end+1} = sprintf('%s  %.2f  %s', prefix, ax.value, ax.interpretation);
    if ~strcmp(ax.status, 'ok') && ~isempty(ax.status_detail)
        lines{end+1} = sprintf('    [%s] %s', upper(ax.status), ax.status_detail);
    end
end
end
