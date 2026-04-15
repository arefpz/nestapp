function summaryText = summarizeReports(reports)
% SUMMARIZEREPORTS  Build a cross-file summary from a cell array of PipelineReport structs.
%
%   summaryText = SUMMARIZEREPORTS(reports)
%
%   reports - 1×N cell array of structs returned by initPipelineReport and
%             populated by runPipeline. N must be >= 2.
%
%   Returns a formatted char suitable for display in the pipeline report dialog
%   above the individual per-file reports.
%
%   See also: initPipelineReport, exportReport, runPipeline

N = numel(reports);
lines = {};
lines{end+1} = sprintf('=== PIPELINE SUMMARY  (%d files) ===', N);
lines{end+1} = '';

%% Channels
origCh  = cellfun(@(r) r.channels.original,      reports);
rejCh   = cellfun(@(r) r.channels.nRejected,     reports);
intpCh  = cellfun(@(r) r.channels.nInterpolated, reports);
finalCh = cellfun(@(r) r.channels.final,          reports);

lines{end+1} = 'CHANNELS';
lines{end+1} = sprintf('  Original:     %s', fmtStat(origCh));
if any(rejCh > 0)
    lines{end+1} = sprintf('  Rejected:     %s', fmtStat(rejCh));
end
if any(intpCh > 0)
    lines{end+1} = sprintf('  Interpolated: %s', fmtStat(intpCh));
end
lines{end+1} = sprintf('  Final:        %s', fmtStat(finalCh));
lines{end+1} = '';

%% Trials
origTr  = cellfun(@(r) r.trials.original, reports);
rejTr   = cellfun(@(r) r.trials.rejected, reports);
finalTr = cellfun(@(r) r.trials.final,    reports);

epoched = origTr > 0;
if any(epoched)
    lines{end+1} = 'TRIALS';
    if ~all(epoched)
        lines{end+1} = sprintf('  %d of %d files epoched', sum(epoched), N);
    end
    ep = find(epoched);
    lines{end+1} = sprintf('  Original:  %s', fmtStat(origTr(ep)));
    if any(rejTr(ep) > 0)
        pctRej = rejTr(ep) ./ origTr(ep) * 100;
        lines{end+1} = sprintf('  Rejected:  %s  (%.1f%%±%.1f%% of trials)', ...
            fmtStat(rejTr(ep)), mean(pctRej), std(pctRej));
    end
    lines{end+1} = sprintf('  Final:     %s', fmtStat(finalTr(ep)));
    lines{end+1} = '';
end

%% ICA
nComp  = cellfun(@(r) r.ica.nComponents, reports);
nRej   = cellfun(@(r) r.ica.nRejected,  reports);
varRem = cellfun(@(r) r.ica.varRemoved,  reports);

hasICA = nComp > 0;
if any(hasICA)
    lines{end+1} = 'ICA';
    icaIdx = find(hasICA);
    lines{end+1} = sprintf('  Identified: %s components', fmtStat(nComp(icaIdx)));
    lines{end+1} = sprintf('  Removed:    %s', fmtStat(nRej(icaIdx)));

    hasVar = ~isnan(varRem) & hasICA;
    if any(hasVar)
        lines{end+1} = sprintf('  ICA var removed: %s%%', fmtStat(varRem(hasVar)));
    end

    % Per-category totals (only if all ICA files share the same category scheme)
    icaReports = reports(hasICA);
    catNames = icaReports{1}.ica.categories.names;
    allSameCats = all(cellfun(@(r) isequal(r.ica.categories.names, catNames), icaReports));
    anyRejected = any(cellfun(@(r) any(r.ica.categories.nRemoved > 0), icaReports));

    if allSameCats && anyRejected
        lines{end+1} = '  By category (mean per file):';
        nCats = numel(catNames);
        for ci = 1:nCats
            counts = cellfun(@(r) r.ica.categories.nRemoved(ci), icaReports);
            if any(counts > 0)
                lines{end+1} = sprintf('    %-12s %s', [catNames{ci} ':'], fmtStat(counts));
            end
        end
    end
    lines{end+1} = '';
end

%% TEP Quality
hasQ = cellfun(@(r) isfield(r.teps,'composite') && ~isnan(r.teps.composite), reports);
if any(hasQ)
    qReports = reports(hasQ);
    tstat = cellfun(@(r) r.teps.tstat,     qReports);
    sh    = cellfun(@(r) r.teps.splitHalf, qReports);
    snr   = cellfun(@(r) r.teps.snr,       qReports);
    comp  = cellfun(@(r) r.teps.composite, qReports);
    lines{end+1} = 'TEP QUALITY';
    if any(~isnan(tstat))
        lines{end+1} = sprintf('  T-statistic:  %s', fmtStat(tstat(~isnan(tstat))));
    end
    if any(~isnan(sh))
        lines{end+1} = sprintf('  Split-half r: %s', fmtStat(sh(~isnan(sh))));
    end
    if any(~isnan(snr))
        lines{end+1} = sprintf('  SNR:          %s', fmtStat(snr(~isnan(snr))));
    end
    lines{end+1} = sprintf('  Composite:    %s', fmtStat(comp));
    lines{end+1} = '';
end

summaryText = strjoin(lines, newline);
end

%% ---- helpers ---------------------------------------------------------------

function s = fmtStat(v)
% Format a numeric vector as "mean ± SD" or just the value when all equal.
v = double(v(:));
if isscalar(v)
    s = sprintf('%.1f', v);
elseif std(v) < 1e-9
    s = sprintf('%.1f (all equal)', mean(v));
else
    s = sprintf('%.1f ± %.1f  [%g – %g]', mean(v), std(v), min(v), max(v));
end
end
