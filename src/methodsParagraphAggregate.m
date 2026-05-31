
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function txt = methodsParagraphAggregate(reports)
% METHODSPARAGRAPHAGGREGATE  Cross-file methods prose for a session.
%   txt = METHODSPARAGRAPHAGGREGATE(reports)
%
%   reports - cell array of PipelineReport structs. Returns a publication-ready
%   paragraph summarising the session as mean +/- SD across files (channels
%   retained, epochs retained, ICA components removed). Falls back to the
%   single-file paragraph when only one report is given. Used by
%   summarizeReports (session summary) and the "Copy Methods" button.
%
%   See also: methodsParagraph, summarizeReports

if iscell(reports) && isscalar(reports)
    txt = methodsParagraph(reports{1});
    return
end
N = numel(reports);

parts = {};

% Channels (over files that recorded an original count)
origCh = cellfun(@(r) r.channels.original, reports);
finCh  = cellfun(@(r) r.channels.final,    reports);
rejCh  = cellfun(@(r) r.channels.nRejected, reports);
intpCh = cellfun(@(r) r.channels.nInterpolated, reports);
hasCh  = origCh > 0;
if any(hasCh)
    chSent = sprintf('%s of %s channels were retained', ...
        meanSd(finCh(hasCh)), meanSd(origCh(hasCh)));
    extra = {};
    if any(rejCh(hasCh) > 0); extra{end+1} = sprintf('%s removed', meanSd(rejCh(hasCh))); end
    if any(intpCh(hasCh) > 0); extra{end+1} = sprintf('%s interpolated', meanSd(intpCh(hasCh))); end
    if ~isempty(extra); chSent = sprintf('%s (%s)', chSent, strjoin(extra, ', ')); end
    parts{end+1} = chSent;
end

% Trials (over epoched files)
origTr = cellfun(@(r) r.trials.original, reports);
finTr  = cellfun(@(r) r.trials.final,    reports);
rejTr  = cellfun(@(r) r.trials.rejected, reports);
hasTr  = origTr > 0;
if any(hasTr)
    parts{end+1} = sprintf('%s of %s epochs were retained (%s rejected)', ...
        meanSd(finTr(hasTr)), meanSd(origTr(hasTr)), meanSd(rejTr(hasTr)));
end

% ICA (over files where ICA ran)
nComp = cellfun(@(r) r.ica.nComponents, reports);
nRej  = cellfun(@(r) r.ica.nRejected,  reports);
hasICA = nComp > 0;
if any(hasICA)
    parts{end+1} = sprintf('ICA identified %s components, of which %s were removed', ...
        meanSd(nComp(hasICA)), meanSd(nRej(hasICA)));
end

if isempty(parts)
    txt = sprintf('Across %d files, TMS-EEG data were preprocessed using nestapp.', N);
else
    txt = sprintf('Across %d files, TMS-EEG data were preprocessed using nestapp. %s. Values are mean +/- SD across files.', ...
        N, strjoin(parts, '; '));
end
end

function s = meanSd(v)
% Compact "mean" or "mean +/- SD" depending on spread.
v = double(v(:));
if isempty(v)
    s = '0';
elseif isscalar(v) || std(v) < 1e-9
    s = sprintf('%.0f', mean(v));
else
    s = sprintf('%.0f +/- %.0f', mean(v), std(v));
end
end
