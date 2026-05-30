
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function txt = methodsParagraph(report)
% METHODSPARAGRAPH  One-paragraph methods prose for a single file's report.
%   txt = METHODSPARAGRAPH(report)
%
%   Returns a concise, copy-paste-ready sentence describing what the pipeline
%   did to one file (channels retained, epochs retained, ICA components
%   removed with a category breakdown when available). Shared by buildReportText
%   (the brief per-file note) and the Reports tab "Copy Methods" button.
%
%   See also: methodsParagraphAggregate, buildReportText, summarizeReports

parts = {};

ch = report.channels;
if ch.original > 0
    if ch.nRejected > 0 && ch.nInterpolated > 0
        parts{end+1} = sprintf('%d of %d channels were retained (%d removed, %d interpolated)', ...
            ch.final, ch.original, ch.nRejected, ch.nInterpolated);
    elseif ch.nRejected > 0
        parts{end+1} = sprintf('%d of %d channels were retained (%d removed)', ...
            ch.final, ch.original, ch.nRejected);
    elseif ch.nInterpolated > 0
        parts{end+1} = sprintf('%d channels were retained (%d interpolated)', ...
            ch.final, ch.nInterpolated);
    else
        parts{end+1} = sprintf('all %d channels were retained', ch.final);
    end
end

if report.trials.original > 0
    tr = report.trials;
    parts{end+1} = sprintf('%d of %d epochs were retained (%d rejected)', ...
        tr.final, tr.original, tr.rejected);
end

if report.ica.nComponents > 0
    parts{end+1} = icaSentence(report.ica);
end

if isempty(parts)
    txt = 'TMS-EEG data were preprocessed using nestapp.';
else
    txt = sprintf('TMS-EEG data were preprocessed using nestapp. %s.', ...
        strjoin(parts, '; '));
end
end

function s = icaSentence(ica)
% Brief ICA sentence with a per-category breakdown when categories exist.
if ica.nRejected == 0
    s = sprintf('ICA identified %d components, none removed', ica.nComponents);
    return
end
catStr = '';
if isfield(ica, 'categories') && any(ica.categories.nRemoved > 0)
    cats = ica.categories;
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
roundStr = '';
if isfield(ica, 'rounds') && numel(ica.rounds) > 1
    roundStr = sprintf(' over %d rounds', numel(ica.rounds));
end
s = sprintf('ICA identified %d components, of which %d were removed%s%s', ...
    ica.nComponents, ica.nRejected, roundStr, catStr);
end
