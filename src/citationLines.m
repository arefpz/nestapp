
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function lines = citationLines(stepNames)
% CITATIONLINES  Format a CITATION block for the methods a pipeline used.
%   lines = CITATIONLINES(stepNames)
%
%   Returns a cellstr CITATION block listing one reference (+ DOI) per method
%   whose trigger steps appear in stepNames, or {} when no cited method was
%   used. Derives the list from the actual steps via stepCitations, so the
%   report cites exactly what ran - no template-name guessing, no conditional
%   "if you also used X" notes.
%
%   See also: stepCitations, buildReportText, summarizeReports

lines = {};
cites = stepCitations(stepNames);
if isempty(cites)
    return
end

lines{end+1} = 'CITATION';
lines{end+1} = '  Methods used in this pipeline - please cite:';
for i = 1:numel(cites)
    lines{end+1} = sprintf('    %s', cites(i).reference); %#ok<AGROW>
    if ~isempty(cites(i).doi)
        lines{end+1} = sprintf('      DOI: %s', cites(i).doi); %#ok<AGROW>
    end
end
end
