
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function lines = citationLines(pipelineName)
% CITATIONLINES  Format a built-in template's citation as report text lines.
%   lines = CITATIONLINES(pipelineName)
%
%   Returns a cellstr CITATION block (a header plus 'Cite as:' / 'DOI:' / notes)
%   for the built-in template named pipelineName, or {} when the name has no
%   citation registered. Wraps templateCitation so buildReportText (per-file
%   report + PDF cover) and summarizeReports (session summary) share one source
%   of truth with the batch-log citation printed by runPipelineCore.
%
%   See also: templateCitation, buildReportText, summarizeReports

lines = {};
if nargin < 1 || isempty(pipelineName)
    return
end

c = templateCitation(pipelineName);
if isempty(c.reference)
    return
end

lines{end+1} = 'CITATION';
lines{end+1} = sprintf('  Cite as: %s', c.reference);
if ~isempty(c.doi)
    lines{end+1} = sprintf('  DOI:     %s', c.doi);
end
if ~isempty(c.notes)
    lines{end+1} = sprintf('  %s', c.notes);
end
end
