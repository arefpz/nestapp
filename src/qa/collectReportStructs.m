
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function reports = collectReportStructs(entries)
% COLLECTREPORTSTRUCTS  Extract the report struct from each non-synthetic entry.
%   entries : cell array of Reports-tab entries
%   reports : cell array of report structs (Summary and Dashboard
%             synthetic entries are skipped).
reports = {};
for k = 1:numel(entries)
    e = entries{k};
    if isfield(e, 'isSummary') && e.isSummary, continue, end
    if isfield(e, 'isDashboard') && e.isDashboard, continue, end
    if ~isfield(e, 'report'), continue, end
    reports{end+1} = e.report; %#ok<AGROW>
end
end
