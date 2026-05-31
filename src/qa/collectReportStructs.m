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
