function tf = anyReportHasGates(entries)
% ANYREPORTHASGATES  True if any Reports-tab entry carries Quality Gate data.
%   entries : cell array of Reports-tab entries (mix of session-report,
%             loaded-report, and synthetic Session Summary entries).
%   Returns true if at least one non-summary entry has a non-empty
%   report.quality.gates cell. Used to decide whether the listbox
%   should append a "Session Quality Dashboard" entry.
tf = false;
for k = 1:numel(entries)
    e = entries{k};
    if isfield(e, 'isSummary') && e.isSummary, continue, end
    if isfield(e, 'isDashboard') && e.isDashboard, continue, end
    if ~isfield(e, 'report') || ~isstruct(e.report), continue, end
    r = e.report;
    if isfield(r, 'quality') && isfield(r.quality, 'gates') ...
            && ~isempty(r.quality.gates)
        tf = true; return
    end
end
end
