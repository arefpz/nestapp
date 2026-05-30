
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function names = reportStepNames(report)
% REPORTSTEPNAMES  Names of the steps that ran for a PipelineReport.
%   names = REPORTSTEPNAMES(report) returns a 1xN cellstr of step names from
%   report.steps (the steps that actually completed), or {} when none ran.
%   Used to derive method citations from the pipeline that was applied.
%
%   See also: stepCitations, citationLines, buildReportText, summarizeReports

names = {};
if isfield(report, 'steps') && ~isempty(report.steps)
    names = cellfun(@(s) s.name, report.steps, 'UniformOutput', false);
end
end
