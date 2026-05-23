function [summaryText, matPath] = exportReport(report, target)
% EXPORTREPORT  Export a PipelineReport to disk and return a formatted summary.
%
%   [summaryText, matPath] = EXPORTREPORT(report, target)
%
%   target may be:
%     - a batchCtx struct (from buildBatchContext) - normal pipeline path
%     - a char folder path - legacy callers and tests
%     - omitted - .mat lands next to the input file, with `pwd` as a last resort
%
%   Saves a .mat file containing the full report struct and returns a
%   multi-line text summary suitable for display in the GUI.
%
%   summaryText - formatted char suitable for uitextarea or uialert
%   matPath     - full path to saved .mat file, or '' if save failed
%
%   See also: initPipelineReport, runPipelineCore, outputPaths

[~, stem] = fileparts(report.inputFile);
if nargin < 2 || isempty(target)
    outputDir = fileparts(report.inputFile);
    if isempty(outputDir), outputDir = pwd; end
elseif isstruct(target)
    outputDir = outputPaths(target, 'reports', stem);
else
    outputDir = target;
end

summaryText = buildReportText(report);

matPath = fullfile(outputDir, reportArtifactName(report, 'mat'));

try
    pipelineReport = report;
    save(matPath, 'pipelineReport');
catch
    matPath = '';
end
end
