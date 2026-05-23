function [summaryText, matPath] = exportReport(report, outputDir)
% EXPORTREPORT  Export a PipelineReport to disk and return a formatted summary.
%
%   [summaryText, matPath] = EXPORTREPORT(report, outputDir)
%
%   Saves a .mat file containing the full report struct into outputDir and
%   returns a multi-line text summary suitable for display in the GUI.
%
%   summaryText - formatted char suitable for uitextarea or uialert
%   matPath     - full path to saved .mat file, or '' if save failed
%
%   See also: initPipelineReport, runPipelineCore

if nargin < 2 || isempty(outputDir)
    % Use the user-specified report folder if set, otherwise data folder
    prefFolder = getpref('nestapp', 'reportFolder', '');
    if ~isempty(prefFolder) && isfolder(prefFolder)
        outputDir = prefFolder;
    else
        outputDir = fileparts(report.inputFile);
        if isempty(outputDir)
            outputDir = pwd;
        end
    end
end

summaryText = buildReportText(report);


%% Save .mat file
[~, baseName] = fileparts(report.inputFile);
if isempty(baseName)
    baseName = 'pipeline';
end
if getpref('nestapp', 'overwriteReports', false)
    matFileName = sprintf('%s_report.mat', baseName);
else
    timestamp   = string(report.processedAt, 'yyyyMMdd_HHmmss');
    matFileName = sprintf('%s_report_%s.mat', baseName, timestamp);
end
matPath = fullfile(outputDir, matFileName);

try
    pipelineReport = report;
    save(matPath, 'pipelineReport');
catch
    matPath = '';
end
end
