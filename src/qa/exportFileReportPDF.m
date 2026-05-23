function pdfPath = exportFileReportPDF(report, target)
% EXPORTFILEREPORTPDF  One-PDF-per-file: text summary + checkpoint PNGs.
%   pdfPath = EXPORTFILEREPORTPDF(report, target)
%
%   target may be:
%     - a batchCtx struct from buildBatchContext (normal pipeline path)
%     - a folder path (legacy callers and tests): PDF named via
%       reportArtifactName(report, 'pdf')
%     - a full *.pdf path: written there verbatim (used by the GUI's
%       Save As dialog so the user's chosen filename is honored).
%     - omitted: written next to the input file (last-resort fallback).
%
%   Page 1   : the formatted text report (from buildReportText) rendered
%              as monospace text in a regular figure.
%   Page 2+  : one page per PNG in report.quality.figures (auto QC
%              checkpoints), rendered via imread + image().
%
%   Uses exportgraphics(fig, file, 'Append', true) for multi-page
%   output. The first call uses 'Append', false which overwrites any
%   existing file - re-running the pipeline on the same file produces
%   a fresh artifact, not an appended one.
%
%   Returns the absolute path to the written PDF. Throws if the text
%   page cannot render. PNG-page errors are caught per-page and a
%   placeholder ("Could not load <path>") is rendered in its place so
%   the rest of the bundle still ships.

[~, stem] = fileparts(report.inputFile);
if nargin < 2 || isempty(target)
    target = fileparts(report.inputFile);
    if isempty(target), target = pwd; end
end

if isstruct(target)
    pdfPath = fullfile(outputPaths(target, 'reports', stem), ...
        reportArtifactName(report, 'pdf'));
elseif endsWith(lower(target), '.pdf')
    pdfPath = target;
else
    pdfPath = fullfile(target, reportArtifactName(report, 'pdf'));
end

% Page 1: text report.
textFig = renderTextPage(report);
try
    exportgraphics(textFig, pdfPath, 'Append', false);
catch err
    closeFigSafely(textFig);
    rethrow(err);
end
closeFigSafely(textFig);

% Subsequent pages: one per QC PNG.
if isfield(report, 'quality') && isfield(report.quality, 'figures')
    for k = 1:numel(report.quality.figures)
        pngPath = report.quality.figures{k};
        imgFig = renderImagePage(pngPath);
        try
            exportgraphics(imgFig, pdfPath, 'Append', true);
        catch err
            warning('exportFileReportPDF:appendFailed', ...
                'Failed to append page for %s: %s', pngPath, err.message);
        end
        closeFigSafely(imgFig);
    end
end
end

% -- page renderers --------------------------------------------------------

function fig = renderTextPage(report)
% Regular figure (not uifigure) so exportgraphics captures the text;
% uifigure components would require exportapp instead.
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'inches', 'Position', [1 1 8.5 11]);
ax = axes(fig, 'Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
xlim(ax, [0 1]); ylim(ax, [0 1]);
text(ax, 0.0, 1.0, buildReportText(report), ...
    'FontName', 'Courier New', 'FontSize', 8, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'Interpreter', 'none');
end

function fig = renderImagePage(pngPath)
% Render one PNG as a full-page image. On read failure, render a
% placeholder so the bundle still has a page in the right slot.
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'inches', 'Position', [1 1 11 8.5]);
ax = axes(fig, 'Position', [0.02 0.02 0.96 0.96]);
try
    img = imread(pngPath);
    image(ax, img);
    axis(ax, 'image');
    axis(ax, 'off');
    [~, base, ext] = fileparts(pngPath);
    title(ax, [base, ext], 'Interpreter', 'none');
catch err
    axis(ax, [0 1 0 1]);
    text(ax, 0.5, 0.5, ...
        sprintf('Could not load:\n%s\n\n%s', pngPath, err.message), ...
        'HorizontalAlignment', 'center', 'FontSize', 14, ...
        'Interpreter', 'none');
    axis(ax, 'off');
end
end

function closeFigSafely(fig)
if isvalid(fig)
    try
        close(fig, 'force');
    catch
    end
end
end
