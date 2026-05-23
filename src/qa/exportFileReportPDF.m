function pdfPath = exportFileReportPDF(report, outputDir)
% EXPORTFILEREPORTPDF  One-PDF-per-file: text summary + checkpoint PNGs.
%   pdfPath = EXPORTFILEREPORTPDF(report, outputDir)
%
%   Page 1   : the formatted text report (from buildReportText) rendered
%              into a uifigure as monospace text.
%   Page 2+  : one page per PNG in report.quality.figures (auto QC
%              checkpoints), rendered via imread + image().
%
%   Uses exportgraphics(fig, file, 'Append', true) for multi-page
%   output. Always overwrites the PDF - re-running the pipeline on the
%   same file should produce a fresh, current artifact, not append to
%   a stale one.
%
%   Returns the absolute path to the written PDF. Throws if the text
%   page cannot render. PNG-page errors are caught per-page and a
%   placeholder ("Could not load <path>") is rendered in place so the
%   rest of the bundle still ships.

if nargin < 2 || isempty(outputDir)
    outputDir = fileparts(report.inputFile);
    if isempty(outputDir), outputDir = pwd; end
end

[~, baseName] = fileparts(report.inputFile);
if isempty(baseName), baseName = 'pipeline'; end

ts = string(report.processedAt, 'yyyyMMdd_HHmmss');
if getpref('nestapp', 'overwriteReports', false)
    pdfName = sprintf('%s_report.pdf', baseName);
else
    pdfName = sprintf('%s_report_%s.pdf', baseName, ts);
end
pdfPath = fullfile(outputDir, pdfName);

if exist(pdfPath, 'file'), delete(pdfPath); end

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
% Render the buildReportText output as monospace text in an invisible
% regular figure. exportgraphics on a regular figure captures the text;
% uifigure components would need exportapp instead.
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'inches', 'Position', [1 1 8.5 11], ...
    'PaperUnits', 'inches', 'PaperSize', [8.5 11], ...
    'PaperPosition', [0 0 8.5 11]);
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
    'Units', 'inches', 'Position', [1 1 11 8.5], ...
    'PaperUnits', 'inches', 'PaperSize', [11 8.5], ...
    'PaperPosition', [0 0 11 8.5]);
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
