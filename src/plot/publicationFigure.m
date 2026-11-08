% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [fig, info] = publicationFigure(drawFcn, opts)
% PUBLICATIONFIGURE  Compose a figure at journal size and export it at print DPI.
%   [fig, info] = PUBLICATIONFIGURE(drawFcn, opts)
%
%   drawFcn  @(parent, axesFcn) that draws into parent, a pixel-sized uipanel
%            filling the figure, minting its axes with axesFcn. Exactly what
%            the Explore tab draws with, so an exported figure is composed by
%            the same code that draws the canvas - not copied off the screen,
%            and not a second layout that can drift from the first.
%
%            axesFcn hands back CLASSIC axes, and that is not a detail. Every
%            export path in MATLAB - exportgraphics, print, saveas - silently
%            omits UI components: a figure of uiaxes exports with the plots
%            missing and only a warning on the command line to say so. A
%            drawing API that takes its axes maker as an argument is what makes
%            one layout serve both the screen and the page.
%
%   opts:
%     .width      'single' (89 mm) | 'double' (183 mm) | width in mm
%     .height     height in mm, default 0.62 x width
%     .dpi        export resolution, default 600
%     .fontSize   point size applied to every axes. Default follows the width -
%                 8 pt at double column, down to a 5 pt floor at single - because
%                 a fixed size that reads well at 183 mm collides with itself at
%                 89 mm, and most journals set their floor at 5 to 7 pt
%     .file       write here and close the figure; omit to leave it open
%     .title      figure name
%     .provenance struct of name/value pairs stamped under the figure and
%                 stored in fig.UserData
%
%   info returns .widthMm, .heightMm, .dpi, .pixels and .file.
%
%   Why compose fresh at the final size rather than resize afterwards: text does
%   not scale with the axes. A figure laid out at 900 px and squeezed to 89 mm
%   keeps its 10 pt labels, which then cover the data; laid out at the size it
%   will be printed, what you see is the proportion you get. This is also why
%   the figure is built at the true physical width - screen pixels per inch is
%   read from the root rather than assumed, so 89 mm is 89 mm on a HiDPI
%   display too.
%
%   Saving goes through print, not exportgraphics, for one reason: exportgraphics
%   crops to the content's bounding box. That is the right default for dropping
%   a plot into a document, and the wrong one here - a figure asked for at 89 mm
%   came out at 80 mm, which is exactly the requirement a journal checks. print
%   with an explicit PaperPosition writes the size that was asked for.
%
%   Provenance is stamped INTO the image, not just into UserData. A figure
%   leaves the app and travels: into a slide, an email, a manuscript draft, and
%   six months later nobody can say which cohort or which pipeline it came
%   from. A footer line costs 3 mm and answers that.
%
%   See also: exploreResults, drawTEPTopo, exportgraphics

if nargin < 2; opts = struct(); end
opts = fillDefaults(opts, struct('width', 'double', 'height', [], ...
    'dpi', 600, 'fontSize', [], 'file', '', 'title', 'nestapp figure', ...
    'provenance', struct()));
if isempty(opts.dpi); opts.dpi = 600; end

% The millimetre lookup, the 0.62 default height and the width-derived type
% size are all pure arithmetic, and they used to be reachable only by composing
% a figure - nine of them, in the test that checked this. publicationFigureSize
% owns them now and is tested without any graphics at all.
sz            = publicationFigureSize(opts);
widthMm       = sz.widthMm;
heightMm      = sz.heightMm;
opts.fontSize = sz.fontSize;

ppi = get(groot, 'ScreenPixelsPerInch');
px  = round([widthMm heightMm] / 25.4 * ppi);

% Scaled with the type: a fixed 11 px strip is 2% of a double-column figure and
% 16% of a single-column one, where it would crowd out the plot it is supposed
% to be captioning. One blank line plus one per line of text.
footerPts = max(5, opts.fontSize - 1);
lineH     = ceil(footerPts * 1.6);
lines     = provenanceLines(opts.provenance);
footerPx  = 0;
if ~isempty(lines)
    footerPx = lineH * (numel(lines) + 1);
end

% The figure IS the drawing parent - no uipanel holder. R2026a's print rejects
% any figure containing a uipanel as "UI components are not supported", so the
% holder that made the on-screen popout convenient is exactly what would stop
% the file from ever being written. A figure has the .Position the drawers read,
% so it substitutes directly.
%
% It is created at the DRAWING height and grown afterwards, which is how the
% footer gets its strip without every drawer having to know about it. Pixel-
% positioned children are anchored bottom-left, so growing the height moves
% nothing; shifting them up by the footer height then clears the strip.
fig = outputFigure('Name', opts.title, 'NumberTitle', 'off', ...
             'Units', 'pixels', 'Position', [80 80 px(1) px(2) - footerPx], ...
             'MenuBar', 'figure', 'ToolBar', 'figure');
fig.UserData = opts.provenance;

drawFcn(fig, @classicAxes);
applyFontSize(fig, opts.fontSize);

if footerPx > 0
    fig.Position(4) = px(2);
    shiftUp(fig, footerPx);
    % annotation, not uilabel: a uilabel is a UI component and would be dropped
    % by every export path, taking the provenance with it - the one thing on the
    % page that most has to survive the trip.
    for k = 1:numel(lines)
        annotation(fig, 'textbox', ...
            [8 / px(1), (footerPx - lineH * k) / px(2), 1 - 16 / px(1), lineH / px(2)], ...
            'String', lines{k}, 'EdgeColor', 'none', 'Margin', 0, ...
            'FontSize', footerPts, 'Color', [0.45 0.48 0.53], ...
            'VerticalAlignment', 'middle', 'Interpreter', 'none');
    end
end

info = struct('widthMm', widthMm, 'heightMm', heightMm, 'dpi', opts.dpi, ...
              'pixels', px, 'file', opts.file);

if ~isempty(opts.file)
    drawnow;
    printAtSize(fig, opts.file, widthMm, heightMm, opts.dpi);
    delete(fig);
    fig = gobjects(0);
end
end

% ── helpers ─────────────────────────────────────────────────────────────────

function ax = classicAxes(parent, pos)
% A classic axes positioned in pixels, matching the uiaxes call it replaces.
ax = axes('Parent', parent, 'Units', 'pixels', 'Position', pos);
end

function shiftUp(fig, dy)
% Move every pixel-positioned axes up, to clear the footer strip. Colorbars and
% legends are attached to their axes and follow on their own.
kids = findall(fig, 'Type', 'axes');
for k = 1:numel(kids)
    if strcmp(kids(k).Units, 'pixels')
        kids(k).Position(2) = kids(k).Position(2) + dy;
    end
end
end

function printAtSize(fig, file, widthMm, heightMm, dpi)
% Exact physical size, and a format chosen from the extension. Vector formats
% ignore the resolution for the drawing itself, which is the point of them.
cm = [widthMm heightMm] / 10;
fig.PaperUnits        = 'centimeters';
fig.PaperSize         = cm;
fig.PaperPosition     = [0 0 cm];
fig.PaperPositionMode = 'manual';
fig.InvertHardcopy    = 'off';   % keep the white background we set

[~, ~, ext] = fileparts(file);
switch lower(ext)
    case '.pdf',            driver = '-dpdf';
    case {'.eps', '.epsc'}, driver = '-depsc';
    case {'.tif', '.tiff'}, driver = '-dtiff';
    case '.svg',            driver = '-dsvg';
    otherwise,              driver = '-dpng';
end
print(fig, file, driver, sprintf('-r%d', round(dpi)));
end

function applyFontSize(root, pts)
% Every axes at one size, so a composite does not print its curve panel in one
% size and its maps in another - which is what happens when each drawer picks
% its own and the whole thing is then scaled.
small = max(pts - 1, 5);

% Free text and legends FIRST, then the axes: an axes Title is itself a text
% object, so doing it the other way round would set every title back to the
% smaller size a moment after setting it to the larger one.
set(findall(root, 'Type', 'text'),   'FontSize', small);
set(findall(root, 'Type', 'legend'), 'FontSize', small);

kids = [findall(root, 'Type', 'axes'); findall(root, 'Type', 'uiaxes')];
for k = 1:numel(kids)
    kids(k).FontSize = pts;
    labels = {kids(k).Title, kids(k).XLabel, kids(k).YLabel};
    for j = 1:numel(labels)
        if ~isempty(labels{j}) && isgraphics(labels{j})
            labels{j}.FontSize = pts;
        end
    end
end
end

function lines = provenanceLines(p)
% One "key: value" line per scalar field, wrapped across at most two lines.
% Anything not renderable as short text is skipped rather than printed as a
% struct - a footer is a caption, not a dump.
lines = {};
if ~isstruct(p) || isempty(fieldnames(p)); return; end
names = fieldnames(p);
bits  = {};
for k = 1:numel(names)
    v = p.(names{k});
    if ischar(v) || (isstring(v) && isscalar(v))
        s = char(v);
    elseif isnumeric(v) && isscalar(v)
        s = num2str(v);
    else
        continue
    end
    if isempty(strtrim(s)); continue; end
    bits{end+1} = sprintf('%s: %s', names{k}, s); %#ok<AGROW>
end
if isempty(bits); return; end

lines = {strjoin(bits, '   |   ')};
if numel(lines{1}) > 150
    half  = ceil(numel(bits) / 2);
    lines = {strjoin(bits(1:half), '   |   '), strjoin(bits(half+1:end), '   |   ')};
end
end
