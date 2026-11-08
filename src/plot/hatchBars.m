% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function h = hatchBars(ax, x, heights, width, patterns, fillColors)
% HATCHBARS  Draw a texture over each bar of a bar chart, so groups survive greyscale.
%   h = HATCHBARS(ax, x, heights, width, patterns, fillColors)
%
%   x, heights  bar centres and signed heights, in data units (bars start at 0)
%   width       bar width in data units, as passed to bar()
%   patterns    cellstr, one hatchSegments pattern per bar
%   fillColors  n x 3, each bar's fill, so the texture is drawn in whichever of
%               black or white stands out against it
%
%   MATLAB (R2026a) has no hatch property on bars or patches, and the File
%   Exchange tools for it carry their own licences, so the texture is drawn
%   here as plain line and marker objects: crisp vectors in PDF/SVG, correct in
%   print, and hidden from legends. Call after the axes limits are final - the
%   geometry is laid out in pixels so 45 degrees stays 45 degrees whatever the
%   panel's aspect ratio.
%
%   See also: hatchSegments, drawWindowBars, plotThemes

h = gobjects(0);
if all(strcmpi(patterns, 'solid')); return; end

keepUnits = ax.Units;
ax.Units = 'pixels';
pos = ax.InnerPosition;            % the plot box, where the data are drawn
ax.Units = keepUnits;
W = max(pos(3), 1); H = max(pos(4), 1);
xl = xlim(ax); yl = ylim(ax);
ylim(ax, yl);                      % the texture must not move the limits it was laid out on
toPx   = @(xd, yd) deal((xd - xl(1)) / diff(xl) * W, (yd - yl(1)) / diff(yl) * H);
toData = @(xp, yp) deal(xl(1) + xp / W * diff(xl), yl(1) + yp / H * diff(yl));
SPACING = 6;   % pixels between hatch lines

wasHeld = ishold(ax);
hold(ax, 'on');
for i = 1:numel(x)
    if strcmpi(patterns{i}, 'solid') || heights(i) == 0; continue; end
    [px0, py0] = toPx(x(i) - width / 2, min(0, heights(i)));
    [px1, py1] = toPx(x(i) + width / 2, max(0, heights(i)));
    % Keep the texture within what is visible of the bar.
    py0 = max(py0, 0); py1 = min(py1, H);
    [segs, dots] = hatchSegments([px0 px1 py0 py1], patterns{i}, SPACING);
    ink = inkFor(fillColors(i, :));
    if ~isempty(segs)
        [xd, yd] = toData(segs(1, :), segs(2, :));
        h(end+1) = plot(ax, xd, yd, '-', 'Color', ink, 'LineWidth', 0.7, ...
                        'HandleVisibility', 'off'); %#ok<AGROW>
    end
    if ~isempty(dots)
        [xd, yd] = toData(dots(1, :), dots(2, :));
        h(end+1) = plot(ax, xd, yd, '.', 'Color', ink, 'MarkerSize', 6, ...
                        'HandleVisibility', 'off'); %#ok<AGROW>
    end
end
if ~wasHeld; hold(ax, 'off'); end
end

function ink = inkFor(rgb)
% Black on light fills, white on dark ones (relative luminance).
lin = ((rgb + 0.055) / 1.055) .^ 2.4;
small = rgb <= 0.04045;
lin(small) = rgb(small) / 12.92;
if lin * [0.2126; 0.7152; 0.0722] > 0.25
    ink = [0 0 0];
else
    ink = [1 1 1];
end
end
