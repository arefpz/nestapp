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
%   print, and hidden from legends. All bars of one ink go into one object, so
%   a panel costs at most four. Call after the axes limits are final - the
%   geometry is laid out in pixels so 45 degrees stays 45 degrees whatever the
%   panel's aspect ratio; the app redraws on resize.
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
toPx = @(xd, yd) [(xd - xl(1)) / diff(xl) * W; (yd - yl(1)) / diff(yl) * H];
SPACING = 6;                       % pixels between hatch lines

% Per ink (1 = black, 2 = white): every bar's segments and dots, concatenated.
segs = {zeros(2, 0), zeros(2, 0)};
dots = {zeros(2, 0), zeros(2, 0)};
dark = relativeLuminance(fillColors) <= 0.25;
for i = 1:numel(x)
    lo = toPx(x(i) - width / 2, min(0, heights(i)));
    hi = toPx(x(i) + width / 2, max(0, heights(i)));
    rect = [lo(1) hi(1) max(lo(2), 0) min(hi(2), H)];   % the visible part of the bar
    [s, d] = hatchSegments(rect, patterns{i}, SPACING);
    ink = dark(i) + 1;
    segs{ink} = [segs{ink}, s];
    dots{ink} = [dots{ink}, d];
end

INK = [0 0 0; 1 1 1];
toData = @(p) [xl(1) + p(1, :) / W * diff(xl); yl(1) + p(2, :) / H * diff(yl)];
wasHeld = ishold(ax);
hold(ax, 'on');
for ink = 1:2
    if ~isempty(segs{ink})
        p = toData(segs{ink});
        h(end+1) = plot(ax, p(1, :), p(2, :), '-', 'Color', INK(ink, :), ...
                        'LineWidth', 0.7, 'HandleVisibility', 'off'); %#ok<AGROW>
    end
    if ~isempty(dots{ink})
        p = toData(dots{ink});
        h(end+1) = plot(ax, p(1, :), p(2, :), '.', 'Color', INK(ink, :), ...
                        'MarkerSize', 6, 'HandleVisibility', 'off'); %#ok<AGROW>
    end
end
if ~wasHeld; hold(ax, 'off'); end
end
