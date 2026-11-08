% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function info = drawTEPTopo(parent, res, opts)
% DRAWTEPTOPO  Curve panel with scalp maps above it, one column per window.
%   info = DRAWTEPTOPO(parent, res, opts) draws the composite figure: the group
%   waveforms across the bottom, and a grid of scalp maps above - one COLUMN
%   per window of interest, one ROW per group.
%
%   info returns .clim (the shared colour limits) and .axes (every axes made).
%
%   Rows are groups because that is what makes it a comparison. The classic
%   timtopo shows one dataset's maps in a strip; with n groups the same strip
%   per group, aligned on the same columns and the same colour scale, is what
%   lets a reader see that a component differs between conditions rather than
%   just that it exists. With one group it degenerates to the classic figure.
%
%   The maps AVERAGE OVER EACH WINDOW rather than sampling a single latency, so
%   a map and the number in the measures table describe the same interval by
%   construction. Change N45 to 40-55 ms and both move together. Sampling at a
%   window's midpoint instead would put a map beside a mean it does not match.
%
%   One colour scale across every map by default - per-map limits would make a
%   1 uV map and a 10 uV map look identical, which is the one thing a grid of
%   maps must not do.
%
%   'per window' scales each COLUMN instead, and is a different claim rather
%   than a weaker one: it keeps the groups comparable within a window, which is
%   the comparison the rows exist for, while letting a late component read at
%   all. P180 against a scale set by an early muscle-adjacent peak is a
%   uniformly neutral row that says nothing. Each column then states its own
%   limit in its title - with a symmetric diverging map that fixes the whole
%   mapping, white at 0 and the extremes at the stated value, so it needs no
%   bar of its own - and info.clim comes back EMPTY, the same signal
%   drawGroupTopo gives, so a caller cannot hang one shared bar over columns
%   that no longer share a scale.
%
%   The curve panel shades each window, so the reader can see which slice of the
%   waveform each column of maps came from - the linkage the dashed leader lines
%   in AARATEP's c_EEG_plotTimtopo draw. Shading is cheaper and does not break
%   when windows overlap.
%
%   opts:
%     .windows    window structs; omit for defaultTEPComponentDefs. Supplying an
%                 EMPTY set means no windows, and draws the curve alone
%     .mapWindows names of the windows to map, a subset of .windows; omit for
%                 all of them. Six columns is what crowds an 89 mm figure
%     .mapScale   'shared' (default) | 'per window'
%     .curve      false to omit the curve panel and give the grid the height
%     .mode       'TEP' | 'GMFP' | 'LMFP' for the curve's axis label
%     .xlim      curve time limits, default [-50 300]
%     .colors    group colours, default groupColors(nGroups)
%     .showBands shade the windows on the curve, default true
%     .axesFcn   @(parent, position) -> axes. Default makes uiaxes, which is
%                what an in-app panel needs; a publication figure passes a
%                classic-axes maker, because every export path drops UI
%                components. The alternative is sniffing the parent's type,
%                which is unreliable - a uifigure and a figure are both
%                matlab.ui.Figure.
%
%   See also: drawTEPOverlay, drawScalpTopo, drawGroupTopo, groupCurves

if nargin < 3; opts = struct(); end
nG = numel(res.groups);
% "Not supplied" and "supplied empty" mean opposite things here: no windows
% named means fall back to the standard components, but a caller who emptied
% its windows table has said there are none, and reviving the defaults there
% would show six windows the user just deleted. So .windows is settled BEFORE
% fillDefaults and left out of its defaults struct, which cannot then overwrite
% a deliberate empty.
if ~isfield(opts, 'windows'); opts.windows = defaultTEPComponentDefs(); end
opts = fillDefaults(opts, struct( ...
    'mode', 'TEP', 'xlim', [-50 300], ...
    'colors', groupColors(max(nG, 1)), 'showBands', true, ...
    'mapWindows', [], 'mapScale', 'shared', 'curve', true, ...
    'axesFcn', @(p, pos) uiaxes(p, 'Position', pos)));

% The subset is applied to opts.windows itself, so the columns and the curve's
% shading are the same set by construction. Shading six windows under three
% columns of maps would break the linkage the shading exists to draw.
opts.windows = keepNamedWindows(opts.windows, opts.mapWindows);

info = struct('clim', [0 0], 'axes', gobjects(0));
if nG == 0 || isempty(res.time); return; end

w  = opts.windows;
nW = numel(w);
if nW == 0
    % No windows means no maps; the curve alone is still worth drawing.
    ax = opts.axesFcn(parent, curveRect(parent, 0));   % scale derived inside
    drawTEPOverlay(ax, res, opts);
    info.axes = ax;
    return
end

% ── window means per group, all of them, before anything is drawn ────────
% The shared scale can only be computed once every value is known.
vals = cell(nG, nW);
for k = 1:nW
    % The sample selection is a property of the WINDOW, so it is computed once
    % per column rather than once per cell.
    sel = res.time >= min(w(k).winStart, w(k).winEnd) & ...
          res.time <= max(w(k).winStart, w(k).winEnd);
    for g = 1:nG
        if ~any(sel)
            vals{g, k} = zeros(numel(res.channelLabels), 1);
        else
            vals{g, k} = mean(res.groups(g).chanMeans(:, sel), 2, 'omitnan');
        end
    end
end
% Shared by default; 'per window' gives each column its own limit across the
% groups in it. topoColourScale decides both, and returns an EMPTY shared limit
% under 'per window' so nothing downstream can hang one colour bar over columns
% that no longer share a scale.
[colClim, info.clim] = topoColourScale(vals, ...
                            matchesChoice(opts.mapScale, 'per window'));
perWindow = isempty(info.clim);

% ── geometry ────────────────────────────────────────────────────────────
% Maps are SQUARE and sized by whichever budget binds - the width available per
% column, or the height available per row. A head is round, so a cell wider than
% it is tall just adds dead space between columns while the heads stay small;
% sizing off the smaller dimension and centring the block keeps the heads as
% large as the panel allows at any window count and any window size.
P = parentRect(parent);
s = chromeScale(P);
PAD     = max(4,  round(10 * s));
TITLE_H = max(9,  round(18 * s));  % the column titles, above the first row
LABEL_W = max(24, round(58 * s));  % room for the row (group) names
% The bar's reserved strip is built from the parts that are actually drawn -
% gap, bar, tick labels - so widening the bar cannot silently overrun the space
% set aside for it sixty lines away. The tick-label floor is what the earlier
% flat 28 px missed: it left room for the bar and none for the numbers.
BAR_GAP = max(6,  round(10 * s));
BAR_W   = max(10, round(12 * s));
TICK_W  = max(24, round(40 * s));
CBAR_W  = BAR_GAP + BAR_W + TICK_W;
availW  = P(3) - 2*PAD - LABEL_W - CBAR_W;
% With no curve under it the grid takes the height the curve was using. The
% strip reserved for the colour bar stays reserved either way: under 'per
% window' each column carries its limit in its own title instead, and
% reclaiming the strip would move every map when the setting changed.
if opts.curve
    availH = P(4) * 0.55 - TITLE_H - PAD;
else
    availH = P(4) - 2*PAD - TITLE_H;
end
side    = max(min(availW / nW, availH / nG), 40);
gridW   = side * nW;
gridH   = side * nG + TITLE_H + PAD;
x0      = PAD + LABEL_W + max((availW - gridW) / 2, 0);
yTop    = P(4) - PAD - TITLE_H;

axesMade = gobjects(1, nG * nW + 2);   % maps, the shared bar, the curve panel
nMade    = 0;

% ── the map grid ────────────────────────────────────────────────────────
for g = 1:nG
    for k = 1:nW
        x = x0 + (k - 1) * side;
        y = yTop - g * side;
        ax = opts.axesFcn(parent, [x, y, side - 6, side - 6]);
        drawScalpTopo(ax, vals{g, k}, res.chanlocs, ...
                      struct('clim', colClim{k}, 'colorbar', false));
        % Only the top row is titled: one label per column, not nG copies.
        if g == 1
            ttl = sprintf('%s  %g-%g', w(k).name, w(k).winStart, w(k).winEnd);
            if perWindow
                % The column's own scale, since there is no bar to read it off.
                ttl = sprintf('%s\n\\pm%.3g \\muV', ttl, colClim{k}(2));
            end
            title(ax, ttl, 'FontSize', 9);
        else
            title(ax, '');
        end
        nMade = nMade + 1;
        axesMade(nMade) = ax;
        if k == 1; firstInRow = ax; end
    end
    % Row label, in the group's own colour so it ties to the curve below.
    % A text object on the row's FIRST map, not a uilabel on the parent: a
    % uilabel cannot live in a classic figure, and every export path -
    % exportgraphics, print, saveas - silently drops UI components, so a
    % publication figure built that way would come out with no group names on
    % it and no warning that they were missing. Normalized units with clipping
    % off put it just outside the axes whatever the map's data limits are.
    text(firstInRow, -0.06, 0.5, res.groups(g).name, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 10, 'Clipping', 'off', ...
        'Color', opts.colors(min(g, size(opts.colors, 1)), :));
end

if ~perWindow
    % One bar for the whole grid, since every map shares the scale.
    cbAx = sharedColorbar(parent, opts.axesFcn, ...
        [x0 + gridW + BAR_GAP, yTop - nG * side + 6, BAR_W, side * nG - 16], ...
        info.clim);
    nMade = nMade + 1;
    axesMade(nMade) = cbAx;
end

% ── the curve panel ─────────────────────────────────────────────────────
% TEP-topo's registry entry deliberately offers FEWER overlay options than
% drawTEPOverlay accepts - no showBand, no legend, no showTraces. The
% pass-through below is not an implicit promise to expose all of them: this
% curve panel is a sixth of the figure under a grid of scalp maps, and it is
% the panel least able to absorb extra ink.
%
% opts carries .showBands and .windows straight through, and drawTEPOverlay
% owns the shading - calling shadeTimeWindows again here painted every band,
% boundary and label twice, stacking the alpha and overprinting the names.
if opts.curve
    cAx = opts.axesFcn(parent, curveRect(parent, gridH));
    drawTEPOverlay(cAx, res, opts);
    nMade = nMade + 1;
    axesMade(nMade) = cAx;
end
info.axes = axesMade(1:nMade);
end

function w = keepNamedWindows(w, names)
% The named subset, in the ORDER THE WINDOWS ARE DEFINED rather than the order
% they were picked: the columns read left to right in time, and letting a
% selection reorder them would put a late component before an early one.
%
% An empty or absent selection means all of them, and a name matching nothing -
% a window renamed or deleted since the setting was saved - is simply not
% matched. Falling back to all when NOTHING matches keeps a stale saved
% selection from silently emptying the figure.
if isempty(names) || isempty(w); return; end
if ~iscell(names); names = cellstr(string(names)); end
keep = ismember(lower({w.name}), lower(names));
if ~any(keep); return; end
w = w(keep);
end

% ── helpers ─────────────────────────────────────────────────────────────────

function r = parentRect(parent)
p = parent.Position;
r = [0 0 p(3) p(4)];
end

function r = curveRect(parent, gridH)
% The axis labels need room, and how much depends on how big the type is - which
% scales with the panel. Fixed margins tuned for a 900 px panel take 20% of an
% 89 mm one and push the curve up into the maps.
%
% The scale is derived here rather than passed in: as a parameter it had to be
% kept in step by every caller, and the no-windows path above had already
% forgotten it and was laying out an 89 mm figure with 900 px margins.
P    = parentRect(parent);
s    = chromeScale(P);
% Floors, because tick labels and an axis label do not shrink below legibility:
% a proportional bottom margin alone left "Time (ms)" sitting on the footer.
PAD  = max(4,  round(10 * s));
left = max(32, round(45 * s));
bot  = max(24, round(34 * s));
top  = max(28, round(44 * s));   % the panel title needs room above the box
rght = max(12, round(20 * s));   % the last x tick label needs room beside it
r    = [PAD + left, PAD + bot, P(3) - PAD - left - rght, ...
        max(P(4) - gridH - PAD - top, 40)];
end

function s = chromeScale(P)
% Margins, labels and the colour bar are chrome: at the size the app draws they
% are the pixel constants below, and on a smaller canvas they shrink with it
% rather than eating the plot. Reference is 900x500, roughly the in-app panel;
% the floor stops a very small figure from losing its labels entirely.
s = min([P(3) / 900, P(4) / 500, 1]);
s = max(s, 0.45);
end
