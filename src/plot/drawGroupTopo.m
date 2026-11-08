% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function clim = drawGroupTopo(axList, res, opts)
% DRAWGROUPTOPO  One scalp map per group, by default on one shared scale.
%   clim = DRAWGROUPTOPO(axList, res, opts) draws the window-averaged scalp
%   distribution for each group into the corresponding axes of axList, and
%   returns the symmetric colour limits used - or [] when each map carries its
%   own scale, because then there is nothing one bar could describe. A caller
%   minting a shared colour bar must check for that.
%
%   axList must have at least one axes per group; the caller mints them, so
%   this works the same into a uifigure panel and into a publication figure.
%
%   opts:
%     .window   [t1 t2] in ms to average over (default the whole epoch)
%     .scale    'shared' (default) | 'per map'
%     .climit   scalar uV pinning the scale to [-climit climit]; overrides
%               .scale, since a stated number is a stated scale. 0 and
%               non-finite are read as no scale at all, i.e. derived
%     .clim     force colour limits outright (2-element); wins over both
%     .titles   cellstr overriding the per-map titles
%     .colorbar false when the caller draws its own shared bar (see
%               sharedColorbar); ignored under 'per map', where each map needs
%               its own
%     .markers  passed to drawScalpTopo: 'off' | 'dots' | 'labels'
%     .contours passed to drawScalpTopo: contour-line count
%
%   THE SHARED SCALE IS STILL THE DEFAULT, and for the usual reason: each map
%   drawn with its own limits makes two groups look alike whatever their
%   amplitudes - a 1 uV map and a 10 uV map both render as a full-range red and
%   blue dipole. Comparability by eye is the only reason to put maps side by
%   side, so nothing about the default view changes.
%
%   'per map' exists for the case the shared scale is bad at: one group an
%   order of magnitude larger flattens every other map to near-neutral, and
%   the TOPOGRAPHY of the smaller groups - which is a separate question from
%   their amplitude - becomes unreadable. Each map then gets its own bar, and
%   clim comes back empty so a caller cannot hang one shared bar over maps that
%   no longer share a scale.
%
%   'climit' answers a third question: holding the scale FIXED across runs, so
%   two exported figures can be compared. It is a scalar rather than a pair
%   because the limits must stay symmetric about zero - with a diverging map an
%   off-centre zero puts the neutral colour somewhere other than
%   no-deflection and misreads the polarity, which is drawScalpTopo's
%   invariant, not a preference.
%
%   See also: drawScalpTopo, groupCurves, divergingColormap, sharedColorbar

if nargin < 3; opts = struct(); end
opts = fillDefaults(opts, struct('scale', 'shared', 'climit', [], ...
                                'markers', [], 'contours', []));
% markers and contours stay EMPTY when unset rather than being defaulted
% here: drawScalpTopo owns those defaults, and fillDefaults reads empty as
% absent, so passing [] through leaves them in exactly one place.
nG = numel(res.groups);
if nG == 0; clim = [0 0]; return; end
if numel(axList) < nG
    error('nestapp:tooFewAxes', ...
          'drawGroupTopo needs one axes per group: %d groups, %d axes.', ...
          nG, numel(axList));
end

if ~isfield(opts, 'window') || isempty(opts.window)
    opts.window = [res.time(1) res.time(end)];
end

% Window average per group, all computed before anything is drawn so the shared
% scale can be derived from the complete set.
vals = cell(1, nG);
sel  = res.time >= opts.window(1) & res.time <= opts.window(2);
if ~any(sel)
    error('nestapp:windowOutsideEpoch', ...
        'The window %g-%g ms lies outside the epoch (%g to %g ms).', ...
        opts.window(1), opts.window(2), res.time(1), res.time(end));
end
for g = 1:nG
    vals{g} = mean(res.groups(g).chanMeans(:, sel), 2, 'omitnan');
end

% One map per group is a single ROW, so this plot's "per map" is the same
% operation as TEP-topo's "per window": one limit per column. topoColourScale
% owns the precedence (an explicit pair beats a pinned magnitude beats the
% mode) and the symmetry, and is unit-tested without drawing anything.
override = fieldOr(opts, 'clim', opts.climit);
[perMapClim, clim] = topoColourScale(vals, matchesChoice(opts.scale, 'per map'), ...
                                     override);
perMap = isempty(clim);

% No bars on the maps. A bar shrinks the axes it attaches to, so putting one on
% the last map alone left that group's head smaller than its siblings - in a
% figure whose whole purpose is that the maps are comparable by eye. The caller
% mints the axes, so the caller reserves a strip and hangs one sharedColorbar
% there; opts.colorbar exists for a caller that has nowhere to put one.
% Under 'per map' the caller's shared strip is empty, so each map has to carry
% its own bar whatever the caller asked for - a map with no scale at all is
% worse than a slightly smaller head.
showBar = ~isfield(opts, 'colorbar') || opts.colorbar;
mapOpts = struct('markers', opts.markers, 'contours', opts.contours);
for g = 1:nG
    mapOpts.clim     = perMapClim{g};
    mapOpts.colorbar = perMap || (showBar && g == nG);
    drawScalpTopo(axList(g), vals{g}, res.chanlocs, mapOpts);
    if isfield(opts, 'titles') && numel(opts.titles) >= g
        title(axList(g), opts.titles{g});
    else
        title(axList(g), sprintf('%s  (%g-%g ms)', res.groups(g).name, ...
                                 opts.window(1), opts.window(2)));
    end
end
end
