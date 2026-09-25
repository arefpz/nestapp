
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function cLim = drawScalpTopo(ax, values, chanlocs, opts)
% DRAWSCALPTOPO  Render an EEGLAB scalp topography into any axes, with a scale.
%   cLim = DRAWSCALPTOPO(ax, values, chanlocs)
%   cLim = DRAWSCALPTOPO(ax, values, chanlocs, opts)
%
%   ax       - destination axes. A classic axes is drawn into directly; a
%              uiaxes goes through a hidden figure (topoplot calls gca/clim/
%              axis, which a uiaxes does not support) and the result is copied
%              back.
%   values   - one scalar per channel (the time-window-averaged voltage, uV).
%   chanlocs - matching EEGLAB chanlocs struct array.
%   opts     - optional struct:
%                .clim      force these colour limits instead of deriving them
%                           from this map alone.
%                .colorbar  false to omit the colour bar. A grid of maps on one
%                           shared scale needs ONE bar, not one per map: twelve
%                           identical bars cost about 40% of the width and
%                           shrink every head to pay for it.
%                .markers   'off' (default) | 'dots' | 'labels'. Whether the
%                           electrode positions are drawn, and named.
%                .contours  number of iso-voltage contour lines, default 5.
%                           0 for a plain interpolated field.
%
%   Returns the symmetric colour limits actually applied, in uV.
%
%   opts.clim exists for comparing maps. Limits derived per map make any two
%   maps look alike whatever their amplitudes - a 1 uV map and a 10 uV map both
%   render as a full-range dipole - so a caller drawing one map per group
%   computes the limits once across all of them and passes them in. Without an
%   override the behaviour is unchanged: symmetric absmax from this map.
%
%   Markers are off by default because a scalp map's job here is the field, and
%   sixty-odd dots over a 40 mm head obscure it. They matter when the question
%   is WHICH electrode - a reader checking that a frontal focus really sits on
%   the channels the ROI names - so 'labels' exists for the figure that has to
%   answer that, and is unreadable at grid size on purpose: it is for one large
%   map, not twelve small ones.
%
%   The colour scale is plain linear microvolts: topoplot's default maplimits
%   is 'absmax', so the limits are symmetric about zero and recomputed for
%   every timepoint. We re-read them off the drawing axes and re-apply them to
%   the destination - copying only the children drops them, leaving the
%   destination silently autoscaled to the copied surface instead of the
%   limits topoplot chose. With a symmetric CLim and a diverging map, white is
%   exactly 0 uV and the colorbar reads directly.
%
%   See also: divergingColormap, drawGroupTopo, topoplot

INTRAD   = 0.55;   % EEGLAB default interpolation radius
CB_LABEL = 'uV';

if nargin < 4 || ~isstruct(opts); opts = struct(); end
opts = fillDefaults(opts, struct('markers', 'off', 'contours', 5));

topoArgs = {'electrodes', topoElectrodes(opts.markers), ...
            'numcontour', max(0, round(opts.contours)), ...
            'intsquare', 'on', 'style', 'map', 'conv', 'on', ...
            'intrad', INTRAD};

% topoplot draws into gca, so it mutates the current figure/axes either way.
% Callers should not have to think about that, so restore it here - and note
% get(groot,'CurrentFigure') rather than gcf, which would CREATE a figure when
% none is open.
prevFig = get(groot, 'CurrentFigure');
restore = onCleanup(@() restoreCurrentFigure(prevFig));

if isa(ax, 'matlab.ui.control.UIAxes')
    cLim = drawViaHiddenFigure(ax, values, chanlocs, topoArgs);
else
    % topoplot calls axis equal/square, which MOVES AND RESIZES the axes it was
    % given. A caller laying out a grid has already decided where this map goes,
    % so the frame is put back afterwards - otherwise the maps grow out of their
    % cells and overlap whatever is below the grid. The uiaxes path never showed
    % this because it draws offscreen and copies only the children back.
    %
    % topoplot also paints the whole FIGURE its own colour: icadefs sets
    % BACKCOLOR = [.93 .96 1], the pale blue EEGLAB uses, and topoplot applies
    % it to the enclosing figure rather than to its own axes. Any figure holding
    % a scalp map therefore turned blue, on screen and in every exported file,
    % including publication figures asked for on a white page. The colour is put
    % back for the same reason the position is.
    keepUnits = ax.Units;
    keepPos   = ax.Position;
    fig       = ancestor(ax, 'figure');
    keepColor = fig.Color;
    cla(ax, 'reset');
    axes(ax);
    topoplot(values, chanlocs, topoArgs{:});
    cLim = symmetricLimits(ax);
    ax.Units    = keepUnits;
    ax.Position = keepPos;
    fig.Color   = keepColor;
end

% A caller comparing several maps supplies one scale for all of them; it wins
% over the per-map limits derived above.
if isfield(opts, 'clim') && ~isempty(opts.clim)
    cLim = opts.clim(:)';
end

ax.CLim = cLim;
colormap(ax, divergingColormap());
if ~isfield(opts, 'colorbar') || opts.colorbar
    cb = colorbar(ax);
    cb.Label.String = CB_LABEL;
end
quietAxes(ax);
end

% ── local helpers ─────────────────────────────────────────────────────────────

function cLim = drawViaHiddenFigure(uiAx, values, chanlocs, topoArgs)
% topoplot needs a classic axes (it calls axes/gca/clim/axis internally), so
% draw into an offscreen one and copy the graphics across.
%
% A fresh figure per call, deliberately, even though it costs about 8.5 ms and
% a TEP-topo grid pays it twelve times a repaint. Keeping one scratch figure
% for the session was tried: topoplot resolves its target through gca, which
% SKIPS handle-invisible figures, so a hidden-handle scratch figure silently
% received nothing and the maps came out blank. Making it handle-visible works
% but leaves a stray figure in the user's session for close all and any figure
% enumeration to trip over, which this function promises not to do. The repaint
% cost is addressed where it actually hurt - the resize handler coalesces a
% drag into one repaint rather than one per pixel.
hiddenFig = figure('Visible', 'off');
closeFig  = onCleanup(@() delete(hiddenFig));

% A blank axes, NOT a clone of the destination: topoplot does not clear what
% it draws into, so cloning a uiaxes that already holds a map would copy the
% previous render back alongside the new one, growing the child list on every
% redraw. Only allchild is copied back, so nothing about the source's own
% properties is needed here.
tmpAx = axes(hiddenFig);
topoplot(values, chanlocs, topoArgs{:});
cLim = symmetricLimits(tmpAx);

cla(uiAx, 'reset');
copyobj(allchild(tmpAx), uiAx);
axis(uiAx, 'equal');
axis(uiAx, 'off');
end

function e = topoElectrodes(markers)
% Our three words for topoplot's option. Named for what the reader sees rather
% than for topoplot's spelling, where 'on' means dots and there is no word for
% the difference between a dot and a named dot.
switch lower(strtrim(char(markers)))
    case {'off', ''};  e = 'off';
    case 'dots';       e = 'on';
    case 'labels';     e = 'labels';
    otherwise
        error('nestapp:badMarkers', ...
              'markers must be off, dots or labels; got "%s".', char(markers));
end
end

function restoreCurrentFigure(prevFig)
if ~isempty(prevFig) && isvalid(prevFig)
    set(groot, 'CurrentFigure', prevFig);
end
end

function cLim = symmetricLimits(ax)
% topoplot's 'absmax' limits are already symmetric bar rounding; force exact
% symmetry so the neutral midpoint of the colormap lands on 0 uV. Falls back
% to a unit range for an all-zero (or degenerate) map.
m = max(abs(ax.CLim));
if ~isfinite(m) || m <= 0
    m = 1;
end
cLim = [-m m];
end
