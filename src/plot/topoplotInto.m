% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function topoplotInto(ax, values, chanlocs, varargin)
% TOPOPLOTINTO  EEGLAB topoplot into one classic axes, leaving the rest alone.
%   TOPOPLOTINTO(ax, values, chanlocs, ...) draws topoplot(values, chanlocs,
%   ...) into ax and undoes what topoplot does to everything AROUND it:
%
%   - It sets the FIGURE colormap (EEGLAB's turbo). That flips every axes in
%     the figure back to inheriting it, so in a grid every map drawn before
%     this one went rainbow, and a parula heatmap beside ICA maps in the QC
%     figure went turbo. Sibling axes that had their own colormap get it back.
%   - It paints the figure its own pale blue (icadefs BACKCOLOR), on screen
%     and in every exported file. The figure colour is put back.
%   - It calls axis equal/square, which moves and resizes ax. A caller laying
%     out a grid has already decided where the map goes, so the frame is put
%     back.
%
%   Every topoplot into a figure shared with other axes goes through here, so
%   a figure that holds a scalp map behaves like any other figure. (The uiaxes
%   path in drawScalpTopo draws into its own throwaway figure and needs none.)
%
%   See also: drawScalpTopo, renderQualityFigure, topoplot

fig       = ancestor(ax, 'figure');
keepUnits = ax.Units;
keepPos   = ax.Position;
keepColor = fig.Color;
siblings  = findall(fig, 'Type', 'axes', 'ColormapMode', 'manual');
siblings  = siblings(siblings ~= ax);
keepMaps  = arrayfun(@(a) a.Colormap, siblings, 'UniformOutput', false);

axes(ax);
topoplot(values, chanlocs, varargin{:});

ax.Units    = keepUnits;
ax.Position = keepPos;
fig.Color   = keepColor;
for k = 1:numel(siblings)
    colormap(siblings(k), keepMaps{k});
end
end
