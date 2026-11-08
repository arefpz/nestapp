% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef TopoColourTest < NestappTestCase
% TOPOCOLOURTEST  Every map in an exported grid keeps the diverging colormap.
%
%   Reported from the field: TEP-topo maps were red-blue in the app and rainbow
%   in the exported PNG - all but the last one. The export draws into CLASSIC
%   axes, where each topoplot call resets the colormap of the maps already
%   drawn, and the per-window scale skipped the pass that re-applies it after
%   the grid. The app's uiaxes did not show it, so only a classic-axes draw
%   catches it. Needs EEGLAB for topoplot.

    methods (TestClassSetup)
        function bootEeglab(tc)
            startEeglab(tc);
        end
    end

    methods (Test)

        function everyTepTopoMapIsDivergingOnClassicAxes(tc, mapScale)
            res = fakeGroupResult('groups', 2, 'subjects', 3);
            fig = tc.classicFigure();
            info = drawTEPTopo(fig, res, struct('mapScale', mapScale, ...
                'curve', false, 'axesFcn', @(p, pos) axes('Parent', p, ...
                'Units', 'pixels', 'Position', pos)));
            maps = tc.mapAxes(info.axes);
            tc.verifyNotEmpty(maps);
            tc.verifyAllDiverging(maps);
        end

        function everyGroupTopoMapIsDivergingOnClassicAxes(tc, groupScale)
            res = fakeGroupResult('groups', 3, 'subjects', 3);
            fig = tc.classicFigure();
            ax = gobjects(1, 3);
            for g = 1:3
                ax(g) = axes('Parent', fig, 'Units', 'pixels', ...
                             'Position', [20 + 140 * (g - 1) 20 120 120]);
            end
            drawGroupTopo(ax, res, struct('scale', groupScale, 'window', [0 50]));
            tc.verifyAllDiverging(ax);
        end
    end

    properties (TestParameter)
        mapScale   = {'shared', 'per window'}
        groupScale = {'shared', 'per map'}
    end

    methods (Access = private)
        function fig = classicFigure(tc)
            fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
                         'Position', [50 50 900 500]);
            tc.addTeardown(@() delete(fig));
        end

        function maps = mapAxes(~, axList)
        % The maps, not the colour bar axes or the curve panel: maps hold the
        % topoplot surface.
            keep = arrayfun(@(a) isgraphics(a) && ...
                ~isempty(findobj(a, 'Type', 'surface')), axList);
            maps = axList(keep);
        end

        function verifyAllDiverging(tc, axList)
            cmap = divergingColormap();
            for k = 1:numel(axList)
                tc.verifyEqual(axList(k).Colormap, cmap, 'AbsTol', 1e-12, ...
                    sprintf('map %d of %d lost the diverging colormap', k, numel(axList)));
                lim = axList(k).CLim;
                tc.verifyEqual(lim(1), -lim(2), 'AbsTol', 1e-9, ...
                    sprintf('map %d is not on a symmetric scale', k));
            end
        end
    end
end
