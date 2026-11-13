% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef PlotThemeTest < NestappTestCase
% PLOTTHEMETEST  Group styles by theme, kept by name, and the bar textures.
%
%   The promise users see is "a group keeps its colour": the same condition
%   must not change colour between two figures because another group was
%   added or removed. And past a theme's colours, groups must still differ -
%   by line style and texture, not by a colour too close to an earlier one.

    methods (Test)

        % ── themes ───────────────────────────────────────────────────────────

        function theDefaultThemeIsOkabeItoWithAmber(tc)
            themes = plotThemes();
            tc.verifyEqual(themes(1).name, 'Okabe-Ito');
            tc.verifyFalse(any(all(abs(themes(1).colors - [0.94 0.89 0.26]) < 0.02, 2)), ...
                'the yellow that vanished on white is gone');
        end

        function everyThemeIsWellFormed(tc)
            for t = plotThemes()'
                tc.verifyGreaterThanOrEqual(size(t.colors, 1), 3, t.name);
                tc.verifyTrue(all(t.colors(:) >= 0 & t.colors(:) <= 1), t.name);
                tc.verifyNotEmpty(t.lineStyles, t.name);
                tc.verifyNotEmpty(t.barPatterns, t.name);
            end
        end

        function groupColorsIsTheDefaultThemeInOrder(tc)
            themes = plotThemes();
            tc.verifyEqual(groupColors(3), themes(1).colors(1:3, :));
            tc.verifySize(groupColors(0), [0 3]);
        end

        % ── stable by name ───────────────────────────────────────────────────

        function removingAGroupDoesNotRecolourTheOthers(tc)
            [before, mem] = groupStyles({'Pre', 'Mid', 'Post'}, '', [], []);
            after = groupStyles({'Pre', 'Post'}, '', [], mem);
            tc.verifyEqual(after.colors(2, :), before.colors(3, :), ...
                'Post kept its colour when Mid was removed');
        end

        function reorderingGroupsDoesNotRecolourThem(tc)
            [before, mem] = groupStyles({'A', 'B', 'C'}, '', [], []);
            after = groupStyles({'C', 'A', 'B'}, '', [], mem);
            tc.verifyEqual(after.colors, before.colors([3 1 2], :));
        end

        function aNewGroupTakesTheFirstFreeColour(tc)
            [~, mem] = groupStyles({'A', 'B', 'C'}, '', [], []);
            [s, ~] = groupStyles({'A', 'C', 'D'}, '', [], mem);
            themes = plotThemes();
            tc.verifyEqual(s.colors(3, :), themes(1).colors(2, :), ...
                'D took the slot B left, not a fourth colour');
        end

        function aNewGroupListedFirstDoesNotTakeAKnownGroupsColour(tc)
            [before, mem] = groupStyles({'Pre', 'Post'}, '', [], []);
            after = groupStyles({'Sham', 'Pre', 'Post'}, '', [], mem);
            tc.verifyEqual(after.colors(2:3, :), before.colors, ...
                'Pre and Post keep their colours when a new group comes first');
        end

        function memoryHoldsEachNameOnce(tc)
            [~, mem] = groupStyles({'A', 'B'}, '', [], []);
            [~, mem] = groupStyles({'B', 'A', 'C'}, '', [], mem);
            tc.verifyNumElements(unique({mem.name}), numel(mem));
        end

        % ── what drawers receive ─────────────────────────────────────────────

        function scriptsGetLineStylesPastThePaletteToo(tc)
        % A script that passes no style must not get identical repeated lines.
            names = arrayfun(@(k) sprintf('g%d', k), 1:9, 'UniformOutput', false);
            sty = resolveGroupStyle(struct(), names);
            tc.verifyNotEqual(sty.lineStyles{end}, sty.lineStyles{1});
        end

        function passedColoursReplaceOnlyThePalette(tc)
            sty = resolveGroupStyle(struct('colors', [1 0 0; 0 0 1]), {'a', 'b', 'c'});
            tc.verifyEqual(sty.colors, [1 0 0; 0 0 1; 1 0 0], 'cycled to the group count');
            tc.verifyNumElements(sty.lineStyles, 3);
        end

        function theDifferenceIsDrawnInTheContrastColourNotAGroups(tc)
        % Group 1's colour made a Pre-minus-Post wave read as the Pre curve.
            over = struct('name', 'Pre', 'rgb', [1 0 0]);
            s = groupStyles({'Pre', 'Post'}, '', over, []);
            tc.verifyNotEqual(s.contrast, s.colors(1, :));
        end

        function aUserColourOverridesTheTheme(tc)
            over = struct('name', 'Post', 'rgb', [0.1 0.2 0.3]);
            s = groupStyles({'Pre', 'Post'}, '', over, []);
            tc.verifyEqual(s.colors(2, :), [0.1 0.2 0.3]);
        end

        function anUnknownThemeFallsBackToTheDefault(tc)
            a = groupStyles({'A', 'B'}, 'No such theme', [], []);
            b = groupStyles({'A', 'B'}, '', [], []);
            tc.verifyEqual(a.colors, b.colors);
        end

        % ── many groups ──────────────────────────────────────────────────────

        function pastTheColoursGroupsDifferByLineStyle(tc)
            themes = plotThemes();
            nC = size(themes(1).colors, 1);
            names = arrayfun(@(k) sprintf('g%d', k), 1:nC + 1, 'UniformOutput', false);
            s = groupStyles(names, '', [], []);
            tc.verifyEqual(s.colors(end, :), s.colors(1, :), 'colours repeat');
            tc.verifyNotEqual(s.lineStyles{end}, s.lineStyles{1}, 'but the style changes');
            tc.verifyNotEqual(s.barPatterns{end}, s.barPatterns{1}, 'and the texture');
            tc.verifyEqual(unique(s.lineStyles(1:nC)), {'-'}, 'first lap is solid');
        end

        function greyscaleVariesStyleFromTheFirstGroup(tc)
            s = groupStyles({'A', 'B', 'C'}, 'Greyscale', [], []);
            tc.verifyNumElements(unique(s.lineStyles), 3);
            tc.verifyNumElements(unique(s.barPatterns), 3);
            tc.verifyGreaterThan(s.lineWidth, 1.75, 'heavier lines for print');
        end
    end
end

