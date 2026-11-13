% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef WindowThemeTest < NestappTestCase
% WINDOWTHEMETEST  Every nestapp window is light, whatever the desktop or OS uses.
%
%   Under OS dark mode (R2025a+, e.g. a Mac set to match the system) MATLAB
%   re-coloured whatever nestapp had not fixed, and windows came out half
%   light, half dark. The fix is made once, where windows are created:
%   appWindow and outputFigure pin the light theme. This pins the helpers,
%   and - because a window created any other way silently brings the mix
%   back - checks that src/ creates no window any other way.
%
%   In tests/gui: it creates real figures, but needs no EEGLAB and no app.

    methods (Test)

        function anAppWindowIsPinnedLight(tc)
            fig = appWindow('Visible', 'off');
            tc.addTeardown(@() delete(fig));
            tc.verifyPinnedLight(fig);
        end

        function anOutputFigureIsWhiteAndPinnedLight(tc)
            fig = outputFigure('Visible', 'off');
            tc.addTeardown(@() delete(fig));
            tc.verifyEqual(fig.Color, [1 1 1]);
            tc.verifyPinnedLight(fig);
        end

        function aCallersColourStillWins(tc)
            fig = outputFigure('Visible', 'off', 'Color', [0.9 0.9 0.8]);
            tc.addTeardown(@() delete(fig));
            tc.verifyEqual(fig.Color, [0.9 0.9 0.8]);
        end

        function aFigureWithoutThemesIsLeftAlone(tc)
        % Releases before R2025a have no ThemeMode; pinning must be a no-op.
            thing = struct('Name', 'no themes here');
            tc.verifyEqual(pinLightTheme(thing), thing);
        end

        function noWindowIsCreatedOutsideTheHelpers(tc)
        % A creation call is figure/uifigure with no arguments or a
        % name/value string first. figure(handle) only focuses an existing
        % window, so it is allowed.
            src = fullfile(addNestappPath(), 'src');
            files = dir(fullfile(src, '**', '*.m'));
            helpers = {'appWindow.m', 'outputFigure.m'};
            bad = {};
            for k = 1:numel(files)
                if ismember(files(k).name, helpers); continue; end
                p = fullfile(files(k).folder, files(k).name);
                lines = strsplit(fileread(p), newline);
                for n = find(cellfun(@createsWindow, lines))
                    bad{end+1} = sprintf('%s:%d  %s', ...
                        erase(p, [src filesep]), n, strtrim(lines{n})); %#ok<AGROW>
                end
            end
            tc.verifyEmpty(bad, sprintf(['Create windows with appWindow (app ' ...
                'windows) or outputFigure (files, scratch) so they stay light ' ...
                'under OS dark mode:\n  %s'], strjoin(bad, '\n  ')));
        end
    end

    methods (Test, TestTags = {'rule'})
        function theRuleCatchesCreationAndNothingElse(tc, line)
        % The rule above is only worth as much as its pattern: these pin what
        % it must catch and what it must leave alone.
            tc.verifyEqual(createsWindow(line.code), line.creates, line.code);
        end
    end

    properties (TestParameter)
        line = struct( ...
            'figureWithArgs',   struct('code', 'f = figure(''Visible'', ''off'');', 'creates', true), ...
            'bareUifigure',     struct('code', 'f = uifigure;', 'creates', true), ...
            'bareStatement',    struct('code', 'figure;', 'creates', true), ...
            'emptyParens',      struct('code', 'h = figure();', 'creates', true), ...
            'focusCall',        struct('code', 'figure(app.UIFigure)', 'creates', false), ...
            'insideAString',    struct('code', 'error(''could not write the QC figure'')', 'creates', false), ...
            'findallType',      struct('code', 'h = findall(groot, ''Type'', ''figure'');', 'creates', false), ...
            'inAComment',       struct('code', 'x = 1;  % figure(''Name'', ''a'')', 'creates', false), ...
            'helperName',       struct('code', 'f = outputFigure(''Visible'', ''off'');', 'creates', false))
    end

    methods (Access = private)
        function verifyPinnedLight(tc, fig)
            if ~isprop(fig, 'ThemeMode'); return; end   % before R2025a: no themes
            tc.verifyEqual(char(fig.ThemeMode), 'manual', 'theme not pinned');
            tc.verifyEqual(char(fig.Theme.BaseColorStyle), 'light');
        end
    end
end

function tf = createsWindow(line)
% True when a line of MATLAB creates a figure or uifigure: figure/uifigure
% with no arguments or a name/value string first. figure(handle) only focuses
% an existing window. String literals are blanked first - "...the QC figure"
% in a message is not a call - then comments; a quote after a word, ) ] } . or
% another quote is a transpose, not the start of a string.
code = regexprep(line, '(?<![\w)\]}.''])''([^'']|'''')*''', '''''');
code = regexprep(code, '"([^"]|"")*"', '""');
code = regexprep(code, '%.*$', '');
tf = ~isempty(regexp(code, '(?<![\w.])(ui)?figure\s*(\(\s*(\)|''|")|;|$)', 'once'));
end
