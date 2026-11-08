% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [state, accepted] = groupColoursDialog(names, state, anchor, styleFcn, onApply)
% GROUPCOLOURSDIALOG  Pick a colour theme, and optionally a colour for single groups.
%   [state, accepted] = GROUPCOLOURSDIALOG(names, state, anchor, styleFcn, onApply)
%
%   names     cellstr of the current group names
%   state     .theme (a plotThemes name, '' = default) and .overrides (struct
%             array .name .rgb - may hold groups not shown here; they are kept)
%   anchor    figure to centre on
%   styleFcn  @(state) -> groupStyles result for names under that state, so the
%             swatches show exactly what the plots will draw
%   onApply   @(state) called on every change, and again on cancel with the
%             original state, so the plot behind previews and then restores
%
%   One row per group: its swatch (click to pick a colour), its line style and
%   bar texture as words - which is what tells groups apart in greyscale - and
%   Reset to drop a picked colour back to the theme's.
%
%   Same modal discipline as plotOptionsDialog: waitfor(fig) and a plain
%   delete, never uiwait/uiresume, because a nested picker inside a uiwait'ed
%   figure can leave it unreleasable.
%
%   See also: plotThemes, groupStyles, plotOptionsDialog

if nargin < 5 || isempty(onApply); onApply = @(~) []; end
accepted = false;
original = state;
themes = plotThemes();
if isempty(state.theme) || ~any(strcmpi({themes.name}, state.theme))
    state.theme = themes(1).name;
end

ROW = 30; W = 480;
listH = min(numel(names), 10) * ROW + 6;
H = 150 + listH;

fig = uifigure('Name', 'Plot colours', 'Position', centreOn(anchor, W, H), 'Resize', 'off');
fig.CloseRequestFcn = @(src, ~) delete(src);   % X == cancel

uilabel(fig, 'Text', 'Theme', 'FontWeight', 'bold', 'Position', [16 H-40 60 22]);
uidropdown(fig, 'Items', {themes.name}, 'Value', themeName(state.theme), ...
                'Position', [80 H-40 W-96 22], 'ValueChangedFcn', @(s, ~) onTheme(s.Value));
strip = uipanel(fig, 'Position', [80 H-64 W-96 16], 'BorderType', 'none');

uilabel(fig, 'Text', 'Groups', 'FontWeight', 'bold', 'Position', [16 H-96 100 22]);
list = uipanel(fig, 'Position', [16 52 W-32 listH], 'BorderType', 'none', ...
               'Scrollable', 'on');

uibutton(fig, 'Text', 'Reset all', 'Position', [16 16 90 26], ...
    'Tooltip', 'Use the theme''s colour for every group.', ...
    'ButtonPushedFcn', @(~,~) onResetAll());
uibutton(fig, 'Text', 'Cancel', 'Position', [W-206 16 90 26], ...
    'ButtonPushedFcn', @(~,~) delete(fig));
uibutton(fig, 'Text', 'OK', 'Position', [W-106 16 90 26], ...
    'ButtonPushedFcn', @(~,~) onOk());

refresh();
waitfor(fig);

if ~accepted
    state = original;
    onApply(state);
end

    function refresh()
        if ~isvalid(fig); return; end
        delete(strip.Children);
        t = themes(strcmpi({themes.name}, state.theme));
        n = size(t.colors, 1);
        sw = (strip.Position(3)) / n;
        for k = 1:n
            uilabel(strip, 'Text', '', 'BackgroundColor', t.colors(k, :), ...
                    'Position', [(k-1)*sw 0 sw-2 16]);
        end

        delete(list.Children);
        s = styleFcn(state);
        nG = numel(names);
        rowsH = nG * ROW;
        for g = 1:nG
            y = max(listH, rowsH) - g * ROW + 4;
            uilabel(list, 'Text', names{g}, 'Position', [0 y 120 22], ...
                    'Tooltip', names{g});
            uibutton(list, 'Text', '', 'BackgroundColor', s.colors(g, :), ...
                     'Position', [126 y 44 22], 'Tooltip', 'Pick a colour for this group', ...
                     'ButtonPushedFcn', @(~,~) onPick(g, s.colors(g, :)));
            uilabel(list, 'Text', sprintf('%s line, %s bars', styleWord(s.lineStyles{g}), ...
                    patternWord(s.barPatterns{g})), 'Position', [178 y 200 22], 'FontColor', [0.4 0.42 0.46]);
            uibutton(list, 'Text', 'Reset', 'Position', [W-32-64 y 60 22], ...
                     'Enable', onOff(isOverridden(names{g})), ...
                     'ButtonPushedFcn', @(~,~) onResetOne(names{g}));
        end
    end

    function onTheme(name)
        state.theme = name;
        changed();
    end

    function onPick(g, current)
        c = uisetcolor(current, sprintf('Colour for %s', names{g}));
        if isvalid(fig); figure(fig); end   % uisetcolor can leave focus elsewhere
        if isequal(c, 0); return; end       % the picker was cancelled
        dropOverride(names{g});
        state.overrides(end+1) = struct('name', names{g}, 'rgb', c);
        changed();
    end

    function onResetOne(name)
        dropOverride(name);
        changed();
    end

    function onResetAll()
        % Only the groups shown: overrides for groups not loaded now are kept.
        for g = 1:numel(names); dropOverride(names{g}); end
        changed();
    end

    function dropOverride(name)
        if isempty(state.overrides); return; end
        state.overrides(strcmp({state.overrides.name}, name)) = [];
    end

    function tf = isOverridden(name)
        tf = ~isempty(state.overrides) && any(strcmp({state.overrides.name}, name));
    end

    function changed()
        refresh();
        onApply(state);
    end

    function onOk()
        accepted = true;
        delete(fig);
    end
end

function v = themeName(name)
themes = plotThemes();
v = themes(find(strcmpi({themes.name}, name), 1)).name;
end

function w = styleWord(ls)
switch ls
    case '--', w = 'dashed';
    case ':',  w = 'dotted';
    case '-.', w = 'dash-dot';
    otherwise, w = 'solid';
end
end


function s = onOff(tf)
if tf; s = 'on'; else; s = 'off'; end
end

function w = patternWord(p)
switch p
    case 'solid',        w = 'plain';
    case 'backdiagonal', w = 'back-diagonal';
    case 'crosshatch',   w = 'crosshatched';
    case 'horizontal',   w = 'striped';
    case 'dots',         w = 'dotted';
    otherwise,           w = p;
end
end
