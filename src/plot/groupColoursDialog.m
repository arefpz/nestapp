% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [state, accepted] = groupColoursDialog(names, state, anchor, memory, onApply)
% GROUPCOLOURSDIALOG  Pick a colour for single groups, over the current theme.
%   [state, accepted] = GROUPCOLOURSDIALOG(names, state, anchor, memory, onApply)
%
%   names    cellstr of the current group names
%   state    .theme and .overrides (struct array .name .rgb - may hold groups
%            not shown here; they are kept)
%   anchor   figure to centre on
%   memory   groupStyles slot memory, so the swatches show exactly what the
%            plots will draw
%   onApply  @(state) called on every change, and again on cancel with the
%            original state, so the plot behind previews and then restores
%
%   Opened from the Colours row of plot Options, which holds the theme. One
%   row per group: its swatch (click to pick), its line style and bar texture
%   in words - which is what tells groups apart in greyscale - and Reset.
%
%   Same modal discipline as plotOptionsDialog: waitfor(fig) and a plain
%   delete, never uiwait/uiresume.
%
%   See also: plotOptionsDialog, groupStyles, plotThemes

if nargin < 5 || isempty(onApply); onApply = @(~) []; end
accepted = false;
original = state;

ROW = 30; W = 480;
nG = numel(names);
listH = min(nG, 10) * ROW + 6;
H = 110 + listH;

fig = appWindow('Name', sprintf('Group colours - %s theme', state.theme), ...
               'Position', centreOn(anchor, W, H), 'Resize', 'off', 'WindowStyle', 'modal');
fig.CloseRequestFcn = @(src, ~) delete(src);   % X == cancel

uilabel(fig, 'Text', 'Groups', 'FontWeight', 'bold', 'Position', [16 H-40 100 22]);
list = uipanel(fig, 'Position', [16 52 W-32 listH], 'BorderType', 'none', ...
               'Scrollable', 'on');

% Built once; refresh() only updates what a change can affect.
swatch = gobjects(1, nG); desc = gobjects(1, nG); reset = gobjects(1, nG);
for g = 1:nG
    y = max(listH, nG * ROW) - g * ROW + 4;
    uilabel(list, 'Text', names{g}, 'Position', [0 y 120 22], 'Tooltip', names{g});
    swatch(g) = uibutton(list, 'Text', '', 'Position', [126 y 44 22], ...
        'Tooltip', 'Pick a colour for this group', 'ButtonPushedFcn', @(~,~) onPick(g));
    desc(g) = uilabel(list, 'Position', [178 y 200 22], 'FontColor', [0.4 0.42 0.46]);
    reset(g) = uibutton(list, 'Text', 'Reset', 'Position', [W-96 y 60 22], ...
        'ButtonPushedFcn', @(~,~) onResetOne(names{g}));
end

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

    function s = styles()
        s = groupStyles(names, state.theme, state.overrides, memory);
    end

    function refresh()
        s = styles();
        for k = 1:nG
            swatch(k).BackgroundColor = s.colors(k, :);
            desc(k).Text = sprintf('%s line, %s bars', styleWord(s.lineStyles{k}), ...
                                   patternWord(s.barPatterns{k}));
            reset(k).Enable = matlab.lang.OnOffSwitchState(isOverridden(names{k}));
        end
    end

    function onPick(g)
        s = styles();
        c = uisetcolor(s.colors(g, :), sprintf('Colour for %s', names{g}));
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
        for k = 1:nG; dropOverride(names{k}); end
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

function w = styleWord(ls)
switch ls
    case '--', w = 'dashed';
    case ':',  w = 'dotted';
    case '-.', w = 'dash-dot';
    otherwise, w = 'solid';
end
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
