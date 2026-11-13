% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [params, accepted, colours] = plotOptionsDialog(entry, params, anchor, context, onApply, colours)
% PLOTOPTIONSDIALOG  Edit one plot's settings, generated from its registry entry.
%   [params, accepted, colours] = PLOTOPTIONSDIALOG(entry, params, anchor, context, onApply, colours)
%
%   entry    one element of plotRegistry(), whose .params drives the form
%   params   the values already set for this plot (struct, possibly empty)
%   anchor   optional figure to centre on
%   context  optional struct supplying choices for params that declare
%            'choicesFrom' - the windows in the table right now, say. Passed
%            straight to paramForm; nothing here reads it.
%   onApply  optional @(params) called on EVERY change, and again on cancel
%            with the original values, so the plot behind the dialog follows
%            what the controls say
%   colours  optional struct for the Colours row at the top: .state (.theme,
%            .overrides), .names (current groups), .memory (groupStyles slot
%            memory) and .onChange(state), called on every change and with the
%            original state on cancel. Colours apply to EVERY plot - they sit
%            here because this is where people look for how a plot looks. The
%            edited state comes back as the third output.
%
%   Returns the edited params and whether the user accepted them. On cancel -
%   including the window's X - the input is handed back untouched.
%
%   EVERY PLOT SETTING IS A DRAW SETTING, which is what makes applying live
%   safe: none of them changes the estimates, so nothing here can re-run a
%   statistic as a side effect of a marker or a colour scale being adjusted.
%   Turning a band on and looking at a static picture is a poor way to choose,
%   and the reason the settings were unreachable in the first place was that
%   they were hard to evaluate.
%
%   If a COMPUTE setting is ever added - one that changes what groupCurves
%   returns rather than how it is drawn - it must not ride this path: mark it
%   in the registry and gate it behind an explicit Apply, or a colour change
%   starts silently recomputing the cohort.
%
%   Nothing here knows which plot it is editing. paramForm reads the makeParam
%   metadata and picks a control per setting, so a plot that declares a new
%   param in plotRegistry gets a usable editor for it with no change to this
%   file or to the tab. That is the whole reason plot params are makeParam
%   -shaped.
%
%   AN UNSET PARAM STAYS UNSET. Every control can say "default", and a setting
%   the user has not touched is absent from the returned struct rather than
%   frozen to a copy of whatever the default is today. So "I never touched this"
%   and "I set this to the default value" stay distinguishable - the first
%   follows the draw function if its default ever changes, the second does not.
%
%   Modal discipline: waitfor(fig) and a plain delete, never uiwait/uiresume.
%   A nested uiconfirm or uialert inside a figure already sitting in uiwait
%   leaves uiresume unable to release it, and the app hangs - the trap
%   selectDataTree documents and exploreFilesTable was caught by.
%
%   See also: paramForm, plotRegistry, makeParam, plotDrawOpts

if nargin < 2 || isempty(params); params = struct(); end
if nargin < 3; anchor = []; end
if nargin < 4; context = struct(); end
if nargin < 5 || isempty(onApply); onApply = @(~) []; end
if nargin < 6; colours = []; end
accepted = false;

meta = entry.params;
hasColours = ~isempty(colours);
if isempty(meta) && ~hasColours
    return   % nothing to edit; the caller should not have offered the button
end
original = params;
if hasColours; originalColours = colours.state; end

W        = 470;
formH    = 0;
if ~isempty(meta); formH = paramFormHeight(meta); end
COLOURS_H = 40 * hasColours;
H        = formH + COLOURS_H + 108;

fig = appWindow('Name', sprintf('%s - options', entry.name), ...
               'Position', centreOn(anchor, W, H), 'Resize', 'off', 'WindowStyle', 'modal');
fig.CloseRequestFcn = @(src, ~) delete(src);   % X == cancel, and nothing else

uilabel(fig, 'Position', [16 H-32 W-32 22], 'Text', entry.name, ...
        'FontWeight', 'bold');

if hasColours
    themes = plotThemes();
    k = find(strcmpi({themes.name}, colours.state.theme), 1);
    if isempty(k); k = 1; end
    colours.state.theme = themes(k).name;
    y = H - 72;
    uilabel(fig, 'Text', 'Colours', 'Position', [26 y 100 22], ...
            'Tooltip', 'Applies to every plot, not just this one.');
    uidropdown(fig, 'Items', {themes.name}, 'Value', colours.state.theme, ...
        'Position', [130 y 170 22], 'Tooltip', 'Colour theme for every plot', ...
        'ValueChangedFcn', @(s, ~) onTheme(s.Value));
    uibutton(fig, 'Text', 'Group colours...', 'Position', [310 y W-326 22], ...
        'Tooltip', 'Pick a colour for a single group', ...
        'ButtonPushedFcn', @(~,~) onGroups());
end

form = uipanel(fig, 'Position', [16 56 W-32 max(formH, 1)], 'BorderType', 'none');
if ~isempty(meta); paramForm(form, meta, params, @onChange, context); end

uibutton(fig, 'Text', 'Reset', 'Position', [16 16 90 26], ...
    'Tooltip', 'Put every setting back to its default.', ...
    'Enable', matlab.lang.OnOffSwitchState(~isempty(meta)), ...
    'ButtonPushedFcn', @(~,~) onReset());
uibutton(fig, 'Text', 'Cancel', 'Position', [W-206 16 90 26], ...
    'ButtonPushedFcn', @(~,~) delete(fig));
uibutton(fig, 'Text', 'OK', 'Position', [W-106 16 90 26], ...
    'ButtonPushedFcn', @(~,~) onOk());

waitfor(fig);

% Cancel, or the X, leaves the caller with exactly what it passed in. Edits are
% applied to the live struct as they are made - that is what keeps the controls
% and the values in step - so the discard has to happen here.
%
% And because the changes were also being drawn as they were made, cancel has
% to put the PICTURE back too, not just the struct. A dialog that discards its
% values while leaving the plot showing them is worse than one that never
% previewed.
if ~accepted
    params = original;
    onApply(params);
    if hasColours
        colours.state = originalColours;
        colours.onChange(colours.state);
    end
end
if hasColours; colours = colours.state; end

    function onTheme(name)
        colours.state.theme = name;
        colours.onChange(colours.state);
    end

    function onGroups()
        % Its own Cancel already restores (and previews) the state it was given.
        [colours.state, ~] = groupColoursDialog(colours.names, colours.state, fig, ...
                                                colours.memory, colours.onChange);
        if isvalid(fig); figure(fig); end
    end

    function onChange(key, value)
        if isempty(value)
            if isfield(params, key); params = rmfield(params, key); end
        else
            params.(key) = value;
        end
        onApply(params);
    end

    function onReset()
        params = struct();
        delete(form.Children);
        paramForm(form, meta, params, @onChange, context);
        onApply(params);
    end

    function onOk()
        accepted = true;
        delete(fig);
    end
end
