% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function selected = roiPicker(current, available, opts)
% ROIPICKER  Modal head-diagram picker for a region of interest.
%   selected = ROIPICKER(current) opens the montage with the electrodes in
%   `current` switched on and returns the chosen cellstr, or [] if cancelled.
%   Returning [] rather than {} for a cancel matters: an empty ROI is a legal
%   thing to choose, so the two cannot share a value.
%
%   selected = ROIPICKER(current, available) greys out electrodes the loaded
%   data does not offer. They stay visible - a missing electrode is information
%   about the data, and hiding it makes the montage look wrong.
%
%   selected = ROIPICKER(current, available, opts):
%     .parent         figure to centre the dialog over
%     .title          window title
%     .availableNote  why the greyed ones are unavailable, e.g.
%                     'not on the modal cap'. The picker renders this reason
%                     rather than asserting one: availability used to mean
%                     "in every loaded file", and groupCurves has since moved
%                     to a modal montage, so a hard-coded sentence would
%                     misstate why an electrode is greyed.
%     .presets        struct array (.name .labels .userDefined) overriding the
%                     stored presets
%     .savePreset     @(name, labels) overriding where a preset is written
%
%   The last two exist so a caller can supply its own preset store. Without
%   them a modal dialog performs global persistent writes with no interception
%   point, which is also what made the picker awkward to test against live
%   user preferences.
%
%   Why this is a popup rather than part of a tab. The montage occupies 350x336
%   of an 867x529 tab - the largest single consumer of space there - to serve a
%   choice most users make once a session.
%
%   It also moves where the ROI lives. The old tab read the selection back out
%   of the widgets, looping over app.(NAMEButton).Value, so the UI *is* the
%   state: an ROI cannot be set programmatically, saved into a session, or
%   applied without the buttons existing. The picker takes a cellstr in and
%   hands one back, with the buttons as a view over it.
%
%   See also: roiSelectionState, applyRoiPreset, roiMontageLayout, roiPresets

if nargin < 1 || isempty(current);   current = {};   end
if nargin < 2;                       available = {}; end
if nargin < 3;                       opts = struct(); end
opts = fillDefaults(opts, struct( ...
    'parent', [], 'title', 'Select ROI electrodes', ...
    'availableNote', 'not in every loaded file', ...
    'optional', {{}}, ...
    'presets', [], 'savePreset', @saveRoiPreset));

if ischar(current) || isstring(current);     current = cellstr(current);     end
if ischar(available) || isstring(available); available = cellstr(available); end

[layout, headSize] = roiMontageLayout();
% Every set-of-labels decision lives in roiSelectionState, which is pure and
% unit-tested; this function is the view over it.
st          = roiSelectionState(current, available, opts.optional);
labels      = st.labels;
enabled     = st.enabled;
partial     = st.partial;
selected    = st.selected;
offLabels   = st.offLabels;
offSelected = st.offSelected;
offEnabled  = st.offEnabled;
includeAll  = false;    % the "not in every file" override, off by default
accepted    = false;

PAD    = 12;
SIDE_W = 210;
BAR_H  = 44;
BTN_W  = 150;
GAP    = 5;
figW   = headSize(1) + SIDE_W + PAD * 3;
figH   = headSize(2) + BAR_H + PAD * 2;

fig = appWindow('Name', char(opts.title), 'Resize', 'off', ...
               'Position', centreOn(opts.parent, figW, figH), ...
               'WindowStyle', 'modal');
% Exit is by DELETING the figure, never by uiresume, and the wait is waitfor
% rather than uiwait. selectDataTree already documents why: uiresume's
% close-on-X path can leave the window up and soft-lock the app. The concrete
% failure is a nested modal - uiconfirm or uialert on this same figure runs its
% own wait, and afterwards a uiresume no longer releases the outer uiwait, so
% the X silently does nothing and the app is stuck behind a window that will
% not close. waitfor returns the moment the figure is destroyed, and a plain
% delete cannot be vetoed or missed.
%
% Not an onCleanup either: the callbacks are nested functions holding this
% workspace, which would hold the onCleanup, which would hold the figure - a
% cycle that leaks a window on every call.
%
% This dialog has both nested modals - uialert when a preset is partly
% unavailable, inputdlg when one is saved - so it had the same latent freeze.
fig.CloseRequestFcn = @(src, ~) delete(src);

try
% ── head diagram with one toggle per electrode ───────────────────────────
headPanel = uipanel(fig, 'BorderType', 'none', ...
    'Position', [PAD, BAR_H + PAD, headSize(1), headSize(2)]);
uiimage(headPanel, 'Position', [0 0 headSize], ...
        'ImageSource', headImagePath(), 'ScaleMethod', 'fit');

btns = gobjects(1, numel(layout));
for i = 1:numel(layout)
    btns(i) = uibutton(headPanel, 'state', ...
        'Text', labels{i}, 'Position', layout(i).pos, ...
        'FontSize', 8, 'FontWeight', 'bold', ...
        'IconAlignment', 'center', 'HorizontalAlignment', 'left', ...
        'Value', selected(i), 'Enable', enabled(i));
    btns(i).ValueChangedFcn = @(src, ~) onToggle(i, src.Value);
end

% ── preset / bulk-action column ──────────────────────────────────────────
sx = PAD * 2 + headSize(1);
y  = figH - PAD;

uilabel(fig, 'Text', 'Presets', 'FontWeight', 'bold', 'Position', row(22, 0));
presetDrop = uidropdown(fig, 'Position', row(22, 4), 'Items', {}, ...
    'Editable', 'off', 'ValueChangedFcn', @(~, ~) updateDeleteEnable());
uibutton(fig, 'Text', 'Apply preset',     'Position', row(24, 4), ...
    'ButtonPushedFcn', @(~, ~) applyPreset());
uibutton(fig, 'Text', 'Save current as…', 'Position', row(24, 4), ...
    'ButtonPushedFcn', @(~, ~) savePreset());
deleteBtn = uibutton(fig, 'Text', 'Delete preset', 'Position', row(24, 4), ...
    'ButtonPushedFcn', @(~, ~) deletePreset());

uilabel(fig, 'Text', 'Selection', 'FontWeight', 'bold', 'Position', row(22, 16));
uibutton(fig, 'Text', 'Select all available', 'Position', row(24, 4), ...
    'ButtonPushedFcn', @(~, ~) setAll(true));
uibutton(fig, 'Text', 'Clear all', 'Position', row(24, 4), ...
    'ButtonPushedFcn', @(~, ~) setAll(false));

partialBox = uicheckbox(fig, 'Position', row(22, 8), ...
    'Text', 'Include partial electrodes', 'Value', false, ...
    'Tooltip', ['Offer electrodes present in some files but not all. The ROI ' ...
                'mean is then taken over different channels in different ' ...
                'files, so this is a methodological choice, not a display one.'], ...
    'ValueChangedFcn', @(src, ~) setIncludeAll(src.Value));
if ~any(partial); partialBox.Enable = 'off'; end

% Electrodes with no spot on the head image. Checkboxes, not a text line:
% they are legal ROI members, and listing them without offering them is how
% accepting the dialog used to delete them.
offBox = gobjects(1, numel(offLabels));
if ~isempty(offLabels)
    uilabel(fig, 'Text', 'Other electrodes', 'FontWeight', 'bold', ...
            'Position', row(20, 10));
    for i = 1:numel(offLabels)
        offBox(i) = uicheckbox(fig, 'Position', row(20, 0), ...
            'Text', offLabels{i}, 'Value', offSelected(i), ...
            'Enable', offEnabled(i), ...
            'ValueChangedFcn', @(src, ~) onOffToggle(i, src.Value));
    end
end

countLabel = uilabel(fig, 'Position', row(46, 12), 'WordWrap', 'on', ...
                     'VerticalAlignment', 'top', 'FontColor', [0.35 0.38 0.43]);

% ── accept / cancel ──────────────────────────────────────────────────────
right = figW - PAD;
uibutton(fig, 'Text', 'Use these electrodes', ...
    'Position', [right - 2 * BTN_W - GAP, PAD, BTN_W, 26], ...
    'ButtonPushedFcn', @(~, ~) accept());
uibutton(fig, 'Text', 'Cancel', ...
    'Position', [right - BTN_W, PAD, BTN_W, 26], ...
    'ButtonPushedFcn', @(~, ~) delete(fig));

presets = struct('name', {}, 'labels', {}, 'userDefined', {});
refreshPresets();
refreshCount();

% Announce that the dialog is about to block. Anything driving it - a test, a
% script - otherwise has to guess when uiwait has started, and uiresume issued
% before that point is silently a no-op, leaving the dialog up for good. The
% flag makes "is it waiting yet" answerable instead of a race against however
% long the 69 buttons took to build.
setappdata(fig, 'nestappModalReady', true);
waitfor(fig);
catch ME
    if isvalid(fig); delete(fig); end
    rethrow(ME);
end

if accepted
    % Diagram picks PLUS the off-diagram ones. Rebuilding from labels(selected)
    % alone is what silently dropped FT9 from an ROI that already held it.
    selected = currentSelection();
else
    selected = [];        % cancelled - distinct from an empty selection
end
if isvalid(fig); delete(fig); end

% ── nested callbacks ─────────────────────────────────────────────────────
    function p = row(h, gap)
    % One cursor down the side column, so a gap is a named number rather than
    % an offset every control below has to be re-derived from by hand.
        y = y - gap - h;
        p = [sx, y, SIDE_W, h];
    end

    function sel = currentSelection()
    % The ROI as it stands: diagram picks plus off-diagram ones. Anything that
    % reads the selection back out must go through here - rebuilding it from
    % labels(selected) is what silently dropped FT9 on accept and on save.
        sel = [labels(selected), offLabels(offSelected)];
    end

    function onOffToggle(i, value)
        offSelected(i) = value;
        refreshCount();
    end

    function setIncludeAll(value)
    % The override only ever ADDS what the data has; it never offers an
    % electrode no loaded file carries.
        includeAll = value;
        for k = 1:numel(btns)
            set(btns(k), 'Enable', enabled(k) || (includeAll && partial(k)));
        end
        refreshCount();
    end

    function onToggle(i, value)
        selected(i) = value;
        refreshCount();
    end

    function setAll(value)
        % "All" means all SELECTABLE: switching on an electrode no loaded file
        % has would put a channel in the ROI that cannot be averaged. What
        % counts as selectable widens when the partial override is on.
        selected    = value & (enabled | (includeAll & partial));
        offSelected = value & offEnabled;
        syncButtons();
    end

    function applyPreset()
        k = find(strcmp({presets.name}, presetDrop.Value), 1);
        if isempty(k); return; end
        want = presets(k).labels;
        pickable = enabled | (includeAll & partial);
        [selected, missing] = applyRoiPreset(labels, pickable, want);
        % A preset may legitimately name an off-diagram electrode - FT9 is in
        % this cohort's own recordings - so those are applied here rather than
        % reported as missing.
        offSelected = ismember(lower(offLabels), lower(want)) & offEnabled;
        missing     = missing(~ismember(lower(missing), lower(offLabels(offSelected))));
        syncButtons();
        if ~isempty(missing)
            % Silently dropping part of a named ROI would change what the
            % preset means without saying so.
            uialert(fig, sprintf(['This preset names %s, which %s %s. The rest ' ...
                'has been applied.'], strjoin(missing, ', '), ...
                isAre(numel(missing)), opts.availableNote), ...
                'Preset partly unavailable', 'Icon', 'warning');
        end
    end

    function savePreset()
        answer = inputdlg('Name for this ROI:', 'Save ROI preset', [1 45]);
        if isempty(answer) || isempty(strtrim(answer{1})); return; end
        opts.savePreset(strtrim(answer{1}), currentSelection());
        refreshPresets(strtrim(answer{1}));
    end

    function deletePreset()
        if isempty(presetDrop.Value); return; end
        opts.savePreset(presetDrop.Value, {});
        refreshPresets();
    end

    function refreshPresets(selectName)
        if isempty(opts.presets)
            presets = roiPresets();
        else
            presets = opts.presets;
        end
        presetDrop.Items = {presets.name};
        if nargin >= 1 && any(strcmp(presetDrop.Items, selectName))
            presetDrop.Value = selectName;
        end
        updateDeleteEnable();
    end

    function updateDeleteEnable()
        % Only a user-defined preset can be deleted. Every name in the list is
        % either a built-in or came from the store, so roiPresets' own
        % userDefined flag settles it - no second look at the preferences.
        k = find(strcmp({presets.name}, presetDrop.Value), 1);
        deleteBtn.Enable = ~isempty(k) && presets(k).userDefined;
    end

    function syncOffBoxes()
        for k = 1:numel(offBox)
            if isgraphics(offBox(k)); offBox(k).Value = offSelected(k); end
        end
    end

    function syncButtons()
        for j = 1:numel(btns)
            btns(j).Value = selected(j);
        end
        syncOffBoxes();
        refreshCount();
    end

    function refreshCount()
        n    = sum(selected) + sum(offSelected);
        nOff = sum(~enabled & ~partial);
        txt  = sprintf('%d electrode%s selected', n, plural(n));
        if nOff > 0
            txt = sprintf('%s\n%d %s', txt, nOff, opts.availableNote);
        end
        if any(partial)
            txt = sprintf('%s\n%d in some files but not all', txt, sum(partial));
        end
        countLabel.Text = txt;
    end

    function accept()
        accepted = true;
        delete(fig);
    end
end

% ── helpers ─────────────────────────────────────────────────────────────────

function p = headImagePath()
% Head.png sits in src/, alongside this file. createComponents reaches the same
% file as fileparts(fileparts(mfilename)) from inside @nestapp.
p = fullfile(nestappRoot(), 'src', 'Head.png');
end
