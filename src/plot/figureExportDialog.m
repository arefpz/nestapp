% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [opts, action] = figureExportDialog(defaults, anchor)
% FIGUREEXPORTDIALOG  Choose the size and resolution of a figure for print.
%   [opts, action] = FIGUREEXPORTDIALOG(defaults, anchor)
%
%   Returns a publicationFigure opts struct and one of:
%     'open'   open the figure in MATLAB, leave it there to be edited
%     'save'   write it to opts.file at the chosen resolution
%     'cancel' do nothing
%
%   Both exits matter, and offering only one is why the old tab confused
%   people. A figure that lands in MATLAB is for adjusting - the MATLAB figure
%   tools are better at that than anything the app could grow. A figure that
%   lands in a file is finished, at the width the journal asked for and at print
%   resolution. The same composition backs both, so what is saved is what was
%   opened.
%
%   EVERY SETTING GETS THE CONTROL ITS VALUES CALL FOR, rather than a row in a
%   table of text cells. Each has a shape a text cell cannot express:
%
%     Width       three choices, two of them named constants - radio buttons,
%                 showing the millimetre figure rather than leaving it to be known
%     Height      derived from the width, or a number - a checkbox that
%                 discloses the field it controls
%     Resolution  a closed list of print standards - a dropdown
%     Font size   derived from the width, or a number in a narrow range - a
%                 checkbox and a bounded spinner
%
%   Width is what prompted this: as a text cell it required the user to know
%   that the words were "single" and "double", said nothing about what they
%   measured, and rejected the millimetre value the setting is documented to
%   accept.
%
%   An "auto" box left ticked leaves that field EMPTY in the returned struct,
%   which is how publicationFigure is told to derive it. The derivation stays in
%   the one function that applies it rather than being copied into this dialog.
%
%   Modal discipline: waitfor + delete, never uiwait/uiresume - a nested
%   uiputfile inside a figure sitting in uiwait is exactly the freeze that
%   selectDataTree documents.
%
%   See also: publicationFigure, plotOptionsDialog

SINGLE_MM = 89;
DOUBLE_MM = 183;

if nargin < 1 || isempty(defaults); defaults = struct(); end
if nargin < 2; anchor = []; end

% Seeded EMPTY. publicationFigure owns these defaults and derives some of them
% from each other - the font size follows the column width and floors at 5 pt
% for a single column. Seeding real values here made that derivation
% unreachable, so single-column figures printed at 8 pt with exactly the
% collision it exists to prevent.
blank  = struct('width', [], 'height', [], 'dpi', [], 'fontSize', [], 'file', '');
opts   = fillDefaults(defaults, blank);
given  = opts;
action = 'cancel';

W = 470; H = 360;
fig = appWindow('Name', 'Publication figure', 'Resize', 'off', ...
               'Position', centreOn(anchor, W, H));
fig.CloseRequestFcn = @(src, ~) delete(src);

uilabel(fig, 'Position', [16 H-38 W-32 22], ...
        'Text', 'Publication figure', 'FontWeight', 'bold');

% ── width ───────────────────────────────────────────────────────────────
uilabel(fig, 'Position', [16 H-64 200 18], 'Text', 'Width', ...
        'FontWeight', 'bold', 'FontSize', 11);
widthGroup = uibuttongroup(fig, 'Position', [16 H-152 W-32 84], ...
    'BorderType', 'line', 'SelectionChangedFcn', @onWidth);
rbSingle = uiradiobutton(widthGroup, 'Position', [12 54 240 22], ...
    'Text', sprintf('Single column  (%g mm)', SINGLE_MM));
rbDouble = uiradiobutton(widthGroup, 'Position', [12 30 240 22], ...
    'Text', sprintf('Double column  (%g mm)', DOUBLE_MM));
rbCustom = uiradiobutton(widthGroup, 'Position', [12 6 76 22], 'Text', 'Custom');
customW  = uieditfield(widthGroup, 'numeric', 'Position', [92 6 70 22], ...
    'Limits', [1 Inf], 'Value', DOUBLE_MM, 'ValueDisplayFormat', '%g', ...
    'ValueChangedFcn', @onWidth);
uilabel(widthGroup, 'Position', [168 6 40 22], 'Text', 'mm');

% ── height ──────────────────────────────────────────────────────────────
uilabel(fig, 'Position', [16 H-180 200 18], 'Text', 'Height', ...
        'FontWeight', 'bold', 'FontSize', 11);
autoH = uicheckbox(fig, 'Position', [16 H-206 210 22], ...
    'Text', 'Scale with width (0.62x)', 'Value', true, ...
    'ValueChangedFcn', @onHeight, ...
    'Tooltip', 'A height proportional to the width, which suits most figures.');
customH = uieditfield(fig, 'numeric', 'Position', [232 H-206 70 22], ...
    'Limits', [1 Inf], 'Value', round(DOUBLE_MM * 0.62), ...
    'ValueDisplayFormat', '%g', 'ValueChangedFcn', @onHeight);
uilabel(fig, 'Position', [308 H-206 40 22], 'Text', 'mm');

% ── resolution ──────────────────────────────────────────────────────────
uilabel(fig, 'Position', [16 H-240 100 22], 'Text', 'Resolution', ...
        'FontWeight', 'bold', 'FontSize', 11);
dpiDrop = uidropdown(fig, 'Position', [120 H-240 96 22], ...
    'Items', {'300', '600', '1200'}, 'ItemsData', {300, 600, 1200}, ...
    'Value', 600, 'ValueChangedFcn', @onDpi, ...
    'Tooltip', ['600 dpi is the usual requirement for a figure containing ' ...
                'line art. 300 is the floor for photographic content; some ' ...
                'journals ask 1200 for pure line drawings.']);
uilabel(fig, 'Position', [222 H-240 40 22], 'Text', 'dpi');

% ── font size ───────────────────────────────────────────────────────────
uilabel(fig, 'Position', [16 H-274 100 22], 'Text', 'Font size', ...
        'FontWeight', 'bold', 'FontSize', 11);
autoF = uicheckbox(fig, 'Position', [120 H-274 130 22], ...
    'Text', 'Match width', 'Value', true, 'ValueChangedFcn', @onFont, ...
    'Tooltip', ['Type shrinks with the page: 8 pt at double column, down to a ' ...
                '5 pt floor at single, because a size that reads well across ' ...
                '183 mm collides with itself at 89 mm.']);
customF = uispinner(fig, 'Position', [256 H-274 70 22], ...
    'Limits', [4 14], 'Step', 0.5, 'Value', 8, ...
    'ValueDisplayFormat', '%g', 'ValueChangedFcn', @onFont);
uilabel(fig, 'Position', [332 H-274 40 22], 'Text', 'pt');

uilabel(fig, 'Position', [16 50 W-32 28], 'WordWrap', 'on', ...
    'FontSize', 11, 'FontColor', [0.45 0.48 0.53], ...
    'Text', ['Saved figures record the groups, sample sizes, design and ' ...
             'region of interest in a caption line.']);

uibutton(fig, 'Text', 'Cancel', 'Position', [16 16 90 26], ...
    'ButtonPushedFcn', @(~,~) delete(fig));
uibutton(fig, 'Text', 'Open in MATLAB', 'Position', [W-236 16 120 26], ...
    'ButtonPushedFcn', @(~,~) finish('open'));
uibutton(fig, 'Text', 'Save as...', 'Position', [W-106 16 90 26], ...
    'ButtonPushedFcn', @(~,~) onSave());

showDefaults();
waitfor(fig);

% Cancel, or the X, hands back exactly what was passed in. The controls write
% into opts as they are touched, which is what keeps the struct and the dialog
% in step, so the discard has to happen here.
if strcmp(action, 'cancel')
    opts = given;
end

    function showDefaults()
    % Reflect whatever came in, then read the controls back into opts so the
    % struct and the dialog agree before anything is touched.
        w = opts.width;
        if isnumeric(w) && isscalar(w) && isfinite(w) && w > 0
            widthGroup.SelectedObject = rbCustom;
            customW.Value = w;
        elseif (ischar(w) || isstring(w)) && strcmpi(char(w), 'single')
            widthGroup.SelectedObject = rbSingle;
        else
            widthGroup.SelectedObject = rbDouble;
        end

        autoH.Value = isempty(opts.height);
        if ~isempty(opts.height); customH.Value = opts.height; end

        if ~isempty(opts.dpi) && ismember(opts.dpi, cell2mat(dpiDrop.ItemsData))
            dpiDrop.Value = opts.dpi;
        end

        autoF.Value = isempty(opts.fontSize);
        if ~isempty(opts.fontSize); customF.Value = opts.fontSize; end

        onWidth(); onHeight(); onFont(); onDpi();
    end

    function onWidth(~, ~)
        isCustom = widthGroup.SelectedObject == rbCustom;
        customW.Enable = matlab.lang.OnOffSwitchState(isCustom);
        if isCustom
            opts.width = customW.Value;
        elseif widthGroup.SelectedObject == rbSingle
            opts.width = 'single';
        else
            opts.width = 'double';
        end
    end

    function onHeight(~, ~)
        % Disabled rather than hidden, so the value it would take is visible
        % while the automatic one is in force.
        customH.Enable = matlab.lang.OnOffSwitchState(~autoH.Value);
        if autoH.Value
            opts.height = [];      % empty means "derive it"
        else
            opts.height = customH.Value;
        end
    end

    function onFont(~, ~)
        customF.Enable = matlab.lang.OnOffSwitchState(~autoF.Value);
        if autoF.Value
            opts.fontSize = [];
        else
            opts.fontSize = customF.Value;
        end
    end

    function onDpi(~, ~)
        opts.dpi = dpiDrop.Value;
    end

    function onSave()
        [f, p] = uiputfile({'*.png', 'PNG image'; '*.pdf', 'PDF (vector)'; ...
                            '*.tiff', 'TIFF image'; '*.eps', 'EPS (vector)'}, ...
                           'Save figure', 'tep_figure.png');
        if isequal(f, 0); return; end
        opts.file = fullfile(p, f);
        finish('save');
    end

    function finish(what)
        action = what;
        delete(fig);
    end
end
