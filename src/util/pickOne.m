% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function idx = pickOne(title, prompt, labels, anchor)
% PICKONE  Ask which one of a list, modally.
%   idx = PICKONE(title, prompt, labels, anchor)
%
%   labels   cellstr of choices, shown in the order given
%   anchor   optional figure to centre on
%
%   Returns the 1-based index chosen, or [] if the user cancelled - including
%   by closing the window. A single-element list still asks, because the
%   caller that wants "don't ask when there is only one" can check numel
%   itself and say so, whereas a dialog that sometimes does not appear is a
%   surprise.
%
%   Modal discipline: waitfor(fig) and a plain delete, never uiwait/uiresume.
%   A nested uiconfirm or uialert inside a figure already sitting in uiwait
%   leaves uiresume unable to release it, and the app hangs - the trap
%   selectDataTree documents and exploreFilesTable was caught by.
%
%   See also: plotOptionsDialog, figureExportDialog, centreOn

if nargin < 4; anchor = []; end
idx = [];

labels = cellstr(labels);
labels = labels(:)';
if isempty(labels); return; end

W = 460;
listH = min(max(numel(labels), 3), 12) * 18 + 8;
H = listH + 108;

fig = appWindow('Name', title, 'Position', centreOn(anchor, W, H), 'Resize', 'off');
fig.CloseRequestFcn = @(src, ~) delete(src);   % X == cancel, and nothing else

uilabel(fig, 'Position', [16 H-34 W-32 22], 'Text', prompt);

list = uilistbox(fig, 'Position', [16 56 W-32 listH], ...
                 'Items', labels, 'ItemsData', 1:numel(labels), ...
                 'Multiselect', 'off', 'Value', 1, ...
                 'DoubleClickedFcn', @(~,~) onOk());

uibutton(fig, 'Text', 'Cancel', 'Position', [W-206 16 90 26], ...
    'ButtonPushedFcn', @(~,~) delete(fig));
uibutton(fig, 'Text', 'Open', 'Position', [W-106 16 90 26], ...
    'ButtonPushedFcn', @(~,~) onOk());

% Announce that the dialog is about to block, the same way roiPicker and
% exploreFilesTable do. Anything driving it otherwise has to guess when the
% wait has started, and a release issued before that point is silently a
% no-op. This is what lets tests use driveModalDialog instead of a fixed
% sleep long enough to hope the window exists.
setappdata(fig, 'nestappModalReady', true);
waitfor(fig);

    function onOk()
        idx = list.Value;
        delete(fig);
    end
end
