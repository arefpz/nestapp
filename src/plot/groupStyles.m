% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [s, memory] = groupStyles(names, themeName, overrides, memory)
% GROUPSTYLES  Colour, line style and bar texture for each group, kept by group NAME.
%   [s, memory] = GROUPSTYLES(names, themeName, overrides, memory)
%
%   names      cellstr of group names, in drawing order
%   themeName  a plotThemes name; '' or unknown gives the default theme
%   overrides  struct array .name .rgb - colours the user picked for a group,
%              which win over the theme. May be empty
%   memory     struct array .name .slot from a previous call, or empty
%
%   s.colors (n x 3), s.lineStyles and s.barPatterns (1 x n cellstr),
%   s.lineWidth, s.contrast (colour for non-group lines, e.g. a difference)
%   memory     updated; pass it back next time
%
%   Each group holds a SLOT, and the slot decides its style. A group keeps its
%   slot for as long as the caller keeps the memory, so adding, removing or
%   reordering groups never restyles the others. A new group takes the lowest
%   slot no current group holds.
%
%   Slot k takes theme colour mod(k, nColours). In a colour theme line style
%   and texture change only when the colours run out - the second lap is
%   dashed and hatched rather than darker or invented, because more hues than a
%   palette was designed for cannot be told apart. With styleFirst
%   (greyscale) they vary from slot 0.
%
%   See also: plotThemes, resolveGroupStyle, hatchBars

if nargin < 3; overrides = []; end
if nargin < 4 || isempty(memory); memory = struct('name', {}, 'slot', {}); end
names = cellstr(names);
t = pickTheme(themeName);

% Remembered slots first, so a new group listed ahead of a known one cannot
% take its slot; then new groups (and any clash) take the lowest free slot.
n = numel(names);
slots = -ones(1, n);
[~, at] = ismember(names, {memory.name});
for i = find(at)
    if ~any(slots == memory(at(i)).slot); slots(i) = memory(at(i)).slot; end
end
for i = find(slots < 0)
    free = 0;
    while any(slots == free); free = free + 1; end
    slots(i) = free;
    if at(i) > 0
        memory(at(i)).slot = free;
    else
        memory(end+1) = struct('name', names{i}, 'slot', free); %#ok<AGROW>
    end
end

nC  = size(t.colors, 1);
lap = slots;
if ~t.styleFirst; lap = floor(slots / nC); end
s = struct();
s.colors      = t.colors(mod(slots, nC) + 1, :);
s.lineStyles  = t.lineStyles(mod(lap, numel(t.lineStyles)) + 1);
s.barPatterns = t.barPatterns(mod(lap, numel(t.barPatterns)) + 1);
s.lineWidth   = t.lineWidth;
s.contrast    = t.contrast;

for i = 1:numel(overrides)
    k = find(strcmp(names, overrides(i).name), 1);
    if ~isempty(k); s.colors(k, :) = overrides(i).rgb; end
end
end

function t = pickTheme(name)
themes = plotThemes();
k = find(strcmpi({themes.name}, char(name)), 1);
if isempty(k); k = 1; end
t = themes(k);
end
