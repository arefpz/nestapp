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
%   s.colors (n x 3), s.lineStyles and s.barPatterns (1 x n cellstr), s.lineWidth
%   memory     updated; pass it back next time
%
%   Each group holds a SLOT, and the slot decides its style. A group keeps its
%   slot for as long as the caller keeps the memory, so adding, removing or
%   reordering groups never recolours the others: "Post-TBS" stays the same
%   orange in every figure of a session. A new group takes the lowest slot no
%   current group holds.
%
%   Slot k maps to theme colour mod(k, nColours). In a colour theme the line
%   style and bar texture change only when the colours run out - the second
%   lap of colours is dashed and diagonally hatched, not darker or invented,
%   because more hues than a palette was designed for cannot be told apart. In
%   a theme with styleFirst (greyscale) style and texture vary from slot 0.
%
%   See also: plotThemes, groupColors, hatchBars

if nargin < 3; overrides = []; end
if nargin < 4 || isempty(memory); memory = struct('name', {}, 'slot', {}); end
names = cellstr(names);
t = pickTheme(themeName);

n = numel(names);
slots = zeros(1, n);
known = {memory.name};
for i = 1:n
    k = find(strcmp(known, names{i}), 1);
    if ~isempty(k); slots(i) = memory(k).slot; else; slots(i) = -1; end
end
% Two current groups cannot share a slot, even if a stale memory says so.
for i = 1:n
    if slots(i) >= 0 && any(slots(1:i-1) == slots(i)); slots(i) = -1; end
end
for i = find(slots < 0)
    free = 0;
    while any(slots == free); free = free + 1; end
    slots(i) = free;
    memory(end+1) = struct('name', names{i}, 'slot', free); %#ok<AGROW>
end
for i = 1:n
    k = find(strcmp({memory.name}, names{i}), 1);
    memory(k).slot = slots(i);
end

nC = size(t.colors, 1);
if t.styleFirst
    lap = slots;
else
    lap = floor(slots / nC);
end
s = struct();
s.colors      = t.colors(mod(slots, nC) + 1, :);
s.lineStyles  = t.lineStyles(mod(lap, numel(t.lineStyles)) + 1);
s.barPatterns = t.barPatterns(mod(lap, numel(t.barPatterns)) + 1);
s.lineWidth   = t.lineWidth;

if ~isempty(overrides)
    for i = 1:n
        k = find(strcmp({overrides.name}, names{i}), 1);
        if ~isempty(k); s.colors(i, :) = overrides(k).rgb; end
    end
end
end

function t = pickTheme(name)
themes = plotThemes();
t = themes(1);
if isempty(name); return; end
k = find(strcmpi({themes.name}, char(name)), 1);
if ~isempty(k); t = themes(k); end
end
