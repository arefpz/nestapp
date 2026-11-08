% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function C = groupColors(n)
% GROUPCOLORS  Default-theme colours for n groups, in order.
%   C = GROUPCOLORS(n) returns an n-by-3 RGB matrix, one row per group.
%
%   The drawers' fallback when a caller passes no colours, and a convenience
%   for scripts. The app instead passes colours from groupStyles, which keys
%   them by group NAME and follows the user's theme; this gives what that
%   would for groups seen for the first time, in the default theme.
%
%   The default is Okabe & Ito's colour-blind-safe set with amber in place of
%   their yellow, which was nearly invisible on white. Its first colour is a
%   dark neutral, so the common single-group case draws in near-black. Past
%   its seven colours, colours repeat - groupStyles is what then tells those
%   groups apart by line style.
%
%   See also: groupStyles, plotThemes, drawTEPOverlay

if nargin < 1 || isempty(n); n = 1; end
n = max(0, round(n));
if n == 0; C = zeros(0, 3); return; end
s = groupStyles(arrayfun(@(k) sprintf('group %d', k), 1:n, 'UniformOutput', false), '');
C = s.colors;
end
