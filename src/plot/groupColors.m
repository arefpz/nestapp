% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function C = groupColors(n)
% GROUPCOLORS  The default theme's colours for n groups, cycled.
%   C = GROUPCOLORS(n) returns an n-by-3 RGB matrix. Kept for scripts; plots
%   take their styles through resolveGroupStyle, which also varies line style
%   and texture past the palette (see groupStyles).
%
%   See also: groupStyles, plotThemes

if nargin < 1 || isempty(n); n = 1; end
n = max(0, round(n));
themes = plotThemes();
P = themes(1).colors;
C = P(mod(0:n - 1, size(P, 1)) + 1, :);
end
