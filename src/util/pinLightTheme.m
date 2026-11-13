% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function fig = pinLightTheme(fig)
% PINLIGHTTHEME  Keep a figure in MATLAB's light theme, whatever the desktop or OS uses.
%   fig = PINLIGHTTHEME(fig)
%
%   The policy is "always light", decided here and nowhere else. nestapp's
%   colours are fixed for light, so following OS dark mode (R2025a+) left
%   windows half themed. Every window is created through appWindow or
%   outputFigure, which call this. Releases before R2025a have no figure
%   themes, so it does nothing there.
%
%   See also: appWindow, outputFigure

if isprop(fig, 'ThemeMode')   % ThemeMode and theme() arrived together (R2025a)
    theme(fig, 'light');
end
end
