% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function fig = pinLightTheme(fig)
% PINLIGHTTHEME  Keep a figure in MATLAB's light theme, whatever the desktop or OS uses.
%   fig = PINLIGHTTHEME(fig)
%
%   From R2025a figures follow the desktop theme (ThemeMode 'auto'), and on a
%   Mac set to match the system that means OS dark mode. nestapp is designed
%   for light - many of its colours are fixed - so under dark mode MATLAB
%   re-coloured only what nestapp had not pinned, and windows came out half
%   light, half dark: a light status bar with grey text, white text on light
%   disabled buttons. Exports could come out with dark axes too.
%
%   The policy is "always light", decided here and nowhere else. Every window
%   is created through appWindow or outputFigure, which call this, and
%   WindowThemeTest fails if one is not - fixing it once where windows are
%   born, rather than chasing colours one at a time afterwards.
%
%   Releases before R2025a have no figure themes and never follow dark mode,
%   so this does nothing there.
%
%   See also: appWindow, outputFigure

if isprop(fig, 'ThemeMode') && exist('theme', 'file')
    theme(fig, 'light');
end
end
