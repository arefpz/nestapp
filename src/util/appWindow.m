% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function fig = appWindow(varargin)
% APPWINDOW  Create a nestapp window: a uifigure in the light theme.
%   fig = APPWINDOW(...) takes exactly the arguments uifigure takes.
%
%   Every app window - the main one and every dialog - is created here, so
%   none can follow OS dark mode into a half-themed mess (see pinLightTheme).
%   The theme is pinned before anything is drawn, and what is created inside
%   the window - uiaxes, uialert, uiconfirm, uiprogressdlg - inherits it.
%
%   See also: outputFigure, pinLightTheme

fig = pinLightTheme(uifigure(varargin{:}));
end
