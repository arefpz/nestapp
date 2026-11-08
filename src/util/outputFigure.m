% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function fig = outputFigure(varargin)
% OUTPUTFIGURE  Create a classic figure for output: white and light-themed.
%   fig = OUTPUTFIGURE(...) takes the arguments figure takes; the background
%   defaults to white and the caller's own 'Color' still wins.
%
%   For figures that become files - publication figures, the QC figure, PDF
%   report pages - and for offscreen scratch figures. A figure exported on a
%   dark-mode machine must look the same as one exported anywhere else, so it
%   is pinned light (see pinLightTheme). Classic rather than uifigure because
%   every export path drops UI components.
%
%   See also: appWindow, pinLightTheme, publicationFigure

fig = pinLightTheme(figure('Color', 'w', varargin{:}));
end
