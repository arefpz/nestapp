% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function sty = resolveGroupStyle(opts, names)
% RESOLVEGROUPSTYLE  The group styles a drawer should use, decided in one place.
%   sty = RESOLVEGROUPSTYLE(opts, names) returns a groupStyles struct for the
%   groups named, from:
%     opts.groupStyle  a groupStyles struct - what the app passes, already
%                      following the user's theme and picked colours
%     otherwise        the default theme for these names, so a script gets the
%                      same line styles and textures past the palette as the app
%     opts.colors      if given, n x 3 colours that replace the palette (kept
%                      for scripts that pass their own colours)
%
%   Drawers used to default and size-check colours, line styles, textures and
%   width separately, so the four could come from different sources. They now
%   read one struct.
%
%   See also: groupStyles, drawTEPOverlay, drawWindowBars, drawTEPTopo

names = cellstr(names);
n = numel(names);
if isfield(opts, 'groupStyle') && ~isempty(opts.groupStyle) && ...
        size(opts.groupStyle.colors, 1) >= n
    sty = opts.groupStyle;
else
    sty = groupStyles(names, '');
end
if isfield(opts, 'colors') && ~isempty(opts.colors)
    c = opts.colors;
    sty.colors = c(mod(0:n - 1, size(c, 1)) + 1, :);
end
end
