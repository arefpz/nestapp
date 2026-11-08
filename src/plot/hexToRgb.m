% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function rgb = hexToRgb(varargin)
% HEXTORGB  '#RRGGBB' strings to an n x 3 RGB matrix in [0, 1].
%   rgb = HEXTORGB('#0077BB', '#EE7733', ...)
%
%   Written out rather than using MATLAB's hex2rgb, which needs R2024a.
%
%   See also: plotThemes, candidatePalettes

rgb = zeros(numel(varargin), 3);
for k = 1:numel(varargin)
    h = varargin{k}(2:end);
    rgb(k, :) = [hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))] / 255;
end
end
