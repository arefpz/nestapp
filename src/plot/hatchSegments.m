% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [lines, dots] = hatchSegments(rect, pattern, spacing)
% HATCHSEGMENTS  The geometry of a fill texture inside one rectangle.
%   [lines, dots] = HATCHSEGMENTS(rect, pattern, spacing)
%
%   rect     [x0 x1 y0 y1] in an ISOTROPIC space (pixels), so angles are true
%   pattern  'solid' | 'diagonal' | 'backdiagonal' | 'crosshatch' |
%            'horizontal' | 'vertical' | 'dots'
%   spacing  distance between neighbouring lines or dots, same units
%
%   lines    2 x M [x; y], segments separated by NaN columns - one plot call
%   dots     2 x K [x; y] dot centres
%
%   Every segment is clipped exactly to the rectangle, and the families are
%   anchored to the rectangle's corner, so neighbouring bars of the same group
%   look alike. Pure geometry, no graphics: hatchBars draws it.
%
%   See also: hatchBars, drawWindowBars

lines = zeros(2, 0);
dots  = zeros(2, 0);
x0 = min(rect(1:2)); x1 = max(rect(1:2));
y0 = min(rect(3:4)); y1 = max(rect(3:4));
if x1 - x0 <= 0 || y1 - y0 <= 0 || spacing <= 0; return; end

switch lower(char(pattern))
    case 'solid'
        return
    case 'diagonal',     angles = 45;
    case 'backdiagonal', angles = -45;
    case 'crosshatch',   angles = [45 -45];
    case 'horizontal',   angles = 0;
    case 'vertical',     angles = 90;
    case 'dots'
        [gx, gy] = meshgrid(x0 + spacing / 2 : spacing : x1, ...
                            y0 + spacing / 2 : spacing : y1);
        dots = [gx(:)'; gy(:)'];
        return
    otherwise
        error('nestapp:badPattern', 'Unknown bar pattern "%s".', char(pattern));
end

corners = [0, x1 - x0, x1 - x0, 0; 0, 0, y1 - y0, y1 - y0];   % relative to (x0, y0)
for a = angles
    d = [cosd(a); sind(a)];           % along the lines
    nrm = [-d(2); d(1)];              % across them
    c = nrm' * corners;               % offsets of the corners across the family
    for off = ceil(min(c) / spacing) * spacing : spacing : max(c)
        p0 = [x0; y0] + off * nrm;
        seg = clipToRect(p0, d, x0, x1, y0, y1);
        if ~isempty(seg); lines = [lines, seg, [NaN; NaN]]; end %#ok<AGROW>
    end
end
end

function seg = clipToRect(p0, d, x0, x1, y0, y1)
% The part of the line p0 + t*d inside the rectangle (Liang-Barsky).
tLo = -Inf; tHi = Inf;
lims = [x0 x1; y0 y1];
for k = 1:2
    if abs(d(k)) < 1e-12
        if p0(k) < lims(k, 1) || p0(k) > lims(k, 2); seg = []; return; end
    else
        t = (lims(k, :) - p0(k)) / d(k);
        tLo = max(tLo, min(t)); tHi = min(tHi, max(t));
    end
end
if tHi - tLo <= 1e-9; seg = []; return; end
seg = [p0 + tLo * d, p0 + tHi * d];
end
