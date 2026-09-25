% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function s = shortPath(p)
% SHORTPATH  Enough of a path's tail to identify the file without a 200-character cell.
%   s = SHORTPATH(p) keeps the file name and the two folders above it,
%   prefixed with '.../' when anything was cut. Two folders because the
%   cohort has the same basename in two pipeline output folders, one level up
%   from a shared 'data' folder.
%
%   See also: pathSegments, uniqueFileLabels

seg = pathSegments(p);
if numel(seg) <= 3
    s = strjoin(seg, '/');
else
    s = ['.../' strjoin(seg(end-2:end), '/')];
end
end
