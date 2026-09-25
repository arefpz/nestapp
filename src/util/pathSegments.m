% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function seg = pathSegments(p)
% PATHSEGMENTS  A path's folder and file names, with separators and roots removed.
%   seg = PATHSEGMENTS(p) splits p at / or \ into a 1xK cellstr, dropping
%   empty pieces and a drive letter, so 'C:\data\S01\rec.set' gives
%   {'data', 'S01', 'rec.set'}. The drive is where the data happens to live,
%   never part of what names it.
%
%   See also: shortPath, filenameParts, uniqueFileLabels

seg = strsplit(strrep(char(p), '\', '/'), '/');
seg = seg(~cellfun(@isempty, seg));
if ~isempty(seg) && ~isempty(regexp(seg{1}, '^[A-Za-z]:$', 'once'))
    seg(1) = [];
end
end
