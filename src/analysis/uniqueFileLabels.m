% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function labels = uniqueFileLabels(paths)
% UNIQUEFILELABELS  One readable, unique label per file: the one-subject-per-file ids.
%   labels = UNIQUEFILELABELS(paths) returns the basename, with as much of the
%   parent path prepended as it takes to be unique. The cohort has the same
%   basename in two pipeline output folders, so basename alone is not an
%   identity - and two files sharing a label would silently merge into one
%   "subject".
%
%   See also: exploreDataset, subjectIdsFromRule, pathSegments

if ischar(paths) || isstring(paths); paths = cellstr(paths); end

n      = numel(paths);
labels = cell(1, n);
parts  = cell(1, n);
for i = 1:n
    seg = pathSegments(paths{i});
    [~, seg{end}] = fileparts(seg{end});
    parts{i}  = seg;
    labels{i} = seg{end};
end

depth = 1;
while numel(unique(labels)) < n && depth < 8
    depth = depth + 1;
    for i = 1:n
        labels{i} = strjoin(parts{i}(max(1, end - depth + 1):end), '/');
    end
end
end
