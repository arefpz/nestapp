% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function use = groupedCacheEntries(cache, entries)
% GROUPEDCACHEENTRIES  The files a comparison uses: grouped, and loaded.
%   use = GROUPEDCACHEENTRIES(cache, entries) joins the files table to the
%   cache from loadReducedSets and keeps the entries that have a group and
%   whose file loaded. use(j) carries .ci (the cache index), .subject, .group,
%   .label (file name) and .nTrials, plus empty .chanAvg and .curve slots that
%   groupCurves fills.
%
%   Public so the app counts exactly the files groupCurves will use - the cap
%   prompt states how many files each rule keeps, and a second copy of this
%   filter could disagree with the result about which files those are.
%
%   Holds a cache INDEX, not the cached arrays. Copy-on-write means copying
%   never memcopied the 63x2000 trial averages, so the win is avoided
%   REALLOCATION: growing a struct array whose elements each carried them, then
%   filtering it again by montage. Measured 0.20 -> 0.09 ms at 12 files and
%   2.9 -> 1.4 ms at 200, on a path that re-runs on every ROI change.
%
%   See also: groupCurves, capSummary, loadReducedSets

use = struct('ci', {}, 'subject', {}, 'group', {}, 'label', {}, ...
             'nTrials', {}, 'chanAvg', {}, 'curve', {});
if isempty(entries) || isempty(cache); return; end

% One lookup for all entries; ismember takes the first matching cache row.
[~, ci] = ismember({entries.path}, {cache.path});
grouped = ~cellfun(@isempty, {entries.group});
loaded  = false(size(ci));
loaded(ci > 0) = [cache(ci(ci > 0)).ok];
idxE = find(grouped & loaded);
if isempty(idxE); return; end

use = repmat(struct('ci', 0, 'subject', '', 'group', '', 'label', '', ...
                    'nTrials', 1, 'chanAvg', [], 'curve', []), 1, numel(idxE));
for j = 1:numel(idxE)
    i = idxE(j);
    [~, nm, ex]    = fileparts(entries(i).path);
    use(j).ci      = ci(i);
    use(j).subject = char(entries(i).subject);
    use(j).group   = char(entries(i).group);
    use(j).label   = [nm ex];
    use(j).nTrials = max(cache(ci(i)).nTrials, 1);
end
end
