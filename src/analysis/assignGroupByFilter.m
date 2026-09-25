% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [entries, nMatched] = assignGroupByFilter(entries, pattern, groupName, rootFolder)
% ASSIGNGROUPBYFILTER  Put every entry whose path matches a filter into a group.
%   [entries, nMatched] = ASSIGNGROUPBYFILTER(entries, pattern, groupName)
%   assigns groupName to each entry whose path contains pattern, and returns how
%   many were matched.
%
%   [...] = ASSIGNGROUPBYFILTER(entries, pattern, groupName, rootFolder) matches
%   against the path relative to rootFolder instead of the absolute path, so a
%   pattern cannot accidentally match something in the directories above the
%   data - typing "pre" should not match every file because the cohort happens
%   to live under C:\preprocessed.
%
%   Matching is a case-insensitive substring test on the path with separators
%   normalised to '/', deliberately identical to the filter in selectDataTree
%   (`contains(allRelLower, lower(query))`). That is the whole point: the filter
%   the user already types to FIND a condition's files is the same expression
%   that NAMES the group, so what you searched is what you grouped. Keeping the
%   two matchers in step is a correctness requirement, not a convenience.
%
%   An empty pattern matches nothing rather than everything - "I have not typed
%   a filter yet" must not silently sweep the entire cohort into one group.
%
%   Inputs:
%     entries    - struct array with at least a .path field; .group is written.
%     pattern    - char/string filter text.
%     groupName  - char/string group to assign to the matches.
%     rootFolder - optional char; paths are made relative to it before matching.
%
%   Outputs:
%     entries  - the input with .group set on matching elements.
%     nMatched - number of entries assigned.
%
%   See also: selectDataTree, subjectIdsFromRule, datasetSummary

if nargin < 4, rootFolder = ''; end
nMatched = 0;

if isempty(entries) || ~isfield(entries, 'path')
    return
end

pattern = lower(strtrim(char(pattern)));
if isempty(pattern)
    return
end

groupName = char(groupName);
root      = normalisePath(char(rootFolder));

for i = 1:numel(entries)
    hay = normalisePath(char(entries(i).path));
    if ~isempty(root)
        hay = makeRelative(hay, root);
    end
    if contains(lower(hay), pattern)
        entries(i).group = groupName;
        nMatched = nMatched + 1;
    end
end
end

% ── helpers ─────────────────────────────────────────────────────────────────

function p = normalisePath(p)
p = strrep(p, '\', '/');
end

function rel = makeRelative(p, root)
% Strip a leading root from p. Case-insensitive because Windows paths are, and
% the cohort is processed on Windows.
root = regexprep(root, '/+$', '');
if strncmpi(p, [root '/'], numel(root) + 1)
    rel = p(numel(root) + 2:end);
else
    rel = p;
end
end
