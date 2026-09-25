% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function parts = filenameParts(path, splitAlnum)
% FILENAMEPARTS  Break a data file's path into the parts a subject rule picks from.
%   parts = FILENAMEPARTS(path) returns a struct with
%     .folders  1x3 cellstr: the file's own folder, the one above it, and the
%               one above that ('' where the path is not that deep). Counted
%               upward so that "1" always means the folder the file sits in.
%     .name     1xK cellstr: the basename cut into parts
%     .isDate   1xK logical: which name parts are a date or a time of day
%
%   The basename is cut at anything that is not a letter or a digit, so
%   _ - . spaces brackets commas and the like all separate, and letters with
%   accents stay whole ('Müller'). The one exception: a date (2024-03-15,
%   20240315, 15.03.2024) or a time (14-31-05) is ONE part, not three, and is
%   flagged in .isDate. A date is never a person, and cutting it up would
%   scatter digits that look like ids.
%
%   parts = FILENAMEPARTS(paths) with a cellstr returns a 1xN struct array.
%
%   parts = FILENAMEPARTS(path, true) also splits each non-date part where
%   letters meet digits, so 'S01pre' gives {'S', '01', 'pre'}, for names
%   written with no separators at all.
%
%   Splitting depends on the name alone, so the same name always gives the
%   same parts - there is no guessing in here. The guessing, such as it is,
%   lives in suggestSubjectRule and is always shown before it applies.
%
%   See also: subjectIdsFromRule, suggestSubjectRule, subjectRuleDialog

if nargin < 2; splitAlnum = false; end
if iscell(path)
    % A cellstr gives a 1xN struct array, one element per path.
    parts = cellfun(@(p) filenameParts(p, splitAlnum), path(:)', 'UniformOutput', false);
    parts = [parts{:}];
    return
end

seg = pathSegments(path);

base = '';
if ~isempty(seg)
    [~, base] = fileparts(seg{end});
    seg(end) = [];
end

parts.folders = {'', '', ''};
for f = 1:min(3, numel(seg))
    parts.folders{f} = seg{end - f + 1};
end

DATE = ['(?:19|20)\d\d[-._]?[01]\d[-._]?[0-3]\d' ...     % 2024-03-15, 20240315
        '|[0-3]\d[-._][01]\d[-._](?:19|20)\d\d' ...      % 15.03.2024
        '|[0-2]\d[-.:][0-5]\d[-.:][0-5]\d'];             % 14-31-05
name   = regexp(base, [DATE '|[^\W_]+'], 'match');
isDate = ~cellfun(@isempty, regexp(name, ['^(?:' DATE ')$'], 'once'));

if splitAlnum
    runs  = cell(1, numel(name));
    flags = cell(1, numel(name));
    for k = 1:numel(name)
        if isDate(k)
            runs{k} = name(k);
        else
            runs{k} = regexp(name{k}, '\d+|[^\W\d_]+', 'match');
        end
        flags{k} = repmat(isDate(k), 1, numel(runs{k}));
    end
    name   = [runs{:}];
    isDate = [flags{:}];
end
parts.name   = reshape(name, 1, []);
parts.isDate = reshape(logical(isDate), 1, []);
end
