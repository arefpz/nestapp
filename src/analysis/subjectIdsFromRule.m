% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [ids, matched, nRespelled] = subjectIdsFromRule(paths, rule, parts, labels)
% SUBJECTIDSFROMRULE  Subject id per file, from the filename parts the user picked.
%   [ids, matched, nRespelled] = SUBJECTIDSFROMRULE(paths, rule) applies
%
%       rule.folderIdx   folders to take: 1 the file's own folder, 2 the one
%                        above, 3 the one above that
%       rule.nameIdx     name parts to take (see filenameParts)
%       rule.fromEnd     count name parts from the END of the name, for names
%                        whose front varies in length (tms_S01 vs
%                        resting_state_S01) - the subject is then "the last part"
%       rule.splitAlnum  split letters from digits before counting
%
%   The chosen parts are joined with '_' in the order they appear in the path,
%   outermost folder first and then left to right along the name, so the id
%   never depends on the order they were clicked in.
%
%   Ids that differ only in case or zero padding (s01, S01, S1) are one person
%   and get one spelling: the alphabetically first. nRespelled counts the
%   files whose id was changed to that spelling, so the caller can say it did.
%
%   A file that does not HAVE one of the chosen parts (its name is shorter, or
%   it is not that deep) keeps its one-subject-per-file label and gets
%   matched=false. It is never folded into a neighbour: a wrong merge corrupts
%   n without looking wrong, so the file stays on its own and the dialog flags
%   it for the user to look at.
%
%   An empty rule picks nothing and gives one subject per file.
%
%   [...] = SUBJECTIDSFROMRULE(paths, rule, parts, labels) skips the parsing:
%   parts is filenameParts of every path at rule.splitAlnum, labels is
%   uniqueFileLabels(paths). Both depend only on the paths, so a caller that
%   applies rule after rule to one cohort - the dialog, on every click -
%   computes them once. Either may be [] to have it computed here.
%
%   Outputs:
%     ids        - 1xN cellstr of subject ids.
%     matched    - 1xN logical, false where the file lacked a chosen part.
%     nRespelled - files whose id was unified with another spelling.
%
%   See also: filenameParts, suggestSubjectRule, subjectKey, subjectRuleDialog

if ischar(paths) || isstring(paths); paths = cellstr(paths); end
rule = fillDefaults(rule, struct('folderIdx', [], 'nameIdx', [], ...
                                 'fromEnd', false, 'splitAlnum', false));
if nargin < 4 || isempty(labels); labels = uniqueFileLabels(paths); end

n          = numel(paths);
ids        = labels;
matched    = true(1, n);
nRespelled = 0;

folderIdx = sort(unique(rule.folderIdx(:)'), 'descend');   % outermost first
nameIdx   = unique(rule.nameIdx(:)');
if isempty(folderIdx) && isempty(nameIdx)
    return
end
if nargin < 3 || isempty(parts)
    parts = filenameParts(paths, rule.splitAlnum);
end

for i = 1:n
    K   = numel(parts(i).name);
    pos = sort(namePartPos(nameIdx, K, rule.fromEnd));     % left to right
    if any(pos < 1 | pos > K)
        matched(i) = false;
        continue
    end
    picked = [parts(i).folders(folderIdx), parts(i).name(pos)];
    if any(cellfun(@isempty, picked))
        matched(i) = false;
        continue
    end
    ids{i} = strjoin(picked, '_');
end

% One spelling per person. Only keys with more than one spelling need work.
got       = ids(matched);
[~, ~, g] = unique(subjectKey(got));
[~, ~, s] = unique(got);
nSpell    = accumarray(g(:), s(:), [], @(x) numel(unique(x)));
for k = find(nSpell > 1)'
    spellings   = unique(got(g == k));
    nRespelled  = nRespelled + sum(g == k & ~strcmp(got(:), spellings{1}));
    got(g == k) = spellings(1);
end
ids(matched) = got;
end
