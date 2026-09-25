% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [rule, counts, parts] = suggestSubjectRule(paths, splitAlnum, parts)
% SUGGESTSUBJECTRULE  A starting rule for which filename part is the subject.
%   [rule, counts] = SUGGESTSUBJECTRULE(paths) returns a rule for
%   subjectIdsFromRule, plus how many distinct values each part position takes
%   across the whole cohort:
%
%       counts.folders  1x3, per folder level (1 = the file's own folder)
%       counts.name     1xK, per name part counted from the start
%       counts.nameEnd  1xK, per name part counted from the end
%       counts.varies   true when files have different numbers of name parts,
%                       which is the only time counting from the end differs
%
%   Values are counted the way ids are compared (subjectKey), so S1 and S01
%   are one value.
%
%   The counts are the real product here. They are what the subject dialog
%   prints under each chip, and they let someone who has never seen a pattern
%   read the naming scheme off the screen: a part with 1 value is the same in
%   every file (a pipeline or study name), a part with as many values as there
%   are files is unique per file, and a subject is something in between.
%
%   THE SUGGESTION reads the cohort as a whole, never one file at a time. A
%   part is a candidate when it is present in every file, is not a date or
%   time, and takes more than one value but fewer than there are files - a
%   part that collapses nothing or everything is not an identity. Among the
%   candidates:
%     1. parts with a digit in (nearly) every file beat parts with one in
%        most, which beat parts without ('S01', '007' over 'pre' - a bare word
%        is far more often a condition than a person);
%     2. then the MOST distinct values wins. Sessions, runs, visits and
%        conditions come in twos and threes; people come in dozens. So 'ses1'
%        or 'run2' in front of the subject no longer steals the suggestion;
%     3. ties go to the name over a folder, counting from the start over from
%        the end, the leftmost name part, and the outermost folder.
%   With splitAlnum omitted, the name is first read as written, and split
%   between letters and digits only if nothing qualifies that way ('S01pre').
%   The returned rule.splitAlnum says which was used.
%
%   If nothing qualifies the rule is empty: one subject per file, the app's
%   safe default. The result is only ever a preselection in a dialog that
%   shows its effect; it is never applied unseen.
%
%   [...] = SUGGESTSUBJECTRULE(paths, splitAlnum, parts) takes filenameParts
%   of the paths at that splitAlnum instead of parsing them, and the third
%   output returns the parts used, so a caller asking again - the dialog, as
%   the split is toggled - parses the cohort once per split mode.
%
%   See also: subjectIdsFromRule, filenameParts, subjectKey, subjectRuleDialog

if ischar(paths) || isstring(paths); paths = cellstr(paths); end

if nargin < 2 || isempty(splitAlnum)
    parts = filenameParts(paths, false);
    [rule, counts] = suggestFor(parts, false);
    % Only the name depends on the split, so only a name suggestion can
    % appear by splitting.
    if isempty(rule.nameIdx) && isempty(rule.folderIdx)
        split = filenameParts(paths, true);
        [r2, c2] = suggestFor(split, true);
        if ~isempty(r2.nameIdx)
            rule = r2; counts = c2; parts = split;
        end
    end
else
    if nargin < 3 || isempty(parts)
        parts = filenameParts(paths, logical(splitAlnum));
    end
    [rule, counts] = suggestFor(parts, logical(splitAlnum));
end
end

% ── helpers ─────────────────────────────────────────────────────────────────

function [rule, counts] = suggestFor(parts, splitAlnum)
n    = numel(parts);
rule = struct('folderIdx', [], 'nameIdx', [], 'fromEnd', false, ...
              'splitAlnum', splitAlnum);
counts = struct('folders', zeros(1, 3), 'name', zeros(1, 0), ...
                'nameEnd', zeros(1, 0), 'varies', false);
if n == 0; return; end

len = arrayfun(@(p) numel(p.name), parts);
K   = max([len 0]);
counts.varies = numel(unique(len)) > 1;

% One row per candidate position.
KIND = 1; IDX = 2; NDIST = 3; DIGIT = 4; ALL = 5; DATE = 6;
NAME = 1; NAME_END = 2; FOLDER = 3;                    % values of KIND
cand = zeros(0, 6);
for f = 1:3
    vals = arrayfun(@(p) p.folders{f}, parts, 'UniformOutput', false);
    row  = describe(vals, false(1, n));
    counts.folders(f) = row(1);
    cand(end + 1, :) = [FOLDER f row]; %#ok<AGROW>
end
for fromEnd = [false true]
    c = zeros(1, K);
    for k = 1:K
        vals = repmat({''}, 1, n); dates = false(1, n);
        for i = 1:n
            pos = namePartPos(k, len(i), fromEnd);
            if pos >= 1 && pos <= len(i)
                vals{i}  = parts(i).name{pos};
                dates(i) = parts(i).isDate(pos);
            end
        end
        row  = describe(vals, dates);
        c(k) = row(1);
        cand(end + 1, :) = [NAME + fromEnd, k, row]; %#ok<AGROW>
    end
    if fromEnd; counts.nameEnd = c; else; counts.name = c; end
end
if ~counts.varies
    cand(cand(:, KIND) == NAME_END, :) = [];   % identical to counting from the start
end
if n < 2; return; end

ok = cand(:, ALL) == 1 & cand(:, DATE) == 0 & cand(:, NDIST) > 1 & cand(:, NDIST) < n;
cand = cand(ok, :);
if isempty(cand); return; end

% Rank: digits first, then most distinct, then name/start/leftmost. "Digits"
% has two levels because a part that is misaligned in a few files (S01_pre
% next to S02_pre_rerun, counted from the end) picks up words there - it
% still has digits in MOST files, and its stray words even inflate its count.
% A real id column has them in essentially all.
% Among folders a tie goes to the OUTER one: a person's folder holds their
% sessions, never the other way round (S01/ses1/, sub-01/ses-1/).
tier  = 3 - (cand(:, DIGIT) > 0.5) - (cand(:, DIGIT) >= 0.9);
order = cand(:, IDX);
isFolder = cand(:, KIND) == FOLDER;
order(isFolder) = -order(isFolder);
[~, o] = sortrows([tier, -cand(:, NDIST), cand(:, KIND), order]);
best = cand(o(1), :);
switch best(KIND)
    case NAME;     rule.nameIdx = best(IDX);
    case NAME_END; rule.nameIdx = best(IDX); rule.fromEnd = true;
    case FOLDER;   rule.folderIdx = best(IDX);
end
end

function row = describe(vals, dates)
% [distinct values, share with a digit, present in all, any is a date]
present = ~cellfun(@isempty, vals);
got     = vals(present);
digit   = sum(~cellfun(@isempty, regexp(got, '\d', 'once'))) / max(numel(vals), 1);
row     = [numel(unique(subjectKey(got))), digit, all(present), any(dates)];
end
