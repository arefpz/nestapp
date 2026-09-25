% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [ids, rule] = subjectRuleDialog(paths, opts)
% SUBJECTRULEDIALOG  Let the user click which part of the filename is the subject.
%   [ids, rule] = SUBJECTRULEDIALOG(paths) opens a modal that shows one example
%   filename broken into clickable parts, and returns a subject id per path
%   plus the rule that produced them - or {} and [] if cancelled.
%
%   Nobody should need to write a pattern to say "the subject is the bit
%   before the first underscore". So the name is shown as chips
%
%       [ rtmsct001 ] [ 1 ] [ pre ] [ SPL ] [ tesa ]
%         8 values     3     2       2       1
%
%   and the user clicks the one that is the person. The number under each chip
%   is how many different values that part takes across ALL the files, which is
%   the whole guide: 1 means it never changes, as many as there are files means
%   it is unique per file, and the subject is usually the one in between.
%
%   Everything updates on every click - the table of file -> subject, sorted so
%   files of one subject sit together, and the count line stating how n
%   changes. The dialog IS the preview; nothing applies until "Apply".
%
%   The rule is by position and nothing else (see subjectIdsFromRule), so it is
%   deterministic, and exploreDataset's opts.subjectRule reproduces it in a
%   script.
%
%   opts:
%     .parent   figure to centre over
%     .rule     a previous rule to start from; otherwise suggestSubjectRule's
%
%   See also: subjectIdsFromRule, suggestSubjectRule, filenameParts, exploreFilesTable

if nargin < 2; opts = struct(); end
opts = fillDefaults(opts, struct('parent', [], 'rule', []));

if ischar(paths) || isstring(paths); paths = cellstr(paths); end
paths = paths(:)';
ids   = {};
rule  = [];
if isempty(paths); return; end

nFiles   = numel(paths);
accepted = false;

% Everything below depends only on the paths, or on the paths and the split
% mode, so it is worked out once here rather than on every click.
labels = uniqueFileLabels(paths);
short  = cellfun(@shortPath, paths, 'UniformOutput', false);
cache  = {[], []};                        % filenameParts per split mode

% The automatic suggestion also decides whether letters and digits need
% splitting; a previous rule, if given, decides instead.
[autoRule, counts, autoParts] = suggestSubjectRule(paths);
cache{autoRule.splitAlnum + 1} = autoParts;
if isstruct(opts.rule)
    work = fillDefaults(opts.rule, struct('folderIdx', [], 'nameIdx', [], ...
                                          'fromEnd', false, 'splitAlnum', false));
    work.fromEnd    = logical(work.fromEnd);
    work.splitAlnum = logical(work.splitAlnum);
    [~, counts] = suggestSubjectRule(paths, work.splitAlnum, partsFor(work.splitAlnum));
else
    work = autoRule;
end
exampleIdx = typicalExample(partsFor(work.splitAlnum));
ex         = [];                          % the example's parts, set by buildChips
curIds     = {};
curMatched = [];

W = 820; H = 620;
MUTED = [0.35 0.38 0.43];

fig = uifigure('Name', 'Subjects from filenames', 'Resize', 'off', ...
               'Position', centreOn(opts.parent, W, H), 'WindowStyle', 'modal');
% Same exit discipline as exploreFilesTable: delete + waitfor, never uiresume.
fig.CloseRequestFcn = @(src, ~) delete(src);

try
g = uigridlayout(fig, [9 1], 'RowHeight', {22, 22, 150, 22, 22, '1x', 22, 30, 30}, ...
                 'RowSpacing', 6, 'Padding', [12 12 12 12]);

exRow = uigridlayout(g, [1 2], 'ColumnWidth', {100, '1x'}, 'Padding', [0 0 0 0]);
uilabel(exRow, 'Text', 'Example file:');
uidropdown(exRow, 'Items', labels, 'ItemsData', 1:nFiles, ...
    'Value', exampleIdx, 'ValueChangedFcn', @(src, ~) onExample(src.Value));

uilabel(g, 'FontWeight', 'bold', ...
    'Text', 'Click the part of the name that identifies the person (click more than one to combine them).');

chipPanel = uipanel(g, 'BorderType', 'line');
chipGrid  = [];

optRow = uigridlayout(g, [1 2], 'ColumnWidth', {'1x', '1x'}, 'Padding', [0 0 0 0]);
splitBox = uicheckbox(optRow, 'Value', work.splitAlnum, ...
    'Text', 'Split letters from numbers (for names like S01pre)', ...
    'ValueChangedFcn', @(src, ~) resetTo(logical(src.Value)));
endBox = uicheckbox(optRow, 'Value', work.fromEnd, ...
    'Text', 'Count name parts from the end', ...
    'Tooltip', ['For names whose front varies in length, like tms_S01 and ' ...
                'resting_state_S01: the subject is then "the last part".'], ...
    'ValueChangedFcn', @(src, ~) onFromEnd(src.Value));

ruleLabel = uilabel(g, 'FontColor', MUTED);

tbl = uitable(g, 'ColumnName', {'File', 'Subject', 'Files for this subject'}, ...
    'ColumnWidth', {'1x', 200, 150}, 'RowName', {}, ...
    'ColumnEditable', [false false false]);

countLabel = uilabel(g, 'FontWeight', 'bold');

uilabel(g, 'FontSize', 11, 'FontColor', MUTED, 'WordWrap', 'on', ...
    'Text', ['The number under each part is how many different values it has across ' ...
             'all files. Highlighted rows lack the chosen part and stay their own ' ...
             'subject - pick another part, or fix them in the table afterwards. ' ...
             'S1, S01 and s01 count as the same subject.']);

btnRow = uigridlayout(g, [1 4], 'ColumnWidth', {120, '1x', 120, 120}, 'Padding', [0 0 0 0]);
uibutton(btnRow, 'Text', 'Suggest', 'Tooltip', 'Go back to the automatic suggestion', ...
    'ButtonPushedFcn', @(~, ~) resetTo(autoRule.splitAlnum));
uilabel(btnRow, 'Text', '');
uibutton(btnRow, 'Text', 'Apply', 'ButtonPushedFcn', @(~, ~) accept());
uibutton(btnRow, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) delete(fig));

syncOptions();
buildChips();
refresh();
setappdata(fig, 'nestappModalReady', true);
waitfor(fig);
catch ME
    if isvalid(fig); delete(fig); end
    rethrow(ME);
end

if accepted
    ids  = curIds;
    rule = work;
end
if isvalid(fig); delete(fig); end

% ── nested ───────────────────────────────────────────────────────────────
    function p = partsFor(split)
        if isempty(cache{split + 1})
            cache{split + 1} = filenameParts(paths, split);
        end
        p = cache{split + 1};
    end

    function buildChips()
        if ~isempty(chipGrid) && isvalid(chipGrid); delete(chipGrid); end
        cohort = partsFor(work.splitAlnum);
        ex = cohort(exampleIdx);
        K  = numel(ex.name);
        nCols = 1 + max(K, 3);
        chipGrid = uigridlayout(chipPanel, [4 nCols], ...
            'RowHeight', {26, 16, 26, 16}, ...
            'ColumnWidth', [{60}, repmat({'fit'}, 1, nCols - 1)], ...
            'RowSpacing', 2, 'ColumnSpacing', 6, 'Scrollable', 'on');

        % Folders shown outermost first, as they read in the path; level 1 is
        % the folder the file sits in.
        rowLabel(1, 'Folder');
        for f = 3:-1:1
            if isempty(ex.folders{f}); continue; end
            chip(1, 4 - f, [ex.folders{f} '/'], 'folder', f, counts.folders(f));
        end
        rowLabel(3, 'Name');
        perPos = counts.name;
        if work.fromEnd; perPos = counts.nameEnd; end
        for k = 1:K
            p = namePartPos(k, K, work.fromEnd);
            chip(3, k, ex.name{k}, 'name', p, perPos(p));
        end
    end

    function rowLabel(r, txt)
        lb = uilabel(chipGrid, 'Text', txt, 'FontColor', MUTED);
        lb.Layout.Row = r; lb.Layout.Column = 1;
    end

    function chip(r, c, txt, kind, idx, nDistinct)
        on = ismember(idx, work.([kind 'Idx']));
        b = uibutton(chipGrid, 'state', 'Text', txt, 'Value', on, ...
            'ValueChangedFcn', @(src, ~) onChip(kind, idx, src.Value));
        b.Layout.Row = r; b.Layout.Column = c + 1;
        lb = uilabel(chipGrid, 'Text', valuesText(nDistinct), ...
            'FontSize', 10, 'FontColor', MUTED, 'HorizontalAlignment', 'center');
        lb.Layout.Row = r + 1; lb.Layout.Column = c + 1;
    end

    function s = valuesText(nDistinct)
        if nDistinct == nFiles && nFiles > 1
            s = sprintf('%d (all differ)', nDistinct);
        elseif nDistinct == 1
            s = '1 (never changes)';
        else
            s = sprintf('%d values', nDistinct);
        end
    end

    function onChip(kind, idx, on)
        field = [kind 'Idx'];
        if on
            work.(field) = union(work.(field), idx);
        else
            work.(field) = setdiff(work.(field), idx);
        end
        refresh();
    end

    function onExample(v)
        exampleIdx = v;
        buildChips();
        refresh();
    end

    function resetTo(split)
        % Back to the suggestion at this split. Positions mean different
        % parts once the split changes, so a selection cannot carry over.
        [work, counts] = suggestSubjectRule(paths, split, partsFor(split));
        syncOptions();
        buildChips();
        refresh();
    end

    function onFromEnd(v)
        % Keep the same chips selected on the example; only how they are
        % counted changes.
        K = numel(ex.name);
        keep = work.nameIdx(work.nameIdx <= K);
        work.nameIdx = namePartPos(keep, K, true);
        work.fromEnd = logical(v);
        buildChips();
        refresh();
    end

    function syncOptions()
        splitBox.Value = work.splitAlnum;
        endBox.Value   = work.fromEnd;
        % Counting from the end only differs when names differ in length.
        endBox.Enable  = counts.varies || work.fromEnd;
    end

    function refresh()
        [curIds, curMatched, nRespelled] = subjectIdsFromRule(paths, work, ...
            partsFor(work.splitAlnum), labels);
        [u, ~, which] = unique(curIds);
        perSubject = accumarray(which(:), 1);

        % Sorted by subject then path, so the files of one person sit together
        % and a wrong merge or split is visible by eye.
        [~, order] = sortrows([curIds(:), paths(:)]);
        tbl.Data = [short(order)', curIds(order)', num2cell(perSubject(which(order)))];
        removeStyle(tbl);
        bad = find(~curMatched(order));
        if ~isempty(bad)
            addStyle(tbl, uistyle('BackgroundColor', [1 0.9 0.8]), 'row', bad);
        end

        ruleLabel.Text  = ruleText();
        countLabel.Text = countText(numel(u), sum(~curMatched), nRespelled);
    end

    function s = ruleText()
        K = numel(ex.name);
        bits = {};
        LEVEL = {'the file''s folder', 'the folder above', 'two folders up'};
        for f = sort(work.folderIdx(:)', 'descend')
            bits{end + 1} = sprintf('%s ("%s" here)', LEVEL{f}, partOr(ex.folders, f)); %#ok<AGROW>
        end
        for p = sort(work.nameIdx(:)', 'ascend')
            if ~work.fromEnd
                where = sprintf('name part %d', p);
            elseif p == 1
                where = 'the last name part';
            else
                where = sprintf('name part %d from the end', p);
            end
            k = namePartPos(p, K, work.fromEnd);
            bits{end + 1} = sprintf('%s ("%s" here)', where, partOr(ex.name, k)); %#ok<AGROW>
        end
        if isempty(bits)
            s = 'Nothing selected: every file is its own subject.';
        elseif isscalar(bits)
            s = ['Subject = ' bits{1} '.'];
        else
            s = ['Subject = ' strjoin(bits, ' + ') ', joined with _.'];
        end
    end

    function s = countText(nSub, nUnmatched, nRespelled)
        s = sprintf('%d files  ->  %d subject%s', nFiles, nSub, plural(nSub));
        if nSub < nFiles
            s = sprintf('%s   |   n falls from %d to %d', s, nFiles, nSub);
        end
        if nUnmatched > 0
            s = sprintf('%s   |   %d file%s without that part', s, nUnmatched, plural(nUnmatched));
        end
        if nRespelled > 0
            s = sprintf('%s   |   %d spelling%s unified', s, nRespelled, plural(nRespelled));
        end
    end

    function accept()
        accepted = true;
        delete(fig);
    end
end

% ── helpers ─────────────────────────────────────────────────────────────────

function idx = typicalExample(parts)
% The first file with the most common number of name parts, so the chips
% show the cohort's usual shape rather than an odd one out.
len = arrayfun(@(p) numel(p.name), parts);
idx = find(len == mode(len), 1);
end

function v = partOr(c, k)
if k >= 1 && k <= numel(c) && ~isempty(c{k}); v = c{k}; else; v = '?'; end
end
