% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [entries, subjectRule] = exploreFilesTable(entries, opts)
% EXPLOREFILESTABLE  Show and edit which file belongs to whom and to which group.
%   entries = EXPLOREFILESTABLE(entries) opens a modal table of
%   file | subject | group, with subject and group editable, and returns the
%   edited entries - or [] if cancelled.
%
%   This is the backstop the whole subject story rests on. Subject identity
%   cannot be derived reliably from arbitrary naming, and the count it produces
%   sets n for every confidence interval, so the app's position is: never guess
%   silently, always show what was assumed, and make it correctable here.
%
%   Two actions fill the subject column wholesale:
%
%     One subject per file        every recording is its own person, so n is the
%                                 number of files. The default, and correct
%                                 whenever each file IS a different participant.
%     Subjects from filenames     opens subjectRuleDialog, where the user
%                                 clicks which part of the filename is the
%                                 person and sees the result live before it
%                                 applies - on a real cohort this turns 148
%                                 files into 95 subjects, and a change that
%                                 large to n should never happen unseen.
%
%   The count line is the point of the whole dialog: it states files, subjects
%   and groups at all times, so "why is n 92 and not 148" is answerable by
%   looking rather than by reading source.
%
%   opts:
%     .parent       figure to centre over
%     .title        window title
%     .subjectRule  rule the filename dialog starts from (last one used)
%
%   [entries, subjectRule] = EXPLOREFILESTABLE(...) also returns the rule the
%   filename dialog should reopen on next time: the last one applied, else
%   the one passed in.
%
%   See also: exploreDataset, subjectRuleDialog, datasetSummary

if nargin < 2; opts = struct(); end
opts = fillDefaults(opts, struct('parent', [], 'title', 'Files, subjects and groups', ...
                                 'subjectRule', []));
subjectRule = [];

if isempty(entries)
    entries = [];
    return
end

work     = entries;
accepted = false;
lastRule = opts.subjectRule;

W = 760; H = 520;
PAD = 12;

fig = uifigure('Name', char(opts.title), 'Resize', 'off', ...
               'Position', centreOn(opts.parent, W, H), 'WindowStyle', 'modal');
% Exit is by DELETING the figure, never by uiresume, and the wait is waitfor
% rather than uiwait. selectDataTree already documents why: uiresume's
% close-on-X path can leave the window up and soft-lock the app. The concrete
% failure is a nested modal - uiconfirm or uialert on this same figure runs its
% own wait, and afterwards a uiresume no longer releases the outer uiwait, so
% the X silently does nothing and the app is stuck behind a window that will
% not close. waitfor returns the moment the figure is destroyed, and a plain
% delete cannot be vetoed or missed.
%
% Not an onCleanup either: the callbacks are nested functions holding this
% workspace, which would hold the onCleanup, which would hold the figure - a
% cycle that leaks a window on every call.
fig.CloseRequestFcn = @(src, ~) delete(src);

try
tbl = uitable(fig, 'Position', [PAD, 96, W - 2*PAD, H - 96 - 56], ...
    'ColumnName', {'File', 'Subject', 'Group'}, ...
    'ColumnWidth', {430, 150, 130}, ...
    'ColumnEditable', [false true true], 'RowName', {});
tbl.CellEditCallback = @(~, ev) onEdit(ev);

countLabel = uilabel(fig, 'Position', [PAD, H - 34, W - 2*PAD, 22], ...
    'FontWeight', 'bold');

uibutton(fig, 'Text', 'One subject per file', 'Position', [PAD, 62, 170, 26], ...
    'ButtonPushedFcn', @(~, ~) subjectsPerFile());
uibutton(fig, 'Text', 'Subjects from filenames...', 'Position', [PAD + 178, 62, 190, 26], ...
    'ButtonPushedFcn', @(~, ~) subjectsFromFilenames());

uilabel(fig, 'Position', [PAD, 34, W - 2*PAD, 22], 'FontSize', 11, ...
    'FontColor', [0.35 0.38 0.43], ...
    'Text', ['Subject identity only matters for two things: collapsing repeat ' ...
             'recordings of one person, and pairing groups. Otherwise leave it.']);

uibutton(fig, 'Text', 'Use these assignments', 'Position', [W - PAD - 320, PAD, 155, 26], ...
    'ButtonPushedFcn', @(~, ~) accept());
uibutton(fig, 'Text', 'Cancel', 'Position', [W - PAD - 155, PAD, 155, 26], ...
    'ButtonPushedFcn', @(~, ~) delete(fig));

refresh();
setappdata(fig, 'nestappModalReady', true);
waitfor(fig);
catch ME
    if isvalid(fig); delete(fig); end
    rethrow(ME);
end

% `accepted` and `work` live in this workspace, which outlives the figure, so
% the answer survives the deletion that released waitfor.
if accepted
    entries     = work;
    subjectRule = lastRule;
else
    entries = [];
end
if isvalid(fig); delete(fig); end

% ── nested ───────────────────────────────────────────────────────────────
    function refresh()
        n = numel(work);
        data = cell(n, 3);
        for i = 1:n
            data{i, 1} = shortPath(work(i).path);
            data{i, 2} = work(i).subject;
            data{i, 3} = work(i).group;
        end
        tbl.Data = data;
        countLabel.Text = summaryText(work);
    end

    function onEdit(ev)
        r = ev.Indices(1);
        c = ev.Indices(2);
        if r < 1 || r > numel(work); return; end
        v = strtrim(char(string(ev.NewData)));
        switch c
            case 2
                if isempty(v)
                    % An empty subject would silently merge rows; refuse it.
                    refresh();
                    return
                end
                work(r).subject = v;
                work(r).subjectConfident = true;   % given, not guessed
            case 3
                work(r).group = v;
        end
        refresh();
    end

    function subjectsPerFile()
        ids = uniqueFileLabels({work.path});
        for i = 1:numel(work)
            work(i).subject          = ids{i};
            work(i).subjectConfident = true;
        end
        refresh();
    end

    function subjectsFromFilenames()
        % The dialog is its own preview: it shows the file -> subject table and
        % the change in n on every click, so nothing here asks again.
        [ids, rule] = subjectRuleDialog({work.path}, ...
            struct('parent', fig, 'rule', lastRule));
        if isempty(ids); return; end
        for i = 1:numel(work)
            work(i).subject          = ids{i};
            work(i).subjectConfident = true;   % the user chose the rule
        end
        lastRule = rule;
        refresh();
    end

    function accept()
        accepted = true;
        delete(fig);
    end
end

% ── helpers ─────────────────────────────────────────────────────────────────

function s = summaryText(work)
[~, overall] = datasetSummary(work);
nSub = numel(unique({work.subject}));
s = sprintf('%d files   |   %d subjects   |   %d group(s)', ...
            numel(work), nSub, overall.nGroups);
if overall.nUngrouped > 0
    s = sprintf('%s   |   %d not in a group', s, overall.nUngrouped);
end
if nSub < numel(work)
    s = sprintf('%s   |   %d files share a subject', s, numel(work) - nSub);
end
end
