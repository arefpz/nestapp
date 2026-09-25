% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef SubjectRuleTest < NestappTestCase
% SUBJECTRULETEST  Picking the subject out of a filename, by position.
%
%   The subject dialog is only a face on three pure functions: filenameParts
%   cuts a name into parts, subjectIdsFromRule takes the chosen positions, and
%   suggestSubjectRule picks a starting position from the whole cohort. What
%   matters is that the rule is deterministic, that it never merges a file it
%   cannot read into someone else, and that the suggestion lands on the person
%   across the naming schemes people actually come up with.
%
%   SCHEMES is a table of those schemes, each generated with a KNOWN subject
%   for 4 people x 2 recordings, from plain to deliberately awkward. Each row
%   asserts the suggestion recovers the true grouping. Comparison is by
%   grouping, not spelling: which files share a subject is what sets n.

    properties (TestParameter)
        scheme = struct( ...
            'plain',             {{'S%s_%s.set',                       'id', 'cond'}}, ...
            'sessionAndSite',    {{'rtmsct0%s_1_%s_SPL_tesa.set',      'id', 'cond'}}, ...
            'bidsNested',        {{'sub-%s/ses-%d/eeg/sub-%s_ses-%d_eeg.set', 'bids'}}, ...
            'conditionFirst',    {{'%s_S%s.set',                       'cond', 'id'}}, ...
            'noSeparators',      {{'S%s%s.set',                        'id', 'cond'}}, ...
            'digitsThenWord',    {{'0%s%s.set',                        'id', 'cond'}}, ...
            'spacesAndParens',   {{'Patient %s (visit %d).set',        'id', 'visit'}}, ...
            'bracketsAndComma',  {{'[S%s],%s.set',                     'id', 'cond'}}, ...
            'dateFirst',         {{'2024031%d_S%s.set',                'day', 'id'}}, ...
            'dottedDateFirst',   {{'2024.03.1%d_P%s_TMS.set',          'day', 'id'}}, ...
            'sessionFirst',      {{'ses%d_S%s.set',                    'visit', 'id'}}, ...
            'initials',          {{'%s_%s.set',                        'initials', 'cond'}}, ...
            'folderPerSubject',  {{'S%s/%s.set',                       'id', 'cond'}}, ...
            'subjectTwoUp',      {{'S%s/session%d/data.set',           'id', 'visit'}}, ...
            'variablePrefix',    {{'%sS%s.set',                        'prefix', 'id'}}, ...
            'mixedCase',         {{'%s%s_%s.set',                      'case', 'id', 'cond'}}, ...
            'zeroPadding',       {{'S%s_%s.set',                       'pad', 'cond'}}, ...
            'accentedNames',     {{'Proband_%s_%s.set',                'names', 'cond'}});
    end

    methods (Test)
        function suggestionFindsTheSubject(tc, scheme)
            [paths, truth] = generate(scheme);
            rule = suggestSubjectRule(paths);
            [ids, matched] = subjectIdsFromRule(paths, rule);
            tc.verifyTrue(all(matched), 'every file should have the suggested part');
            tc.verifyTrue(samePartition(ids, truth), ...
                sprintf('grouped as %s', strjoin(ids, ', ')));
        end

        function countsReadTheNamingScheme(tc)
            paths = {'D:/d/rtmsct001_1_pre_SPL_tesa.set', 'D:/d/rtmsct001_2_post_SPL_tesa.set', ...
                     'D:/d/rtmsct003_1_pre_SPL_tesa.set', 'D:/d/rtmsct003_2_post_SPL_tesa.set'};
            [~, counts] = suggestSubjectRule(paths);
            tc.verifyEqual(counts.name, [2 2 2 1 1]);
            tc.verifyFalse(counts.varies);
        end

        function folderTieGoesToTheOuterFolder(tc)
            % 2 people x 2 sessions: both folder levels take 2 values, and the
            % person's folder is the one that holds the sessions.
            paths = {'D:/c/S01/ses1/eeg/data.set', 'D:/c/S01/ses2/eeg/data.set', ...
                     'D:/c/S02/ses1/eeg/data.set', 'D:/c/S02/ses2/eeg/data.set'};
            rule = suggestSubjectRule(paths);
            tc.verifyEqual(rule.folderIdx, 3);
            tc.verifyEqual(subjectIdsFromRule(paths, rule), {'S01', 'S01', 'S02', 'S02'});
        end

        function driveLetterIsNotAFolder(tc)
            parts = filenameParts('D:\rec.set');
            tc.verifyEqual(parts.folders, {'', '', ''});
        end

        function datesAndTimesAreOnePart(tc)
            parts = filenameParts('x/S01_2024-03-15_14-31-05.set');
            tc.verifyEqual(parts.name, {'S01', '2024-03-15', '14-31-05'});
            tc.verifyEqual(parts.isDate, [false true true]);
        end

        function everyNonAlphanumericSeparates(tc)
            parts = filenameParts('x/[S01] pre,(v2)+final.set');
            tc.verifyEqual(parts.name, {'S01', 'pre', 'v2', 'final'});
        end

        function foldersCountUpwardFromTheFile(tc)
            parts = filenameParts('C:/a/b/c/rec.set');
            tc.verifyEqual(parts.folders, {'c', 'b', 'a'});
        end

        function splitLettersFromNumbers(tc)
            parts = filenameParts('x/S01pre.set', true);
            tc.verifyEqual(parts.name, {'S', '01', 'pre'});
        end

        function combinedPartsJoinInPathOrder(tc)
            rule = struct('nameIdx', [3 1], 'folderIdx', 1);
            ids  = subjectIdsFromRule({'siteA/x_pre_07.set'}, rule);
            tc.verifyEqual(ids, {'siteA_x_07'});
        end

        function countingFromTheEnd(tc)
            paths = {'tms_S01.set', 'resting_state_S01.set', 'tms_S02.set'};
            ids = subjectIdsFromRule(paths, struct('nameIdx', 1, 'fromEnd', true));
            tc.verifyEqual(ids, {'S01', 'S01', 'S02'});
        end

        function spellingVariantsAreOnePerson(tc)
            [ids, ~, nRespelled] = subjectIdsFromRule( ...
                {'s01_pre.set', 'S1_post.set', 'S02_pre.set'}, struct('nameIdx', 1));
            tc.verifyEqual(ids{1}, ids{2});
            tc.verifyNotEqual(ids{1}, ids{3});
            tc.verifyEqual(nRespelled, 1);
        end

        function fileWithoutThePartStaysItsOwnSubject(tc)
            paths = {'a/P1_pre.set', 'a/P1_post.set', 'a/odd.set'};
            [ids, matched] = subjectIdsFromRule(paths, struct('nameIdx', 2));
            tc.verifyEqual(matched, [true true false]);
            tc.verifyEqual(ids{3}, 'odd');
            tc.verifyNotEqual(ids{3}, ids{1});
        end

        function emptyRuleIsOneSubjectPerFile(tc)
            paths = {'a/x.set', 'b/x.set'};
            [ids, matched] = subjectIdsFromRule(paths, struct());
            tc.verifyEqual(ids, uniqueFileLabels(paths));
            tc.verifyTrue(all(matched));
        end

        function nothingPlausibleSuggestsNothing(tc)
            % Every part is either constant or unique per file.
            rule = suggestSubjectRule({'study_A1.set', 'study_B2.set', 'study_C3.set'});
            tc.verifyEmpty(rule.nameIdx);
            tc.verifyEmpty(rule.folderIdx);
        end

        function exploreDatasetReproducesTheRule(tc)
            paths = {'P1_pre.set', 'P1_post.set', 'P2_pre.set'};
            rule  = struct('nameIdx', 1);
            entries = exploreDataset(paths, {}, struct('subjectRule', rule));
            tc.verifyEqual({entries.subject}, subjectIdsFromRule(paths, rule));
            guessed = exploreDataset(paths, {}, struct('subjectMode', 'guess'));
            tc.verifyEqual({guessed.subject}, {'P1', 'P1', 'P2'});
        end
    end
end

% ── helpers ─────────────────────────────────────────────────────────────────

function [paths, truth] = generate(scheme)
% 4 people x 2 recordings under D:/cohort, from a format and its fields.
fmt  = scheme{1};
IDS  = {'01', '02', '03', '04'};
COND = {'pre', 'post'};
paths = {}; truth = {};
for s = 1:4
    for c = 1:2
        args = {};
        for f = scheme(2:end)
            switch f{1}
                case 'id';       args{end + 1} = IDS{s};                              %#ok<AGROW>
                case 'cond';     args{end + 1} = COND{c};                             %#ok<AGROW>
                case 'visit';    args{end + 1} = c;                                   %#ok<AGROW>
                case 'day';      args{end + 1} = c + s;                               %#ok<AGROW>
                case 'bids';     args = {IDS{s}, c, IDS{s}, c};
                case 'initials'; args{end + 1} = pick({'JD', 'AB', 'KL', 'MN'}, s);   %#ok<AGROW>
                case 'names';    args{end + 1} = pick({'Müller', 'Øberg', 'Núñez', 'Zoë'}, s); %#ok<AGROW>
                case 'prefix';   args{end + 1} = pick({'tms_', 'resting_state_'}, c); %#ok<AGROW>
                case 'case';     args{end + 1} = pick({'s', 'S'}, c);                 %#ok<AGROW>
                case 'pad';      args{end + 1} = pick({sprintf('%d', s), IDS{s}}, c); %#ok<AGROW>
            end
        end
        paths{end + 1} = ['D:/cohort/' sprintf(fmt, args{:})]; %#ok<AGROW>
        truth{end + 1} = sprintf('person%d', s);               %#ok<AGROW>
    end
end
end

function tf = samePartition(a, b)
[~, ~, ia] = unique(a);
[~, ~, ib] = unique(b);
nPairs = size(unique([ia(:) ib(:)], 'rows'), 1);
tf = nPairs == max(ia) && nPairs == max(ib);
end

function v = pick(c, i)
v = c{i};
end
