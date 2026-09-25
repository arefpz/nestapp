% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [groups, overall] = datasetSummary(entries)
% DATASETSUMMARY  Files and SUBJECTS per group, plus complete-case counts.
%   [groups, overall] = DATASETSUMMARY(entries) describes a grouped dataset.
%
%   groups - 1xG struct array, one per group present, in first-seen order:
%     .name         group name
%     .nFiles       recordings assigned to it
%     .nSubjects    distinct subjects contributing - THE n for any estimate
%     .subjects     their ids, sorted
%
%   overall - struct describing the dataset as a whole:
%     .nFiles           total assigned recordings
%     .nGroups          number of groups
%     .nSubjects        distinct subjects anywhere
%     .completeSubjects subjects present in EVERY group
%     .nComplete        numel(completeSubjects)
%     .nUngrouped       entries with no group yet
%
%   Why nSubjects and nComplete are computed here rather than at each call site:
%   they are the two numbers most easily got wrong, and getting them wrong is
%   invisible. Averaging across FILES when a subject contributed four recordings
%   is pseudo-replication - it shrinks the confidence interval by up to a factor
%   of two on this cohort (8 participants, 32 recordings) while looking entirely
%   normal on screen. And a paired contrast is only defined on subjects present
%   in every group, so a figure claiming n=8 when one subject is missing their
%   post session is simply wrong. Both numbers belong on the figure.
%
%   Entries with an empty group are counted in overall.nUngrouped and excluded
%   from every group - they are not yet part of any comparison.
%
%   See also: assignGroupByFilter, subjectIdsFromRule, groupCurves

groups  = struct('name', {}, 'nFiles', {}, 'nSubjects', {}, 'subjects', {});
overall = struct('nFiles', 0, 'nGroups', 0, 'nSubjects', 0, ...
                 'completeSubjects', {{}}, 'nComplete', 0, 'nUngrouped', 0);

if isempty(entries) || ~isfield(entries, 'group') || ~isfield(entries, 'subject')
    return
end

allGroups   = cellfun(@char, {entries.group},   'UniformOutput', false);
allSubjects = cellfun(@char, {entries.subject}, 'UniformOutput', false);

assigned            = ~cellfun(@isempty, allGroups);
overall.nUngrouped  = sum(~assigned);
overall.nFiles      = sum(assigned);
overall.nSubjects   = numel(unique(allSubjects(assigned)));

names = unique(allGroups(assigned), 'stable');   % first-seen order
if isempty(names)
    % Files present but none grouped yet - a normal state while the user is
    % still defining groups, and groups(0) = ... is an error.
    return
end
groups(numel(names)) = struct('name', '', 'nFiles', 0, 'nSubjects', 0, 'subjects', {{}});
for g = 1:numel(names)
    inGroup = assigned & strcmp(allGroups, names{g});
    subs    = unique(allSubjects(inGroup));
    groups(g).name      = names{g};
    groups(g).nFiles    = sum(inGroup);
    groups(g).nSubjects = numel(subs);
    groups(g).subjects  = subs;
end
overall.nGroups = numel(groups);

% Complete cases: present in every group. With one group that is just its own
% subjects; with none it is empty.
if ~isempty(groups)
    complete = groups(1).subjects;
    for g = 2:numel(groups)
        complete = intersect(complete, groups(g).subjects);
    end
    overall.completeSubjects = complete;
    overall.nComplete        = numel(complete);
end
end
