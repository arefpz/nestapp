% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function res = groupCurves(cache, entries, opts)
% GROUPCURVES  Per-group subject curves and intervals, from cached trial averages.
%   res = GROUPCURVES(cache, entries, opts) turns the reduced datasets from
%   loadReducedSets into everything a comparison view needs. It touches no
%   files, so it is cheap enough to re-run on every ROI, mode or window change.
%
%   Inputs:
%     cache   - struct array from loadReducedSets, keyed by .path.
%     entries - struct array with .path .subject .group (see assignGroupByFilter).
%     opts    - .roi        cellstr of ROI electrode labels (required for TEP/LMFP)
%               .mode       'TEP' | 'GMFP' | 'LMFP'   (default 'TEP')
%               .design     'paired' | 'unpaired'     (default 'unpaired')
%               .smoothWin  moving-average samples    (default 5, as the app uses)
%               .level      confidence level          (default 0.95)
%               .capRule    'exclude' | 'shared'      (default 'exclude')
%                           what to do when files were recorded on different
%                           caps - see "Mixed caps" below
%
%   Output res:
%     .time           1xT common time base (ms)
%     .channelLabels  electrodes common to EVERY file in EVERY group
%     .chanlocs       locations for those electrodes, for scalp maps
%     .groups         1xG struct: .name .subjects .curves (nSubj x T)
%                     .chanMeans (C x T, subject-averaged) .nFiles .nSubjects
%                     .files - 1xnFiles struct (.name .subject .curve), the
%                     individual recordings BEFORE the subject collapse
%     .est            1xG from curveInterval, aligned with .groups
%     .contrast       differenceInterval for groups(2) minus groups(1) when
%                     there are exactly two groups; struct([]) otherwise
%     .design         the design actually used
%     .complete       subjects present in every group
%     .dropped        subjects excluded from a paired estimate, and why
%     .info           time-base info from commonTimeBase, plus .roi and .level
%                     (.requested .matched .missing - a partial ROI is
%                     averaged, not refused, so the difference is recorded)
%                     and .montage:
%                       .uniform    true when every file has the same channels
%                       .rule       the capRule applied
%                       .nChannels  size of the montage the curves use
%                       .excluded   files left out, by name ('exclude' only)
%                       .excludedByGroup  struct array .group .n - how the
%                                   exclusions fall across groups
%                       .dropped    electrodes of the modal cap left out
%                                   because some file lacks them ('shared' only)
%                       .variants   each distinct montage with its files and
%                                   what it lacks relative to the modal one
%
%   Three decisions here are the ones that were previously wrong or absent:
%
%   Subjects, then groups. Files are averaged within a subject first - weighted
%   by trial count, so a 40-trial and a 120-trial session do not count equally -
%   and only then across subjects. The old code averaged files directly and
%   divided by sqrt(nFiles), which on this cohort (8 people, 32 recordings)
%   treats one person's four sessions as four independent observations.
%
%   Mixed caps: the analyst's call, and whichever rule applies is recorded.
%   'exclude' keeps the montage most files carry and leaves the other files
%   out, named in res.info.montage.excluded. 'shared' keeps every file and
%   uses only the electrodes all of them carry, naming the electrodes lost in
%   res.info.montage.dropped. Neither is neutral - see chooseMontage.
%
%   Complete cases for paired designs. A paired interval is only defined on
%   subjects present in every group; the rest are dropped and NAMED in
%   res.dropped so the figure can say "n = 7 of 8 complete pairs" instead of
%   quietly reporting 8.
%
%   See also: loadReducedSets, curveInterval, commonTimeBase, tepFieldCurve

if nargin < 3; opts = struct(); end
opts = fillDefaults(opts, struct('roi', {{}}, 'mode', 'TEP', ...
    'design', 'unpaired', 'smoothWin', 5, 'level', 0.95, 'capRule', 'exclude'));
opts.capRule = lower(char(string(opts.capRule)));
if ~any(strcmp(opts.capRule, {'exclude', 'shared'}))
    error('nestapp:badCapRule', ...
        'capRule must be ''exclude'' or ''shared''; got "%s".', opts.capRule);
end

res = struct('time', [], 'channelLabels', {{}}, 'chanlocs', [], ...
             'groups', struct('name', {}, 'subjects', {}, 'curves', {}, ...
                              'chanMeans', {}, 'nFiles', {}, 'nSubjects', {}, ...
                              'files', {}), ...
             'est', struct([]), 'contrast', struct([]), ...
             'design', opts.design, 'complete', {{}}, ...
             'dropped', {{}}, 'info', struct());

use = groupedCacheEntries(cache, entries);
if isempty(use); return; end

% ── one montage, then one time base for the files that carry it ──────────
% Montage first: under 'exclude' files on a different cap are left out
% outright, so they must not get a say in the time base either.
[res.channelLabels, conforms, chanIdx, montage] = chooseMontage(cache, use, ...
                                                                opts.capRule);
use     = use(conforms);
chanIdx = chanIdx(conforms);
if isempty(res.channelLabels)
    error('nestapp:noCommonChannels', ...
        'The selected files have no channel labels, so they cannot be compared.');
end

[res.time, timeIdx, res.info] = commonTimeBase({cache([use.ci]).time}, {use.label});
res.info.montage = montage;
% The level the estimates were computed at, recorded so nothing downstream has
% to restate it. A status line or figure footer that hardcoded 95 would keep
% saying 95 after this default changed, or after a draw option overrode it -
% and a wrong coverage on an exported figure is the one label that does real
% damage, because the figure outlives the session that made it.
res.info.level   = opts.level;
res.chanlocs = cache(use(1).ci).chanlocs(chanIdx{1});

% ── per file: crop to the common montage and time base, then reduce ──────
roiIdx = roiChannelIndex(res.channelLabels, opts.roi);
if isempty(roiIdx) && ~strcmpi(opts.mode, 'GMFP')
    error('nestapp:emptyROI', ...
        ['The ROI matches none of the %d electrodes common to these files. ' ...
         'Pick electrodes that exist in every selected file.'], ...
        numel(res.channelLabels));
end

% What was ASKED FOR versus what was AVERAGED. Only a total miss errored
% above; a partial match - three of five electrodes present - was averaged
% silently while the status line and the exported figure's footer went on
% naming all five. A caption that overstates the ROI is the one place this
% does real damage, so the difference is recorded and the footer stamps
% .matched rather than the request.
res.info.roi = struct( ...
    'requested', {opts.roi(:)'}, ...
    'matched',   {res.channelLabels(roiIdx)}, ...
    'missing',   {opts.roi(~ismember(lower(opts.roi), lower(res.channelLabels)))});

for i = 1:numel(use)
    avg            = cache(use(i).ci).trialAvg(chanIdx{i}, timeIdx{i});
    use(i).chanAvg = avg;
    % tepFieldCurve takes channels x time x trials; a 2-D trial average is the
    % degenerate case its mean(...,3) leaves untouched, so the one definition of
    % each mode is reused rather than restated here.
    use(i).curve  = smoothCurve(tepFieldCurve(avg, roiIdx, opts.mode), opts.smoothWin);
end

% ── collapse to subjects, then to groups ─────────────────────────────────
groupNames = unique({use.group}, 'stable');
G          = numel(groupNames);
groups     = res.groups;
for g = 1:G
    inG   = strcmp({use.group}, groupNames{g});
    block = use(inG);
    subs  = unique({block.subject}, 'stable');

    curves    = zeros(numel(subs), numel(res.time));
    chanMeans = zeros(numel(res.channelLabels), numel(res.time));
    for s = 1:numel(subs)
        rows            = strcmp({block.subject}, subs{s});
        [curves(s, :), sChan] = collapseSubject(block(rows));
        chanMeans       = chanMeans + sChan;
    end

    groups(g).name      = groupNames{g};
    groups(g).subjects  = subs;
    groups(g).curves    = curves;
    groups(g).chanMeans = chanMeans / max(numel(subs), 1);
    groups(g).nFiles    = numel(block);
    groups(g).nSubjects = numel(subs);

    % The per-FILE records, kept rather than discarded. They were already
    % computed above and thrown away with `use`, so retaining them costs one
    % projection and nothing in the reduction. They answer "is that one
    % recording an outlier" - a question a group mean cannot be asked, and
    % which previously needed the old Add-to-current-Figure gesture.
    %
    % A struct array, not parallel arrays: it is the shape `cache`, `entries`
    % and `block` itself already use, and it is what keeps the paired-design
    % filter below a single indexing operation however many per-file fields
    % this grows. .name carries the identity a curve alone does not have, so a
    % saved result can say WHICH recording an outlying row was.
    groups(g).files = struct('name',    {block.label}, ...
                             'subject', {block.subject}, ...
                             'curve',   {block.curve});
end

% ── intervals ────────────────────────────────────────────────────────────
res.complete = completeSubjects(groups);
if strcmpi(opts.design, 'paired')
    [groups, res.dropped] = restrictToComplete(groups, res.complete);
end
res.groups = groups;
res.est    = curveInterval({groups.curves}, opts.design, opts.level);

% The contrast is an estimate, so it belongs here with the others rather than
% at draw time. drawDifferenceWave was deriving it with no level argument,
% which drew a 95% difference band beside 99% group bands whenever the run's
% level was changed.
if numel(groups) == 2
    res.contrast = differenceInterval(groups(1).curves, groups(2).curves, ...
                                      opts.design, opts.level);
end
end

% ── helpers ─────────────────────────────────────────────────────────────────

function [labels, conforms, idx, montage] = chooseMontage(cache, use, rule)
% The channel set the curves are computed on, and which files carry it.
%
% On the development cohort most files carry AFz/FCz/FT9/FT10 and a few carry
% CPz/FPz/Iz instead - two cap configurations, not bad channels awaiting
% interpolation. There is no neutral way to combine them, so the rule is the
% caller's:
%
%   'exclude'  The montage most files carry wins and the rest are left out and
%              NAMED. Every included file is exact and GMFP stays comparable to
%              a full-cap analysis - but there are fewer files, and if the cap
%              changed with site, period or condition the exclusions are not
%              random and can unbalance the groups.
%   'shared'   Every file stays, on the electrodes all of them carry. ROI
%              curves are unchanged when the ROI is on every cap, but maps lose
%              those electrodes and GMFP - a standard deviation ACROSS channels
%              - is over fewer of them than a full-cap analysis.
%
% Either way what was given up is recorded, so the status line and the figure
% footer can say it.
%
% Ties (two montages with equally many files) go to the one seen first, which
% is stable for a given file order and reported either way.
labelSets = {cache([use.ci]).labels};
caps      = capSummary(labelSets);

if strcmp(rule, 'shared')
    labels   = caps.sharedLabels;
    conforms = true(1, numel(use));
    dropped  = caps.dropped;
else
    labels   = caps.modalLabels;
    conforms = caps.modal;
    dropped  = {};
end

idx = cell(1, numel(use));
lowLabels = lower(labels);
for i = find(conforms)
    [~, idx{i}] = ismember(lowLabels, lower(labelSets{i}));
end

montage = struct();
montage.uniform   = (caps.nVariants == 1);
montage.rule      = rule;
montage.nChannels = numel(labels);
montage.excluded  = {use(~conforms).label};
montage.dropped   = dropped;

% How the exclusions fall across groups. Lopsided exclusions are how 'exclude'
% biases a comparison, so they are counted for the app to flag.
groupOf = {use.group};
gNames  = unique(groupOf, 'stable');
montage.excludedByGroup = struct('group', gNames, 'n', 0);
for g = 1:numel(gNames)
    montage.excludedByGroup(g).n = sum(~conforms & strcmp(groupOf, gNames{g}));
end

montage.variants = struct('nChannels', {}, 'nFiles', {}, 'files', {}, 'missing', {});
lowModal = lower(caps.modalLabels);
for v = 1:caps.nVariants
    inV = (caps.which == v);
    vLabels = labelSets{find(inV, 1)};
    montage.variants(v).nChannels = numel(vLabels);
    montage.variants(v).nFiles    = sum(inV);
    montage.variants(v).files     = {use(inV).label};
    % What this variant lacks relative to the modal one - the actionable bit.
    montage.variants(v).missing   = caps.modalLabels(~ismember(lowModal, lower(vLabels)));
end
end

function [curve, chanAvg] = collapseSubject(files)
% One subject's recordings, weighted by trial count: a 120-trial session is a
% better estimate of that person than a 40-trial one and should count more.
w       = [files.nTrials];
w       = w(:) / sum(w);
curve   = sum(cell2mat({files.curve}') .* w, 1);
chanAvg = zeros(size(files(1).chanAvg));
for i = 1:numel(files)
    chanAvg = chanAvg + files(i).chanAvg * w(i);
end
end

function subs = completeSubjects(groups)
if isempty(groups); subs = {}; return; end
subs = groups(1).subjects;
for g = 2:numel(groups)
    subs = intersect(subs, groups(g).subjects);
end
end

function [groups, dropped] = restrictToComplete(groups, complete)
% Keep only complete cases, in the SAME subject order in every group - the
% paired interval pairs by row position, so a mismatch here would silently
% compare different people.
dropped = {};
for g = 1:numel(groups)
    [tf, loc] = ismember(complete, groups(g).subjects);
    dropped   = [dropped, setdiff(groups(g).subjects, complete)]; %#ok<AGROW>
    keep      = loc(tf);
    groups(g).subjects  = groups(g).subjects(keep);
    groups(g).curves    = groups(g).curves(keep, :);
    groups(g).nSubjects = numel(keep);

    % The per-file records have to follow, or a traces overlay would draw
    % recordings belonging to subjects this estimate just excluded - the
    % picture would disagree with the interval drawn on top of it. One
    % indexing operation, so a future per-file field cannot fall out of step.
    keepFile        = ismember({groups(g).files.subject}, groups(g).subjects);
    groups(g).files = groups(g).files(keepFile);
    groups(g).nFiles = numel(groups(g).files);
end
dropped = unique(dropped);
end

function y = smoothCurve(x, win)
if isempty(win) || win <= 1; y = x; return; end
y = smoothdata(x, 'movmean', win);
end
