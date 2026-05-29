
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function ctx = buildBatchContext(filePaths, pipelineName, layout, outputRootOverride)
% BUILDBATCHCONTEXT  Build the shared batch-output struct for one pipeline run.
%   ctx = BUILDBATCHCONTEXT(filePaths, pipelineName, layout)
%   ctx = BUILDBATCHCONTEXT(filePaths, pipelineName, layout, outputRoot)
%
%   Resolves the output root in this priority:
%     1. outputRootOverride argument (when non-empty) - used by tests
%        and programmatic callers that need to bypass user prefs.
%     2. 'outputRoot' preference under the 'nestapp' group.
%     3. commonResultsRoot(filePaths) as a last-resort fallback.
%
%   Names a per-run batch folder yyyyMMdd_HHmmss_<pipelineSlug> and
%   bundles the fields every downstream writer needs:
%     .outputRoot     absolute path to the resolved output root
%     .batchRoot      absolute path to this run's batch folder
%     .batchId        the folder name (without parent dirs)
%     .pipelineSlug   sanitized pipeline name used in batchId
%     .layout         'typeBased' | 'perInput'
%     .startedAt      datetime stamp for this run
%     .nFiles         numel(filePaths)
%
%   Creates the batch folder on disk. Worker processes inherit ctx
%   verbatim (no graphics handles, all serializable).

if nargin < 2, pipelineName = ''; end
if nargin < 3 || isempty(layout), layout = 'typeBased'; end
if nargin < 4, outputRootOverride = ''; end
layout = validateLayout(layout);

ctx.outputRoot   = resolveOutputRoot(filePaths, outputRootOverride);
ctx.pipelineSlug = pipelineSlug(pipelineName);
ctx.startedAt    = datetime('now');
ctx.batchId      = sprintf('%s_%s', ...
    char(ctx.startedAt, 'yyyyMMdd_HHmmss'), ctx.pipelineSlug);
ctx.batchRoot    = fullfile(ctx.outputRoot, ctx.batchId);
ctx.layout       = layout;
ctx.nFiles       = numel(filePaths);

if ~exist(ctx.batchRoot, 'dir')
    mkdir(ctx.batchRoot);
end
end

% -- helpers ---------------------------------------------------------------

function root = resolveOutputRoot(filePaths, override)
% Override beats pref beats common-parent fallback. Empty override
% values are skipped so callers can pass '' to opt into the normal
% pref behaviour.
if ~isempty(override)
    root = override;
    if ~isfolder(root), mkdir(root); end
    return
end
pref = getpref('nestapp', 'outputRoot', '');
if ~isempty(pref) && isfolder(pref)
    root = pref;
else
    root = commonResultsRoot(filePaths);
end
end

function slug = pipelineSlug(name)
% Lowercase + [a-z0-9_]+, runs collapsed, max 40 chars. Empty -> 'pipeline'.
if isempty(name)
    slug = 'pipeline';
    return
end
slug = lower(char(name));
slug = regexprep(slug, '[^a-z0-9]+', '_');
slug = regexprep(slug, '_+', '_');
slug = regexprep(slug, '^_|_$', '');
if isempty(slug)
    slug = 'pipeline';
end
if numel(slug) > 40
    slug = slug(1:40);
end
end

function layout = validateLayout(layout)
allowed = {'typeBased', 'perInput'};
if ~any(strcmp(layout, allowed))
    layout = 'typeBased';
end
end
