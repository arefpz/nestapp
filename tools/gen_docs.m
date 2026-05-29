
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function gen_docs(outDir)
% GEN_DOCS  Generate the nestapp API reference from source help blocks.
%
%   Syntax
%     gen_docs              % write Markdown into <repo>/site
%     gen_docs(outDir)      % write into outDir
%
%   Inputs
%     outDir  char - output directory (default <repo>/site). Created if absent.
%
%   Walks src/, extracts each function's H1 line and help block, and emits a
%   browsable Markdown site: an index grouped by module folder, one page per
%   file, and a generated reference of every pipeline step (name + info) from
%   stepRegistry. The function headers are the single source of truth, so the
%   site never drifts from the code (run this in CI; do not commit the output).
%
%   Side effects   Writes Markdown files under outDir.
%   Requires       nestapp src/ on the path (for the steps reference).
%   See also: run_lint, stepRegistry

toolsRoot = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(toolsRoot);
if nargin < 1 || isempty(outDir)
    outDir = fullfile(repoRoot, 'site');
end
srcRoot = fullfile(repoRoot, 'src');
if ~exist(outDir, 'dir'); mkdir(outDir); end
apiDir = fullfile(outDir, 'api');
if ~exist(apiDir, 'dir'); mkdir(apiDir); end

addpath(genpath(srcRoot));

files = dir(fullfile(srcRoot, '**', '*.m'));
[modules, items] = groupByModule(repoRoot, srcRoot, files);

writeIndex(outDir, modules, items);
for i = 1:numel(items)
    writeFilePage(apiDir, items(i));
end
writeStepsReference(outDir, repoRoot);

fprintf('gen_docs: wrote %d API pages + index + steps reference to %s\n', ...
    numel(items), outDir);
end

% ── extraction ──────────────────────────────────────────────────────────────

function [modules, items] = groupByModule(repoRoot, srcRoot, files)
items = struct('name', {}, 'module', {}, 'relpath', {}, 'slug', {}, ...
               'h1', {}, 'help', {});
for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    [h1, helpText] = extractHelp(fpath);
    rel    = strrep(fpath, [repoRoot filesep], '');
    modRel = strrep(files(i).folder, [srcRoot filesep], '');
    if strcmp(files(i).folder, srcRoot); modRel = 'src (top level)'; else; modRel = ['src/' strrep(modRel, filesep, '/')]; end
    [~, base] = fileparts(files(i).name);
    items(end+1) = struct( ...
        'name', files(i).name, 'module', modRel, ...
        'relpath', strrep(rel, filesep, '/'), 'slug', base, ...
        'h1', h1, 'help', helpText); %#ok<AGROW>
end
modules = unique({items.module});
end

function [h1, helpText] = extractHelp(fpath)
% Return the H1 summary and the full help block of a .m file. The help block
% is the first contiguous run of comment lines that is NOT the SPDX license
% header (which sits at the top of every file) - i.e. the function's own
% docstring after its signature. Robust to BOM and to license/no-license files.
raw = fileread(fpath);
if ~isempty(raw) && double(raw(1)) == 65279   % strip a leading UTF-8 BOM
    raw(1) = [];
end
lines = regexp(raw, '\r\n|\n|\r', 'split');
isComment = ~cellfun(@(L) isempty(regexp(L, '^\s*%', 'once')), lines);

% Walk every maximal comment run; return the first that is not a license block.
i = 1; n = numel(lines);
while i <= n
    if ~isComment(i); i = i + 1; continue; end
    runStart = i;
    while i < n && isComment(i + 1); i = i + 1; end
    runStop = i;
    block = regexprep(lines(runStart:runStop), '^\s*%\s?', '');
    if ~any(contains(block, 'SPDX-License-Identifier'))
        h1 = strtrim(block{1});
        helpText = strjoin(block, newline);
        return
    end
    i = i + 1;
end

h1 = '(no description)';
helpText = '_No header comment found._';
end

% ── writers ───────────────────────────────────────────────────────────────

function writeIndex(outDir, modules, items)
L = {'# nestapp API reference', '', ...
     'Auto-generated from source header comments by `tools/gen_docs.m`. ', ...
     'Do not edit by hand — update the function headers instead.', '', ...
     'See [architecture](../docs/architecture.md) for the high-level map and', ...
     'where to make common changes.', ''};
for m = 1:numel(modules)
    L{end+1} = sprintf('## %s', modules{m}); %#ok<AGROW>
    L{end+1} = ''; %#ok<AGROW>
    sel = find(strcmp({items.module}, modules{m}));
    for k = sel
        L{end+1} = sprintf('- [`%s`](api/%s.md) — %s', ...
            items(k).name, items(k).slug, items(k).h1); %#ok<AGROW>
    end
    L{end+1} = ''; %#ok<AGROW>
end
writeFile(fullfile(outDir, 'index.md'), strjoin(L, newline));
end

function writeFilePage(apiDir, item)
L = {sprintf('# `%s`', item.name), '', ...
     sprintf('_Module: %s — source: `%s`_', item.module, item.relpath), '', ...
     '```text', item.help, '```', ''};
writeFile(fullfile(apiDir, [item.slug '.md']), strjoin(L, newline));
end

function writeStepsReference(outDir, repoRoot) %#ok<INUSD>
% Build a table of every pipeline step (name + info) from the registry.
L = {'# Pipeline step reference', '', ...
     'Every step available in the pipeline builder, generated from', ...
     '`stepRegistry()`. To add a step see .github/CONTRIBUTING.md.', '', ...
     '| Step | Description |', '|---|---|'};
try
    steps = stepRegistry();
    for i = 1:numel(steps)
        info = steps(i).info;
        if iscell(info); info = strjoin(info, ' '); end
        info = regexprep(char(info), '\s+', ' ');
        if numel(info) > 200; info = [info(1:197) '...']; end
        L{end+1} = sprintf('| `%s` | %s |', steps(i).name, info); %#ok<AGROW>
    end
catch ME
    L{end+1} = sprintf('_Could not load stepRegistry: %s_', ME.message);
end
writeFile(fullfile(outDir, 'steps.md'), strjoin(L, newline));
end

function writeFile(fpath, text)
fid = fopen(fpath, 'w');
fwrite(fid, text);
fclose(fid);
end
