function changed = apply_headers(mode)
% APPLY_HEADERS  Add a short GPL-3.0 license header to source files.
%
%   Syntax
%     apply_headers              % DRY RUN: list files that lack a header
%     apply_headers('apply')     % insert the header into files that lack one
%     changed = apply_headers(...)% return the list of affected files
%
%   Inputs
%     mode  char - 'dryrun' (default) or 'apply'.
%
%   Outputs
%     changed  cellstr - repo-relative paths of files that were (or would be)
%              modified.
%
%   The header is inserted after the function's H1/help comment block so the
%   help text and the generated API docs are unaffected. Idempotent: files
%   that already contain the SPDX tag are skipped. GPL-3.0 does not strictly
%   require per-file headers (the repo LICENSE governs the whole project), so
%   this is a convenience for consistency, not a compliance gate.
%
%   See also: run_lint

if nargin < 1 || isempty(mode); mode = 'dryrun'; end
doApply = strcmpi(mode, 'apply');

toolsRoot = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(toolsRoot);

header = strjoin({
    '% SPDX-License-Identifier: GPL-3.0-or-later'
    '% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.'
    '% Part of nestapp; see the LICENSE file for full terms.'
    ''}, newline);

roots   = {fullfile(repoRoot,'src'), fullfile(repoRoot,'tests'), fullfile(repoRoot,'tools')};
changed = {};
for r = 1:numel(roots)
    files = dir(fullfile(roots{r}, '**', '*.m'));
    for i = 1:numel(files)
        fpath = fullfile(files(i).folder, files(i).name);
        txt   = fileread(fpath);
        if contains(txt, 'SPDX-License-Identifier')
            continue   % already has a header
        end
        changed{end+1} = strrep(fpath, [repoRoot filesep], ''); %#ok<AGROW>
        if doApply
            insertHeaderAfterHelp(fpath, txt, header);
        end
    end
end

verb = 'need';
if doApply; verb = 'updated'; end
fprintf('apply_headers (%s): %d file(s) %s a header.\n', mode, numel(changed), verb);
for i = 1:numel(changed)
    fprintf('  %s\n', changed{i});
end
if ~doApply && ~isempty(changed)
    fprintf('Run apply_headers(''apply'') to insert headers.\n');
end
end

% ── helpers ───────────────────────────────────────────────────────────────

function insertHeaderAfterHelp(fpath, txt, header)
% Insert the header block immediately after the leading function line + its
% contiguous comment (help) block, so help text stays first for `help`/`doc`.
if ~isempty(txt) && double(txt(1)) == 65279   % ignore a leading UTF-8 BOM
    txt(1) = [];
end
lines = regexp(txt, '\r\n|\n|\r', 'split');
n = numel(lines);
% Find the function signature line (first non-empty, non-comment line).
sigIdx = find(~cellfun(@isempty, regexp(lines, '^\s*function\b', 'once')), 1);
if isempty(sigIdx); sigIdx = 0; end
% Help block = contiguous comment lines immediately after the signature.
j = sigIdx + 1;
while j <= n && ~isempty(regexp(lines{j}, '^\s*%', 'once'))
    j = j + 1;
end
insertAt = j - 1;   % after the last help-comment line
out = [lines(1:insertAt), {''}, strsplit(strip(header, 'right', newline), newline), lines(insertAt+1:end)];
fid = fopen(fpath, 'w');
fwrite(fid, strjoin(out, newline));
fclose(fid);
end
