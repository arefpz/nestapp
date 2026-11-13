% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function result = installAaratep(varargin)
% INSTALLAARATEP  Download and install the AARATEP helper functions.
%   result = INSTALLAARATEP()
%   result = INSTALLAARATEP('Progress', fcn, 'Force', tf, 'TargetDir', d)
%
%   Fetches the pinned AARATEP release as a zip, extracts it, and puts it where
%   ensureAaratepOnPath will find it. Replaces the documented manual route -
%   open a terminal, cd third_party, git clone - which needs git installed and
%   a shell, neither of which the researchers this app is for necessarily have.
%
%   Name-value options:
%     'Progress'   fcn(fraction, message) called as work proceeds, so a dialog
%                  can show it. Default: no output.
%     'Force'      re-download even if a tree is already installed (false).
%     'TargetDir'  install here instead of the first writable default.
%
%   result fields:
%     .installed   true if a usable tree is in place when this returns
%     .dir         where it is
%     .action      'already-present' | 'installed' | 'failed'
%     .message     one line suitable for a dialog
%     .newerTag    a newer upstream release, if one exists ('' otherwise)
%
%   NO GIT, DELIBERATELY. websave + unzip are built into MATLAB, so this works
%   on a machine with nothing but MATLAB installed. Measured: 0.6 MB, about a
%   second, 324 files.
%
%   EXTRACTED TO A TEMP DIR AND MOVED IN ONE STEP. A download that fails
%   halfway must not leave a half-populated aaratep/ folder behind, because
%   ensureAaratepOnPath decides by isfolder and would then put a broken tree on
%   the path and fail deep inside a pipeline run. So nothing appears at the
%   destination until the extracted copy has been verified to contain the
%   sentinel function.
%
%   See also: aaratepRelease, ensureAaratepOnPath, nestappDoctor

p = inputParser;
p.addParameter('Progress',  @(~,~) [], @(f) isa(f, 'function_handle'));
p.addParameter('Force',     false, @(x) islogical(x) && isscalar(x));
p.addParameter('TargetDir', '',    @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
report    = p.Results.Progress;
force     = p.Results.Force;
targetDir = char(p.Results.TargetDir);

rel    = aaratepRelease();
result = struct('installed', false, 'dir', '', 'action', 'failed', ...
                'message', '', 'newerTag', '');

% ── already there? ────────────────────────────────────────────────────────────
% Scoped to TargetDir when the caller named one - aaratepStatus owns that
% distinction, and conflating it was a real bug (see its header).
existing = aaratepStatus(targetDir);
if existing.installed && ~force
    result.installed = true;
    result.dir       = existing.dir;
    result.action    = 'already-present';
    result.message   = sprintf('%s is already installed at %s', ...
                               existing.label, existing.dir);
    result.newerTag  = newerReleaseOrEmpty(rel);
    return
end

% ── where to put it ───────────────────────────────────────────────────────────
if isempty(targetDir)
    targetDir = firstWritable(rel.writeDirs);
end
if isempty(targetDir)
    result.message = sprintf( ...
        ['Could not find a writable place to install AARATEP. Tried:\n  %s\n' ...
         'Set a folder with:  setpref(''nestapp'', ''aaratepPath'', yourFolder)'], ...
        strjoin(rel.writeDirs, sprintf('\n  ')));
    return
end

% ── fetch, verify, then move ──────────────────────────────────────────────────
staging = tempname;
cleanup = onCleanup(@() removeQuietly(staging));
mkdir(staging);

try
    report(0.05, sprintf('Downloading AARATEP %s...', rel.tag));
    zipFile = fullfile(staging, 'aaratep.zip');
    websave(zipFile, rel.zipUrl);

    report(0.55, 'Extracting...');
    unzip(zipFile, staging);

    % GitHub archives extract to <name>-<version>/, so the real tree is one
    % level down and its name is not something to hardcode.
    extracted = soleSubfolder(staging);

    report(0.80, 'Verifying...');
    if ~isfile(fullfile(extracted, [rel.sentinel '.m']))
        error('nestapp:aaratepDownloadIncomplete', ...
            ['The downloaded archive does not contain %s.m, so it is not the ' ...
             'AARATEP pipeline. Nothing has been installed.'], rel.sentinel);
    end

    report(0.90, 'Installing...');
    if isfolder(targetDir)
        removeQuietly(targetDir);   % a Force re-install replaces cleanly
    end
    parent = fileparts(targetDir);
    if ~isfolder(parent); mkdir(parent); end
    [ok, msg] = movefile(extracted, targetDir);
    if ~ok
        error('nestapp:aaratepInstallFailed', ...
            'Could not move the extracted files to %s: %s', targetDir, msg);
    end

    stampVersion(targetDir, rel.tag);

    % The path memo may hold a "missing" answer from an earlier attempt.
    ensureAaratepOnPath('reset');

    report(1.0, 'Done.');
    result.installed = true;
    result.dir       = targetDir;
    result.action    = 'installed';
    result.message   = sprintf('AARATEP %s installed to %s', rel.tag, targetDir);
    result.newerTag  = newerReleaseOrEmpty(rel);

catch err
    result.message = sprintf('AARATEP install failed: %s', err.message);
end
end

% ── helpers ───────────────────────────────────────────────────────────────────

function d = firstWritable(candidates)
% Probe by actually creating a file. isfolder plus a permissions guess is wrong
% often enough - a read-only network share, a managed institutional machine, a
% toolbox folder owned by the installer - and the failure would otherwise
% surface after a download rather than before one.
d = '';
for k = 1:numel(candidates)
    parent = fileparts(candidates{k});
    if ~isfolder(parent)
        [ok, ~] = mkdir(parent);
        if ~ok; continue; end
    end
    probe = fullfile(parent, ['.nestapp_write_probe_' char(matlab.lang.internal.uuid())]);
    fid = fopen(probe, 'w');
    if fid < 0; continue; end
    fclose(fid);
    delete(probe);
    d = candidates{k};
    return
end
end

function sub = soleSubfolder(dir_)
% GitHub archives extract to a single <name>-<version>/ folder, so the tree is
% one level down under a name that must not be hardcoded.
e = dir(dir_);
e = e([e.isdir] & ~ismember({e.name}, {'.', '..'}));
if numel(e) ~= 1
    error('nestapp:aaratepUnexpectedArchive', ...
        ['Expected the archive to hold one top-level folder, found %d. ' ...
         'Nothing has been installed.'], numel(e));
end
sub = fullfile(dir_, e(1).name);
end

function tag = newerReleaseOrEmpty(rel)
% Report a newer upstream release WITHOUT installing it - see aaratepRelease
% for why the pin is the pin. Best effort: no network, a proxy, or a GitHub
% outage must never turn into a failed install.
tag = '';
try
    url  = sprintf('https://api.github.com/repos/%s/releases/latest', rel.repo);
    opts = weboptions('Timeout', 5, 'ContentType', 'json');
    latest = webread(url, opts);
    if isfield(latest, 'tag_name') && ~strcmp(latest.tag_name, rel.tag)
        tag = latest.tag_name;
    end
catch
    % Silent on purpose.
end
end

function stampVersion(dir_, tag)
% A one-line marker naming what was installed. Written last, so it is present
% only if everything before it succeeded.
fid = fopen(fullfile(dir_, '.nestapp-aaratep-version'), 'w');
if fid < 0; return; end   % a missing stamp degrades to "unknown", not an error
fprintf(fid, '%s\n', tag);
fclose(fid);
end

function removeQuietly(d)
if ~isfolder(d); return; end
try
    rmdir(d, 's');
catch
    % A locked file costs disk, not correctness.
end
end
