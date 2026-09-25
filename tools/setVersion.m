% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function newVersion = setVersion(target, varargin)
% SETVERSION  Cut a release number: set it once, and stamp it everywhere else.
%   newVersion = SETVERSION('minor')     2.1.0 -> 2.2.0   (also 'major', 'patch')
%   newVersion = SETVERSION('2.2.0')     an explicit version
%   newVersion = SETVERSION(..., 'Date', 'YYYY-MM-DD')   release date (default today)
%   newVersion = SETVERSION(..., 'Root', dir)            repository root (default:
%                                                        the one this file is in)
%
%   src/nestappVersion.m holds THE version; every running consumer (About,
%   report headers, error bundles, saved results, package_toolbox) calls it.
%   Three static documents also state it, and cannot compute it, because
%   GitHub and citation tools read them as text:
%
%     CHANGELOG.md   [Unreleased] becomes "## [x.y.z] - date", a fresh empty
%                    [Unreleased] opens above it, and the compare/tag links at
%                    the foot move to the new version
%     CITATION.cff   version: and date-released:
%     README.md      the version badge
%
%   This function is the only thing meant to write any of them. It builds all
%   four new texts first and writes nothing unless every expected pattern was
%   found exactly once, so a file whose layout changed stops the release with
%   a message rather than being skipped and left behind. VersionTest then
%   fails if anything is edited by hand out of step.
%
%   It refuses a version that is not higher than the current one, and a
%   release whose [Unreleased] section is empty: a release with no notes is
%   usually a release cut by mistake.
%
%   After it runs: review the diff, commit, and tag vX.Y.Z on the release
%   commit.
%
%   See also: nestappVersion, package_toolbox

opts = struct('Date', char(datetime('today', 'Format', 'yyyy-MM-dd')), ...
              'Root', fileparts(fileparts(mfilename('fullpath'))));
for k = 1:2:numel(varargin)
    opts.(varargin{k}) = char(varargin{k + 1});
end
root = opts.Root;
if isempty(regexp(opts.Date, '^\d{4}-\d{2}-\d{2}$', 'once'))
    error('nestapp:setVersion:badDate', 'Date must be YYYY-MM-DD, got "%s".', opts.Date);
end

files.version   = fullfile(root, 'src', 'nestappVersion.m');
files.changelog = fullfile(root, 'CHANGELOG.md');
files.citation  = fullfile(root, 'CITATION.cff');
files.readme    = fullfile(root, 'README.md');

txt = structfun(@readText, files, 'UniformOutput', false);

tok = regexp(txt.version, '(?m)^v = ''(\d+\.\d+\.\d+)'';', 'tokens');
if ~isscalar(tok)
    error('nestapp:setVersion:noConstant', ...
        'nestappVersion.m must contain exactly one line  v = ''x.y.z'';');
end
oldVersion = tok{1}{1};
newVersion = nextVersion(oldVersion, char(target));

% ── nestappVersion.m ──────────────────────────────────────────────────────
out.version = replaceOnce(txt.version, '(?m)^v = ''\d+\.\d+\.\d+'';', ...
    sprintf('v = ''%s'';', newVersion), 'the version constant');

% ── CHANGELOG.md ──────────────────────────────────────────────────────────
notes = regexp(txt.changelog, '(?s)## \[Unreleased\]\s*(.*?)\n## ', 'tokens', 'once');
if isempty(notes) || isempty(strtrim(notes{1}))
    error('nestapp:setVersion:noNotes', ...
        'CHANGELOG.md''s [Unreleased] section is empty - write the release notes first.');
end
nl = newlineOf(txt.changelog);
c = replaceOnce(txt.changelog, '(?m)^## \[Unreleased\][ \t]*\r?$', ...
    sprintf('## [Unreleased]%s%s## [%s] - %s', nl, nl, newVersion, opts.Date), ...
    'the ## [Unreleased] heading');
link = regexp(c, '(?m)^\[Unreleased\]: (\S+)/compare/v\d+\.\d+\.\d+\.\.\.HEAD', 'tokens');
if ~isscalar(link)
    error('nestapp:setVersion:noLink', ...
        'CHANGELOG.md must end with one  [Unreleased]: <repo>/compare/vX.Y.Z...HEAD  link.');
end
repo = link{1}{1};
out.changelog = replaceOnce(c, ...
    '(?m)^\[Unreleased\]: \S+/compare/v\d+\.\d+\.\d+\.\.\.HEAD', ...
    sprintf('[Unreleased]: %s/compare/v%s...HEAD%s[%s]: %s/releases/tag/v%s', ...
        repo, newVersion, nl, newVersion, repo, newVersion), 'the [Unreleased] link');

% ── CITATION.cff ──────────────────────────────────────────────────────────
cc = replaceOnce(txt.citation, '(?m)^version:[ \t]*\S+', ...
    sprintf('version: %s', newVersion), 'version:');
out.citation = replaceOnce(cc, '(?m)^date-released:[ \t]*\S+', ...
    sprintf('date-released: "%s"', opts.Date), 'date-released:');

% ── README.md ─────────────────────────────────────────────────────────────
out.readme = replaceOnce(txt.readme, 'badge/version-\d+\.\d+\.\d+-', ...
    sprintf('badge/version-%s-', newVersion), 'the version badge');

% Every pattern was found; only now touch the disk.
for f = fieldnames(files)'
    writeText(files.(f{1}), out.(f{1}));
end
fprintf('nestapp %s -> %s (%s). Review the diff, commit, then tag v%s.\n', ...
    oldVersion, newVersion, opts.Date, newVersion);
end

% ── helpers ─────────────────────────────────────────────────────────────────

function v = nextVersion(old, target)
    cur = sscanf(old, '%d.%d.%d')';
    switch lower(target)
        case 'major'; nxt = [cur(1) + 1, 0, 0];
        case 'minor'; nxt = [cur(1), cur(2) + 1, 0];
        case 'patch'; nxt = [cur(1), cur(2), cur(3) + 1];
        otherwise
            if isempty(regexp(target, '^\d+\.\d+\.\d+$', 'once'))
                error('nestapp:setVersion:badTarget', ...
                    'Target must be major, minor, patch or X.Y.Z; got "%s".', target);
            end
            nxt = sscanf(target, '%d.%d.%d')';
    end
    v = sprintf('%d.%d.%d', nxt);
    if ~isHigher(nxt, cur)
        error('nestapp:setVersion:notHigher', ...
            '%s is not higher than the current %s.', v, old);
    end
end

function tf = isHigher(a, b)
    d = find(a ~= b, 1);
    tf = ~isempty(d) && a(d) > b(d);
end

function s = replaceOnce(s, pattern, replacement, what)
% Exactly one match, or the release stops: zero means the file's layout
% changed, two means the pattern is no longer specific enough to trust.
    n = numel(regexp(s, pattern, 'start'));
    if n ~= 1
        error('nestapp:setVersion:pattern', ...
            'Expected exactly one %s, found %d. Nothing was written.', what, n);
    end
    s = regexprep(s, pattern, regexptranslate('escape', replacement));
end

function nl = newlineOf(s)
    if contains(s, sprintf('\r\n')); nl = sprintf('\r\n'); else; nl = newline; end
end

function s = readText(p)
    if ~isfile(p)
        error('nestapp:setVersion:missing', 'Missing %s.', p);
    end
    fid = fopen(p, 'r', 'n', 'UTF-8');
    s = fread(fid, [1 Inf], '*char');
    fclose(fid);
end

function writeText(p, s)
    fid = fopen(p, 'w', 'n', 'UTF-8');
    fwrite(fid, s, 'char');
    fclose(fid);
end
