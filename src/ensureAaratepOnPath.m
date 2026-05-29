function ensureAaratepOnPath()
% ENSUREAARATEPONPATH  Idempotent addpath for the vendored AARATEP tree.
%   Adds third_party/aaratep and its Common subtree to the MATLAB path on
%   first call. Subsequent calls are cheap no-ops.
%
%   Two bundled-fork subtrees are kept OFF the path so they cannot shadow
%   the user's own installs:
%     - Common/ThirdParty/FromEEGLab  forked EEGLAB functions (topoplot,
%       epoch, pop_loadbv, pop_resample_mod, ...).
%     - Common/ThirdParty/FastICA     a bundled FastICA copy. The user
%       should run the FastICA they think they have (normally EEGLAB's).
%       If NO other fastica is on the path, the bundled copy is added back
%       as a fallback so the AARATEP pipeline still runs. When the user's
%       FastICA differs from the bundled (tested) one, a one-time warning is
%       printed - the pipeline was validated against the bundled version.
%
%   Vendored from chriscline/AARATEPPipeline @ be75262 under MIT license.
%   See THIRD_PARTY_NOTICES.md.

persistent ready
if ~isempty(ready) && ready
    return
end

repoRoot   = fileparts(fileparts(mfilename('fullpath')));   % src/ -> repo root
aaratepDir = fullfile(repoRoot, 'third_party', 'aaratep');

if ~isfolder(aaratepDir)
    error('ensureAaratepOnPath:Missing', ...
        ['AARATEP vendored tree not found at %s. ' ...
         'Re-run: cd third_party && git clone --depth 1 ' ...
         'https://github.com/chriscline/AARATEPPipeline.git aaratep'], ...
        aaratepDir);
end

% Drop the bundled-fork subtrees from genpath so they do not shadow the
% user's EEGLAB / FastICA installs.
thirdParty = fullfile(aaratepDir, 'Common', 'ThirdParty');
shadowDirs = {fullfile(thirdParty, 'FromEEGLab'), fullfile(thirdParty, 'FastICA')};
allPaths   = strsplit(genpath(aaratepDir), pathsep);
allPaths   = allPaths(~cellfun(@isempty, allPaths));
keep       = true(1, numel(allPaths));
for i = 1:numel(shadowDirs)
    keep = keep & ~startsWith(allPaths, shadowDirs{i});
end
addpath(strjoin(allPaths(keep), pathsep));

resolveFastICA(fullfile(thirdParty, 'FastICA'));

ready = true;
end

% ── helpers ───────────────────────────────────────────────────────────────────

function resolveFastICA(bundledFasticaDir)
% Prefer the user's fastica; fall back to the bundled copy only if none
% exists; warn once when the user's version differs from the bundled one.
bundledFile = fullfile(bundledFasticaDir, 'fastica.m');
userFile    = which('fastica');

if isempty(userFile)
    % No FastICA anywhere - add the bundled copy so AARATEP can still run.
    if isfile(bundledFile)
        addpath(bundledFasticaDir);
        warning('nestapp:aaratepFastICAFallback', ...
            ['No FastICA found on the path; using the FastICA bundled with ' ...
             'AARATEP. Install FastICA (or EEGLAB''s) for your expected version.']);
    end
    return
end

if ~isfile(bundledFile)
    return   % nothing to compare against
end

vUser    = fasticaVersion(userFile);
vBundled = fasticaVersion(bundledFile);
if fasticaDiffers(userFile, bundledFile, vUser, vBundled)
    warning('nestapp:aaratepFastICAMismatch', ...
        ['Using FastICA %s (%s), but AARATEP was tested with the bundled ' ...
         'FastICA %s. Results may differ slightly from the published pipeline.'], ...
        verLabel(vUser), userFile, verLabel(vBundled));
end
end

function tf = fasticaDiffers(userFile, bundledFile, vUser, vBundled)
% Prefer a version-string comparison; fall back to byte comparison.
if ~isempty(vUser) && ~isempty(vBundled)
    tf = ~strcmp(vUser, vBundled);
else
    tf = ~isequal(fileread(userFile), fileread(bundledFile));
end
end

function v = fasticaVersion(file)
% Extract a version token (e.g. '2.5') from a fastica.m, '' if not found.
v = '';
try
    txt = fileread(file);
catch
    return
end
tok = regexp(txt, 'version\s*[:=]?\s*v?(\d+\.\d+)', 'tokens', 'once', 'ignorecase');
if ~isempty(tok)
    v = tok{1};
end
end

function s = verLabel(v)
if isempty(v); s = '(unknown version)'; else; s = ['v' v]; end
end
