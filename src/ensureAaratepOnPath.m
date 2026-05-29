function ensureAaratepOnPath()
% ENSUREAARATEPONPATH  Idempotent addpath for the vendored AARATEP tree.
%   Adds third_party/aaratep and its Common subtree to the MATLAB path on
%   first call. Excludes the bundled Common/ThirdParty/FromEEGLab subtree
%   because it ships forked copies of EEGLAB functions (topoplot, epoch,
%   pop_loadbv, pop_resample_mod, ...) that would otherwise shadow the
%   user's real EEGLAB install. Subsequent calls are cheap no-ops.
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

% Drop the FromEEGLab subtree from genpath - those forks would shadow
% the user's EEGLAB install.
shadowDir = fullfile(aaratepDir, 'Common', 'ThirdParty', 'FromEEGLab');
allPaths  = strsplit(genpath(aaratepDir), pathsep);
keep      = ~startsWith(allPaths, shadowDir);
safePath  = strjoin(allPaths(keep & ~cellfun(@isempty, allPaths)), pathsep);

addpath(safePath);
ready = true;
end
