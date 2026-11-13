% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function outFile = package_toolbox(varargin)
% PACKAGE_TOOLBOX  Build nestapp as an installable MATLAB toolbox (.mltbx).
%   outFile = PACKAGE_TOOLBOX()
%   outFile = PACKAGE_TOOLBOX('OutputDir', d)
%
%   Produces nestapp-<version>.mltbx. A user double-clicks it and MATLAB
%   installs nestapp, puts it on the path permanently, and lists it under
%   Add-Ons where it can be updated or uninstalled. That replaces the current
%   install instructions - download a zip, unzip it somewhere, cd there, run a
%   script, and re-do the cd every session - which is the step novice users
%   most often get wrong.
%
%   WHAT GOES IN IS WHAT GIT TRACKS. The file list comes from `git ls-files`,
%   not from a glob of the working tree, and that is deliberate: this
%   repository keeps a full EEGLAB install under eeglab2026.0.0/ (2,699
%   untracked files) and a vendored AARATEP tree under third_party/, neither of
%   which nestapp may redistribute. A glob would quietly ship both - a
%   licensing problem and a ~200 MB package. Asking git means the package
%   contains exactly what the repository does.
%
%   EEGLAB IS DECLARED, NOT BUNDLED, via RequiredAdditionalSoftware, so the
%   installer tells the user what else they need and where to get it rather
%   than failing later with an undefined-function error.
%
%   Tests are excluded. They are for developers, who have the repository.
%
%   See also: matlab.addons.toolbox.ToolboxOptions, nestappVersion, nestappDoctor

p = inputParser;
p.addParameter('OutputDir', pwd, @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
outDir = char(p.Results.OutputDir);

root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));   % genpath: src/ has subfolders
version = nestappVersion();

% ── the file list, from git ───────────────────────────────────────────────────
tracked = gitTracked(root);
keep    = startsWith(tracked, 'src/') | ...
          ismember(tracked, {'run_nestapp.m', 'LICENSE', 'README.md', ...
                             'THIRD_PARTY_NOTICES.md', 'CITATION.cff', ...
                             'CHANGELOG.md'}) | ...
          startsWith(tracked, 'docs/');
files = fullfile(root, strrep(tracked(keep), '/', filesep));

missing = files(~isfile(files) & ~isfolder(files));
if ~isempty(missing)
    error('nestapp:packageMissingFile', ...
        'git tracks %d file(s) that are not on disk: %s', ...
        numel(missing), strjoin(missing(1:min(3,end)), ', '));
end

% ── the toolbox ───────────────────────────────────────────────────────────────
% A FIXED identifier. MATLAB uses it to recognise an upgrade of an already
% installed toolbox rather than a second copy, so it must never change between
% releases - which is why it is a literal here and not a fresh uuid per build.
IDENTIFIER = 'dd153273-6fb0-425f-bf4c-b0cb2a3ed799';

opts = matlab.addons.toolbox.ToolboxOptions(root, IDENTIFIER, ...
    'ToolboxFiles', files);

opts.ToolboxName    = 'nestapp';
opts.ToolboxVersion = version;
opts.AuthorName     = 'Aref Pariz and Wesley Dunne';
opts.Summary        = 'Point-and-click TMS-EEG cleaning and TEP analysis, built on EEGLAB.';
opts.Description    = [ ...
    'nestapp is a MATLAB app for cleaning and analysing TMS-EEG recordings. ' ...
    'It wraps EEGLAB, TESA and the surrounding plugin stack in a pipeline ' ...
    'builder, so a researcher can preprocess a cohort, inspect what happened ' ...
    'to every file, and produce publication-ready TEP figures and peak ' ...
    'measurements without writing code. Every processed file carries its full ' ...
    'pipeline - steps, parameters, timestamp - inside EEG.history.' newline newline ...
    'Launch with:  nestapp' newline newline ...
    'Requires a separate EEGLAB installation; run nestappDoctor to check ' ...
    'what is missing.'];

opts.MinimumMatlabRelease = 'R2023b';

% Only src/ and its subfolders go on the user's path - not docs/, and not the
% repository root, which would put README.md's folder on the path for nothing.
opts.ToolboxMatlabPath = srcPathEntries(root);

logo = fullfile(root, 'src', 'LogoNest.jpg');
if isfile(logo); opts.ToolboxImageFile = logo; end

% Declared rather than shipped: nestapp cannot redistribute EEGLAB, and the
% installer surfacing this is far better than an undefined-function error on
% first launch.
% One entry PER PLATFORM: the packager rejects 'all' and wants a concrete
% platform per record, so EEGLAB is declared four times over the same URL.
% The URL must be the archive itself: the installer downloads and unzips it
% into <add-ons>/Toolboxes/AdditionalSoftware, where findEeglab looks. The
% 2.1.0 package pointed at eeglab.org/download/, a web page (since a 404),
% so the offered download could never have worked.
platforms = {'win64', 'glnxa64', 'maca64', 'mac'};
req = struct('Name', {}, 'Platform', {}, 'DownloadURL', {}, 'LicenseURL', {});
for k = 1:numel(platforms)
    req(k) = struct( ...
        'Name',        'EEGLAB', ...
        'Platform',    platforms{k}, ...
        'DownloadURL', 'https://sccn.ucsd.edu/eeglab/currentversion/eeglab_current.zip', ...
        'LicenseURL',  'https://github.com/sccn/eeglab/blob/develop/LICENSE');
end
opts.RequiredAdditionalSoftware = req;

if ~isfolder(outDir); mkdir(outDir); end
outFile = fullfile(outDir, sprintf('nestapp-%s.mltbx', version));
opts.OutputFile = outFile;

matlab.addons.toolbox.packageToolbox(opts);
fprintf('packaged %d files -> %s (%.1f MB)\n', numel(files), outFile, ...
        dir(outFile).bytes / 1e6);
end

% ── helpers ───────────────────────────────────────────────────────────────────

function files = gitTracked(root)
[status, out] = system(sprintf('git -C "%s" ls-files', root));
if status ~= 0
    error('nestapp:packageNeedsGit', ...
        ['package_toolbox asks git what to ship, so it must run inside a ' ...
         'checkout. `git ls-files` failed: %s'], strtrim(out));
end
files = strtrim(strsplit(strtrim(out), newline));
end

function entries = srcPathEntries(root)
% src/ plus every subfolder that holds code, minus class folders - MATLAB
% resolves @nestapp from src/ being on the path and rejects a class folder
% added directly.
d = dir(fullfile(root, 'src', '**', '*.m'));
folders = unique({d.folder});
folders = folders(~contains(folders, [filesep '@']));
entries = folders(:)';
end
