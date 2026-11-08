% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function root = findEeglab(searchRoots)
% FINDEEGLAB  Look for an EEGLAB folder in the places a new user is likely to have one.
%   root = FINDEEGLAB() returns the folder holding eeglab.m, or '' if none
%   was found. Where several are found, the highest-numbered one wins
%   (eeglab2026.1.0 over eeglab2025.0.0).
%
%   root = FINDEEGLAB(searchRoots) searches these folders instead of the
%   defaults (for tests).
%
%   Why this exists: "set the EEGLAB folder in Preferences" was the step a new
%   user was most likely to miss, and the two usual places EEGLAB ends up are
%   predictable:
%     - <add-ons>/Toolboxes/AdditionalSoftware, where the .mltbx installer puts
%       the EEGLAB it offers to download;
%     - userpath (Documents/MATLAB), where people unzip things by hand.
%   Only a few folder levels are searched, never the whole tree, so a miss
%   costs milliseconds.
%
%   See also: ensureEeglabReady

if nargin < 1
    searchRoots = defaultRoots();
end

hits = {};
for r = searchRoots(:)'
    if isempty(r{1}) || ~isfolder(r{1}); continue; end
    for depth = 1:3
        pattern = fullfile(r{1}, strjoin(repmat({'*'}, 1, depth), filesep), 'eeglab.m');
        d = dir(pattern);
        hits = [hits, {d.folder}]; %#ok<AGROW>
    end
end
% The folder name carries the version (eeglab2026.1.0), so the newest sorts last.
if isempty(hits)
    root = '';
    return
end
hits = unique(hits);
[~, order] = sort(cellfun(@leaf, hits, 'UniformOutput', false));
root = hits{order(end)};
end

% ── helpers ───────────────────────────────────────────────────────────────────

function roots = defaultRoots()
roots = {};
addons = addonsFolder();
if ~isempty(addons)
    roots{end+1} = fullfile(addons, 'Toolboxes', 'AdditionalSoftware');
end
up = strtok(userpath, pathsep);
if ~isempty(up)
    roots{end+1} = up;
end
end

function f = addonsFolder()
% The add-ons folder setting arrived after R2023b on some platforms, so fall
% back to inferring it from where this installed toolbox lives
% (<add-ons>/Toolboxes/nestapp/...).
f = '';
try
    s = settings;
    f = char(s.matlab.addons.InstallationFolder.ActiveValue);
catch
end
if isempty(f)
    tok = regexp(nestappRoot(), ['^(.*)' regexptranslate('escape', [filesep 'Toolboxes' filesep])], ...
                 'tokens', 'once');
    if ~isempty(tok); f = tok{1}; end
end
end

function n = leaf(p)
[~, n, ext] = fileparts(p);
n = [n ext];   % eeglab2026.1.0 would otherwise lose its ".0"
end
