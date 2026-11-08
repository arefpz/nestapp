% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function info = aaratepRelease()
% AARATEPRELEASE  The AARATEP release nestapp is pinned to, and where it goes.
%   info = AARATEPRELEASE()
%
%   Fields:
%     .tag          the pinned upstream release, e.g. 'v2.1.1'
%     .commit       the commit that tag resolves to
%     .repo         'owner/name' on GitHub
%     .zipUrl       archive URL for .tag
%     .sentinel     a function that must exist once installed
%     .searchDirs   where to look for an installed tree, in order
%     .writeDirs    where installAaratep may put a new tree, in order
%
%   A PINNED RELEASE, NOT "LATEST", and this is the important decision in the
%   whole feature. AARATEP defines the behaviour of a published pipeline, so
%   two users running the same nestapp template must get the same AARATEP -
%   otherwise a result depends on the day the helpers were downloaded, which is
%   exactly the kind of silent divergence a preprocessing tool exists to
%   prevent. installAaratep therefore fetches .tag; it can REPORT that a newer
%   release exists, and will install one only when explicitly asked.
%
%   The pin is the same commit THIRD_PARTY_NOTICES.md records, and it happens
%   to be exactly upstream's v2.1.1 (the annotated tag dereferences to it), so
%   the tag is a stable, human-checkable name for what we already vendored
%   rather than a second source of truth.
%
%   See also: installAaratep, ensureAaratepOnPath, THIRD_PARTY_NOTICES.md

info.repo     = 'chriscline/AARATEPPipeline';
info.tag      = 'v2.1.1';
info.commit   = 'be75262af689d4e8e5053c05aaa4ed3be258350a';
info.zipUrl   = sprintf('https://github.com/%s/archive/refs/tags/%s.zip', ...
                        info.repo, info.tag);
info.sentinel = 'c_TMSEEG_Preprocess_AARATEPPipeline';
info.searchDirs = searchDirs();
% Writing goes to prefdir BEFORE third_party/. An installed toolbox's own
% folder is writable, so trying third_party/ first put the download inside the
% toolbox, and the next nestapp update deleted it (seen installing 2.1.0 and
% reinstalling it). A developer's own clone in third_party/ is still found
% first by searchDirs; only where a NEW download lands changes.
toolboxDir = fullfile(nestappRoot(), 'third_party', 'aaratep');
info.writeDirs = [setdiff(info.searchDirs, {toolboxDir}, 'stable'), {toolboxDir}];
end

function dirs = searchDirs()
% Where an installed tree may live, in the order ensureAaratepOnPath should
% prefer, and the order installAaratep should try to write.
%
% THREE LOCATIONS, because nestapp is used three ways and only the first two
% existed before it could be installed as a toolbox:
%
%   1. A path the user set explicitly - someone with their own clone, or a
%      shared lab copy on a network drive, should not have a second download
%      forced on them.
%   2. third_party/ next to the source. This is the developer checkout, the
%      layout THIRD_PARTY_NOTICES and CONTRIBUTING both describe, and it must
%      keep working exactly as before.
%   3. A folder under prefdir. This is the one that makes the wizard possible
%      for an installed toolbox: writing into an installed toolbox's own
%      folder would work today and then be deleted by the next upgrade or
%      uninstall, so a ~300-file download would silently vanish.
dirs = {};

% GUARDED BY ispref, NOT getpref-with-a-default. getpref(group, pref, default)
% CREATES the preference when it is absent, so the three-argument form would
% make this metadata function write to the user's settings merely by being
% called - and it did: reading the search order once left an empty
% 'aaratepPath' behind. A function that reports where things are has no
% business changing anything.
if ispref('nestapp', 'aaratepPath')
    userSet = getpref('nestapp', 'aaratepPath');
    if ~isempty(userSet)
        dirs{end+1} = char(userSet);
    end
end

dirs{end+1} = fullfile(nestappRoot(), 'third_party', 'aaratep');
dirs{end+1} = fullfile(prefdir, 'nestapp', 'aaratep');
end
