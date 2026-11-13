% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef AaratepInstallTest < NestappTestCase
% AARATEPINSTALLTEST  The AARATEP installer's decisions, without the download.
%
%   installAaratep fetches a ~300-file tree from GitHub, so the download itself
%   is not testable here - a suite that reaches the network fails on a train,
%   behind a proxy, and on a runner with no egress, and a test that skips in
%   those conditions is a test that never runs.
%
%   The decisions AROUND the download are where the bugs were, and they are all
%   in aaratepStatus, which is a pure query with no network and no side effects.
%   That split exists for this file: an earlier version tested through
%   installAaratep, and one case reached GitHub and downloaded 0.6 MB whenever
%   the network happened to be up - putting real external I/O in the suite whose
%   entire contract is that it needs nothing.
%
%   `pure/` because nothing here needs EEGLAB or a display. The fake trees are
%   a scratch directory plus one file named after the sentinel, which is all
%   the "is it installed" check looks at.

    properties (Constant)
        % installAaratep decides a tree is present by looking for this. Read
        % from aaratepRelease rather than restated, so the fake trees below
        % cannot drift from what the real check wants.
        Sentinel = [aaratepRelease().sentinel '.m']
    end

    methods (Test)

        function thePinIsAReleaseAndACommitThatAgree(tc)
        % The pin is the contract: it is what makes a template reproduce for
        % two different people. THIRD_PARTY_NOTICES.md records the same commit,
        % and the tag is a human-checkable name for it.
            rel = aaratepRelease();
            tc.verifyMatches(rel.tag, '^v\d+\.\d+', 'the pin should name a release');
            tc.verifyMatches(rel.commit, '^[0-9a-f]{40}$', 'a full commit sha');
            tc.verifySubstring(rel.zipUrl, rel.tag, ...
                'the download URL must fetch the pinned tag, not a branch');
            tc.verifyNotEmpty(rel.searchDirs);

            notices = fileread(fullfile(addNestappPath(), 'THIRD_PARTY_NOTICES.md'));
            tc.verifySubstring(notices, rel.commit, ...
                ['the licensing document records a different commit than the ' ...
                 'installer fetches - one of them is wrong']);
        end

        function aNewDownloadDoesNotGoInsideTheToolbox(tc)
        % An installed toolbox's folder is writable, so when third_party/ was
        % the first writable place tried, the helpers landed inside the
        % toolbox and the next nestapp update deleted them.
            rel = aaratepRelease();
            inToolbox = fullfile(nestappRoot(), 'third_party', 'aaratep');
            inPrefdir = fullfile(prefdir, 'nestapp', 'aaratep');
            tc.verifyLessThan(find(strcmp(rel.writeDirs, inPrefdir)), ...
                find(strcmp(rel.writeDirs, inToolbox)), ...
                'prefdir must be tried before the toolbox''s own folder');
            tc.verifyEqual(sort(rel.writeDirs), sort(rel.searchDirs), ...
                'every place written to must also be a place searched');
        end

        function askingOneFolderIsNotAskingTheSearchOrder(tc)
        % THE BUG THIS FILE EXISTS FOR. The "already installed?" check consulted
        % the global search order even when the caller named a destination, so
        % installAaratep('TargetDir', X) found the developer checkout at
        % third_party/aaratep, reported already-present, and never touched X -
        % an explicit destination silently ignored.
        %
        % The search order is CONTROLLED rather than observed: aaratepPath is
        % the first place aaratepStatus looks, so pointing it at a fake tree
        % makes the broad query succeed regardless of what this checkout has
        % vendored. An earlier version read the ambient state and guarded it
        % with an assume, which this suite forbids - a skip is a failure, and a
        % test that quietly opts out on a fresh clone is one that never runs.
            root  = scratchDir(tc);
            tree  = fullfile(root, 'planted');
            empty = fullfile(root, 'nothing_here');
            mkdir(tree); mkdir(empty);
            tc.writeFile(fullfile(tree, tc.Sentinel), '% fake');

            isolatePrefs(tc, 'aaratepPath');
            setpref('nestapp', 'aaratepPath', tree);

            broad = aaratepStatus();
            tc.assertTrue(broad.installed, ...
                'the planted tree was not found, so the pair below proves nothing');
            tc.verifyEqual(broad.dir, tree, 'aaratepPath is searched first');

            scoped = aaratepStatus(empty);
            tc.verifyFalse(scoped.installed, ...
                ['reported installed for a folder that is empty, so it is ' ...
                 'answering about somewhere else']);
            tc.verifyEmpty(scoped.dir);
        end

        function anInstalledTreeIsRecognisedAtTheFolderItIsIn(tc)
            d = fullfile(scratchDir(tc), 'aaratep');
            mkdir(d);
            tc.writeFile(fullfile(d, tc.Sentinel), '% fake');

            st = aaratepStatus(d);
            tc.verifyTrue(st.installed);
            tc.verifyEqual(st.dir, d);
        end

        function anUnstampedTreeIsReportedAsUnknownRatherThanAsThePin(tc)
        % The second fault found by hand, and the one that matters for a
        % reproducibility claim. The check finds a tree by its sentinel file,
        % which says nothing about its version - so the first version happily
        % announced "AARATEP v2.1.1 is already installed" over a years-old git
        % clone or a working copy on another branch. Telling someone they have
        % the pinned release when they do not is worse than saying nothing.
            d = fullfile(scratchDir(tc), 'aaratep');
            mkdir(d);
            tc.writeFile(fullfile(d, tc.Sentinel), '% a clone, not our install');

            st = aaratepStatus(d);
            tc.verifyTrue(st.installed, 'the tree is there');
            tc.verifyEmpty(st.version, 'but its version is not knowable');
            tc.verifyFalse(st.isPinned);
            tc.verifySubstring(lower(st.label), 'unknown', ...
                'an unstamped tree must not be reported as the pinned release');
            tc.verifyEmpty(strfind(st.label, aaratepRelease().tag), ...
                'and must not name the pinned tag at all');
        end

        function aStampedTreeIsNamedExactly(tc)
        % The other side of the same rule: a tree this app installed carries a
        % marker naming what it wrote, so it can be identified precisely.
            rel = aaratepRelease();
            d   = fullfile(scratchDir(tc), 'aaratep');
            mkdir(d);
            tc.writeFile(fullfile(d, tc.Sentinel), '% fake');
            tc.writeFile(fullfile(d, '.nestapp-aaratep-version'), rel.tag);

            st = aaratepStatus(d);
            tc.verifyEqual(st.version, rel.tag);
            tc.verifyTrue(st.isPinned);
            tc.verifySubstring(st.label, rel.tag);
            tc.verifyEmpty(strfind(lower(st.label), 'unknown'));
        end

        function aStampFromADifferentReleaseSaysSo(tc)
        % A stale install nestapp made itself - after the pin moves - must be
        % distinguishable from a current one, or upgrading is invisible.
            d = fullfile(scratchDir(tc), 'aaratep');
            mkdir(d);
            tc.writeFile(fullfile(d, tc.Sentinel), '% fake');
            tc.writeFile(fullfile(d, '.nestapp-aaratep-version'), 'v1.0.0');

            st = aaratepStatus(d);
            tc.verifyEqual(st.version, 'v1.0.0', 'what is actually there');
            tc.verifyFalse(st.isPinned);
            tc.verifySubstring(st.label, 'v1.0.0');
            tc.verifySubstring(st.label, aaratepRelease().tag, 'and what is wanted');
        end

    end

    methods (Access = private)
        function writeFile(~, path, text)
            fid = fopen(path, 'w');
            fprintf(fid, '%s\n', text);
            fclose(fid);
        end
    end
end
