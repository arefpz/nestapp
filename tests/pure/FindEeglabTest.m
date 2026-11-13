% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef FindEeglabTest < NestappTestCase
% FINDEEGLABTEST  Finding an EEGLAB the user never pointed nestapp at.
%
%   The search roots are passed in, so these use fake trees (a folder holding
%   an eeglab.m) and never look at what this machine really has installed.

    methods (Test)

        function findsTheInstallerLayout(tc)
        % The .mltbx installer unzips EEGLAB under AdditionalSoftware, a level
        % or two below the root that is searched.
            root = scratchDir(tc);
            want = fakeEeglab(fullfile(root, 'EEGLAB', 'eeglab2026.1.0'));
            tc.verifyEqual(findEeglab({root}), want);
        end

        function prefersTheNewestVersion(tc)
            root = scratchDir(tc);
            fakeEeglab(fullfile(root, 'eeglab2025.0.0'));
            want = fakeEeglab(fullfile(root, 'eeglab2026.1.0'));
            fakeEeglab(fullfile(root, 'eeglab2026.0.0'));
            tc.verifyEqual(findEeglab({root}), want);
        end

        function returnsEmptyWhenThereIsNone(tc)
            root = scratchDir(tc);
            mkdir(fullfile(root, 'something', 'else'));
            tc.verifyEmpty(findEeglab({root}));
            tc.verifyEmpty(findEeglab({fullfile(root, 'missing')}), ...
                'a search root that does not exist is skipped, not an error');
        end

        function doesNotDigThroughWholeTrees(tc)
        % Searching is limited to a few levels so a miss stays cheap. An
        % eeglab.m buried deeper than that is not an install anyone unzipped.
            root = scratchDir(tc);
            fakeEeglab(fullfile(root, 'a', 'b', 'c', 'eeglab'));
            tc.verifyEmpty(findEeglab({root}));
        end

    end
end

function folder = fakeEeglab(folder)
mkdir(folder);
fclose(fopen(fullfile(folder, 'eeglab.m'), 'w'));
end
