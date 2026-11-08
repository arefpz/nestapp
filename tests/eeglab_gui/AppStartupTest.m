% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef AppStartupTest < NestappTestCase
% APPSTARTUPTEST  What must be true of a launched app: its steps, its geometry, its dirty flag.
%
%   The only file in eeglab_gui/, and the folder exists for this file's first
%   test. A2.1 is the regression the plan named when it chose folders over
%   tags: the old suite split it in two because its taxonomy had nowhere to put
%   a test that needs BOTH EEGLAB and a display -
%   regression/test_startupStepAvailability.m grepped the ordering out of the
%   source, ui/test_startupPicker.m checked the resulting tree, 151 lines and
%   two cross-referencing headers for one fact. Here it is one test.
%
%   EVERY TEST LAUNCHES ITS OWN APP. 3.5 s each, which is the most expensive
%   thing in the suite, and shared deliberately nowhere: two of these four
%   mutate app state (a resize, a save) and one needs a launch under conditions
%   it sets up itself. A shared app would make the results order-dependent,
%   which is the failure the old ui/ suite had - a leaked eeglab('nogui') and a
%   leaked path memo made it fragile on a second run in one session.
%
%   Ledger rows: A2.1 (cold start hid 32 of 54 steps), A3.9 (resize re-entrancy),
%   A3.10 (tabs excluded by type, not by name), A1.7 (a saved pipeline stayed
%   marked unsaved). A3.11 was dropped - it guarded a refactor, not a defect.

    properties (Constant)
        % nestapp/enforceMinWindowSize's MIN_SIZE. Stated once here; the
        % arithmetic itself is WindowClampTest's subject, and this file only
        % needs to drive the window below the minimum.
        MinSize = [650 420]
    end

    methods (Test)

        % ── A2.1 ─────────────────────────────────────────────────────────────

        function aColdStartOffersEveryInstalledStep(tc)
        % THE ROW THIS FOLDER EXISTS FOR. populateStepsTree asks availableSteps
        % what is installed; availableSteps probes which() per requirement; and
        % nothing EEGLAB provides - including every plugin - resolves until
        % eeglab() has run its plugin scan. startupFcn did not do that first,
        % so on a stock cold MATLAB the picker withheld 32 of 54 steps.
        %
        % The old test could not check this and said so in its own header: "a
        % suite running in a session that already has EEGLAB loaded cannot
        % reproduce the cold start at all", so it grepped startupFcn for the
        % order of three calls instead. That is a source check on a sequence,
        % and it passes or fails on how the method is written.
        %
        % It CAN be reproduced: take every EEGLAB folder off the path first.
        % startupFcn's own loadPrefs/initEeglab are then what has to put them
        % back, which is exactly the behaviour under test.
            hidden = tc.hideEeglabEntirely();

            % POSITIVE CONTROL, and the whole reason this test means anything:
            % prove the cold state was actually reached before launching. If
            % EEGLAB were still resolving, the app would trivially pass and the
            % test would be measuring nothing. 23 of 55 is the documented
            % cold-start number.
            coldCount = numel(availableSteps());
            tc.assertLessThan(coldCount, 30, sprintf( ...
                ['the cold state was not reached - %d steps still resolve, so ' ...
                 'this test cannot detect the regression it exists for'], coldCount));

            app = launchApp(tc);
            clear hidden   %#ok<NASGU> % path restored; the app has done its work

            offered = tc.stepLeafNames(app);
            % BY NAME, not by count. Agreement by name is strictly stronger and
            % subsumes what used to be a second test with a second 3.5 s app
            % launch: it catches the regression (a short tree) and also a
            % picker that offers a step the pre-flight would refuse, which is
            % the same contract seen from the other side.
            tc.verifyEqual(sort(offered), sort({availableSteps().name}), ...
                ['the picker and the pre-flight disagree about what is ' ...
                 'runnable. startupFcn must bring EEGLAB up (loadPrefs, then ' ...
                 'initEeglab) BEFORE populateStepsTree asks which() what is ' ...
                 'available.']);
            tc.verifyGreaterThan(numel(offered), coldCount, sprintf( ...
                ['the tree was built while EEGLAB was still hidden: %d steps ' ...
                 'offered, the same as the cold count'], numel(offered)));
            tc.verifyGreaterThan(numel(offered), 0.9 * numel(stepRegistry()), ...
                'almost the whole registry is installed here, so almost all of it should be offered');
        end

        % ── A3.9 / A3.10 ─────────────────────────────────────────────────────

        function aResizeBelowTheMinimumSettlesInsteadOfRecursing(tc)
        % Ledger A3.9. enforceMinWindowSize writes Position, which re-fires this
        % very callback, and the `drawnow limitrate` inside lets the re-entry
        % actually run - so unguarded the two feed each other for as long as
        % resize events keep arriving.
        %
        % Behavioural, where the old check was a source grep for the absence of
        % one exact assignment string. If the guard were gone this errors on a
        % recursion limit rather than failing an assertion, which is still red
        % and still names the right test.
            app = launchApp(tc);
            pos = app.UIFigure.Position;

            % Drag the bottom edge up and the right edge in, both past the
            % minimum: the shape that made the window ratchet up the screen.
            requested = [pos(1), pos(2) + 200, 500, 300];
            app.UIFigure.Position = requested;
            % Waited for, not slept on - and the predicate is keyed on the
            % position having CHANGED, not on it being big enough. "Big enough"
            % is already true of the window we started from, so it would wait
            % for nothing and then read back the pre-resize geometry. That was
            % the second version of this test, and it failed one run in four
            % for a different reason than the first.
            waitFor(tc, @() ~isequal(app.UIFigure.Position, pos) && ...
                          all(app.UIFigure.Position(3:4) >= tc.MinSize), ...
                    'the clamp to run after the resize');

            got = app.UIFigure.Position;
            tc.verifyGreaterThanOrEqual(got(3:4), tc.MinSize, ...
                'the window was left smaller than its minimum');
            tc.verifyEqual(got(2) + got(4), requested(2) + requested(4), ...
                'the top edge moved - this is the creep the clamp exists to stop');
            tc.verifyFalse(app.isResizing, ...
                ['the re-entrancy guard was left set, so every later resize is ' ...
                 'now a no-op and the window can never be clamped again']);
        end

        function aSecondResizeIsStillHonouredAfterTheFirst(tc)
        % The guard's other failure mode, and the one a single resize cannot
        % see: a guard that is set and never cleared makes the FIRST resize
        % work and every one after it silently do nothing.
            app = launchApp(tc);
            pos = app.UIFigure.Position;

            % WIDTH ONLY, and the two resizes differ in `left`. Two reasons,
            % both about making the test decide one thing: shrinking only the
            % width means the clamp does not move the window vertically, so
            % nothing here depends on the top-edge arithmetic that
            % WindowClampTest already owns and nothing risks being pushed
            % off-screen by the window manager. And moving `left` between them
            % means the second expected position differs from the first, so a
            % dropped second resize cannot satisfy the wait with a value left
            % over from the first.
            wantH = pos(4);
            app.UIFigure.Position = [pos(1), pos(2), 500, wantH];
            waitFor(tc, @() isequal(app.UIFigure.Position, ...
                                    [pos(1), pos(2), tc.MinSize(1), wantH]), ...
                    'the first clamp');
            first = app.UIFigure.Position;

            app.UIFigure.Position = [pos(1) + 20, pos(2), 400, wantH];
            waitFor(tc, @() isequal(app.UIFigure.Position, ...
                                    [pos(1) + 20, pos(2), tc.MinSize(1), wantH]), ...
                    'the second clamp - a guard that is never cleared drops it');
            second = app.UIFigure.Position;

            tc.verifyEqual(first,  [pos(1),      pos(2), tc.MinSize(1), wantH]);
            tc.verifyEqual(second, [pos(1) + 20, pos(2), tc.MinSize(1), wantH], ...
                'the second resize was dropped - the guard is not being cleared');
        end

        function noTabIsCapturedForRescaling(tc)
        % Ledger A3.10. A Tab is sized by its TabGroup and its Position is
        % READ-ONLY, so capturing one makes the NEXT resize throw. The old
        % implementation excluded them by listing the four tab names that
        % existed, and adding a fifth broke resizing silently until someone
        % resized - and Stage 7 did change the tab set.
            app    = launchApp(tc);
            layout = app.baseLayout;
            names  = fieldnames(layout);

            % Positive control first: a captured layout that is empty, or an app
            % with no tabs, would make the exclusion below vacuously true.
            tc.assertNotEmpty(names, 'nothing was captured, so nothing is excluded');
            tabs = findall(app.TabGroup, 'Type', 'uitab');
            tc.assertNotEmpty(tabs, 'the app has no tabs, so this rule proves nothing');

            isTab = cellfun(@(n) isa(app.(n), 'matlab.ui.container.Tab'), names);
            tc.verifyFalse(any(isTab), sprintf( ...
                'captured %d Tab(s) whose Position is read-only: %s', ...
                sum(isTab), strjoin(names(isTab), ', ')));

            % And the rule stated as what it protects: every captured Position
            % must be writable, because that is what the next resize does.
            for k = 1:numel(names)
                h = app.(names{k});
                tc.verifyWarningFree(@() set(h, 'Position', layout.(names{k}).pos), ...
                    sprintf('%s was captured but its Position cannot be written', names{k}));
            end
        end

        function theMainWindowStaysLightUnderAnyDesktopTheme(tc)
        % Under OS dark mode the window came out half light, half dark: MATLAB
        % re-coloured what nestapp had not fixed. createComponents goes
        % through appWindow, which pins the light theme.
            app = launchApp(tc);
            fig = app.UIFigure;
            if ~isprop(fig, 'ThemeMode'); return; end   % before R2025a: no themes
            tc.verifyEqual(char(fig.ThemeMode), 'manual');
            tc.verifyEqual(char(fig.Theme.BaseColorStyle), 'light');
        end

        % ── A1.7 ─────────────────────────────────────────────────────────────

        function savingAPipelineClearsTheUnsavedMarker(tc)
        % Ledger A1.7. The original defect was reading a path out of `uisave`,
        % which does not return one, so the save appeared to fail and
        % pipelineDirty was never cleared - the app kept warning about unsaved
        % changes to a pipeline that was on disk.
            isolatePrefs(tc, {'lastPipelineFolder', 'recentPipelines'});
            outDir = scratchDir(tc);
            shadowFunction(tc, 'uiputfile', {'saved_pipeline.mat', outDir});

            app = launchApp(tc);
            tc.makePipelineDirty(app);
            tc.assertTrue(app.pipelineDirty, ...
                'could not get the app into an unsaved state, so nothing is being cleared');

            tc.invokeMenu(app, 'Save Pipeline');

            tc.verifyTrue(isfile(fullfile(outDir, 'saved_pipeline.mat')), ...
                'the pipeline was not written where the dialog said');
            tc.verifyFalse(app.pipelineDirty, ...
                'the pipeline is on disk but still marked unsaved');
            tc.verifyEqual(app.pipelineName, 'saved_pipeline', ...
                'the app should adopt the name it was saved under');
        end

        function cancellingTheSaveDialogLeavesThePipelineUnsaved(tc)
        % The branch that must NOT clear the flag, and the reason the fix is
        % about reading the dialog's answer rather than about clearing the flag
        % unconditionally. uiputfile returns 0 for a cancel.
            isolatePrefs(tc, {'lastPipelineFolder', 'recentPipelines'});
            shadowFunction(tc, 'uiputfile', {0, 0});

            app = launchApp(tc);
            tc.makePipelineDirty(app);

            tc.invokeMenu(app, 'Save Pipeline');

            tc.verifyTrue(app.pipelineDirty, ...
                'a cancelled save cleared the unsaved marker');
        end
    end

    % ── driving the real app ──────────────────────────────────────────────────
    methods (Access = private)

        function names = stepLeafNames(~, app)
        % The step names the picker actually offers. The tree is stage-grouped,
        % so the grouping nodes have to be dropped: a step leaf is the one that
        % carries NodeData, which is what StepsTreeSelectionChanged reads.
            nodes = findall(app.StepsTree, 'Type', 'uitreenode');
            keep  = arrayfun(@(n) ~isempty(n.NodeData), nodes);
            names = arrayfun(@(n) string(n.NodeData), nodes(keep));
            names = cellstr(names(:)');
        end

        function cleanup = hideEeglabEntirely(tc)
        % Every path folder under the EEGLAB root, removed. hideFromPath takes
        % one function name and EEGLAB spreads over ~48 folders whose plugin
        % functions the steps probe individually, so hiding 'eeglab' alone
        % leaves most of them resolving and the cold start unreproduced.
            root = fileparts(which('eeglab'));
            tc.assertNotEmpty(root, 'EEGLAB is not on the path to begin with');
            saved = path();
            cleanup = onCleanup(@() restorePath(saved));

            folders = strsplit(path, pathsep);
            under   = folders(startsWith(lower(folders), lower(root)));
            cellfun(@rmpath, under);
            rehash;
            tc.assertEmpty(which('pop_saveset'), 'could not hide EEGLAB');
        end

        function makePipelineDirty(tc, app)
        % Through the real mutation path rather than by writing the flag: the
        % test should not be able to pass against an app that never sets it.
            tree = app.StepsTree;
            nodes = findall(tree, 'Type', 'uitreenode');
            leaf  = nodes(arrayfun(@(n) ~isempty(n.NodeData), nodes));
            tc.assertNotEmpty(leaf, 'no step to add, so the pipeline cannot be made dirty');
            tree.SelectedNodes = leaf(1);
            invokeCallback(tc, tree, 'SelectionChangedFcn');

            addBtn = findall(app.UIFigure, 'Type', 'uibutton', 'Text', 'Add');
            tc.assertNotEmpty(addBtn, 'no Add button');
            invokeCallback(tc, addBtn(1), 'ButtonPushedFcn');
        end

        function invokeMenu(tc, app, label)
            m = findall(app.UIFigure, 'Type', 'uimenu', 'Text', label);
            tc.assertNotEmpty(m, sprintf('no "%s" menu item', label));
            invokeCallback(tc, m(1), 'MenuSelectedFcn');
        end
    end
end

% ── local helpers ─────────────────────────────────────────────────────────────

function invokeCallback(tc, h, prop)
% Fire a component's real callback. The app's methods are all private, so this
% is the only route in - and it is the honest one: it is what the click does.
cb = h.(prop);
tc.assertNotEmpty(cb, sprintf('%s has no %s', class(h), prop));
cb(h, []);
drawnow;
end

function restorePath(saved)
path(saved);
rehash;
end
