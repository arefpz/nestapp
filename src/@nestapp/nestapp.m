% WARNING: Do not open nestapp_designer.mlapp and save - App Designer will
% regenerate this file and overwrite startupFcn and other hand-edited methods.
% All edits must be made directly to nestapp.m.

% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef nestapp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        CleaningTab                     matlab.ui.container.Tab
        ReStartStepsButton              matlab.ui.control.Button
        NESTAPPLabel                    matlab.ui.control.Label
        Image                           matlab.ui.control.Image
        RunAnalysisButton               matlab.ui.control.Button
        SelectDatatoPerformAnalysisPanel  matlab.ui.container.Panel
        SelectDataButton                matlab.ui.control.Button
        SelectedFilesListBox            matlab.ui.control.ListBox
        SelectedListBoxLabel_2          matlab.ui.control.Label
        TextArea                        matlab.ui.control.TextArea
        DefaultValueButton              matlab.ui.control.Button
        UITable                         matlab.ui.control.Table
        SelectedListBoxLabel            matlab.ui.control.Label
        RemoveButton                    matlab.ui.control.Button
        AddButton                       matlab.ui.control.Button
        MoveDownButton                  matlab.ui.control.Button
        MoveUpButton                    matlab.ui.control.Button
        SelectedListBox                 matlab.ui.control.ListBox
        StepsListBoxLabel               matlab.ui.control.Label
        InfoTextArea                    matlab.ui.control.TextArea
        CommandDescriptionLabel         matlab.ui.control.Label
        StepsTree                       matlab.ui.container.Tree
        StepsTipPanel                   matlab.ui.container.Panel
        StepsTipLabel                   matlab.ui.control.Label
        UIAxes                          matlab.ui.control.UIAxes
        ReportsTab                      matlab.ui.container.Tab
        ReportsListBox                  matlab.ui.control.ListBox
        ReportsListBoxLabel             matlab.ui.control.Label
        LoadReportsButton               matlab.ui.control.Button
        ClearReportsButton              matlab.ui.control.Button
        ReportsFolderLabel              matlab.ui.control.Label
        ReportsStatusLabel              matlab.ui.control.Label
        ReportsTextArea                 matlab.ui.control.TextArea
        ReportsDashboardPanel           matlab.ui.container.Panel
        ReportsImagePanel               matlab.ui.container.Panel
        OpenReportSetButton             matlab.ui.control.Button
        ReportsViewGroup                matlab.ui.container.ButtonGroup
        ReportsTextViewButton           matlab.ui.control.ToggleButton
        ReportsImageViewButton          matlab.ui.control.ToggleButton
        ExportReportsCSVButton          matlab.ui.control.Button
        ExportPDFButton                 matlab.ui.control.Button
        CopyMethodsButton               matlab.ui.control.Button
    end

    properties (Access = private)
        ItemNum % Index for selected Item
        path         % File Path (single-folder selection; '' when files span folders)
        file         % File Name(s) shown in the listbox (basenames, or folder/name when multi-folder)
        filePaths    % Full path of every queued data file - source of truth for Run Analysis
        spec         % PipelineStep struct array (name + typed params)
        NSelecFiles  % Number of selected files for EEG preprocessing
    end
    properties (Access = public)
        % Tab Cleaning
        selectedItem % Selected Table Item Values
        info % containers.Map: step name -> info/description text (shown in InfoTextArea)
        hoverTimer    % singleShot timer: fires once the pointer has rested on the Steps tree
        treeRect      % cached figure-relative rect of StepsTree (getpixelposition is slow)
        treeRectStamp % UIFigure.Position the cached rect was computed for
        stepBlocks    % containers.Map: step name -> true if it ALWAYS waits for a human
        stepProviders % containers.Map: step name -> provider (TESA, EEGLAB, ...)
        % Canonical pipeline state - single source of truth for steps and params.
        % appendStep/removeStep/moveStep/clearSteps/loadPipelineData all write here.
        currentParamKey  = ''  % param key selected in UITable (transient)
        currentParamType = ''  % type of selected param (transient)
        originalSize     % [w h] of UIFigure at creation - used by UIFigureSizeChanged
        isResizing = false % guards UIFigureSizeChanged against its own clamp write
        baseLayout = []  % per-component base geometry captured at startup - drives rescaleComponents

        % Tab Explore - the grouped-comparison workspace
        ExploreTab                      matlab.ui.container.Tab
        ExploreGroupsLabel              matlab.ui.control.Label
        ExploreGroupsListBox            matlab.ui.control.ListBox
        ExploreAddGroupButton           matlab.ui.control.Button
        ExploreRemoveGroupButton        matlab.ui.control.Button
        ExploreRoiLabel                 matlab.ui.control.Label
        ExploreRoiDropDown              matlab.ui.control.DropDown
        ExploreRoiEditButton            matlab.ui.control.Button
        ExploreRoiSummaryLabel          matlab.ui.control.Label
        ExploreWindowsLabel             matlab.ui.control.Label
        ExploreWindowsTable             matlab.ui.control.Table
        ExploreFilesButton              matlab.ui.control.Button
        ExploreDesignLabel              matlab.ui.control.Label
        ExploreDesignGroup              matlab.ui.container.ButtonGroup
        ExploreUnpairedButton           matlab.ui.control.RadioButton
        ExplorePairedButton             matlab.ui.control.RadioButton
        ExploreDesignNoteLabel          matlab.ui.control.Label
        ExploreWindowsModeDropDown      matlab.ui.control.DropDown
        ExploreWindowsAddButton         matlab.ui.control.Button
        ExploreWindowsRemoveButton      matlab.ui.control.Button
        ExploreWindowsResetButton       matlab.ui.control.Button
        ExplorePlotLabel                matlab.ui.control.Label
        ExplorePlotDropDown             matlab.ui.control.DropDown
        ExplorePlotOptionsButton        matlab.ui.control.Button
        ExplorePlotInfoLabel            matlab.ui.control.Label
        ExploreCanvas                   matlab.ui.container.Panel
        ExploreEmptyLabel               matlab.ui.control.Label
        ExploreFigureButton             matlab.ui.control.Button
        ExploreCsvButton                matlab.ui.control.Button
        ExploreResultsButton            matlab.ui.control.Button
        ExploreStatusLabel              matlab.ui.control.Label

        % Menus
        MenuRecentFiles     % Handle to 'Recent Files' submenu - rebuilt on open
        MenuRecentPipelines % Handle to 'Recent Pipelines' submenu - rebuilt on open
        StatusBar           % uilabel pinned to bottom of UIFigure - visible on both tabs
        pipelineDirty   = false    % true when pipeline has unsaved changes
        pipelineName    = ''       % filename of last saved/loaded pipeline
        allPipelineReports = {}    % cell array of report entry structs from current session
        loadedReports      = {}    % cell array of report entry structs loaded from disk
        lastFailed         = struct([]) % failure records from the most recent run (for the dashboard + summary)
        preSelectedChanFile = ''   % channel location file selected once before a run
        ParallelCheckBox           % uicheckbox - enable parallel participant processing

        % Explore state. entries and cache are the dataset; everything else is
        % derived on demand by groupCurves, which is cheap enough to re-run on
        % every change (measured ~30 ms for 35 files), so nothing here caches a
        % result that could go stale against the controls.
        exploreEntries = struct('path', {}, 'subject', {}, 'group', {}, ...
                                'subjectConfident', {})
        exploreCache   = struct('path', {}, 'trialAvg', {}, 'labels', {}, ...
                                'chanlocs', {}, 'time', {}, 'nTrials', {}, 'ok', {})
        exploreRoi     = {}      % ROI electrode labels (canonical spelling)
        exploreSubjectRule = []  % last filename rule applied, to reopen the dialog on it
        exploreWindows = struct([])
        exploreRes     = struct([])   % last groupCurves result, for the exits
        % Per-plot settings, as a name/params struct array rather than a
        % struct keyed by plot name - "TEP (ROI mean)" is not a valid field
        % name, and renaming plots to suit the storage would be backwards.
        % Only params the user actually SET are held here; the rest stay with
        % the draw function's own defaults.
        explorePlotParams = struct('name', {}, 'params', {})
        exploreFigureOpts = struct()   % remembered publicationFigure settings
        exploreResizeTimer             % coalesces a resize drag into one repaint
        reportsResizeTimer             % same, for the Reports dashboard / QC images
        reportsQcIndex = 1             % which QC checkpoint the image pane shows
        exploreAvailablePlots = struct([])  % registry entries + availability

        % Tab Analysis
    end

    methods (Access = private)
        % -- Pipeline state mutation methods ---------------------------------
        % app.spec, SelectedListBox.Items, and SelectedListBox.ItemsData must
        % stay in sync. These methods are the ONLY permitted way to add,
        % remove, move, or clear steps -- callbacks delegate here.

        function name = selectedStepName(app)
        % SELECTEDSTEPNAME  Registry step name of the selected tree node, or ''
        % when a category/operation header (empty NodeData) is selected.
            name = '';
            n = app.StepsTree.SelectedNodes;
            if isempty(n); return; end
            d = n(1).NodeData;
            if ischar(d) || isstring(d); name = char(d); end
        end

        function UIFigureMouseMoved(app, ~)
        % Any pointer movement hides the tip and restarts its clock, so it
        % only ever appears once the pointer has been still on the tree for
        % HOVER_DELAY seconds - and gets out of the way the instant you move.
        % Fires on every mouse move over the window, so it stays cheap: no
        % graphics write unless the tip is actually up.
            if isempty(app.hoverTimer) || ~isvalid(app.hoverTimer); return; end
            stop(app.hoverTimer);
            hideStepsTip(app);
            if ~pointerOverStepsTree(app); return; end
            start(app.hoverTimer);
        end

        function tf = pointerOverStepsTree(app)
        % True when the pointer is inside the Steps tree, which only exists
        % while the Cleaning tab is showing.
            tf = false;
            if app.TabGroup.SelectedTab ~= app.CleaningTab; return; end
            r = stepsTreeRect(app);
            p = app.UIFigure.CurrentPoint;
            tf = p(1) >= r(1) && p(1) <= r(1) + r(3) ...
              && p(2) >= r(2) && p(2) <= r(2) + r(4);
        end

        function r = stepsTreeRect(app)
        % Figure-relative rect of the tree. getpixelposition costs ~150 us, far
        % too much for a callback that fires on every mouse move, so cache it
        % and recompute only when the window has changed size or moved.
            pos = app.UIFigure.Position;
            if isempty(app.treeRect) || ~isequal(app.treeRectStamp, pos)
                app.treeRect      = getpixelposition(app.StepsTree, true);
                app.treeRectStamp = pos;
            end
            r = app.treeRect;
        end

        function showStepsTip(app)
        % Fired by the dwell timer. Places the tip near the pointer, nudged so
        % it never sits under the cursor and never leaves the window.
            if ~isvalid(app.UIFigure) || ~pointerOverStepsTree(app); return; end
            GAP = 14;
            tip = app.StepsTipPanel.Position(3:4);
            fig = app.UIFigure.Position(3:4);
            p   = app.UIFigure.CurrentPoint;
            x   = min(max(1, p(1) + GAP), fig(1) - tip(1) - 1);
            y   = min(max(1, p(2) - tip(2) - GAP), fig(2) - tip(2) - 1);
            app.StepsTipPanel.Position = [x, y, tip];
            app.StepsTipPanel.Visible = 'on';
        end

        function hideStepsTip(app)
        % Guarded on the current state: this runs on every mouse move, and
        % re-asserting 'off' on an already-hidden panel is a graphics write
        % for nothing.
            if ~isempty(app.StepsTipPanel) && isvalid(app.StepsTipPanel) ...
                    && app.StepsTipPanel.Visible
                app.StepsTipPanel.Visible = 'off';
            end
        end

        function txt = stepsTreeLegend(app)
        % The text of the tree's hover tip, and the one place the design is
        % explained.
        %
        % uitreenode has no Tooltip property, so a per-node hover tip is not
        % possible - this tree-level key carries the dot's meaning instead, and
        % the Info panel below names the provider and the exact wait behaviour
        % for whichever step is selected.
        %
        % It is shown by a dwell timer rather than the native Tooltip: the
        % native one fires on MATLAB's own schedule and cannot be delayed, and
        % assigning Tooltip mid-hover does not reliably re-trigger it.
            n = 0;
            if ~isempty(app.stepBlocks); n = app.stepBlocks.Count; end
            txt = sprintf([ ...
                'Amber dot: the step waits for you - it opens a window you ' ...
                'must close, and cannot run with Parallel Processing on ' ...
                '(%d steps).\n\nSelect a step for details.'], n);
        end

        function lines = stepInfoLines(app, name)
        % What the Info panel shows for one step: its description, then the
        % facts the tree can only hint at - who supplies it (every step has a
        % provider, so that is text rather than a dot) and, for the few that
        % carry the amber dot, what the dot means for this step.
            lines = string.empty(0, 1);
            % The dot's meaning goes first: the Info panel is small, and a
            % note under a long description is a note nobody scrolls to.
            if ~isempty(app.stepBlocks) && isKey(app.stepBlocks, name)
                if app.stepBlocks(name)
                    lines(end+1) = "* Waits for you - opens a window you must " + ...
                        "close before the run continues.";
                else
                    lines(end+1) = "* Waits for you in some modes - opens a " + ...
                        "window you must close when its review option is on.";
                end
                lines(end+1) = "  Cannot run with Parallel Processing on.";
                lines(end+1) = "";
            end
            lines(end+1) = string(app.info(name));
            if ~isempty(app.stepProviders) && isKey(app.stepProviders, name)
                lines(end+1) = "";
                lines(end+1) = "Provided by: " + string(app.stepProviders(name));
            end
        end

        function populateStepsTree(app)
        % POPULATESTEPSTREE  (Re)build the stage-grouped picker tree and the
        % name->info map from the steps this machine can run (availableSteps).
        % Leaf NodeData is the exact registry step name. Any available step the
        % taxonomy forgot is collected under "Other" so nothing silently
        % vanishes from the picker.
            delete(app.StepsTree.Children);
            reg             = stepRegistry();
            [steps, hidden] = availableSteps(reg);
            availNames      = {steps.name};
            % A withheld step leaves no trace in the tree, so "where did my
            % steps go" has to be answerable somewhere. Nearly always a missing
            % or too-old plugin - or EEGLAB not initialised, which used to hide
            % most of the registry here (see startupFcn).
            if ~isempty(hidden)
                nestLog('CFG', 'Step picker: %d of %d steps not offered: %s', ...
                    numel(hidden), numel(reg), strjoin({hidden.name}, ', '));
            end
            app.info        = containers.Map(availNames, {steps.info});

            % uitree node styling (uistyle/addStyle on a tree) is R2023b+; the
            % app's floor is R2023a, so probe once - the picker still builds on
            % older releases, just without bold headers / muted operation rows.
            canStyle = treeStylingSupported(app, app.StepsTree);

            TAX = stepTaxonomy();
            % Steps that can stop and wait for a human get the amber dot. The
            % set is derived from the registry via canStepBlock, so adding an
            % interactive step marks itself in the picker.
            % TODO: runPipelineCore's local findInteractiveSteps is still a
            % hardcoded name list and does not agree with the registry - it is
            % what actually gates the parallel-processing warning.
            app.stepBlocks    = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            app.stepProviders = containers.Map('KeyType', 'char', 'ValueType', 'char');
            for r = 1:numel(reg)
                if ~ismember(reg(r).name, availNames); continue; end
                [tf, always] = canStepBlock(reg(r));
                if tf; app.stepBlocks(reg(r).name) = always; end
            end
            flag = stepInteractiveIcon();
            sCat = uistyle('FontWeight', 'bold');
            sOp  = uistyle('FontColor', [0.42 0.47 0.53]);

            placed = {};
            for c = 1:numel(TAX)
                cat = TAX(c);
                % keep only the operations/variants available on this machine
                ops = struct('name', {}, 'variants', {});
                cnt = 0;
                for o = 1:numel(cat.ops)
                    v = cat.ops(o).variants;
                    v = v(ismember({v.step}, availNames));
                    if isempty(v); continue; end
                    cnt = cnt + 1;
                    ops(cnt).name     = cat.ops(o).name;
                    ops(cnt).variants = v;
                end
                if cnt == 0; continue; end
                nShown = sum(arrayfun(@(o) numel(o.variants), ops));
                cNode = uitreenode(app.StepsTree, ...
                    'Text', sprintf('%s   (%d)', upper(cat.name), nShown), 'NodeData', '');
                if canStyle; addStyle(app.StepsTree, sCat, 'node', cNode); end
                for o = 1:numel(ops)
                    if numel(ops(o).variants) == 1
                        addStepLeaf(app, cNode, ops(o).variants(1), flag);
                    else
                        oNode = uitreenode(cNode, 'Text', ops(o).name, 'NodeData', '');
                        if canStyle; addStyle(app.StepsTree, sOp, 'node', oNode); end
                        for k = 1:numel(ops(o).variants)
                            addStepLeaf(app, oNode, ops(o).variants(k), flag);
                        end
                    end
                    placed = [placed, {ops(o).variants.step}]; %#ok<AGROW>
                end
            end

            % Safety net: surface any available step the taxonomy did not place.
            orphan = setdiff(availNames, placed);
            if ~isempty(orphan)
                oNode = uitreenode(app.StepsTree, ...
                    'Text', sprintf('OTHER   (%d)', numel(orphan)), 'NodeData', '');
                if canStyle; addStyle(app.StepsTree, sCat, 'node', oNode); end
                for i = 1:numel(orphan)
                    uitreenode(oNode, 'Text', orphan{i}, 'NodeData', orphan{i});
                end
            end

            expand(app.StepsTree);
            app.StepsTipLabel.Text = stepsTreeLegend(app);
        end

        function addStepLeaf(app, parent, v, flagIcon)
        % ADDSTEPLEAF  One leaf whose NodeData is the registry step name; an
        % amber dot marks a step that can stop and wait for a human.
            node = uitreenode(parent, 'Text', v.step, 'NodeData', v.step);
            app.stepProviders(v.step) = v.provider;
            if isKey(app.stepBlocks, v.step)
                % Set the flag after creation and tolerate its absence: the leaf
                % must always appear even on a release without node icons.
                try, node.Icon = flagIcon; catch, end %#ok<CTCH>
            end
        end

        function tf = treeStylingSupported(~, tree)
        % TREESTYLINGSUPPORTED  True when uistyle/addStyle work on a uitree node
        % (R2023b+). Probed at runtime rather than gated on a version string, so
        % the picker degrades gracefully without depending on release numbers.
            probe = uitreenode(tree, 'Text', '', 'NodeData', '');
            try
                addStyle(tree, uistyle('FontWeight', 'bold'), 'node', probe);
                tf = true;
            catch
                tf = false;
            end
            delete(probe);
        end

        function appendStep(app, stepName)
        % APPENDSTEP  Append stepName to the pipeline using its default params.
            reg    = stepRegistry();
            regIdx = find(strcmp({reg.name}, stepName), 1);
            if isempty(regIdx); return; end
            n = numel(app.SelectedListBox.Items);
            % Treat a single empty-string sentinel as an empty list
            if n == 1 && isempty(app.SelectedListBox.Items{1})
                n = 0;
                app.SelectedListBox.Items(:)     = [];
                app.SelectedListBox.ItemsData(:) = [];
                app.spec = repmat(struct('name','','params',struct()), 0, 1);
            end
            pos = n + 1;
            app.SelectedListBox.Items{pos}     = stepName;
            app.SelectedListBox.ItemsData{pos} = ['Item' num2str(pos)];
            app.spec(pos)                      = makePipelineStep(stepName, reg);
            app.pipelineDirty = true;
            updateStatusBar(app);
        end

        function removeStep(app, idx)
        % REMOVESTEP  Remove the step at index idx and renumber ItemsData.
            app.SelectedListBox.Items(idx)     = [];
            app.SelectedListBox.ItemsData(idx) = [];
            app.spec(idx)                      = [];
            for i = idx : numel(app.SelectedListBox.ItemsData)
                app.SelectedListBox.ItemsData{i} = ['Item' num2str(i)];
            end
            app.pipelineDirty = true;
            updateStatusBar(app);
        end

        function moveStep(app, idx, direction)
        % MOVESTEP  Swap step at idx with its neighbour in the given direction
        %   (+1 = move down, -1 = move up). No-op at boundaries.
            n    = numel(app.SelectedListBox.Items);
            idx2 = idx + direction;
            if idx2 < 1 || idx2 > n; return; end
            % Swap step names and spec entries
            [app.SelectedListBox.Items{idx}, app.SelectedListBox.Items{idx2}] = ...
                deal(app.SelectedListBox.Items{idx2}, app.SelectedListBox.Items{idx});
            [app.spec(idx), app.spec(idx2)] = deal(app.spec(idx2), app.spec(idx));
            % ItemsData stays in positional order - just renumber both slots
            app.SelectedListBox.ItemsData{idx}  = ['Item' num2str(idx)];
            app.SelectedListBox.ItemsData{idx2} = ['Item' num2str(idx2)];
            app.SelectedListBox.Value = app.SelectedListBox.ItemsData{idx2};
            app.pipelineDirty = true;
            updateStatusBar(app);
        end

        function clearSteps(app)
        % CLEARSTEPS  Remove all pipeline steps and reset state.
            app.spec = repmat(struct('name','','params',struct()), 0, 1);
            app.SelectedListBox.Items(:)     = [];
            app.SelectedListBox.ItemsData(:) = [];
            app.UITable.Data    = [];
            app.ItemNum          = 0;
            app.currentParamKey  = '';
            app.currentParamType = '';
            app.pipelineDirty   = true;
            updateStatusBar(app);
        end

        function idx = selectedStepIndex(app)
        % SELECTEDSTEPINDEX  Decode the current SelectedListBox selection to a 1-based index.
            idx = str2double(strrep(app.SelectedListBox.Value, 'Item', ''));
        end

        function refreshParamTable(app, stepIdx)
        % REFRESHPARAMTABLE  Update UITable from app.spec(stepIdx).
            reg    = stepRegistry();
            step   = app.spec(stepIdx);
            regIdx = find(strcmp({reg.name}, step.name), 1);
            if isempty(regIdx)
                app.UITable.Data = [];
                return
            end
            data = buildParamTableData(step, reg(regIdx));
            % Show params that are disabled for this step - overridden by a
            % mutually-exclusive sibling, or gated off by another param's value
            % - greyed out. The '(' prefix makes styleParamTable grey the cell.
            disabled = disabledParamKeys(reg(regIdx), step.params);
            if ~isempty(disabled)
                keys = {reg(regIdx).params.key};
                for k = 1:numel(disabled)
                    r = find(strcmp(keys, disabled{k}), 1);
                    if ~isempty(r)
                        data{r,2} = '(disabled)';
                    end
                end
            end
            app.UITable.Data = data;
            styleParamTable(app);
        end

        % -----------------------------------------------------------------

        function updateStatusBar(app)
        % UPDATESTATUSBAR  Refresh the status bar text from current app state.
        %   Called after any change to the pipeline list, data selection,
        %   or save/load operations.
            % Pipeline segment
            items = app.SelectedListBox.Items;
            nSteps = numel(items);
            if nSteps == 0 || (nSteps == 1 && isempty(items{1}))
                pipelineStr = 'Pipeline: (empty)';
            else
                name = app.pipelineName;
                if isempty(name); name = 'unsaved'; end
                if app.pipelineDirty
                    pipelineStr = sprintf('Pipeline: %s*  (%d steps)', name, nSteps);
                else
                    pipelineStr = sprintf('Pipeline: %s  (%d steps)', name, nSteps);
                end
            end
            % Data segment
            n = app.NSelecFiles;
            if isempty(n); n = 0; end
            fileWord = 'files'; if n == 1; fileWord = 'file'; end
            if n == 0
                dataStr = 'Data: (none)';
            elseif (ischar(app.path) || isstring(app.path)) && ~isempty(app.path)
                % Single-folder selection - show the folder name.
                parts = strsplit(strtrim(char(app.path)), {'\','/'});
                parts(cellfun(@isempty, parts)) = [];
                folder = parts{end};
                dataStr = sprintf('Data: %s/  (%d %s)', folder, n, fileWord);
            else
                % Multi-folder selection - count the distinct parent folders.
                parents  = cellfun(@fileparts, app.filePaths, 'UniformOutput', false);
                nFolders = numel(unique(parents));
                folderWord = 'folders'; if nFolders == 1; folderWord = 'folder'; end
                dataStr = sprintf('Data: %d %s in %d %s', n, fileWord, nFolders, folderWord);
            end
            app.StatusBar.Text = sprintf('  %s          %s', pipelineStr, dataStr);
        end

        % MenuOpening function: MenuRecentFiles
        function buildRecentFilesMenu(app)
            delete(app.MenuRecentFiles.Children);
            list = getpref('nestapp', 'recentFiles', {});
            if isempty(list)
                uimenu(app.MenuRecentFiles, 'Text', '(none)', 'Enable', 'off');
                return
            end
            for i = 1:numel(list)
                folder = list{i};
                uimenu(app.MenuRecentFiles, 'Text', folder, ...
                    'MenuSelectedFcn', @(~,~) openRecentData(app, folder));
            end
        end

        function buildRecentPipelinesMenu(app)
            delete(app.MenuRecentPipelines.Children);
            list = getpref('nestapp', 'recentPipelines', {});
            if isempty(list)
                uimenu(app.MenuRecentPipelines, 'Text', '(none)', 'Enable', 'off');
                return
            end
            for i = 1:numel(list)
                pPath = list{i};
                [~,nm,ex] = fileparts(pPath);
                uimenu(app.MenuRecentPipelines, 'Text', [nm ex], ...
                    'MenuSelectedFcn', @(~,~) openRecentPipeline(app, pPath));
            end
        end

        function openRecentData(app, folder)
        % Open the data browser straight into a recently used folder.
            if ~isfolder(folder)
                uialert(app.UIFigure, 'Folder no longer exists.', 'Not Found');
                return
            end
            paths = selectDataTree(folder, app.dataFileExts());
            if isempty(paths); return; end
            setFileQueue(app, paths);
            if isvalid(app.UIFigure); figure(app.UIFigure); end
        end

        function openRecentPipeline(app, pPath)
        % Load a pipeline from a recently used full file path.
            if ~isfile(pPath)
                uialert(app.UIFigure, 'Pipeline file no longer exists.', 'Not Found');
                return
            end
            [pFolder, ~, ~] = fileparts(pPath);
            try
                loadPipelineData(app, pPath);
                setpref('nestapp', 'lastPipelineFolder', pFolder);
                pushRecent(app, 'recentPipelines', pPath);
                buildRecentPipelinesMenu(app);
            catch err
                uialert(app.UIFigure, err.message, 'Load Error', 'Icon', 'error');
            end
        end

        function loadPipelineData(app, fullPath)
        % LOADPIPELINEDATA  Load pipeline state from a .mat into app.spec.
        %   Unknown steps produce a warning dialog.
            data = load(fullPath, '-mat');
            reg  = stepRegistry();
            [app.spec, warns] = specFromSaved(data, reg);

            if ~isempty(warns)
                uialert(app.UIFigure, strjoin(warns, newline), ...
                    'Pipeline Warning', 'Icon', 'warning');
            end

            n = numel(app.spec);
            items     = cell(1, n);
            itemsData = cell(1, n);
            for k = 1:n
                items{k}     = app.spec(k).name;
                itemsData{k} = ['Item' num2str(k)];
            end
            app.SelectedListBox.Items     = items;
            app.SelectedListBox.ItemsData = itemsData;

            app.currentParamKey  = '';
            app.currentParamType = '';
            app.UITable.Data     = [];

            if n > 0
                app.SelectedListBox.Value = itemsData{1};
                refreshParamTable(app, 1);
            end
        end

        % MenuSelected callback wrappers (thin shims so uimenu can call private methods)
        function openPreferencesMenu(app, ~)
            openPreferences(app);
        end

        function showAboutMenu(app, ~)
            showAbout(app);
        end

        function copyDiagnosticsMenu(app, ~)
        % COPYDIAGNOSTICSMENU  Help-menu action: copy environment diagnostics.
        %   Runs nestappDoctor, copies the Markdown report to the clipboard,
        %   and shows the problems summary so the user can paste it into a
        %   bug report (see .github/ISSUE_TEMPLATE/bug_report.yml).
            try
                [~, diagInfo] = nestappDoctor('Copy', true, 'Quiet', true);
            catch ME
                uialert(app.UIFigure, ...
                    sprintf('Could not collect diagnostics:\n%s', ME.message), ...
                    'Diagnostics Failed', 'Icon', 'error');
                return
            end
            if isempty(diagInfo.problems)
                msg  = 'Diagnostics copied to the clipboard. No problems detected.';
                icon = 'success';
            else
                msg  = sprintf(['Diagnostics copied to the clipboard.\n\n' ...
                    '%d problem(s) detected:\n  - %s'], ...
                    numel(diagInfo.problems), strjoin(diagInfo.problems, sprintf('\n  - ')));
                icon = 'warning';
            end
            uialert(app.UIFigure, msg, 'nestapp Diagnostics', 'Icon', icon);
        end

        function copyPipelineDescriptionMenu(app, ~)
        % COPYPIPELINEDESCRIPTIONMENU  File-menu action: copy a readable
        %   description of the current pipeline (steps + customised params)
        %   to the clipboard, for methods sections and bug reports.
            if isempty(app.spec)
                uialert(app.UIFigure, 'The pipeline is empty - add steps first.', ...
                    'No Pipeline', 'Icon', 'warning');
                return
            end
            try
                describePipeline(app.spec, 'Copy', true, 'Quiet', true);
            catch ME
                uialert(app.UIFigure, ...
                    sprintf('Could not describe the pipeline:\n%s', ME.message), ...
                    'Export Failed', 'Icon', 'error');
                return
            end
            uialert(app.UIFigure, sprintf( ...
                'Pipeline description (%d steps) copied to the clipboard.', ...
                numel(app.spec)), 'Pipeline Copied', 'Icon', 'success');
        end

        function revealFolder(~, folder)
        % REVEALFOLDER  Open a folder in the OS file browser (best-effort).
            try
                if ispc
                    winopen(folder);
                elseif ismac
                    system(sprintf('open "%s" &', folder));
                else
                    system(sprintf('xdg-open "%s" &', folder));
                end
            catch
                % Non-fatal: the path is shown in the dialog regardless.
            end
        end

        function collectSupportBundleMenu(app, ~)
        % COLLECTSUPPORTBUNDLEMENU  Help action: write a metadata-only support
        %   bundle (environment + current pipeline) and reveal the folder.
            outRoot = getpref('nestapp', 'outputRoot', '');
            if isempty(outRoot) || ~isfolder(outRoot)
                outRoot = tempdir;
            end
            try
                bundleDir = collectSupportBundle(outRoot, app.spec);
            catch ME
                uialert(app.UIFigure, ...
                    sprintf('Could not collect support bundle:\n%s', ME.message), ...
                    'Support Bundle Failed', 'Icon', 'error');
                return
            end
            revealFolder(app, bundleDir);
            uialert(app.UIFigure, sprintf(['Support bundle written to:\n%s\n\n' ...
                'It contains environment + pipeline details only (no recordings). ' ...
                'Attach the folder to your bug report.'], bundleDir), ...
                'Support Bundle', 'Icon', 'success');
        end

        function installAaratepMenu(app, ~)
        % INSTALLAARATEPMENU  Help action: download the AARATEP helper functions.
        %   The AARATEP template needs a ~300-file tree that nestapp may not
        %   redistribute, and the documented alternative is a git clone from a
        %   terminal - which assumes git, a shell, and knowing where to put the
        %   result. This does it in about a second with no tools installed.
            rel = aaratepRelease();

            dlg = uiprogressdlg(app.UIFigure, 'Title', 'AARATEP', ...
                'Message', 'Starting...', 'Indeterminate', 'off');
            closeDlg = onCleanup(@() close(dlg));
            % dlg is a handle, so the closure updates the live dialog.
            onProgress = @(frac, msg) set(dlg, 'Value', frac, 'Message', msg);

            result = installAaratep('Progress', onProgress);
            clear closeDlg   % the alerts below must not sit behind the dialog

            if ~result.installed
                uialert(app.UIFigure, result.message, 'AARATEP Install Failed', ...
                    'Icon', 'error');
                return
            end

            % A newer upstream release is reported, never installed - the pin
            % is what keeps two people running one template on one pipeline.
            % See aaratepRelease.
            extra = '';
            if ~isempty(result.newerTag)
                extra = sprintf(['\n\nNote: AARATEP %s has since been released. ' ...
                    'nestapp pins %s so that a template produces the same result ' ...
                    'for everyone; upgrading is a deliberate change.'], ...
                    result.newerTag, rel.tag);
            end

            uialert(app.UIFigure, [result.message extra], 'AARATEP', ...
                'Icon', 'success');
        end

        function selfTestMenu(app, ~)
        % SELFTESTMENU  Help action: run the fast test suite to verify the
        %   install, reporting pass/fail. Best-effort: needs tests/ present.
            runner = fullfile(nestappRoot(), 'tests', 'run_tests.m');
            if ~isfile(runner)
                uialert(app.UIFigure, ['The test suite (tests/) is not present ' ...
                    'in this installation, so the self-test cannot run.'], ...
                    'Self-test Unavailable', 'Icon', 'warning');
                return
            end
            dlg = uiprogressdlg(app.UIFigure, 'Title', 'Checking install', ...
                'Message', 'Running the fast test suite...', 'Indeterminate', 'on');
            closeDlg = onCleanup(@() close(dlg));
            try
                addpath(fullfile(repo, 'tests'));
                addpath(fullfile(repo, 'tests', 'helpers'));
                results = [];
                evalc('results = run_tests(''fast'')');   % capture verbose output
            catch ME
                clear closeDlg;
                uialert(app.UIFigure, sprintf('Self-test could not run:\n%s', ...
                    ME.message), 'Self-test Error', 'Icon', 'error');
                return
            end
            clear closeDlg;
            nPass = sum([results.Passed]);
            nFail = sum([results.Failed]);
            nInc  = sum([results.Incomplete]);
            if nFail == 0
                uialert(app.UIFigure, sprintf(['Install looks healthy.\n' ...
                    '%d passed, %d skipped (optional plugins).'], nPass, nInc), ...
                    'Self-test Passed', 'Icon', 'success');
            else
                uialert(app.UIFigure, sprintf(['%d test(s) FAILED (%d passed).\n' ...
                    'Run Help > Copy Diagnostics and check your setup.'], nFail, nPass), ...
                    'Self-test Failed', 'Icon', 'error');
            end
        end

        function loadPrefs(~)
        % LOADPREFS  Read persistent preferences and apply to app state.
        %   Called from startupFcn. Uses MATLAB getpref with 'nestapp' group.
        %   The app handle is accepted but not used - prefs apply globally
        %   (addpath) rather than writing to removed UI components.
            eeglabPath = getpref('nestapp', 'eeglabPath', '');
            if ~isempty(eeglabPath) && isfolder(eeglabPath)
                addpath(eeglabPath);
            end

            % One-shot migration: legacy 'reportFolder' pref folds into
            % the new unified 'outputRoot' pref.
            if ispref('nestapp', 'reportFolder')
                rf = getpref('nestapp', 'reportFolder', '');
                if ~ispref('nestapp', 'outputRoot') && ~isempty(rf) && isfolder(rf)
                    setpref('nestapp', 'outputRoot', rf);
                    nestLog('CFG', 'Migrated reportFolder pref -> outputRoot: %s', rf);
                end
                rmpref('nestapp', 'reportFolder');
            end
        end

        function pushRecent(app, prefKey, newEntry) %#ok<INUSL>
        % PUSHRECENT  Prepend newEntry to a 5-item MRU list stored in prefs.
            list = getpref('nestapp', prefKey, {});
            list = [{newEntry}, list(~strcmp(list, newEntry))];
            list = list(1:min(end, 5));
            setpref('nestapp', prefKey, list);
        end

        function openPreferences(~)
        % OPENPREFERENCES  Show a modal Preferences dialog.
        %   Lets users set the EEGLAB path, default data/pipeline folders,
        %   and behavioural options. Changes are written to getpref/setpref
        %   under the 'nestapp' group and applied immediately on Save.
            dlg = uifigure('Name', 'nestapp Preferences', ...
                'Position', [200 200 420 664], ...
                'WindowStyle', 'modal', 'Resize', 'off');

            % --- Quality Screening section (new, at top) ---
            uilabel(dlg, 'Text', 'Quality Screening', 'FontWeight', 'bold', ...
                'Position', [15 629 200 20]);
            cbAutoQC = uicheckbox(dlg, 'Text', 'Auto-generate QC images at each Quality Gate', ...
                'Position', [15 604 380 22], ...
                'Value', getpref('nestapp', 'autoQualityReport', false));
            cbTmsAuto = uicheckbox(dlg, 'Text', 'Auto-detect TMS pulse window from EEG events', ...
                'Position', [15 582 380 22], ...
                'Value', getpref('nestapp', 'qualityTmsAutoDetect', true));
            cbSkipFail = uicheckbox(dlg, 'Text', 'Skip remaining pipeline steps when Quality Gate fails', ...
                'Position', [15 560 380 22], ...
                'Value', getpref('nestapp', 'skipOnQualityFail', false));
            cbAutoPDF = uicheckbox(dlg, 'Text', 'Auto-save PDF report per file (text + checkpoint images)', ...
                'Position', [15 538 380 22], ...
                'Value', getpref('nestapp', 'autoExportPDF', false));
            cbAaratepKeep = uicheckbox(dlg, 'Text', 'Keep AARATEP intermediate datasets (~3x the result on disk)', ...
                'Position', [15 516 380 22], ...
                'Value', getpref('nestapp', 'aaratepKeepIntermediates', true));
            cbAaratepKeep.Tooltip = ['AARATEP writes the dataset before SOUND, before decay ' ...
                'removal and before ICA rejection, as well as the result - four full ' ...
                'datasets per file. Untick to delete the three intermediates once the ' ...
                'step finishes; the result and the QC images are always kept.'];
            uilabel(dlg, 'Text', 'Attribute mode:', ...
                'Position', [15 487 95 22], 'HorizontalAlignment', 'right');
            qcModes = qualityAttributeModes();
            ddAttr = uidropdown(dlg, ...
                'Position', [115 487 150 22], ...
                'Items', qcModes, ...
                'Value', resolveAttributePref());
            uilabel(dlg, 'Text', 'TMS window (ms):', ...
                'Position', [15 460 105 22], 'HorizontalAlignment', 'right');
            qcWin = readTmsWindowPref();
            nfTmsStart = uieditfield(dlg, 'numeric', ...
                'Position', [125 460 55 22], 'Value', qcWin(1));
            uilabel(dlg, 'Text', 'to', ...
                'Position', [185 460 15 22], 'HorizontalAlignment', 'center');
            nfTmsEnd = uieditfield(dlg, 'numeric', ...
                'Position', [205 460 55 22], 'Value', qcWin(2));

            % --- EEGLAB section ---
            uilabel(dlg, 'Text', 'EEGLAB', 'FontWeight', 'bold', ...
                'Position', [15 420 200 20]);
            uilabel(dlg, 'Text', 'Path:', ...
                'Position', [15 395 35 22], 'HorizontalAlignment', 'right');
            fEeglab = uieditfield(dlg, 'text', ...
                'Position', [55 395 275 22], 'Editable', 'on', ...
                'Value', getpref('nestapp','eeglabPath',''));
            uibutton(dlg, 'Text', 'Browse...', 'Position', [335 395 70 22], ...
                'ButtonPushedFcn', @(~,~) browseEeglab());

            % --- Default Locations section ---
            uilabel(dlg, 'Text', 'Default Locations', 'FontWeight', 'bold', ...
                'Position', [15 365 200 20]);
            uilabel(dlg, 'Text', 'Data folder:', ...
                'Position', [15 340 65 22], 'HorizontalAlignment', 'right');
            fData = uieditfield(dlg, 'text', ...
                'Position', [85 340 245 22], 'Editable', 'on', ...
                'Value', getpref('nestapp','lastDataFolder',''));
            uibutton(dlg, 'Text', 'Browse...', 'Position', [335 340 70 22], ...
                'ButtonPushedFcn', @(~,~) browseFolder(fData));
            uilabel(dlg, 'Text', 'Pipeline folder:', ...
                'Position', [15 312 80 22], 'HorizontalAlignment', 'right');
            fPipeline = uieditfield(dlg, 'text', ...
                'Position', [100 312 230 22], 'Editable', 'on', ...
                'Value', getpref('nestapp','lastPipelineFolder',''));
            uibutton(dlg, 'Text', 'Browse...', 'Position', [335 312 70 22], ...
                'ButtonPushedFcn', @(~,~) browseFolder(fPipeline));
            uilabel(dlg, 'Text', 'Output root:', ...
                'Position', [15 284 75 22], 'HorizontalAlignment', 'right');
            fOutputRoot = uieditfield(dlg, 'text', ...
                'Position', [95 284 235 22], 'Editable', 'on', ...
                'Value', getpref('nestapp','outputRoot',''));
            uibutton(dlg, 'Text', 'Browse...', 'Position', [335 284 70 22], ...
                'ButtonPushedFcn', @(~,~) browseFolder(fOutputRoot));
            uilabel(dlg, 'Text', 'Layout:', ...
                'Position', [15 256 75 22], 'HorizontalAlignment', 'right');
            ddLayout = uidropdown(dlg, ...
                'Position', [95 256 175 22], ...
                'Items',     {'By type (data/reports/qc)', 'Per input file'}, ...
                'ItemsData', {'typeBased',                  'perInput'}, ...
                'Value',     getpref('nestapp','outputLayout','typeBased'));
            uilabel(dlg, 'Text', '(blank root = next to inputs)', ...
                'Position', [275 256 140 22], 'FontColor', [0.4 0.4 0.4], 'FontSize', 10);

            % --- Behaviour section ---
            uilabel(dlg, 'Text', 'Behaviour', 'FontWeight', 'bold', ...
                'Position', [15 223 200 20]);
            cbReport = uicheckbox(dlg, 'Text', 'Switch to Reports tab after each run', ...
                'Position', [15 199 300 22], ...
                'Value', getpref('nestapp','showReport',true));
            cbConfirm = uicheckbox(dlg, 'Text', 'Confirm before clearing pipeline', ...
                'Position', [15 175 300 22], ...
                'Value', getpref('nestapp','confirmClear',true));
            cbOverwrite = uicheckbox(dlg, 'Text', 'Overwrite existing report files (no timestamp)', ...
                'Position', [15 151 320 22], ...
                'Value', getpref('nestapp','overwriteReports',false));
            cbSuppressDialogs = uicheckbox(dlg, ...
                'Text', 'Suppress EEGLAB processing dialogs (warn about overwrites before run)', ...
                'Position', [15 127 390 22], ...
                'Value', getpref('nestapp','suppressEEGLABDialogs',true));
            cbHideEEGLAB = uicheckbox(dlg, ...
                'Text', 'Hide EEGLAB window during processing', ...
                'Position', [15 103 300 22], ...
                'Value', getpref('nestapp','hideEEGLABWindow',true));

            % --- Parallel Processing section ---
            uilabel(dlg, 'Text', 'Parallel Processing', 'FontWeight', 'bold', ...
                'Position', [15 72 200 20]);
            spnWorkers = [];
            if license('test', 'Distrib_Computing_Toolbox')
                uilabel(dlg, 'Text', 'Max workers:', ...
                    'Position', [15 48 85 22], 'HorizontalAlignment', 'right');
                spnWorkers = uispinner(dlg, ...
                    'Position', [105 48 60 22], 'Limits', [1 32], 'Step', 1, ...
                    'Value', getpref('nestapp', 'maxParallelWorkers', 4));
                uilabel(dlg, 'Text', 'cap on simultaneous files when Parallel is on', ...
                    'Position', [172 48 240 22], 'FontColor', [0.4 0.4 0.4]);
            else
                uilabel(dlg, 'Text', 'Not available - Parallel Computing Toolbox not licensed.', ...
                    'Position', [15 48 385 22], 'FontColor', [0.5 0.5 0.5]);
            end

            % --- Buttons ---
            uibutton(dlg, 'Text', 'Cancel', 'Position', [220 15 85 28], ...
                'ButtonPushedFcn', @(~,~) close(dlg));
            uibutton(dlg, 'Text', 'Save', 'Position', [315 15 85 28], ...
                'BackgroundColor', [0.20 0.55 0.20], 'FontColor', [1 1 1], ...
                'ButtonPushedFcn', @(~,~) savePrefs());

            uiwait(dlg);

            %% Nested helpers
            function browseEeglab()
                p = uigetdir('', 'Select EEGLAB Folder');
                if ~isequal(p, 0); fEeglab.Value = p; end
            end
            function browseFolder(field)
                p = uigetdir('', 'Select Folder');
                if ~isequal(p, 0); field.Value = p; end
            end
            function savePrefs()
                % EEGLAB path - warn on invalid value but don't abort
                % the rest of the save (otherwise the other preferences
                % the user just toggled would silently be discarded).
                ep = strtrim(fEeglab.Value);
                eeglabPathValid = isempty(ep) || isfolder(ep);
                if ~eeglabPathValid
                    uialert(dlg, ['EEGLAB path does not exist: ' ep ...
                        '. Other preferences were still saved.'], ...
                        'Invalid EEGLAB Path', 'Icon', 'warning');
                elseif ~isempty(ep)
                    addpath(ep);
                end
                setpref('nestapp', 'eeglabPath',          ep);
                setpref('nestapp', 'lastDataFolder',      strtrim(fData.Value));
                setpref('nestapp', 'lastPipelineFolder',  strtrim(fPipeline.Value));
                setpref('nestapp', 'outputRoot',          strtrim(fOutputRoot.Value));
                setpref('nestapp', 'outputLayout',        ddLayout.Value);
                setpref('nestapp', 'showReport',             cbReport.Value);
                setpref('nestapp', 'confirmClear',           cbConfirm.Value);
                setpref('nestapp', 'overwriteReports',       cbOverwrite.Value);
                setpref('nestapp', 'suppressEEGLABDialogs',  cbSuppressDialogs.Value);
                setpref('nestapp', 'hideEEGLABWindow',       cbHideEEGLAB.Value);
                if ~isempty(spnWorkers)
                    setpref('nestapp', 'maxParallelWorkers', round(spnWorkers.Value));
                end

                % Quality Screening prefs - validation mirrors
                % runPipelineCore: invalid mode -> minmax_no_tms,
                % inverted window -> [0 25].
                setpref('nestapp', 'autoQualityReport',    cbAutoQC.Value);
                setpref('nestapp', 'qualityTmsAutoDetect', cbTmsAuto.Value);
                setpref('nestapp', 'skipOnQualityFail',    cbSkipFail.Value);
                setpref('nestapp', 'autoExportPDF',        cbAutoPDF.Value);
                setpref('nestapp', 'aaratepKeepIntermediates', cbAaratepKeep.Value);

                attr = ddAttr.Value;
                if ~any(strcmp(attr, qualityAttributeModes()))
                    attr = 'minmax_no_tms';
                end
                setpref('nestapp', 'qualityAttribute', attr);

                w = [nfTmsStart.Value, nfTmsEnd.Value];
                if ~(isnumeric(w) && numel(w) == 2 && w(2) > w(1))
                    w = [0 25];
                end
                setpref('nestapp', 'qualityTmsWindow', w);

                close(dlg);
            end
            function v = resolveAttributePref()
                v = getpref('nestapp', 'qualityAttribute', 'minmax_no_tms');
                if ~any(strcmp(v, qualityAttributeModes()))
                    v = 'minmax_no_tms';
                end
            end
            function w = readTmsWindowPref()
                w = getpref('nestapp', 'qualityTmsWindow', [0 25]);
                if ~(isnumeric(w) && numel(w) == 2 && w(2) > w(1))
                    w = [0 25];
                end
            end
        end

        function updateReportsTabImpl(app)
        % UPDATEREPORTSTABIMPL  Refresh the Reports tab listbox from session and loaded reports.
        %   Combines app.allPipelineReports (from current run) with app.loadedReports
        %   (loaded from disk). Updates listbox labels and status text. When
        %   any report has Quality Gate data, appends a synthetic "Session
        %   Quality Dashboard" entry that swaps the right-side area to the
        %   dashboard panel on selection.
            allEntries = [app.allPipelineReports, app.loadedReports];
            n = numel(allEntries);
            if n == 0
                app.ReportsListBox.Items = {};
                app.ReportsListBox.ItemsData = {};
                app.ReportsStatusLabel.Text = 'No reports loaded.';
                app.ReportsTextArea.Value = '';
                app.ReportsDashboardPanel.Visible = 'off';
                app.ReportsTextArea.Visible = 'on';
                app.ExportReportsCSVButton.Enable = 'off';
                app.ExportPDFButton.Enable = 'off';
                app.CopyMethodsButton.Enable = 'off';
                app.OpenReportSetButton.Enable = 'off';
                return
            end

            % Append the Dashboard synthetic entry when any report has
            % gates - keeps the listbox tidy when nothing was screened.
            if app.dashboardHasContent(allEntries)
                dashEntry = struct('isDashboard', true, ...
                    'text', '', 'report', struct());
                allEntries{end+1} = dashEntry;
                n = numel(allEntries);
            end

            labels = cell(1, n);
            for i = 1:n
                e = allEntries{i};
                if isfield(e, 'isDashboard') && e.isDashboard
                    labels{i} = 'Session Quality Dashboard';
                elseif isfield(e, 'isSummary') && e.isSummary
                    % Extract file count from the summary header line
                    tok = regexp(e.text, 'PIPELINE SUMMARY\s+\((\d+) files\)', 'tokens', 'once');
                    if ~isempty(tok)
                        labels{i} = sprintf('Session Summary (%s files)', tok{1});
                    else
                        labels{i} = 'Session Summary';
                    end
                else
                    [~, baseName] = fileparts(e.report.inputFile);
                    try
                        dateLabel = string(e.report.processedAt, 'yyyy-MM-dd HH:mm');
                    catch
                        dateLabel = '?';
                    end
                    prefix = '';
                    if isfield(e.report, 'quality') ...
                            && isfield(e.report.quality, 'worstVerdict')
                        switch e.report.quality.worstVerdict
                            case 'Fail',     prefix = '[FAIL] ';
                            case 'Marginal', prefix = '[MARG] ';
                            case 'Pass',     prefix = '[PASS] ';
                        end
                    end
                    labels{i} = sprintf('%s%s (%s)', prefix, baseName, dateLabel);
                end
            end

            % Preserve selection index across refresh if still valid
            prevIdx = app.ReportsListBox.Value;
            app.ReportsListBox.Items = labels;
            app.ReportsListBox.ItemsData = num2cell(1:n);

            if isnumeric(prevIdx) && ~isempty(prevIdx) && prevIdx >= 1 && prevIdx <= n
                app.ReportsListBox.Value = prevIdx;
            else
                app.ReportsListBox.Value = n;
            end
            renderReportsRightPane(app, allEntries);

            nSess   = numel(app.allPipelineReports);
            nLoaded = numel(app.loadedReports);
            parts   = {};
            if nSess   > 0; parts{end+1} = sprintf('%d from session', nSess);   end
            if nLoaded > 0; parts{end+1} = sprintf('%d from disk', nLoaded); end
            app.ReportsStatusLabel.Text = strjoin(parts, ', ');
            app.ExportReportsCSVButton.Enable = 'on';
            app.ExportPDFButton.Enable = 'on';
            app.CopyMethodsButton.Enable = 'on';
        end

        function ReportsListBoxValueChanged(app, ~)
        % Callback - swap the right-side pane based on the selected entry.
            allEntries = [app.allPipelineReports, app.loadedReports];
            if app.dashboardHasContent(allEntries)
                allEntries{end+1} = struct('isDashboard', true, ...
                    'text', '', 'report', struct());
            end
            renderReportsRightPane(app, allEntries);
        end

        function renderReportsRightPane(app, allEntries)
        % Show the dashboard panel for the Dashboard entry; otherwise show
        % the report text in the text area. Called by both
        % updateReportsTabImpl (after a refresh) and the listbox callback.
            idx = app.ReportsListBox.Value;
            if isempty(idx) || ~isnumeric(idx) ...
                    || idx < 1 || idx > numel(allEntries)
                return
            end
            e      = allEntries{idx};
            isFile = isFileEntry(app, e);
            if isfield(e, 'isDashboard') && e.isDashboard
                showReportsPane(app, 'dashboard');
                renderDashboardPanel(app.ReportsDashboardPanel, ...
                    collectReportStructs(allEntries), ...
                    struct( ...
                        'onRefresh',        @() updateReportsTabImpl(app), ...
                        'onExport',         @() exportDashboardPNG(app, allEntries), ...
                        'onFailedRowClick', @(name) jumpToFileEntry(app, allEntries, name), ...
                        'failed',           app.lastFailed));
            elseif app.ReportsImageViewButton.Value && isFile
                % Keep the chosen checkpoint across a re-render: this runs on
                % every resize as well as on a selection change, and silently
                % snapping back to the first image mid-drag is disorienting.
                showReportsPane(app, 'images');
                renderQcImages(app.ReportsImagePanel, qcFiguresOf(app, e), ...
                    struct('selected', app.reportsQcIndex, ...
                           'onSelect', @(k) setQcIndex(app, k), ...
                           'onOpen',   @(d) revealFolder(app, d)));
            else
                showReportsPane(app, 'text');
                if isfield(e, 'text')
                    app.ReportsTextArea.Value = e.text;
                end
            end

            % The switch is meaningless on a synthetic row, which has no file
            % behind it and therefore no images.
            app.ReportsViewGroup.Visible = onOffState(isFile);

            % Open... is offered only when there is something to open: the
            % pipeline saved, and the file is still where it was written.
            app.OpenReportSetButton.Enable = ...
                onOffState(~isempty(selectedReportOutput(app, allEntries)));
        end

        function setQcIndex(app, k)
        % Remember the chosen checkpoint here rather than reading it back out
        % of the rendered panel: the pane is cleared and rebuilt on every
        % resize and selection change, and a second dropdown in that panel
        % would make a findall-based recovery silently pick the wrong one.
            app.reportsQcIndex = k;
        end

        function showReportsPane(app, which)
        % Exactly one of the three right-hand panes is visible. Kept in one
        % place so a new pane cannot leave an old one showing underneath.
            app.ReportsTextArea.Visible       = onOffState(strcmp(which, 'text'));
            app.ReportsDashboardPanel.Visible = onOffState(strcmp(which, 'dashboard'));
            app.ReportsImagePanel.Visible     = onOffState(strcmp(which, 'images'));
        end

        function ReportsViewChanged(app, ~)
        % Routed through the listbox callback because which pane to show is a
        % function of the selected ROW as well as the toggle - a synthetic row
        % has no images whatever the toggle says.
            ReportsListBoxValueChanged(app);
        end

        function OpenReportSetButtonPushed(app, ~)
        % The fast path to the file a report describes: skip the "which one?"
        % question Browse EEG asks and open this report's cleaned output.
            p = selectedReportOutput(app);
            if isempty(p)
                uialert(app.UIFigure, ...
                    ['This report does not name a saved recording. A pipeline ' ...
                     'only produces one when it ends with a Save New Set step, ' ...
                     'and reports written before this version did not record ' ...
                     'where it went.'], 'Nothing to open');
                return
            end
            openInEegplot(app, p);
        end


        function reRenderReportsOnResize(app)
        % Coalesce a drag into ONE repaint, the same restart-a-timer pattern
        % reRenderExploreOnResize uses and for the same reason: a resize event
        % arrives per pixel, and UIFigureSizeChanged's re-entrancy guard is a
        % throttle, not a debounce - it drops events arriving DURING a repaint,
        % not the ones queued behind it.
        %
        % This became necessary when the QC image pane joined the dashboard on
        % this path. Repainting it re-decodes a 1600x1200 PNG, so a drag was
        % hundreds of image decodes on the UI thread.
            if isempty(app.reportsResizeTimer) || ~isvalid(app.reportsResizeTimer)
                app.reportsResizeTimer = timer( ...
                    'StartDelay', 0.20, 'ExecutionMode', 'singleShot', ...
                    'Name', 'nestappReportsResize', ...
                    'TimerFcn', @(~, ~) repaintReportsNow(app));
            end
            stop(app.reportsResizeTimer);
            start(app.reportsResizeTimer);
        end

        function repaintReportsNow(app)
        % Reflow whichever absolutely-positioned pane is showing. Both
        % renderers clear and re-lay-out from the parent's size, so
        % re-rendering IS the reflow. The text area is a single control and
        % rescaleComponents already moves it.
            if ~isvalid(app) || ~isvalid(app.UIFigure); return; end
            panes = {app.ReportsDashboardPanel, app.ReportsImagePanel};
            visible = @(p) ~isempty(p) && isvalid(p) && strcmp(p.Visible, 'on');
            if ~any(cellfun(visible, panes)); return; end
            allEntries = [app.allPipelineReports, app.loadedReports];
            if app.dashboardHasContent(allEntries)
                allEntries{end+1} = struct('isDashboard', true, ...
                    'text', '', 'report', struct());
            end
            renderReportsRightPane(app, allEntries);
        end

        function exportDashboardPNG(app, allEntries)
        % Render the dashboard into an offscreen uifigure and save as PNG.
            [fname, fpath] = uiputfile('*.png', 'Export Quality Dashboard', ...
                'quality_dashboard.png');
            if isequal(fname, 0); return; end
            outPath = fullfile(fpath, fname);
            fig = uifigure('Visible', 'off', 'Position', [100 100 1200 800]);
            cleanup = onCleanup(@() close(fig, 'force'));
            renderDashboardPanel(fig, collectReportStructs(allEntries), ...
                struct('failed', app.lastFailed));
            try
                exportgraphics(fig, outPath, 'Resolution', 150);
                app.ReportsStatusLabel.Text = sprintf('Dashboard saved: %s', fname);
            catch err
                uialert(app.UIFigure, ...
                    sprintf('Export failed: %s', err.message), ...
                    'Export Dashboard PNG', 'Icon', 'error');
            end
        end

        function jumpToFileEntry(app, allEntries, fileName)
        % Failed-files table row click handler - select the listbox
        % entry for the given file basename so the user sees its text.
            for i = 1:numel(allEntries)
                e = allEntries{i};
                if isfield(e, 'isSummary') || isfield(e, 'isDashboard'), continue, end
                if ~isfield(e, 'report') || ~isfield(e.report, 'inputFile'), continue, end
                [~, name] = fileparts(e.report.inputFile);
                if strcmp(name, fileName)
                    app.ReportsListBox.Value = i;
                    renderReportsRightPane(app, allEntries);
                    return
                end
            end
        end

        function tf = dashboardHasContent(app, allEntries)
        % Whether to append the synthetic "Session Quality Dashboard" entry.
        % Worth showing when any report carries Quality Gate data OR when the
        % last run left files that did not complete - so failures are visible
        % even for a gate-less pipeline. Must be consulted at every site that
        % builds allEntries, or the dashboard row's index desyncs.
            tf = anyReportHasGates(allEntries) || ~isempty(app.lastFailed);
        end

        function LoadReportsButtonPushed(app, ~)
        % Browse for a folder of pipeline report .mat files and load them.
            folder = uigetdir(reportsBrowseStart(app), ...
                'Select Folder with Pipeline Reports');
            if isequal(folder, 0); return; end
            setpref('nestapp', 'lastReportsFolder', folder);

            matFiles = findReportMatFiles(folder);
            if isempty(matFiles)
                uialert(app.UIFigure, ...
                    ['No report .mat files (*_report.mat) found in that ' ...
                     'folder or its subfolders.'], 'No Reports Found');
                return
            end
            loaded = 0;
            for k = 1:numel(matFiles)
                fpath = fullfile(matFiles(k).folder, matFiles(k).name);
                try
                    S = load(fpath, 'pipelineReport');
                    if ~isfield(S, 'pipelineReport'); continue; end
                    % buildReportText, not exportReport: the latter is that
                    % plus a save() of the whole struct, so loading N reports
                    % wrote N throwaway .mat files to tempdir.
                    txt = buildReportText(S.pipelineReport);
                    entry.text   = txt;
                    entry.report = S.pipelineReport;
                    app.loadedReports{end+1} = entry;
                    loaded = loaded + 1;
                catch ME
                    warning('nestapp:loadReport', 'Could not load %s: %s', matFiles(k).name, ME.message);
                end
            end

            folderParts = strsplit(folder, {'\','/'});
            folderParts(cellfun(@isempty, folderParts)) = [];
            app.ReportsFolderLabel.Text = folderParts{end};
            updateReportsTabImpl(app);

            if loaded > 0
                app.TabGroup.SelectedTab = app.ReportsTab;
            end
        end

        function ClearReportsButtonPushed(app, ~)
        % Empty the report list, which otherwise accumulates across every run
        % in a session. Discards only the in-memory list - the reports are all
        % on disk - but there is no in-session undo, hence the confirm.
            allEntries = [app.allPipelineReports, app.loadedReports];
            if isempty(allEntries)
                app.ReportsStatusLabel.Text = 'No reports to clear.';
                return
            end

            answer = uiconfirm(app.UIFigure, sprintf( ...
                ['Remove all %d report(s) from the list?\n\n' ...
                 'This only empties the list - nothing on disk is deleted, and ' ...
                 '"Load from Folder" can bring saved reports back.'], numel(allEntries)), ...
                'Clear Report List', ...
                'Options', {'Clear', 'Cancel'}, ...
                'DefaultOption', 2, 'CancelOption', 2);
            if ~strcmp(answer, 'Clear'); return; end

            app.allPipelineReports = {};
            app.loadedReports      = {};
            % Stale failures would otherwise keep the Quality Dashboard entry
            % alive for reports that are no longer listed.
            app.lastFailed         = struct([]);
            app.ReportsFolderLabel.Text = '';
            updateReportsTabImpl(app);
            app.ReportsStatusLabel.Text = 'Report list cleared.';
        end

        function startFolder = reportsBrowseStart(~)
        % Where the folder picker should open. Reports live under the output
        % tree, so prefer the folder reports were last loaded from, then the
        % output root - NOT lastDataFolder, which points at raw input data and
        % is typically a network share. Handing uigetdir a remote or stale
        % path makes the dialog stall while the shell enumerates it, so every
        % candidate is checked with isfolder first.
            candidates = { ...
                getpref('nestapp', 'lastReportsFolder', ''), ...
                getpref('nestapp', 'outputRoot',        ''), ...
                getpref('nestapp', 'lastDataFolder',    '')};
            startFolder = '';
            for k = 1:numel(candidates)
                c = candidates{k};
                if ~isempty(c) && isfolder(c)
                    startFolder = c;
                    return
                end
            end
        end

        function ExportReportsCSVButtonPushed(app, ~)
        % Export a CSV table of key metrics for all visible reports.
            allEntries = [app.allPipelineReports, app.loadedReports];
            if isempty(allEntries)
                uialert(app.UIFigure, 'No reports to export.', 'Export CSV');
                return
            end

            [fname, fpath] = uiputfile('*.csv', 'Export Reports as CSV', 'nestapp_reports.csv');
            if isequal(fname, 0); return; end

            fid = fopen(fullfile(fpath, fname), 'w');
            if fid == -1
                uialert(app.UIFigure, 'Could not open file for writing.', 'Export CSV');
                return
            end

            % Header
            fprintf(fid, ['File,Processed,Channels (orig),Channels (final),' ...
                'Trials (orig),Trials (final),ICA removed,' ...
                'Quality_Verdict,Quality_Reasons\n']);

            for i = 1:numel(allEntries)
                e = allEntries{i};
                if isfield(e, 'isSummary') && e.isSummary; continue; end
                r = e.report;
                [~, baseName] = fileparts(r.inputFile);
                try
                    dStr = string(r.processedAt, 'yyyy-MM-dd HH:mm:ss');
                catch
                    dStr = '?';
                end

                verdict = 'NotChecked';
                reasons = '';
                if isfield(r, 'quality')
                    if isfield(r.quality, 'worstVerdict') ...
                            && ~isempty(r.quality.worstVerdict)
                        verdict = r.quality.worstVerdict;
                    end
                    if isfield(r.quality, 'gates') && ~isempty(r.quality.gates)
                        allReasons = {};
                        for gi = 1:numel(r.quality.gates)
                            g = r.quality.gates{gi};
                            if isfield(g, 'reasons') && ~isempty(g.reasons)
                                allReasons = [allReasons, g.reasons]; %#ok<AGROW>
                            end
                        end
                        reasons = strjoin(allReasons, '; ');
                        reasons = strrep(reasons, ',', ';'); % keep CSV-safe
                    end
                end

                fprintf(fid, '%s,%s,%d,%d,%d,%d,%d,%s,%s\n', ...
                    baseName, dStr, ...
                    r.channels.original, r.channels.final, ...
                    r.trials.original, r.trials.final, ...
                    r.ica.nRejected, verdict, reasons);
            end
            fclose(fid);
            app.ReportsStatusLabel.Text = sprintf('CSV saved: %s', fname);
        end

        function ExportPDFButtonPushed(app, ~)
        % Export report text + checkpoint PNGs as PDF, for the selected file
        % or for every listed report. One file at a time used to be the only
        % option, which made a whole batch tedious whenever Auto-export PDF
        % was off - and it is off by default.
            allEntries  = [app.allPipelineReports, app.loadedReports];
            fileReports = renderableReports(app, allEntries);
            if isempty(fileReports)
                uialert(app.UIFigure, 'No file reports to export.', 'Export PDF');
                return
            end

            selected = selectedFileReport(app, allEntries);
            switch askPdfScope(app, selected, numel(fileReports))
                case 'this'
                    exportOnePDF(app, selected);
                case 'all'
                    exportAllPDFs(app, fileReports);
            end
        end

        function reports = renderableReports(~, allEntries)
        % The report structs that exportFileReportPDF can actually render.
        % collectReportStructs already drops the synthetic Summary/Dashboard
        % rows; inputFile is what names the output file.
            reports = collectReportStructs(allEntries);
            keep    = cellfun(@(r) isstruct(r) && isfield(r, 'inputFile'), reports);
            reports = reports(keep);
        end

        function tf = isFileEntry(~, e)
        % A listbox row backed by a real per-file report, as opposed to the
        % synthetic Session Summary / Quality Dashboard rows. Those have no
        % file, so no images and nothing for the view switch to switch.
        %
        % Those two are already excluded by the test below rather than by name:
        % the summary sets report = [] and the dashboard report = struct(), so
        % neither survives isstruct + isfield(...,'inputFile').
            tf = isstruct(e) && isfield(e, 'report') && isstruct(e.report) ...
                 && isfield(e.report, 'inputFile');
        end

        function figs = qcFiguresOf(~, e)
        % The QC image paths recorded on a report, or {} for a report that
        % predates the field or ran without auto quality reporting.
            figs = {};
            if ~isstruct(e) || ~isfield(e, 'report') || ~isstruct(e.report); return; end
            r = e.report;
            if isfield(r, 'quality') && isstruct(r.quality) ...
                    && isfield(r.quality, 'figures') && iscell(r.quality.figures)
                figs = r.quality.figures;
            end
        end

        function r = selectedFileReport(app, allEntries)
        % The report struct behind the selected listbox row, or empty when the
        % selection is a synthetic row (Session Summary / Quality Dashboard),
        % which has nothing to render. Those rows sort after the real ones, so
        % the plain bounds check rejects them.
            r   = [];
            idx = app.ReportsListBox.Value;
            if isempty(idx) || ~isnumeric(idx) || idx < 1 || idx > numel(allEntries)
                return
            end
            candidate = allEntries{idx};
            if ~isfield(candidate, 'report') || ~isstruct(candidate.report) ...
                    || ~isfield(candidate.report, 'inputFile')
                return
            end
            r = candidate.report;
        end

        function choice = askPdfScope(app, selected, nAll)
        % Ask what to export. With no single report selected the only real
        % option is "all", so do not offer a "this report" button that would
        % have nothing to act on.
            if isempty(selected)
                msg = sprintf(['No single file report is selected - the Summary and ' ...
                    'Dashboard rows cannot be rendered on their own.\n\n' ...
                    'Export all %d file report(s)?'], nAll);
                options = {'All reports', 'Cancel'};
            else
                [~, baseName] = fileparts(selected.inputFile);
                msg = sprintf(['Export a PDF for "%s" only, or for all %d listed ' ...
                    'file report(s)?'], baseName, nAll);
                options = {'This report', 'All reports', 'Cancel'};
            end
            answer = uiconfirm(app.UIFigure, msg, 'Export PDF', ...
                'Options', options, ...
                'DefaultOption', 1, 'CancelOption', numel(options));
            switch answer
                case 'This report', choice = 'this';
                case 'All reports', choice = 'all';
                otherwise,          choice = 'cancel';
            end
        end

        function exportOnePDF(app, report)
        % One report, to a filename the user picks.
            [~, baseName] = fileparts(report.inputFile);
            [fname, fpath] = uiputfile('*.pdf', 'Export Report as PDF', ...
                [baseName, '_report.pdf']);
            if isequal(fname, 0); return; end
            try
                exportFileReportPDF(report, fullfile(fpath, fname));
                app.ReportsStatusLabel.Text = sprintf('PDF saved: %s', fname);
            catch err
                uialert(app.UIFigure, ...
                    sprintf('PDF export failed: %s', err.message), ...
                    'Export PDF', 'Icon', 'error');
            end
        end

        function exportAllPDFs(app, reports)
        % Every listed report into one folder, each named after its input
        % file. One report failing must not cost the user the rest, so
        % failures are collected and reported together at the end.
            folder = uigetdir('', 'Choose a folder for the PDF reports');
            if isequal(folder, 0); return; end

            nOk         = 0;
            failedNames = {};
            dlg = uiprogressdlg(app.UIFigure, 'Title', 'Export PDF', ...
                'Message', 'Rendering reports...', 'Indeterminate', 'off');
            cleanup = onCleanup(@() delete(dlg));
            for k = 1:numel(reports)
                [~, baseName] = fileparts(reports{k}.inputFile);
                dlg.Value   = k / numel(reports);
                dlg.Message = sprintf('Rendering %d of %d: %s', ...
                    k, numel(reports), baseName);
                try
                    exportFileReportPDF(reports{k}, ...
                        fullfile(folder, [baseName, '_report.pdf']));
                    nOk = nOk + 1;
                catch err
                    failedNames{end+1} = sprintf('%s (%s)', baseName, err.message); %#ok<AGROW>
                end
            end

            app.ReportsStatusLabel.Text = sprintf('%d PDF(s) saved to %s', nOk, folder);
            if ~isempty(failedNames)
                uialert(app.UIFigure, sprintf( ...
                    '%d report(s) could not be rendered:\n\n%s', ...
                    numel(failedNames), strjoin(failedNames, newline)), ...
                    'Export PDF', 'Icon', 'warning');
            end
        end

        function CopyMethodsButtonPushed(app, ~)
        % Copy a methods paragraph to the clipboard. A single file report copies
        % that file's full narrative; the Session Summary / Dashboard copies the
        % cross-file aggregate (mean +/- SD). Shares the same builders the text
        % reports use (methodsNarrative / methodsParagraphAggregate).
            idx = app.ReportsListBox.Value;
            if isempty(idx); return; end
            allEntries = [app.allPipelineReports, app.loadedReports];
            if ~isnumeric(idx) || idx < 1 || idx > numel(allEntries); return; end
            e = allEntries{idx};

            if isfield(e, 'report') && isstruct(e.report) && ~isempty(e.report)
                methodsText = methodsNarrative(e.report);
            else
                % Summary / Dashboard entry: aggregate over every per-file report.
                reportStructs = {};
                for k = 1:numel(allEntries)
                    ek = allEntries{k};
                    if (isfield(ek, 'isSummary')   && ek.isSummary)   ...
                            || (isfield(ek, 'isDashboard') && ek.isDashboard)
                        continue
                    end
                    if isfield(ek, 'report') && isstruct(ek.report) && ~isempty(ek.report)
                        reportStructs{end+1} = ek.report; %#ok<AGROW>
                    end
                end
                if isempty(reportStructs)
                    methodsText = 'TMS-EEG data were preprocessed using nestapp.';
                else
                    methodsText = methodsParagraphAggregate(reportStructs);
                end
            end

            clipboard('copy', methodsText);
            app.ReportsStatusLabel.Text = 'Methods text copied to clipboard.';
        end

        function showAbout(app)
        % SHOWABOUT  Display version and citation information.
            eeglabVer = '';
            if ~isempty(which('eeg_getversion'))
                try
                    eeglabVer = eeg_getversion();
                catch
                end
            end
            msg = sprintf([ ...
                'nestapp - TMS-EEG Processing\n\n' ...
                'nestapp: %s\n' ...
                'EEGLAB:  %s\n' ...
                'MATLAB:  %s\n\n' ...
                'Please cite:\n' ...
                'Rogasch et al. (2017) NeuroImage - TESA toolbox\n' ...
                'Delorme & Makeig (2004) J Neurosci Methods - EEGLAB'], ...
                nestappVersion(), eeglabVer, version);
            uialert(app.UIFigure, msg, 'About nestapp', 'Icon', 'info');
        end

        rescaleComponents(app, sX, sY)

        function styleParamTable(app)
        % Grey the value cells that hold a placeholder rather than a value.
        % The rule itself lives in greyPlaceholderCells, so the step table, the
        % plot options dialog and the figure dialog cannot drift apart on what
        % counts as "not set" - they already had, over the literal '[]'.
            greyPlaceholderCells(app.UITable);
        end

        function s = formatParamForDisplay(~, val, paramMeta)
        % FORMATPARAMFORDISPLAY  Render a typed param value as a display string.
        %   val: typed value from spec.params.(key)
        %   paramMeta: params struct entry from stepRegistry (has .placeholder, .type)
            if isempty(val)
                if ~isempty(paramMeta.placeholder)
                    s = paramMeta.placeholder;
                else
                    s = '(not set)';
                end
            elseif ischar(val) && ~isrow(val)
                s = [deblank(val(1,:)), ' ...'];   % char matrix - show first line
            elseif ischar(val) || isstring(val)
                s = char(val);
            elseif iscell(val)
                s = strjoin(cellfun(@char, val, 'UniformOutput', false), ', ');
            else
                s = mat2str(val);
            end
        end

        
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Apply tooltips to controls that are not self-explanatory.
        % Called once from startupFcn. Add new tooltips here as needed.
        function applyTooltips(app)
            % Pipeline builder
            app.AddButton.Tooltip           = 'Add the selected step to the pipeline';
            app.RemoveButton.Tooltip        = 'Remove the selected step from the pipeline';
            app.MoveUpButton.Tooltip        = 'Move the selected step earlier in the pipeline';
            app.MoveDownButton.Tooltip      = 'Move the selected step later in the pipeline';
            app.DefaultValueButton.Tooltip  = 'Reset the selected parameter to its default value';
            app.ReStartStepsButton.Tooltip  = ['Resume processing from the current step index. ' ...
                'Increment the step counter manually to skip steps, or reset to 1 to reprocess from the start.'];
            app.RunAnalysisButton.Tooltip   = 'Run the pipeline on the selected data file(s)';
            app.SelectDataButton.Tooltip    = 'Select a folder or individual files to process';
        end

        % Code that executes after component creation
        function startupFcn(app)
            clc
            % Order matters. The picker asks which() what is installed, and
            % nothing EEGLAB provides - including every plugin - resolves until
            % eeglab() has run its path setup and plugin scan. So prefs (which
            % is where the EEGLAB folder is recorded) and EEGLAB itself have to
            % come up BEFORE the tree is built; otherwise a cold MATLAB hides
            % most of the registry as "not installed". See ensureEeglabReady.
            loadPrefs(app);
            initEeglab(app);
            % Fill the stage-grouped step-picker tree and the name->info map it
            % reads on selection. Both are built from availableSteps, so what is
            % offered here is exactly what the pre-flight will accept.
            populateStepsTree(app);

            app.spec = repmat(struct('name','','params',struct()), 0, 1);
            app.SelectedListBox.Items(:) = [];
            app.SelectedListBox.ItemsData(:) = [];
            app.UITable.Data = [];
            app.ItemNum = 1;
            % Snapshot base geometry before originalSize so the resize handler
            % (guarded on originalSize) always has it. Components are still at
            % their createComponents positions here (no resize has scaled them).
            app.baseLayout        = captureBaseLayout(app);
            app.originalSize      = app.UIFigure.Position(3:4);
            initExploreTab(app);
            applyTooltips(app);
            % Dwell timer for the Steps tree legend. StartDelay is restarted by
            % every mouse move, so the tip appears only after the pointer has
            % been still this long.
            HOVER_DELAY = 3;
            app.hoverTimer = timer( ...
                'StartDelay',    HOVER_DELAY, ...
                'ExecutionMode', 'singleShot', ...
                'BusyMode',      'drop', ...
                'Name',          'nestappStepsHover', ...
                'TimerFcn',      @(~,~) showStepsTip(app));
            buildRecentFilesMenu(app);
            buildRecentPipelinesMenu(app);
            updateStatusBar(app);
            clc
        end


        % ── Explore tab ──────────────────────────────────────────────────
        function initExploreTab(app)
        % Fill the controls that do not depend on data: the plot catalogue, the
        % ROI presets, and the default windows.
            app.exploreWindows = defaultTEPComponentDefs();
            app.exploreRoi     = firstRoiPreset();
            refreshExploreRoi(app);
            refreshExploreWindows(app);
            refreshExplorePlots(app);
            refreshExploreGroups(app);
            refreshExploreDesign(app);
        end

        function refreshExplorePlots(app)
        % The picker is generated from the registry, so a new plot appears here
        % without touching the tab. Unavailable entries are LISTED with their
        % reason rather than hidden - the user is usually one group away from
        % them, and a feature that vanishes looks like a feature that is missing.
            ctx = struct('nGroups', numel(exploreGroupNames(app)), ...
                         'hasWindows', ~isempty(app.exploreWindows), ...
                         'hasChanlocs', exploreHasChanlocs(app));
            [entries, cats] = availablePlots(ctx);
            items = {}; data = {};
            for c = 1:numel(cats)
                for k = 1:numel(cats(c).entries)
                    e = cats(c).entries(k);
                    label = sprintf('%s  -  %s', cats(c).name, e.name);
                    if ~e.available
                        label = sprintf('%s   (%s)', label, e.reason);
                    end
                    items{end+1} = label; %#ok<AGROW>
                    data{end+1}  = e.name; %#ok<AGROW>
                end
            end
            keep = app.ExplorePlotDropDown.Value;
            app.ExplorePlotDropDown.Items     = items;
            app.ExplorePlotDropDown.ItemsData = data;
            if ~isempty(data)
                if ismember(keep, data)
                    app.ExplorePlotDropDown.Value = keep;
                else
                    app.ExplorePlotDropDown.Value = data{1};
                end
            end
            app.exploreAvailablePlots = entries;
        end

        function names = exploreGroupNames(app)
            names = {};
            if isempty(app.exploreEntries); return; end
            g = {app.exploreEntries.group};
            names = unique(g(~cellfun(@isempty, g)), 'stable');
        end

        function tf = exploreHasChanlocs(app)
            tf = false;
            ok = app.exploreCache;
            if isempty(ok); return; end
            for i = 1:numel(ok)
                if ok(i).ok && ~isempty(ok(i).chanlocs) && ...
                        isfield(ok(i).chanlocs, 'theta')
                    tf = true; return
                end
            end
        end

        function refreshExploreGroups(app)
            % datasetSummary returns the per-group array FIRST and the overall
            % totals second; reading groups off the totals silently produced an
            % empty list.
            groups = datasetSummary(app.exploreEntries);
            items = cell(1, numel(groups));
            data  = cell(1, numel(groups));
            for g = 1:numel(groups)
                % Show the n the ESTIMATE rests on, not the n assigned. Files on
                % a different cap are excluded by groupCurves, so the assigned
                % count can be higher - and the rail showing 16 beside a legend
                % showing 15 is how a user stops trusting either number.
                used = usedSubjectCount(app, groups(g).name);
                if isnan(used) || used == groups(g).nSubjects
                    items{g} = sprintf('%s   (n=%d subj, %d files)', ...
                        groups(g).name, groups(g).nSubjects, groups(g).nFiles);
                else
                    items{g} = sprintf('%s   (n=%d of %d subj, %d files)', ...
                        groups(g).name, used, groups(g).nSubjects, groups(g).nFiles);
                end
                data{g}  = groups(g).name;
            end
            app.ExploreGroupsListBox.Items     = items;
            app.ExploreGroupsListBox.ItemsData = data;
            app.ExploreRemoveGroupButton.Enable = onOffState(~isempty(data));
            hasData = ~isempty(data);
            app.ExploreFilesButton.Enable   = onOffState(~isempty(app.exploreEntries));
            app.ExploreFigureButton.Enable  = onOffState(hasData);
            app.ExploreCsvButton.Enable     = onOffState(hasData);
            app.ExploreResultsButton.Enable = onOffState(hasData);
            app.ExploreEmptyLabel.Visible   = onOffState(~hasData);
        end

        function n = usedSubjectCount(app, name)
        % Subjects behind the drawn estimate for this group, or NaN when there
        % is no result yet to ask.
            n = NaN;
            if isempty(app.exploreRes) || isempty(app.exploreRes.groups); return; end
            k = find(strcmp({app.exploreRes.groups.name}, name), 1);
            if ~isempty(k); n = app.exploreRes.est(k).n; end
        end

        function refreshExploreRoi(app)
            p = roiPresets();
            app.ExploreRoiDropDown.Items     = [{'(custom)'}, {p.name}];
            app.ExploreRoiDropDown.ItemsData = [{''}, {p.name}];
            match = '';
            for i = 1:numel(p)
                if isequal(sort(lower(p(i).labels)), sort(lower(app.exploreRoi)))
                    match = p(i).name; break
                end
            end
            app.ExploreRoiDropDown.Value = match;
            n = numel(app.exploreRoi);
            if n == 0
                app.ExploreRoiSummaryLabel.Text = 'No electrodes selected';
            else
                app.ExploreRoiSummaryLabel.Text = sprintf('%d electrode%s: %s', ...
                    n, plural(n), strjoin(app.exploreRoi, ' '));
            end
        end

        function w = exploreWindowColWidths(~, mode)
        % Widths per view. Both views show the same window list in the same
        % place, so the window column is 'Name' in both rather than 'Name'
        % here and 'Win' there - but results carries a Group column that
        % define does not, so one set of widths cannot serve both.
        %
        % Both sum to 232, which is what the 250 px table holds once the
        % vertical scrollbar is allowed for. Anything wider brings a
        % horizontal scrollbar, and then data is hidden rather than merely
        % scrolled to. Measured against the widest data the table holds:
        % N100/P180, a three-digit latency, a negative amplitude.
            if nargin > 1 && strcmp(mode, 'results')
                w = {50, 48, 46, 44, 44};   % Group Name Mean Peak-ms Peak-uV
            else
                w = {56, 46, 65, 65};       % Name T1 T2 Peak
            end
        end

        function refreshExploreWindows(app)
        % Two views of one list in the same space. The Analysis tab showed
        % bounds AND measures in six columns; the rail holds four, so it
        % switches rather than dropping the measures as the first version did.
            w = app.exploreWindows;
            if strcmp(app.ExploreWindowsModeDropDown.Value, 'results')
                showExploreWindowResults(app, w);
            else
                app.ExploreWindowsTable.ColumnName     = {'Name'; 'T1'; 'T2'; 'Peak'};
                app.ExploreWindowsTable.ColumnWidth    = exploreWindowColWidths(app);
                app.ExploreWindowsTable.ColumnEditable = [true true true true];
                app.ExploreWindowsTable.ColumnFormat   = ...
                    {'char', 'numeric', 'numeric', {'auto', 'pos', 'neg'}};
                data = cell(numel(w), 4);
                for i = 1:numel(w)
                    data(i, :) = {w(i).name, w(i).winStart, w(i).winEnd, ...
                                  windowPolarity(w(i))};
                end
                app.ExploreWindowsTable.Data = data;
            end
        end

        function showExploreWindowResults(app, w)
        % Every group's measures, the groups paired under each window.
        %
        % This showed one group's numbers with nothing on screen saying whose:
        % the group went into the section label, which is allotted 70 px and
        % truncated to 'WINDOWS: ...', and an unselected groups list fell back
        % to the first group rather than showing nothing. In a tab whose
        % purpose is comparing groups, one group's measures presented as "the"
        % measures is the worst available answer. The group is a column now,
        % so nothing depends on a selection and nothing can be hidden by a
        % label that does not fit.
            app.ExploreWindowsTable.ColumnName     = {'Group'; 'Name'; 'Mean'; ...
                sprintf('Peak\nms'); sprintf('Peak\nuV')};
            app.ExploreWindowsTable.ColumnWidth    = exploreWindowColWidths(app, 'results');
            app.ExploreWindowsTable.ColumnEditable = false(1, 5);
            app.ExploreWindowsTable.ColumnFormat   = repmat({'char'}, 1, 5);
            app.ExploreWindowsLabel.Text = 'WINDOWS';

            if isempty(app.exploreRes) || isempty(app.exploreRes.groups)
                app.ExploreWindowsTable.Data = [repmat({'-'}, numel(w), 1), ...
                    {w.name}', repmat({'-'}, numel(w), 3)];
                return
            end
            names = {app.exploreRes.groups.name};
            nG    = numel(names);

            % TESA's own detector, once per group, so this table and any peak
            % overlay drawn from the same curve cannot disagree.
            peaks    = cell(1, nG);
            haveTesa = false(1, nG);
            if strcmpi(currentMode(app), 'TEP')
                for k = 1:nG
                    % Read by the evalc below, which checkcode cannot see into.
                    curve = app.exploreRes.est(k).mean; %#ok<NASGU>
                    try
                        % evalc: TESA prints a line per component per call, and
                        % this runs on every window edit - nG times that would
                        % bury the app's own logging.
                        pk = [];
                        evalc('pk = tepPeakFinder(curve, app.exploreRes.time, w);');
                        peaks{k}    = pk;
                        haveTesa(k) = ~isempty(pk);
                    catch
                        peaks{k} = [];   % TESA absent
                    end
                end
            end

            % Windows outermost, groups within: the comparison this tab exists
            % for is between groups at one window, so those rows sit adjacent.
            data = cell(numel(w) * nG, 5);
            r = 0;
            for i = 1:numel(w)
                for k = 1:nG
                    r = r + 1;
                    m = computeWindowMeasures(app.exploreRes.est(k).mean, ...
                            app.exploreRes.time, w(i).winStart, w(i).winEnd, ...
                            windowPolarity(w(i)));
                    data{r, 1} = names{k};
                    data{r, 2} = w(i).name;
                    data{r, 3} = num2str(m.mean, '%.2f');
                    % tepPeakFinder reports latencyMs/amplitudeUV with a
                    % `found` flag; computeWindowMeasures reports
                    % peakLatency/peakAmp. Prefer TESA's detection so the
                    % table agrees with the overlay, and show '-' rather than
                    % a number where no peak was found.
                    if haveTesa(k) && i <= numel(peaks{k})
                        % When TESA ran, its verdict stands - including "no
                        % peak here", shown as '-' rather than quietly
                        % substituting the window extremum: a number where
                        % there is no peak is worse than a dash.
                        if peaks{k}(i).found
                            data{r, 4} = num2str(peaks{k}(i).latencyMs, '%.0f');
                            data{r, 5} = num2str(peaks{k}(i).amplitudeUV, '%.2f');
                        else
                            data{r, 4} = '-';
                            data{r, 5} = '-';
                        end
                    elseif m.found
                        data{r, 4} = num2str(m.peakLatency, '%.0f');
                        data{r, 5} = num2str(m.peakAmp, '%.2f');
                    else
                        data{r, 4} = '-';
                        data{r, 5} = '-';
                    end
                end
            end
            app.ExploreWindowsTable.Data = data;
        end

        function name = selectedExploreGroup(app)
        % '' when nothing is selected. A listbox with ItemsData and no selection
        % returns {}, and strcmp(cellArray, {}) is a size-mismatch error rather
        % than a miss - so every read of the selection goes through here.
            name = '';
            v = app.ExploreGroupsListBox.Value;
            if ischar(v)
                name = v;
            elseif iscell(v) && ~isempty(v) && ischar(v{1})
                name = v{1};
            elseif isstring(v) && isscalar(v)
                name = char(v);
            end
        end

        function recomputeExplore(app)
        % One path from state to picture. Everything the rail changes lands
        % here, and it re-runs groupCurves rather than patching a cached
        % result - the whole point of caching trial averages is that this is
        % arithmetic on a few MB, so there is no stale-state class of bug.
            app.exploreRes = struct([]);
            if isempty(exploreGroupNames(app))
                refreshExploreGroups(app);
                refreshExploreDesign(app);
                refreshExploreWindows(app);
                renderExplorePlot(app);
                return
            end
            % Before reading the design: paired may have become undefined as
            % groups changed, in which case this flips the control back to
            % unpaired. Reading first would use a design the data cannot support.
            refreshExploreDesign(app);
            % Availability depends on the group count, so the catalogue has to
            % be re-evaluated whenever the group set changes - otherwise every
            % plot stays marked with the count it had when the tab was built and
            % rendering refuses to draw.
            refreshExplorePlots(app);
            entry = currentPlotEntry(app);
            mode = 'TEP';
            if ~isempty(entry) && ~isempty(entry.mode); mode = entry.mode; end
            try
                % No 'level' here: groupCurves owns that default and records
                % what it used in res.info.level, which is what the status
                % line and the figure footer read. Naming 95 here made this a
                % second place the number lived.
                app.exploreRes = groupCurves(app.exploreCache, app.exploreEntries, ...
                    struct('roi', {app.exploreRoi}, 'mode', mode, ...
                           'design', exploreDesign(app)));
            catch ME
                app.ExploreStatusLabel.Text = ME.message;
                app.exploreRes = struct([]);
            end
            refreshExploreGroups(app);
            refreshExploreWindows(app);
            renderExplorePlot(app);
        end

        function d = exploreDesign(app)
        % Read from the control. This used to be inferred from subject ids -
        % which are a guess - so a naming coincidence could switch the design to
        % paired and narrow every confidence interval without saying so. The
        % default is unpaired because it is the conservative interval.
            if app.ExplorePairedButton.Value && strcmp(app.ExplorePairedButton.Enable, 'on')
                d = 'paired';
            else
                d = 'unpaired';
            end
        end

        function refreshExploreDesign(app)
        % Offer paired whenever it is defined: two or more groups, and at least
        % two subjects present in every one of them - a paired interval is the
        % spread of within-subject differences, which one subject cannot give.
        %
        % Not "every subject in every group". Real cohorts are rarely that
        % symmetric - one participant misses a post session - and groupCurves
        % already handles it: a paired estimate keeps the complete subjects,
        % in the same order in every group, and names the rest in res.dropped,
        % which the status line reports. Demanding perfect symmetry here locked
        % paired off for the ordinary cohort while the analysis supported it.
            [~, overall] = datasetSummary(app.exploreEntries);
            names   = exploreGroupNames(app);
            canPair = numel(names) >= 2 && overall.nComplete >= 2;

            app.ExplorePairedButton.Enable = onOffState(canPair);
            if ~canPair && app.ExplorePairedButton.Value
                app.ExploreUnpairedButton.Value = true;   % never leave it on a
            end                                           % design that is undefined

            if numel(names) < 2
                app.ExploreDesignNoteLabel.Text = 'paired needs two or more groups';
            elseif canPair && overall.nComplete == overall.nSubjects
                app.ExploreDesignNoteLabel.Text = sprintf( ...
                    'paired available: %d complete sets', overall.nComplete);
            elseif canPair
                app.ExploreDesignNoteLabel.Text = sprintf( ...
                    'paired uses the %d of %d subjects in every group', ...
                    overall.nComplete, overall.nSubjects);
            elseif overall.nComplete == 1
                app.ExploreDesignNoteLabel.Text = ...
                    'paired needs 2+ subjects in every group (only 1)';
            else
                app.ExploreDesignNoteLabel.Text = 'no subject is in every group';
            end
        end

        function ExploreDesignChanged(app, ~)
            recomputeExplore(app);
        end

        function mode = currentMode(app)
        % The curve mode the selected plot reduces with. Read from the registry
        % entry rather than kept as separate state - the mode IS a property of
        % the chosen plot, which is why the old Plot Type radios became registry
        % entries.
            mode  = 'TEP';
            entry = currentPlotEntry(app);
            if ~isempty(entry) && ~isempty(entry.mode)
                mode = entry.mode;
            end
        end

        function entry = currentPlotEntry(app)
            entry = [];
            if isempty(app.exploreAvailablePlots); return; end
            k = find(strcmp({app.exploreAvailablePlots.name}, ...
                            app.ExplorePlotDropDown.Value), 1);
            if ~isempty(k); entry = app.exploreAvailablePlots(k); end
        end

        function renderExplorePlot(app)
        % Clear the canvas and build exactly the axes the chosen plot needs.
            delete(app.ExploreCanvas.Children);
            app.ExploreEmptyLabel = uilabel(app.ExploreCanvas, ...
                'Position', [20 180 app.ExploreCanvas.Position(3) - 40 40], ...
                'HorizontalAlignment', 'center', 'FontSize', 13, ...
                'FontColor', [0.45 0.48 0.53], 'WordWrap', 'on', ...
                'Text', 'Add a group to begin. A group is a set of recordings compared as one condition - pre and post, or one cohort against another.', ...
                'Visible', 'off');

            entry = currentPlotEntry(app);
            refreshExploreOptionsButton(app, entry);
            if isempty(entry); return; end
            app.ExplorePlotInfoLabel.Text = '';

            if ~entry.available
                app.ExplorePlotInfoLabel.Text = entry.reason;
                showExploreEmpty(app, entry.reason);
                return
            end
            if isempty(app.exploreRes) || isempty(app.exploreRes.groups)
                showExploreEmpty(app, '');
                return
            end

            try
                drawExploreInto(app, app.ExploreCanvas, entry);
                app.ExploreStatusLabel.Text = exploreStatusText(app);
            catch ME
                showExploreEmpty(app, ME.message);
                app.ExploreStatusLabel.Text = ME.message;
            end
        end

        function refreshExploreOptionsButton(app, entry)
        % Disabled rather than hidden when a plot has no settings: a button
        % that comes and goes as the picker changes is harder to find than one
        % that is always in the same place and sometimes greyed.
            has = ~isempty(entry) && ~isempty(entry.params);
            app.ExplorePlotOptionsButton.Enable = matlab.lang.OnOffSwitchState(has);
            if has
                n = numel(entry.params);
                app.ExplorePlotOptionsButton.Tooltip = sprintf( ...
                    '%d setting%s for %s', n, plural(n), entry.name);
            else
                app.ExplorePlotOptionsButton.Tooltip = 'This plot has no settings.';
            end
        end

        function reRenderExploreOnResize(app)
        % Coalesce a drag into ONE repaint. A resize event arrives per pixel of
        % a drag, and a TEP-topo grid is twelve topoplots - close to a second
        % of work each time. Redrawing per event queues that up and the window
        % stops following the mouse.
        %
        % Same restart-a-timer pattern the steps-tree hover tip already uses:
        % every event pushes the repaint further out, so it happens once the
        % drag stops. The re-entrancy guard alone could not do this - it drops
        % events that arrive DURING a repaint, not the ones queued behind it.
            if isempty(app.exploreResizeTimer) || ~isvalid(app.exploreResizeTimer)
                app.exploreResizeTimer = timer( ...
                    'StartDelay', 0.20, 'ExecutionMode', 'singleShot', ...
                    'Name', 'nestappExploreResize', ...
                    'TimerFcn', @(~, ~) repaintExploreNow(app));
            end
            stop(app.exploreResizeTimer);
            start(app.exploreResizeTimer);
        end

        function repaintExploreNow(app)
        % No-op unless Explore is the selected tab: repainting a hidden panel
        % is work nobody sees.
            if ~isvalid(app) || ~isvalid(app.UIFigure); return; end
            if app.TabGroup.SelectedTab ~= app.ExploreTab; return; end
            if isempty(app.exploreRes); return; end
            renderExplorePlot(app);
        end

        function drawExploreInto(app, parent, entry, axesFcn)
        % Mint the axes the plot wants inside `parent`, then hand off to the
        % registry's draw function. Shared by the in-app canvas and the
        % exported figure, so what is exported is what was on screen.
        %
        % axesFcn decides what KIND of axes: uiaxes for the canvas, classic
        % axes for anything headed for a file. It is a parameter rather than
        % something sniffed from the parent because every export path -
        % exportgraphics, print, saveas - silently omits UI components, so a
        % figure of uiaxes exports as a blank page with only a warning to say
        % so. The choice belongs to the caller who knows where it is going.
            if nargin < 4 || isempty(axesFcn)
                axesFcn = @(p, pos) uiaxes(p, 'Position', pos);
            end
            res = app.exploreRes;
            pos = parent.Position;
            fn  = str2func(entry.draw);

            % The user's settings go in FIRST and the tab's context fills what
            % is left, so a param the user set always wins and one they never
            % touched is simply absent - leaving the draw function's own
            % default to apply. Every drawer ignores the context fields it does
            % not declare, so the fill is uniform and each branch below states
            % only how many axes to mint.
            opts = fillDefaults( ...
                plotDrawOpts(entry, explorePlotParamsFor(app, entry.name)), ...
                struct( ...
                'mode',    entry.mode, ...
                'windows', app.exploreWindows, ...
                'window',  exploreTopoWindow(app), ...
                'axesFcn', axesFcn));

            switch exploreLayoutOf(app, entry)
                case 'panel'
                    % Owns the canvas: it decides how many axes it needs from
                    % the group and window counts, so it mints them itself.
                    fn(parent, res, opts);

                case 'per-group'
                    n      = numel(res.groups);
                    axList = gobjects(1, n);
                    BAR_W  = 62;    % the strip the shared colour bar sits in
                    w      = (pos(3) - 20 - BAR_W) / max(n, 1);
                    for k = 1:n
                        axList(k) = axesFcn(parent, ...
                            [10 + (k-1)*w, 30, w - 10, pos(4) - 60]);
                    end
                    % One bar for the row, hung beside it rather than on the
                    % last map - attaching one shrinks its host axes, which
                    % would leave that group's head smaller than its siblings
                    % in a figure whose point is that they are comparable.
                    clim = fn(axList, res, ...
                              fillDefaults(struct('colorbar', false), opts));
                    if ~isempty(clim)
                        sharedColorbar(parent, axesFcn, ...
                            [pos(3) - BAR_W + 8, 40, 12, pos(4) - 80], clim);
                    end
                    % Empty clim means the maps no longer share a scale, so
                    % each carries its own bar and one hung here would be a
                    % lie. The strip stays reserved and empty rather than
                    % reclaimed: the axes are minted before the drawer runs,
                    % so the layout cannot know, and asking the tab to read
                    % the plot's own params to find out is exactly the
                    % plot-specific knowledge the layout indirection exists
                    % to keep out of here.

                otherwise   % 'single'
                    ax = axesFcn(parent, [45 45 pos(3) - 70 pos(4) - 70]);
                    fn(ax, res, opts);
            end
        end

        function layout = exploreLayoutOf(~, entry)
        % Registry entries predating the .layout field, and the literal structs
        % the tests build, default to a single axes.
            layout = 'single';
            if isfield(entry, 'layout') && ~isempty(entry.layout)
                layout = entry.layout;
            end
        end

        function params = explorePlotParamsFor(app, name)
            params = struct();
            if isempty(app.explorePlotParams); return; end
            k = find(strcmp({app.explorePlotParams.name}, name), 1);
            if ~isempty(k); params = app.explorePlotParams(k).params; end
        end

        function win = exploreTopoWindow(app)
        % The scalp map averages over the first window of interest, so the map
        % and the measures describe the same interval by construction.
            win = [app.exploreRes.time(1) app.exploreRes.time(end)];
            if ~isempty(app.exploreWindows)
                win = [app.exploreWindows(1).winStart app.exploreWindows(1).winEnd];
            end
        end

        function showExploreEmpty(app, msg)
            if ~isempty(msg)
                app.ExploreEmptyLabel.Text = msg;
            end
            app.ExploreEmptyLabel.Visible = 'on';
        end

        function txt = exploreStatusText(app)
            res = app.exploreRes;
            [groups, design] = exploreCohortText(app);
            txt = sprintf('%s  |  %s  |  %d channels', ...
                groups, design, numel(res.channelLabels));
            m = res.info.montage;
            if ~isempty(m.excluded)
                txt = sprintf('%s  |  %d file%s excluded (different cap)', ...
                    txt, numel(m.excluded), plural(numel(m.excluded)));
            end
            if ~isempty(res.dropped)
                txt = sprintf('%s  |  %d subject%s without a complete set', ...
                    txt, numel(res.dropped), plural(numel(res.dropped)));
            end
            % A partial ROI is averaged rather than refused, so it has to be
            % said - otherwise the curve is over fewer electrodes than the
            % rail claims and nothing on screen differs.
            roi = exploreRoiInfo(app);
            if ~isempty(roi) && ~isempty(roi.missing)
                txt = sprintf('%s  |  ROI missing %s', txt, strjoin(roi.missing, ' '));
            end
        end

        function info = exploreRoiInfo(app)
        % res.info.roi, or [] for a result computed before the field existed.
        % One guard, because both readers below need the same three checks.
            info = [];
            r = app.exploreRes;
            if isempty(r) || ~isfield(r, 'info') || ~isfield(r.info, 'roi'); return; end
            info = r.info.roi;
        end

        % -- callbacks ------------------------------------------------------
        function ExploreAddGroupButtonPushed(app, ~)
            startFolder = getpref('nestapp', 'lastDataFolder', '');
            paths = selectDataTree(startFolder, {'*.set'});
            if isempty(paths); return; end
            if isvalid(app.UIFigure); figure(app.UIFigure); end

            answer = inputdlg('Name for this group:', 'Add group', [1 40], ...
                              {sprintf('group %d', numel(exploreGroupNames(app)) + 1)});
            if isempty(answer) || isempty(strtrim(answer{1})); return; end
            name = strtrim(answer{1});

            % loadReducedSets reads through pop_loadset, so this is an EEGLAB
            % entry point like Run Analysis and Browse EEG, and has to ask.
            [ok, msg] = ensureEeglabReady();
            if ~ok
                uialert(app.UIFigure, msg, 'EEGLAB Not Ready', 'Icon', 'error');
                return
            end

            dlg = uiprogressdlg(app.UIFigure, 'Title', 'Adding group', ...
                'Message', sprintf('Loading %d file(s)...', numel(paths)), ...
                'Indeterminate', 'on');
            closeDlg = onCleanup(@() close(dlg));
            try
                added = exploreDataset(paths, {}, struct());
                for i = 1:numel(added); added(i).group = name; end
                [newCache, warns] = loadReducedSets(paths, struct( ...
                    'progressFcn', @(i, n, ~) 0));
                app.exploreEntries = [app.exploreEntries, added];
                app.exploreCache   = [app.exploreCache, newCache];
                if isempty(app.exploreRoi)
                    app.exploreRoi = firstRoiPreset();
                    refreshExploreRoi(app);
                end
                clear closeDlg
                if ~isempty(warns)
                    uialert(app.UIFigure, strjoin(warns, newline), ...
                        'Some files could not be loaded', 'Icon', 'warning');
                end
            catch ME
                clear closeDlg
                uialert(app.UIFigure, ME.message, 'Could not add group', ...
                        'Icon', 'error');
                return
            end
            refreshExplorePlots(app);
            recomputeExplore(app);
        end

        function ExploreRemoveGroupButtonPushed(app, ~)
            name = selectedExploreGroup(app);
            if isempty(name); return; end
            drop = strcmp({app.exploreEntries.group}, name);
            gone = {app.exploreEntries(drop).path};
            app.exploreEntries(drop) = [];
            % Drop the cached datasets no remaining group refers to.
            keepCache = ~ismember({app.exploreCache.path}, gone);
            app.exploreCache = app.exploreCache(keepCache);
            refreshExplorePlots(app);
            recomputeExplore(app);
        end

        function ExploreFilesButtonPushed(app, ~)
        % The whole subject story is settled here: what n is, and why.
            [edited, rule] = exploreFilesTable(app.exploreEntries, ...
                struct('parent', app.UIFigure, 'subjectRule', app.exploreSubjectRule));
            if isempty(edited); return; end     % cancelled
            app.exploreEntries = edited;
            app.exploreSubjectRule = rule;
            refreshExplorePlots(app);
            recomputeExplore(app);
        end

        function ExploreWindowsModeChanged(app, ~)
            refreshExploreWindows(app);
        end

        function ExploreWindowsAddButtonPushed(app, ~)
            n = numel(app.exploreWindows);
            % 'auto' - the largest absolute deflection - rather than guessing a
            % sign. A new window used to be born 'pos' with no way to change it,
            % so a user-added negative component had its peak read upside down
            % and nothing said so.
            app.exploreWindows(n + 1) = struct('name', sprintf('W%d', n + 1), ...
                'polarity', 'auto', 'nomLatency', 100, 'winStart', 100, 'winEnd', 150);
            app.ExploreWindowsModeDropDown.Value = 'define';
            refreshExploreWindows(app);
            renderExplorePlot(app);
        end

        function ExploreWindowsRemoveButtonPushed(app, ~)
            % Results mode lists group x window, so a row index is not a
            % window index there. Switch to the view where it is rather than
            % deleting whichever window the arithmetic happens to land on.
            if strcmp(app.ExploreWindowsModeDropDown.Value, 'results')
                app.ExploreWindowsModeDropDown.Value = 'define';
                refreshExploreWindows(app);
                return
            end
            sel = app.ExploreWindowsTable.Selection;
            if isempty(sel) || isempty(app.exploreWindows); return; end
            r = sel(1, 1);
            if r < 1 || r > numel(app.exploreWindows); return; end
            app.exploreWindows(r) = [];
            refreshExploreWindows(app);
            renderExplorePlot(app);
        end

        function ExploreRoiEditButtonPushed(app, ~)
            picked = roiPicker(app.exploreRoi, exploreAvailableElectrodes(app), ...
                struct('parent', app.UIFigure, ...
                       'availableNote', 'not on the modal cap', ...
                       'optional', {partialElectrodes(app.exploreCache, ...
                                      exploreAvailableElectrodes(app))}));
            if isempty(picked) && ~iscell(picked)
                return          % cancelled - [] rather than {}
            end
            app.exploreRoi = picked;
            refreshExploreRoi(app);
            recomputeExplore(app);
        end

        function labels = exploreAvailableElectrodes(app)
        % The modal montage of what is loaded, which is what groupCurves will
        % actually compute on - so the picker greys exactly what cannot be used.
            labels = {};
            if isempty(app.exploreRes) || ~isfield(app.exploreRes, 'channelLabels')
                ok = app.exploreCache([app.exploreCache.ok]);
                if ~isempty(ok); labels = ok(1).labels; end
                return
            end
            labels = app.exploreRes.channelLabels;
        end

        function ExploreRoiDropDownValueChanged(app, ~)
            name = app.ExploreRoiDropDown.Value;
            if isempty(name); return; end
            p = roiPresets();
            k = find(strcmp({p.name}, name), 1);
            if isempty(k); return; end
            app.exploreRoi = p(k).labels;
            refreshExploreRoi(app);
            recomputeExplore(app);
        end

        function ctx = explorePlotContext(app)
        % Choices for params that declare 'choicesFrom' - a list the registry
        % cannot state because it is whatever the user has in the table right
        % now. One key, because one param needs one; a context field with no
        % consumer is surface to keep in step for nothing.
            names = {};
            if ~isempty(app.exploreWindows); names = {app.exploreWindows.name}; end
            ctx = struct('windows', {names});
        end

        function ExplorePlotOptionsButtonPushed(app, ~)
            entry = currentPlotEntry(app);
            if isempty(entry) || isempty(entry.params); return; end

            % The plot follows the controls while the dialog is open, and the
            % dialog's cancel calls back once more with the original values,
            % so the picture is restored along with the struct.
            apply = @(p) previewExplorePlotParams(app, entry.name, p);
            [params, accepted] = plotOptionsDialog(entry, ...
                explorePlotParamsFor(app, entry.name), app.UIFigure, ...
                explorePlotContext(app), apply);
            if ~accepted; return; end

            % Store again rather than trusting the preview to have done it:
            % plotOptionsDialog defaults onApply to a no-op, so a caller that
            % passes none must still end up with the accepted values.
            storeExplorePlotParams(app, entry.name, params);

            % The last live edit may have left a debounced repaint pending.
            % Cancel it - the immediate repaint below supersedes it, and
            % letting it fire would redraw the whole grid a second time from
            % identical params, which on TEP-topo is another twelve topoplots.
            cancelExploreRepaint(app);
            repaintExploreNow(app);
        end

        function storeExplorePlotParams(app, name, params)
            app.explorePlotParams = upsertByName(app.explorePlotParams, ...
                struct('name', name, 'params', params));
        end

        function previewExplorePlotParams(app, name, params)
        % Store, then schedule the DEBOUNCED repaint - the same coalescing the
        % resize handler uses, so a run of quick clicks in a listbox does not
        % queue a second of work behind each one.
            storeExplorePlotParams(app, name, params);
            reRenderExploreOnResize(app);
        end

        function cancelExploreRepaint(app)
            if ~isempty(app.exploreResizeTimer) && isvalid(app.exploreResizeTimer)
                stop(app.exploreResizeTimer);
            end
        end

        function ExplorePlotDropDownValueChanged(app, ~)
            recomputeExplore(app);
        end

        function ExploreWindowsTableCellEdit(app, event)
            r = event.Indices(1);
            c = event.Indices(2);
            if r > numel(app.exploreWindows); return; end
            switch c
                case 1, app.exploreWindows(r).name     = char(event.NewData);
                case 2, app.exploreWindows(r).winStart = event.NewData;
                case 3, app.exploreWindows(r).winEnd   = event.NewData;
                case 4, app.exploreWindows(r).polarity = char(event.NewData);
            end
            refreshExploreWindows(app);
            renderExplorePlot(app);
        end

        function ExploreWindowsResetButtonPushed(app, ~)
            app.exploreWindows = defaultTEPComponentDefs();
            refreshExploreWindows(app);
            renderExplorePlot(app);
        end

        function ExploreFigureButtonPushed(app, ~)
        % A real MATLAB figure, drawn fresh rather than copied: the in-app
        % canvas is for looking, and a figure someone will publish wants its own
        % axes at the size it will be printed. Composed by the same
        % drawExploreInto the canvas uses, so what leaves is what was on screen.
            if isempty(app.exploreRes); return; end
            entry = currentPlotEntry(app);
            if isempty(entry) || ~entry.available; return; end

            [opts, action] = figureExportDialog(app.exploreFigureOpts, app.UIFigure);
            if strcmp(action, 'cancel'); return; end
            % Remembered for the session: someone exporting a set of figures
            % for one manuscript wants the same width every time, and re-picking
            % it per figure is how they end up inconsistent.
            app.exploreFigureOpts = rmfield(opts, 'file');

            opts.title      = sprintf('nestapp - %s', entry.name);
            opts.provenance = exploreProvenance(app, entry);
            if strcmp(action, 'open'); opts.file = ''; end

            try
                publicationFigure(@(parent, axesFcn) ...
                    drawExploreInto(app, parent, entry, axesFcn), opts);
            catch ME
                uialert(app.UIFigure, ME.message, 'Could not make the figure');
                return
            end
            if strcmp(action, 'save')
                app.ExploreStatusLabel.Text = sprintf('Saved %s (%s mm, %d dpi)', ...
                    opts.file, num2str(opts.width), opts.dpi);
            end
        end

        function s = exploreRoiText(app)
        % What the curves were ACTUALLY averaged over, which is not always what
        % was asked for. Stamping the request would let a figure caption name
        % five electrodes over a mean of three.
            s = strjoin(app.exploreRoi, ' ');
            roi = exploreRoiInfo(app);
            if isempty(roi); return; end
            s = strjoin(roi.matched, ' ');
            if ~isempty(roi.missing)
                s = sprintf('%s  (asked for %s)', s, strjoin(roi.requested, ' '));
            end
        end

        function p = exploreProvenance(app, entry)
        % What the footer of an exported figure states. Enough to identify the
        % analysis six months later from the image alone: which plot, which
        % groups and how many subjects, which design, which electrodes.
        %
        % The group and design strings come from the same place the status bar
        % reads, so the line under the plot and the line stamped into a
        % published figure cannot disagree about n.
            [groups, design] = exploreCohortText(app);
            p = struct( ...
                'plot',      entry.name, ...
                'groups',    groups, ...
                'design',    design, ...
                'ROI',       exploreRoiText(app), ...
                'nestapp',   nestappVersion(), ...
                'exported',  char(datetime('now', 'Format', 'yyyy-MM-dd')));
        end

        function [groups, design] = exploreCohortText(app)
        % "odd n=5, even n=5" and "unpaired, 90% CI" - the two facts that say
        % what a picture is of. Written once because they appear both on screen
        % and in every exported figure's footer.
        %
        % The level is read, not asserted. It used to be a literal 95 here,
        % which was already only true by coincidence and would have become a
        % wrong number on a published figure the moment the level was
        % selectable.
            res  = app.exploreRes;
            bits = cell(1, numel(res.groups));
            for g = 1:numel(res.groups)
                bits{g} = sprintf('%s n=%d', res.groups(g).name, res.est(g).n);
            end
            groups = strjoin(bits, ', ');
            design = char(res.design);
            lvl    = exploreDrawnLevel(app);
            if ~isempty(lvl)
                design = sprintf('%s, %s', design, ciLabel(lvl));
            end
        end

        function lvl = exploreDrawnLevel(app)
        % The level the CURRENT plot draws its interval at, or [] for a plot
        % that draws no interval at all.
        %
        % A plot declares a 'level' param exactly when it has an interval, so
        % the registry already answers "is there one" and there is no second
        % list to keep in step. This also fixes a standing inaccuracy: a scalp
        % map's status line read "unpaired, 95% CI" while the picture contained
        % no interval whatsoever.
            lvl   = [];
            entry = currentPlotEntry(app);
            if isempty(entry) || isempty(entry.params); return; end
            if ~any(strcmp({entry.params.key}, 'level')); return; end

            % Unset falls back to what groupCurves recorded, not to a
            % literal. A result computed before res.info.level existed simply
            % has none, and then the phrase is left off rather than guessed.
            lvl = fieldOr(fieldOr(app.exploreRes, 'info', struct()), 'level', []);
            lvl = fieldOr(explorePlotParamsFor(app, entry.name), 'level', lvl);
        end

        function ExploreCsvButtonPushed(app, ~)
            if isempty(app.exploreRes); return; end
            T = exploreMeasures(withMode(app.exploreRes, currentMode(app)), ...
                                app.exploreWindows);
            if isempty(T) || height(T) == 0
                uialert(app.UIFigure, 'No measures to export yet.', 'Nothing to save');
                return
            end
            [f, p] = uiputfile({'*.csv', 'Comma-separated values'}, ...
                               'Save measures', 'tep_measures.csv');
            if isequal(f, 0); return; end
            writetable(T, fullfile(p, f));
            app.ExploreStatusLabel.Text = sprintf('Wrote %d rows to %s', height(T), f);
        end

        function ExploreResultsButtonPushed(app, ~)
            if isempty(app.exploreRes); return; end
            entry = currentPlotEntry(app);
            plotName = '';
            if ~isempty(entry); plotName = entry.name; end
            out = exploreResults(app.exploreRes, app.exploreEntries, struct( ...
                'roi', {app.exploreRoi}, 'windows', app.exploreWindows, ...
                'mode', currentMode(app), 'plot', plotName, ...
                'plotParams', app.explorePlotParams));

            choice = uiconfirm(app.UIFigure, ...
                ['The full result - curves at sampling rate, intervals and ' ...
                 'provenance - and everything needed to reopen this analysis ' ...
                 'later with File > Load Analysis. Save it as a .mat, or put ' ...
                 'it in the base workspace to carry on in MATLAB?'], 'Results', ...
                'Options', {'Save as .mat', 'To workspace', 'Cancel'}, ...
                'DefaultOption', 1, 'CancelOption', 3);
            switch choice
                case 'Save as .mat'
                    [f, p] = uiputfile({'*.mat', 'MATLAB data'}, ...
                                       'Save results', 'tep_results.mat');
                    if isequal(f, 0); return; end
                    save(fullfile(p, f), 'out');
                    app.ExploreStatusLabel.Text = sprintf('Saved results to %s', f);
                case 'To workspace'
                    assignin('base', 'tepResults', out);
                    app.ExploreStatusLabel.Text = ...
                        'Results assigned to "tepResults" in the base workspace.';
            end
        end

        function LoadAnalysisMenuSelected(app, ~)
        % Reopen a saved analysis. The Results .mat is the session format - it
        % already carried the files, ROI, windows, design and plot - so this is
        % the read side of an export that existed, not a new artifact.
            startFolder = getpref('nestapp', 'lastDataFolder', '');
            [f, p] = uigetfile({'*.mat', 'nestapp analysis (.mat)'}, ...
                               'Load analysis', startFolder);
            if isequal(f, 0); return; end
            applyExploreState(app, fullfile(p, f));
        end

        function applyExploreState(app, file)
        % Restore the tab from a saved analysis, then recompute from the files.
        %
        % Recomputed, not restored from the stored curves: a saved result holds
        % GROUP averages, which can redraw the figure and nothing else. Change
        % the ROI or move one recording between groups and those averages are
        % wrong. Resuming work means the per-file trial averages, so the files
        % are reloaded and everything downstream follows as if the groups had
        % just been assigned by hand.
            try
                loaded = load(file);
            catch ME
                uialert(app.UIFigure, ME.message, 'Could not read that file');
                return
            end
            [state, report] = exploreStateFromResults(loaded);
            if ~report.ok
                uialert(app.UIFigure, strjoin(report.notes, ' '), 'Not an analysis');
                return
            end

            paths = {state.entries.path};
            if isempty(paths)
                uialert(app.UIFigure, strjoin(report.notes, ' '), 'Nothing to load');
                return
            end

            % Same EEGLAB entry point as Add group: the recordings are read
            % back through pop_loadset.
            [ok, msg] = ensureEeglabReady();
            if ~ok
                uialert(app.UIFigure, msg, 'EEGLAB Not Ready', 'Icon', 'error');
                return
            end

            d = uiprogressdlg(app.UIFigure, 'Title', 'Loading analysis', ...
                'Message', sprintf('Reading %d recording%s...', ...
                                   numel(paths), plural(numel(paths))), ...
                'Indeterminate', 'on');
            closeDlg = onCleanup(@() delete(d));
            try
                cache = loadReducedSets(paths);
            catch ME
                uialert(app.UIFigure, ME.message, 'Could not load the recordings');
                return
            end

            app.exploreEntries     = state.entries;
            app.exploreCache       = cache;
            app.exploreRoi         = state.roi;
            app.exploreWindows     = state.windows;
            app.explorePlotParams  = state.plotParams;

            % Set the design before recomputing. refreshExploreDesign still has
            % the last word - a restored 'paired' that the files no longer
            % support flips back to unpaired rather than being taken on trust.
            if strcmp(state.design, 'paired')
                app.ExploreDesignGroup.SelectedObject = app.ExplorePairedButton;
            else
                app.ExploreDesignGroup.SelectedObject = app.ExploreUnpairedButton;
            end

            refreshExploreRoi(app);
            refreshExploreWindows(app);
            refreshExplorePlots(app);
            if ~isempty(state.plot) && ismember(state.plot, app.ExplorePlotDropDown.ItemsData)
                app.ExplorePlotDropDown.Value = state.plot;
            end
            recomputeExplore(app);

            app.TabGroup.SelectedTab = app.ExploreTab;
            [~, nm, ext] = fileparts(file);
            nGroups = numel(exploreGroupNames(app));
            msg = sprintf('Loaded %s - %d recording%s in %d group%s', [nm ext], ...
                numel(app.exploreEntries), plural(numel(app.exploreEntries)), ...
                nGroups, plural(nGroups));
            if ~isempty(report.notes)
                msg = [msg '  |  ' strjoin(report.notes, '  ')];
            end
            app.ExploreStatusLabel.Text = msg;
            if ~isempty(report.missing)
                uialert(app.UIFigure, strjoin(report.notes, newline), ...
                        'Analysis loaded', 'Icon', 'warning');
            end
        end

        function initEeglab(app)
        % INITEEGLAB  Bring EEGLAB up at launch, with feedback and a verdict.
        %   The plugin scan takes a couple of seconds and the window is already
        %   on screen by now, so it gets a progress dialog rather than a
        %   silent freeze. A failure is worth interrupting for: without EEGLAB
        %   the picker can only offer the handful of steps that need nothing,
        %   and the user would otherwise have to guess why.
            dlg = uiprogressdlg(app.UIFigure, 'Title', 'Starting nestapp', ...
                'Message', 'Initialising EEGLAB and its plugins...', ...
                'Indeterminate', 'on');
            closeDlg = onCleanup(@() close(dlg));
            [ok, msg] = ensureEeglabReady();
            if ~ok
                clear closeDlg
                uialert(app.UIFigure, msg, 'EEGLAB Not Ready', 'Icon', 'warning');
            end
        end

        % Selection changed function: StepsTree
        function StepsTreeSelectionChanged(app, ~)
            % Category and operation nodes have empty NodeData - they group, they
            % are not steps, so selecting one just clears the Info panel.
            name = selectedStepName(app);
            if isempty(name)
                app.InfoTextArea.Value = '';
            elseif isKey(app.info, name)
                app.InfoTextArea.Value = stepInfoLines(app, name);
            end
            app.selectedItem = [];
        end

        % Double-clicked callback: StepsTree
        function StepsTreeDoubleClicked(app, ~)
            % Double-clicking a step adds it; double-clicking a category/
            % operation header (empty NodeData) does nothing here and lets the
            % tree's native expand/collapse take over.
            name = selectedStepName(app);
            if isempty(name)
                return
            end
            appendStep(app, name);
        end

        % Button pushed function: AddButton
        function AddButtonPushed(app, ~)
            name = selectedStepName(app);
            if isempty(name)
                return   % a category/operation header is selected - nothing to add
            end
            appendStep(app, name);
        end

        % Button pushed function: MoveUpButton
        function MoveUpButtonPushed(app, ~)
            ind = selectedStepIndex(app);
            moveStep(app, ind, -1);
        end

        % Button pushed function: SavePipelineButton
        function SavePipelineButtonPushed(app, ~)
            startFolder = getpref('nestapp', 'lastPipelineFolder', '');
            [fName, fPath] = uiputfile('*.mat', 'Save Pipeline', ...
                fullfile(startFolder, 'pipeline.mat'));
            if isequal(fName, 0); return; end   % user cancelled
            spec         = app.spec;
            pipelineName = app.pipelineName;
            version      = '3';
            save(fullfile(fPath, fName), 'spec', 'pipelineName', 'version');
            setpref('nestapp', 'lastPipelineFolder', fPath);
            pushRecent(app, 'recentPipelines', fullfile(fPath, fName));
            [~, baseName, ~] = fileparts(fName);
            app.pipelineName  = baseName;
            app.pipelineDirty = false;
            updateStatusBar(app);
        end

        % Button pushed function: RemoveButton
        function RemoveButtonPushed(app, ~)
            ind = selectedStepIndex(app);
            removeStep(app, ind);
        end

        % Button pushed function: MoveDownButton
        function MoveDownButtonPushed(app, ~)
            ind = selectedStepIndex(app);
            moveStep(app, ind, +1);
        end

        % Button pushed function: LoadPipelineButton
        function LoadPipelineButtonPushed(app, ~)
            startFolder = getpref('nestapp', 'lastPipelineFolder', '');
            [pName, pPath] = uigetfile('*.mat', 'Load Pipeline', startFolder);
            if isequal(pName, 0); return; end
            fullPath = fullfile(pPath, pName);
            loadPipelineData(app, fullPath);
            setpref('nestapp', 'lastPipelineFolder', pPath);
            pushRecent(app, 'recentPipelines', fullPath);
            buildRecentPipelinesMenu(app);
            [~, nm, ~] = fileparts(pName);
            app.pipelineName  = nm;
            app.pipelineDirty = false;
            updateStatusBar(app);
        end

        % Button pushed function: SelectDataButton
        % Button pushed function: SelectDataButton ("Browse...")
        function SelectDataButtonPushed(app, ~)
        % SELECTDATABUTTONPUSHED  Open the data browser and queue checked files.
        %   selectDataTree shows the folder tree beneath a chosen parent as a
        %   checkbox tree with a path filter, so recordings spread across many
        %   subject folders (e.g. just the SPL files, or just session 1) can
        %   be picked in one place. Returns the full paths of every checked
        %   file; setFileQueue then adopts them as the run queue.
            startFolder = getpref('nestapp', 'lastDataFolder', '');
            paths = selectDataTree(startFolder, app.dataFileExts());
            if isempty(paths); return; end
            setFileQueue(app, paths);
            if isvalid(app.UIFigure); figure(app.UIFigure); end
        end

        function setFileQueue(app, fullPaths)
        % SETFILEQUEUE  Adopt a flat list of full file paths as the run queue.
        %   The single source of truth Run Analysis reads from. When every
        %   file shares one folder we keep app.path set and show basenames
        %   (so the TEP "use cleaned data" reuse and single-folder status
        %   bar keep working); when files span folders app.path is cleared
        %   and labels carry the parent-folder name to disambiguate.
            fullPaths = fullPaths(:)';
            if isempty(fullPaths); return; end

            [labels, uniqueParents, parents] = buildFileLabels(app, fullPaths);
            if isscalar(uniqueParents)
                app.path = uniqueParents{1};
            else
                app.path = '';   % files span folders - no single path
            end

            app.filePaths   = fullPaths;
            app.file        = labels;
            app.NSelecFiles = numel(fullPaths);
            app.SelectedFilesListBox.Items = labels;

            % Note: the "Data folder" pref (lastDataFolder) is updated by
            % selectDataTree to the browse root, not clobbered to a subfolder.
            pushRecent(app, 'recentFiles', parents{1});
            buildRecentFilesMenu(app);
            updateStatusBar(app);

            fprintf('[nestapp] Queued %d file(s) from %d folder(s).\n', ...
                numel(fullPaths), numel(uniqueParents));
        end

        function exts = dataFileExts(~)
        % DATAFILEEXTS  Glob patterns for the data formats nestapp loads.
            exts = {'*.set', '*.vhdr', '*.cdt', '*.cnt'};
        end

        function [labels, uniqueParents, parents] = buildFileLabels(~, fullPaths)
        % BUILDFILELABELS  Display labels for a flat list of full file paths.
        %   Basenames when every file shares one folder; parentFolder/basename
        %   when the list spans folders, so same-named files from different
        %   subjects stay distinguishable. Also returns the unique parent
        %   folders (callers use the count / common path) and each file's
        %   parent. Shared by the Cleaning (setFileQueue) and Visualize
        %   (setTEPFileList) tabs so the labelling rule lives in one place.
            [parents, names, exts] = cellfun(@fileparts, fullPaths, ...
                'UniformOutput', false);
            uniqueParents = unique(parents);
            labels = cell(1, numel(fullPaths));
            if isscalar(uniqueParents)
                for i = 1:numel(fullPaths)
                    labels{i} = [names{i} exts{i}];
                end
            else
                for i = 1:numel(fullPaths)
                    [~, folderName] = fileparts(parents{i});
                    labels{i}       = [folderName '/' names{i} exts{i}];
                end
            end
        end

        % Button pushed function: ReStartStepsButton
        function ReStartStepsButtonPushed(app, ~)
            confirmClear = getpref('nestapp', 'confirmClear', true);
            if confirmClear
                answer = uiconfirm(app.UIFigure, ...
                    'Clear all pipeline steps? This cannot be undone.', ...
                    'Clear Pipeline', ...
                    'Options', {'Clear', 'Cancel'}, ...
                    'DefaultOption', 2, 'CancelOption', 2);
                if strcmp(answer, 'Cancel'); return; end
            end
            clc
            clearSteps(app);
        end

        % Menu selected function: Load Template...
        function LoadTemplateMenuSelected(app, ~)
        % LOADTEMPLATEMENUSELECTED  Show a template picker and load the chosen template.
        %   Reads template .mat files from src/templates/ - the same format
        %   as user-saved pipelines.  No override logic runs at runtime.
            % which('nestapp') points at the class folder
            % (src/@nestapp/nestapp.m); templates live one directory up
            % under src/templates/.
            classDir    = fileparts(which('nestapp'));
            templateDir = fullfile(fileparts(classDir), 'templates');
            files = dir(fullfile(templateDir, '*.mat'));
            if isempty(files)
                uialert(app.UIFigure, ...
                    'No template files found in src/templates/.  Run buildTemplates() to generate them.', ...
                    'Templates');
                return
            end

            % Read pipelineName from each file for the picker list.
            n     = numel(files);
            names = cell(n, 1);
            paths = cell(n, 1);
            for i = 1:n
                paths{i} = fullfile(files(i).folder, files(i).name);
                try
                    tmp = load(paths{i});
                    if isfield(tmp, 'pipelineName') && ~isempty(tmp.pipelineName)
                        names{i} = tmp.pipelineName;
                    else
                        [~, names{i}] = fileparts(files(i).name);
                    end
                catch
                    [~, names{i}] = fileparts(files(i).name);
                end
            end

            % Modal picker
            dlg = uifigure('Name', 'Load Template', ...
                'Position', [300 300 320 200], ...
                'WindowStyle', 'modal', 'Resize', 'off');
            uilabel(dlg, 'Text', 'Select a pipeline template:', ...
                'Position', [15 165 290 22]);
            lb = uilistbox(dlg, 'Items', names, ...
                'Position', [15 60 290 100], 'Value', names{1});
            uibutton(dlg, 'Text', 'Cancel', 'Position', [120 15 85 30], ...
                'ButtonPushedFcn', @(~,~) close(dlg));
            uibutton(dlg, 'Text', 'Load', 'Position', [215 15 90 30], ...
                'BackgroundColor', [0.20 0.55 0.20], 'FontColor', [1 1 1], ...
                'ButtonPushedFcn', @(~,~) doLoad());
            uiwait(dlg);

            % Closing the modal picker can leave the main window buried behind
            % other apps (Chrome, etc.); raise nestapp so it isn't lost.
            if isvalid(app.UIFigure)
                figure(app.UIFigure);
            end

            function doLoad()
                idx = find(strcmp(names, lb.Value), 1);
                close(dlg);
                if isempty(idx); return; end
                clearSteps(app);
                clc
                loadPipelineData(app, paths{idx});
                app.pipelineName  = names{idx};
                app.pipelineDirty = true;
                updateStatusBar(app);
                if ~isempty(app.SelectedListBox.Items)
                    app.SelectedListBox.Value = app.SelectedListBox.ItemsData{1};
                    refreshParamTable(app, 1);
                end
                promptForRequiredParams(app);
            end
        end

        function promptForRequiredParams(app)
        % Some templates cannot carry a working default for every setting -
        % AARATEP needs an output folder and there is no reasonable guess. Ask
        % on load, with that step's parameter table already in front of the
        % user, rather than letting them start a run that fails on the first
        % file for a reason they then have to go and find.
        %
        % Which parameters are required is declared on the parameter itself
        % (makeParam(..., 'required', true)) and found by
        % unsetRequiredParams, so this does not need updating when a step
        % gains one.
            [stepIdx, ~, labels] = unsetRequiredParams(app.spec);
            if isempty(stepIdx); return; end

            % Put the offending step's parameters on screen first, so the
            % dialog is pointing at something the user can actually edit.
            if numel(app.SelectedListBox.ItemsData) >= stepIdx
                app.SelectedListBox.Value = app.SelectedListBox.ItemsData{stepIdx};
                refreshParamTable(app, stepIdx);
            end

            uialert(app.UIFigure, sprintf( ...
                ['"%s" needs these before the pipeline can run:\n\n    %s\n\n' ...
                 'Its parameters are shown below - fill them in, then Run.'], ...
                app.spec(stepIdx).name, strjoin(labels, sprintf('\n    '))), ...
                'Settings needed', 'Icon', 'info');
        end

        % Button pushed function: RunAnalysisButton
        function RunAnalysisButtonPushed(app, ~)
            app.RunAnalysisButton.Text = {'Run';'Analysis'};
            if isempty(app.file)
                uialert(app.UIFigure, 'Please select at least one data file.', '');
                return
            end
            if isempty(app.spec)
                uialert(app.UIFigure, 'Please add at least one pipeline step.', '');
                return
            end

            % Normally a no-op - startupFcn has already done this - but the
            % run is the point at which EEGLAB stops being optional, so check
            % rather than assume the session still has it.
            [ok, msg] = ensureEeglabReady();
            if ~ok
                uialert(app.UIFigure, msg, 'EEGLAB Init Failed', 'Icon', 'error');
                return
            end

            % One definition of what the queue means, shared with Browse EEG -
            % otherwise the two could disagree about which files are selected.
            queuePaths = cleaningQueuePaths(app);

            % Pre-select channel location file once if the pipeline needs it.
            app.preSelectedChanFile = '';
            for psi = 1:numel(app.spec)
                if strcmp(app.spec(psi).name, 'Load Channel Location')
                    p = app.spec(psi).params;
                    needChan     = isfield(p, 'needchanloc') && strcmp(p.needchanloc, 'yes');
                    eachFileDiff = isfield(p, 'eachFilediffPath') && strcmp(p.eachFilediffPath, 'yes');
                    if needChan && ~eachFileDiff
                        [chName, chPath] = uigetfile('*.*', 'Select channel location file');
                        if isequal(chName, 0); return; end
                        app.preSelectedChanFile = fullfile(chPath, chName);
                    end
                    break
                end
            end

            % Steps that wait for a human cannot run on a parallel worker:
            % workers get uiFigure = [] and no display, so the run errors on
            % every file or blocks with nothing to click - and only after the
            % batch has started. Warn and offer to turn parallel off first.
            if app.ParallelCheckBox.Value
                blocking = interactivePipelineSteps(app.spec);
                if ~isempty(blocking)
                    choice = uiconfirm(app.UIFigure, sprintf([ ...
                        'This pipeline contains %d step(s) that wait for you:\n\n' ...
                        '    %s\n\n' ...
                        'They cannot run on a parallel worker, so the run would ' ...
                        'fail or hang. Turn parallel processing off and run the ' ...
                        'files one at a time?'], ...
                        numel(blocking), strjoin(blocking, sprintf('\n    '))), ...
                        'Interactive steps in pipeline', ...
                        'Options', {'Turn off parallel', 'Cancel'}, ...
                        'DefaultOption', 1, 'CancelOption', 2, 'Icon', 'warning');
                    if ~strcmp(choice, 'Turn off parallel'); return; end
                    app.ParallelCheckBox.Value = false;
                end
            end

            opts.uiFigure     = app.UIFigure;
            opts.pipelineName = app.pipelineName;
            opts.statusBar    = app.StatusBar;
            opts.parallel     = app.ParallelCheckBox.Value;
            opts.chanLocFile  = app.preSelectedChanFile;

            try
                [allReports, allSummaries, failed] = runPipelineCore(app.spec, queuePaths, opts);
            catch err
                if strcmp(err.identifier, 'nestapp:cancelled')
                    return   % backed out before anything ran - nothing to say
                end
                if strcmp(err.identifier, 'nestapp:cancelledPartial')
                    % The run got far enough to write per-file reports and a
                    % partial session summary; the message names the folder.
                    uialert(app.UIFigure, err.message, 'Run Cancelled', 'Icon', 'warning');
                    return
                end
                uialert(app.UIFigure, err.message, 'Pipeline Error', 'Icon', 'error');
                return
            end

            % Keep this run's failures so the dashboard + summary can show
            % them (they have no report of their own).
            app.lastFailed = failed;

            % A single-file run gets a session summary too - the methods
            % paragraph and citations are just as useful for one file.
            if ~isempty(allReports)
                summEntry.text      = summarizeReports(allReports, failed);
                summEntry.report    = [];
                summEntry.isSummary = true;
                app.allPipelineReports{end+1} = summEntry;
            end
            for ri = 1:numel(allSummaries)
                entry.text      = allSummaries{ri};
                entry.report    = allReports{ri};
                entry.isSummary = false;
                app.allPipelineReports{end+1} = entry;
            end
            updateReportsTab(app);
            if getpref('nestapp', 'showReport', true) && ~isempty(allSummaries)
                app.TabGroup.SelectedTab = app.ReportsTab;
            end
        end

        % Value changed function: TextArea
        function TextAreaValueChanged(app, ~)
            if isempty(app.currentParamKey); return; end
            stepIdx = selectedStepIndex(app);
            if isempty(stepIdx) || stepIdx > numel(app.spec); return; end

            raw = app.TextArea.Value;
            if iscell(raw) || isstring(raw)
                raw = strjoin(raw, ' ');
            end

            val = convertParam(raw, app.currentParamType);
            app.spec(stepIdx).params.(app.currentParamKey) = val;
            refreshParamTable(app, stepIdx);
            app.pipelineDirty = true;
            updateStatusBar(app);
        end

        % Cell selection callback: UITable
        function UITableCellSelection(app, event)
            if isempty(event.Indices); return; end
            row     = event.Indices(1);
            stepIdx = selectedStepIndex(app);
            if isempty(stepIdx) || stepIdx > numel(app.spec); return; end

            reg    = stepRegistry();
            step   = app.spec(stepIdx);
            regIdx = find(strcmp({reg.name}, step.name), 1);
            if isempty(regIdx); return; end
            params = reg(regIdx).params;
            if row > numel(params); return; end

            paramMeta = params(row);
            app.currentParamKey  = paramMeta.key;
            app.currentParamType = paramMeta.type;

            val = step.params.(paramMeta.key);
            if isnumeric(val)
                app.TextArea.Value = num2str(val);
            elseif iscell(val)
                % Cell params may mix char and numeric entries (e.g. SSP-SIR
                % PC = {'data', 90}); strjoin needs char/string, so stringify
                % each element first instead of erroring on the numeric one.
                parts = cell(1, numel(val));
                for ii = 1:numel(val)
                    v = val{ii};
                    if ischar(v)
                        parts{ii} = v;
                    elseif isstring(v)
                        parts{ii} = char(v);
                    elseif isnumeric(v) || islogical(v)
                        parts{ii} = num2str(v);
                    else
                        parts{ii} = char(string(v));
                    end
                end
                app.TextArea.Value = strjoin(parts, ', ');
            elseif ischar(val) && ~isrow(val) && ~isempty(val)
                app.TextArea.Value = cellstr(val);   % char matrix -> cell for TextArea
            else
                app.TextArea.Value = char(val);
            end
        end

        % Button pushed function: DefaultValueButton
        function DefaultValueButtonPushed(app, ~)
            stepIdx = selectedStepIndex(app);
            if isempty(stepIdx) || stepIdx > numel(app.spec); return; end
            reg  = stepRegistry();
            name = app.spec(stepIdx).name;
            app.spec(stepIdx) = makePipelineStep(name, reg);
            app.currentParamKey  = '';
            app.currentParamType = '';
            app.TextArea.Value   = '';
            refreshParamTable(app, stepIdx);
            app.pipelineDirty = true;
            updateStatusBar(app);
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, ~)
            pool = gcp('nocreate');
            if ~isempty(pool)
                nestLog('PAR', 'nestapp closing - shutting down parallel pool (%d workers)', ...
                    pool.NumWorkers);
                delete(pool);
            end
            delete(app);
        end

        % Size changed function: UIFigure
        function UIFigureSizeChanged(app, ~)
            if isempty(app.originalSize); return; end
            % enforceMinWindowSize writes Position, which fires this callback
            % again - and drawnow below lets that re-entry actually run. Left
            % unguarded the two feed each other for as long as resize events
            % keep arriving.
            if app.isResizing; return; end
            app.isResizing = true;
            done = onCleanup(@() endResize(app));

            drawnow limitrate  % throttle: skip redraws that arrive faster than screen refresh
            newSize = enforceMinWindowSize(app);
            sX = newSize(1) / app.originalSize(1);
            sY = newSize(2) / app.originalSize(2);
            rescaleComponents(app, sX, sY);
            reRenderReportsOnResize(app);  % reflow the Reports pane if one is showing
            reRenderExploreOnResize(app);  % and the Explore plot, for the same reason
        end

        function endResize(app)
            app.isResizing = false;
        end

        function newSize = enforceMinWindowSize(app)
        % Grow the window back to the minimum WITHOUT moving it on screen.
        % clampWindowPosition owns the arithmetic - and the explanation of why
        % the top edge has to be pinned - so it can be tested without launching
        % the app. This method is the part that genuinely needs a figure.
            MIN_SIZE = [650 420];
            [newPos, changed] = clampWindowPosition(app.UIFigure.Position, MIN_SIZE);
            newSize = newPos(3:4);
            if changed
                app.UIFigure.Position = newPos;
            end
        end

        % Cell edit callback: UITable
        function UITableCellEdit(app, event)
            stepIdx = selectedStepIndex(app);
            if isempty(stepIdx) || stepIdx > numel(app.spec); return; end

            reg    = stepRegistry();
            step   = app.spec(stepIdx);
            regIdx = find(strcmp({reg.name}, step.name), 1);
            if isempty(regIdx); return; end

            row    = event.Indices(1);
            params = reg(regIdx).params;
            % Reject edits to a parameter that is disabled for this step -
            % overridden by a mutually-exclusive sibling, or gated off by
            % another param's value. Revert the cell and explain why.
            if row <= numel(params)
                [disabled, reasons] = disabledParamKeys(reg(regIdx), step.params);
                key = params(row).key;
                if any(strcmp(key, disabled))
                    uialert(app.UIFigure, reasons.(key), 'Parameter disabled');
                    refreshParamTable(app, stepIdx);   % revert the edited cell
                    return
                end
            end

            app.spec = applyParamEdit(app.spec, stepIdx, row, event.NewData, reg(regIdx));
            app.pipelineDirty = true;
            updateStatusBar(app);
            % Re-render so setting/clearing one mutually-exclusive param updates
            % the disabled state of its sibling.
            refreshParamTable(app, stepIdx);
        end

        % Value changed function: SelectedListBox
        function SelectedListBoxValueChanged(app, ~)
            stepIdx = selectedStepIndex(app);
            if isempty(stepIdx) || stepIdx > numel(app.spec); return; end
            app.currentParamKey  = '';
            app.currentParamType = '';
            app.TextArea.Value   = '';
            refreshParamTable(app, stepIdx);
        end


        % Menu selected function: Tools > Browse Raw EEG...
        function browseRawEegMenu(app, ~)
        % Open EEGLAB's scrolling viewer on a recording - raw or cleaned.
        %
        % It offers what the session already knows about (the Cleaning queue,
        % and the cleaned output of the report selected on the Reports tab)
        % and always offers the file browser as well, because the batch you
        % want to look at is not always the batch this session ran. Nothing
        % here depends on a cohort someone else happened to load first.
            [paths, labels] = browsableRecordings(app);

            k = numel(paths) + 1;          % nothing queued -> straight to the browser
            if ~isempty(paths)
                k = pickOne('Browse EEG', 'Which recording?', ...
                            [labels, {'Browse for a file...'}], app.UIFigure);
                if isempty(k); return; end
            end
            if k > numel(paths)
                chosen = pickSetFromDisk(app);
            else
                chosen = paths{k};
            end
            if isempty(chosen); return; end
            openInEegplot(app, chosen);
        end

        function [paths, labels] = browsableRecordings(app)
        % Everything this session can name: the Cleaning queue, plus the
        % cleaned output of the report selected on Reports. Labelled by where
        % they came from, because "sub-01_t1" as an input and as a cleaned
        % output are different files with nearly the same name.
            paths  = cleaningQueuePaths(app);
            labels = {};
            if ~isempty(paths)
                labels = buildFileLabels(app, paths);
                labels = cellfun(@(s) ['queued:  ' s], labels, 'UniformOutput', false);
            end

            out = selectedReportOutput(app);
            if ~isempty(out) && ~any(strcmpi(paths, out))
                [~, b, e] = fileparts(out);
                paths{end+1}  = out;
                labels{end+1} = ['cleaned: ' b e];
            end
        end

        function p = selectedReportOutput(app, allEntries)
        % Where the report selected on the Reports tab wrote its cleaned .set,
        % or '' when nothing is selected, the pipeline never saved, or the
        % report predates the outputFile field.
        %
        % Takes the entry list rather than rebuilding it: renderReportsRightPane
        % already holds one, and this runs on every selection change.
            if nargin < 2; allEntries = [app.allPipelineReports, app.loadedReports]; end
            p = '';
            r = selectedFileReport(app, allEntries);
            if isempty(r) || ~isfield(r, 'outputFile'); return; end
            if ~isempty(r.outputFile) && isfile(r.outputFile)
                p = r.outputFile;
            end
        end

        function p = pickSetFromDisk(app)
        % A plain file browser, opening where the user last worked.
            p = '';
            start = getpref('nestapp', 'lastDataFolder', '');
            if isempty(start) || ~isfolder(start); start = pwd; end
            exts = dataFileExts(app);
            [f, d] = uigetfile( ...
                {strjoin(exts, ';'), sprintf('EEG recordings (%s)', strjoin(exts, ', ')); ...
                 '*.*', 'All files'}, ...
                'Open EEG recording', start);
            if isvalid(app.UIFigure); figure(app.UIFigure); end
            if isequal(f, 0); return; end
            p = fullfile(d, f);
        end

        function openInEegplot(app, chosen)
        % Load one recording and hand it to EEGLAB's scrolling viewer.
            % The viewer is EEGLAB's, so this is the point at which EEGLAB
            % stops being optional.
            [ok, msg] = ensureEeglabReady();
            if ~ok
                uialert(app.UIFigure, msg, 'EEGLAB Init Failed', 'Icon', 'error');
                return
            end

            try
                EEG = loadEegFile(chosen);
            catch err
                uialert(app.UIFigure, ...
                    sprintf('Could not read %s\n\n%s', chosen, err.message), ...
                    'Load failed', 'Icon', 'error');
                return
            end
            pop_eegplot(EEG, 1, 1, 1);
        end

        function paths = cleaningQueuePaths(app)
        % The Cleaning tab's queue as full paths - the one definition of what
        % "the selected files" means, used by both Run Analysis and Browse EEG.
        %
        % app.filePaths is the source of truth and spans folders when the
        % "Folders..." picker was used; the app.path join is the fallback for a
        % selection made before filePaths was populated.
            if ~isempty(app.filePaths)
                paths = app.filePaths(:)';
            elseif ~isempty(app.file) && ~isempty(app.path)
                paths = cellfun(@(f) fullfile(app.path, f), app.file(:)', ...
                                'UniformOutput', false);
            else
                paths = {};
            end
        end

        % -- Analysis Tab callbacks ----------------------------------------

        % Selection changed function: TabGroup
        function TabGroupSelectionChanged(app, event)
        % Refresh a tab's contents when it becomes active.
            if event.NewValue == app.ExploreTab
                % A resize while another tab was showing is skipped as wasted
                % work, so the canvas may be sized for the old window. Repaint
                % on the way in.
                renderExplorePlot(app);
            end
        end

    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components - body in @nestapp/createComponents.m
        createComponents(app)
    end

    % Public methods callable from external functions (e.g. runPipelineCore)
    methods (Access = public)

        function updateReportsTab(app)
        % UPDATEREPORTSTAB  Public entry point - refreshes the Reports tab.
        %   Delegates to the private implementation. Exposed as public so
        %   runPipelineCore.m can call it after each processing run.
            updateReportsTabImpl(app);
        end

        function loadAnalysis(app, file)
        % LOADANALYSIS  Public entry point - reopen a saved analysis .mat.
        %   loadAnalysis(app, file) is File > Load Analysis without the file
        %   picker. Delegates to the private implementation, the same pattern
        %   updateReportsTab follows.
        %
        %   Public because reopening an analysis is a scriptable operation, not
        %   only a menu click: a batch that regenerates every figure for a
        %   manuscript wants to load each saved .mat in turn without a dialog.
            applyExploreState(app, file);
        end

    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = nestapp

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.UIFigure)

                % Execute the startup function
                runStartupFcn(app, @startupFcn)
            else

                % Focus the running singleton app
                figure(runningApp.UIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)
            % Stop the timers first: a pending callback would otherwise fire
            % against a half-deleted app, and an undeleted timer outlives the
            % app in the MATLAB session, firing into nothing forever.
            for t = [app.hoverTimer, app.exploreResizeTimer, app.reportsResizeTimer]
                if isvalid(t); stop(t); delete(t); end
            end

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end

function s = onOffState(tf)
if tf; s = 'on'; else; s = 'off'; end
end

function labels = firstRoiPreset()
p = roiPresets();
labels = {};
if ~isempty(p); labels = p(1).labels; end
end

function res = withMode(res, mode)
res.mode = mode;
end
