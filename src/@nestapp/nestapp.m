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
        StepsListBox                    matlab.ui.control.ListBox
        VisualizingTab                  matlab.ui.container.Tab
        TEPvarNameEditField             matlab.ui.control.EditField
        TEPvarNameEditFieldLabel        matlab.ui.control.Label
        ExportTEPDataButton             matlab.ui.control.Button
        PlotEEGdataButton               matlab.ui.control.Button
        EEGDatasetDropDown              matlab.ui.control.DropDown
        EEGDatasetDropDownLabel         matlab.ui.control.Label
        TopoplottimeSpinner             matlab.ui.control.Spinner
        TopoplottimeSpinnerLabel        matlab.ui.control.Label
        TEPWindowSlider                 matlab.ui.control.RangeSlider
        TEPWindowSliderLabel            matlab.ui.control.Label
        ReLoadAvailableElectrodesButton  matlab.ui.control.Button
        PO6Button                       matlab.ui.control.StateButton
        PO1Button                       matlab.ui.control.StateButton
        DontfindcommonelectrodesCheckBox  matlab.ui.control.CheckBox
        SelectAllCheckBox               matlab.ui.control.CheckBox
        AF8Button                       matlab.ui.control.StateButton
        AF7Button                       matlab.ui.control.StateButton
        AFZButton                       matlab.ui.control.StateButton
        TP9Button                       matlab.ui.control.StateButton
        TP10Button                      matlab.ui.control.StateButton
        CB2Button                       matlab.ui.control.StateButton
        O2Button                        matlab.ui.control.StateButton
        OZButton                        matlab.ui.control.StateButton
        CB1Button                       matlab.ui.control.StateButton
        PO8Button                       matlab.ui.control.StateButton
        PO2Button                       matlab.ui.control.StateButton
        PO5Button                       matlab.ui.control.StateButton
        PO7Button                       matlab.ui.control.StateButton
        PO4Button                       matlab.ui.control.StateButton
        POZButton                       matlab.ui.control.StateButton
        PO3Button                       matlab.ui.control.StateButton
        O1Button                        matlab.ui.control.StateButton
        P6Button                        matlab.ui.control.StateButton
        P4Button                        matlab.ui.control.StateButton
        PZButton                        matlab.ui.control.StateButton
        P5Button                        matlab.ui.control.StateButton
        P7Button                        matlab.ui.control.StateButton
        P2Button                        matlab.ui.control.StateButton
        P1Button                        matlab.ui.control.StateButton
        P3Button                        matlab.ui.control.StateButton
        P8Button                        matlab.ui.control.StateButton
        TP8Button                       matlab.ui.control.StateButton
        CP6Button                       matlab.ui.control.StateButton
        CP4Button                       matlab.ui.control.StateButton
        CPZButton                       matlab.ui.control.StateButton
        CP5Button                       matlab.ui.control.StateButton
        TP7Button                       matlab.ui.control.StateButton
        T7Button                        matlab.ui.control.StateButton
        C3Button                        matlab.ui.control.StateButton
        C5Button                        matlab.ui.control.StateButton
        FC5Button                       matlab.ui.control.StateButton
        FT7Button                       matlab.ui.control.StateButton
        T8Button                        matlab.ui.control.StateButton
        CP2Button                       matlab.ui.control.StateButton
        CP1Button                       matlab.ui.control.StateButton
        CP3Button                       matlab.ui.control.StateButton
        C2Button                        matlab.ui.control.StateButton
        CZButton                        matlab.ui.control.StateButton
        C1Button                        matlab.ui.control.StateButton
        FC3Button                       matlab.ui.control.StateButton
        FCZButton                       matlab.ui.control.StateButton
        FC1Button                       matlab.ui.control.StateButton
        F7Button                        matlab.ui.control.StateButton
        FT8Button                       matlab.ui.control.StateButton
        C6Button                        matlab.ui.control.StateButton
        C4Button                        matlab.ui.control.StateButton
        F1Button                        matlab.ui.control.StateButton
        FC6Button                       matlab.ui.control.StateButton
        FC4Button                       matlab.ui.control.StateButton
        FC2Button                       matlab.ui.control.StateButton
        FZButton                        matlab.ui.control.StateButton
        F3Button                        matlab.ui.control.StateButton
        F5Button                        matlab.ui.control.StateButton
        F2Button                        matlab.ui.control.StateButton
        F4Button                        matlab.ui.control.StateButton
        F6Button                        matlab.ui.control.StateButton
        F8Button                        matlab.ui.control.StateButton
        AF4Button                       matlab.ui.control.StateButton
        FP2Button                       matlab.ui.control.StateButton
        FPZButton                       matlab.ui.control.StateButton
        FP1Button                       matlab.ui.control.StateButton
        AF3Button                       matlab.ui.control.StateButton
        PlottingModeButtonGroup         matlab.ui.container.ButtonGroup
        AddtocurrentFigureButton        matlab.ui.control.RadioButton
        NewFigureButton                 matlab.ui.control.RadioButton
        ExportTEPFigureButton           matlab.ui.control.Button
        TOPOPLOTButton                  matlab.ui.control.Button
        WindowsizefortimeaveragedTopoplotEditField  matlab.ui.control.NumericEditField
        WindowsizeforTopoplotLabel      matlab.ui.control.Label
        Slider                          matlab.ui.control.Slider
        Image2                          matlab.ui.control.Image
        FilesListBox                    matlab.ui.control.ListBox
        FilesListBoxLabel               matlab.ui.control.Label
        UseCurrentlyCleanedDataCheckBox  matlab.ui.control.CheckBox
        SelectDatatoVisulaizeTEPsPanel  matlab.ui.container.Panel
        SelectDataButton_2              matlab.ui.control.Button
        FolderEditField_2               matlab.ui.control.EditField
        FolderEditField_2Label          matlab.ui.control.Label
        PLOTTEPButton                   matlab.ui.control.Button
        ShowComponentsButton            matlab.ui.control.StateButton
        EditComponentWindowsButton      matlab.ui.control.Button
        TEPComponentTable               matlab.ui.control.Table
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes                          matlab.ui.control.UIAxes
        ReportsTab                      matlab.ui.container.Tab
        ReportsListBox                  matlab.ui.control.ListBox
        ReportsListBoxLabel             matlab.ui.control.Label
        LoadReportsButton               matlab.ui.control.Button
        RefreshReportsButton            matlab.ui.control.Button
        ReportsFolderLabel              matlab.ui.control.Label
        ReportsStatusLabel              matlab.ui.control.Label
        ReportsTextArea                 matlab.ui.control.TextArea
        ReportsDashboardPanel           matlab.ui.container.Panel
        ExportReportsCSVButton          matlab.ui.control.Button
        ExportPDFButton                 matlab.ui.control.Button
        CopyMethodsButton               matlab.ui.control.Button
        % Analysis tab - static elements not auto-resized by MATLAB
        AnalysisSelPanel                matlab.ui.container.Panel
        AnalysisCompWindowsLabel        matlab.ui.control.Label
        AnalysisWorkspaceLabel          matlab.ui.control.Label
        AnalysisBatchLabel              matlab.ui.control.Label
        AnalysisBatchDescLabel          matlab.ui.control.Label
    end

    properties (Access = private)
        ItemNum % Index for selected Item
        elecList = {'FPz','FP1','FP2','AF7','AF3','AFz','AF4','AF8','F7','F5','F3',...
                'F1','F2','F4','F6','F8','Fz','FT7','FT8','FC5','FC3',...
                'FC1','FCz','FC2','FC4','FC6','T7','T8','C5','C3','C1','Cz',...
                'C2','C4','C6','TP7','TP8','CP5','CP3','CP1','CPz',...
                'CP2','CP4','CP6','P7','P5','P3','P1','Pz',...
                'P2','P4','P6','P8','PO7','PO5','PO3','PO1','POz','PO2','PO4','PO6','PO8',...
                'CB1','O1','Oz','O2','CB2','TP9','TP10'}; % All Listed Electrodes
        path         % File Path
        file         % File Name
        spec         % PipelineStep struct array (name + typed params)
        NSelecFiles  % Number of selected files for EEG preprocessing
        cleanedName  % Name used to rename the saved cleaned EEG data
    end
    properties (Access = public)
        % Tab Cleaning
        selectedItem % Selected Table Item Values
        info % Command Information and description
        % Canonical pipeline state - single source of truth for steps and params.
        % appendStep/removeStep/moveStep/clearSteps/loadPipelineData all write here.
        currentParamKey  = ''  % param key selected in UITable (transient)
        currentParamType = ''  % type of selected param (transient)
        originalSize     % [w h] of UIFigure at creation - used by UIFigureSizeChanged
        clickedItem = [];
        doubleClicked = 0;
        TEPfiles % File list for calculating TEPs
        
        % Tab Visualizing
        PathofSelectedFilesforTEP
        NumberOfSelecFilesforTEP
        SelectedFilesforTEP % Selected files to plot the TEP
        Common_Labels % Commong electrod name among files
        ROIelecsLabels % Selected electrodes as Region of Interest
        TEPCreated = false; % true once the TEP plot has been rendered at least once
        EEG_SelectedTEPFiles_Loaded = false;
        EEGofAllSelectedFiles = [];
        DefaulTEPxLim = [-50 300]; % Default xLim for time in TEP
        EEGtime
        TEP2Export
        MenuRecentFiles     % Handle to 'Recent Files' submenu - rebuilt on open
        MenuRecentPipelines % Handle to 'Recent Pipelines' submenu - rebuilt on open
        StatusBar           % uilabel pinned to bottom of UIFigure - visible on both tabs
        pipelineDirty   = false    % true when pipeline has unsaved changes
        pipelineName    = ''       % filename of last saved/loaded pipeline
        lastStepClick   = NaT     % datetime of last StepsListBox click (double-click detection)
        tepPeaks        = struct([]) % struct array from tepPeakFinder; cached after each PLOT TEP
        tepComponentDefs = struct([]) % component window definitions used by tepPeakFinder
        allPipelineReports = {}    % cell array of report entry structs from current session
        loadedReports      = {}    % cell array of report entry structs loaded from disk
        preSelectedChanFile = ''   % channel location file selected once before a run
        ParallelCheckBox           % uicheckbox - enable parallel participant processing

        % Tab Analysis
        AnalysisTab
        ExtractPeaksCSVButton
        AnalysisStatusLabel
        AnalysisSelectionLabel
    end

    methods (Access = private)
        % -- Pipeline state mutation methods ---------------------------------
        % app.spec, SelectedListBox.Items, and SelectedListBox.ItemsData must
        % stay in sync. These methods are the ONLY permitted way to add,
        % remove, move, or clear steps -- callbacks delegate here.

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
            app.UITable.Data = buildParamTableData(step, reg(regIdx));
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
            if ~ischar(app.path) && ~isstring(app.path) || isempty(app.path)
                dataStr = 'Data: (none)';
            else
                parts = strsplit(strtrim(app.path), {'\','/'});
                parts(cellfun(@isempty, parts)) = [];
                folder = parts{end};
                n = app.NSelecFiles;
                if isempty(n) || n == 0; n = 0; end
                fileWord = 'files'; if n == 1; fileWord = 'file'; end
                dataStr = sprintf('Data: %s/  (%d %s)', folder, n, fileWord);
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
        % Load data files from a recently used folder.
            if ~isfolder(folder)
                uialert(app.UIFigure, 'Folder no longer exists.', 'Not Found');
                return
            end
            try
                [app.file, app.path] = uigetfile( ...
                    {'*.set;*.vhdr;*.cdt;*.cnt','Data Files'}, ...
                    'Select File(s)', folder, 'multiSelect', 'on');
                if isequal(app.file, 0); return; end
                if ~iscell(app.file)
                    app.NSelecFiles = 1;
                    app.file = {app.file};
                else
                    app.NSelecFiles = numel(app.file);
                end
                app.SelectedFilesListBox.Items = app.file;
                setpref('nestapp', 'lastDataFolder', app.path);
                pushRecent(app, 'recentFiles', app.path);
                buildRecentFilesMenu(app);
            catch
                warning('nestapp: could not open data from recent folder.');
            end
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

        function selfTestMenu(app, ~)
        % SELFTESTMENU  Help action: run the fast test suite to verify the
        %   install, reporting pass/fail. Best-effort: needs tests/ present.
            repo = fileparts(fileparts(which('nestappVersion')));
            runner = fullfile(repo, 'tests', 'run_tests.m');
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
                'Position', [200 200 420 640], ...
                'WindowStyle', 'modal', 'Resize', 'off');

            % --- Quality Screening section (new, at top) ---
            uilabel(dlg, 'Text', 'Quality Screening', 'FontWeight', 'bold', ...
                'Position', [15 605 200 20]);
            cbAutoQC = uicheckbox(dlg, 'Text', 'Auto-generate QC images at each Quality Gate', ...
                'Position', [15 580 380 22], ...
                'Value', getpref('nestapp', 'autoQualityReport', false));
            cbTmsAuto = uicheckbox(dlg, 'Text', 'Auto-detect TMS pulse window from EEG events', ...
                'Position', [15 558 380 22], ...
                'Value', getpref('nestapp', 'qualityTmsAutoDetect', true));
            cbSkipFail = uicheckbox(dlg, 'Text', 'Skip remaining pipeline steps when Quality Gate fails', ...
                'Position', [15 536 380 22], ...
                'Value', getpref('nestapp', 'skipOnQualityFail', false));
            cbAutoPDF = uicheckbox(dlg, 'Text', 'Auto-save PDF report per file (text + checkpoint images)', ...
                'Position', [15 514 380 22], ...
                'Value', getpref('nestapp', 'autoExportPDF', false));
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
                return
            end

            % Append the Dashboard synthetic entry when any report has
            % gates - keeps the listbox tidy when nothing was screened.
            if anyReportHasGates(allEntries)
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
            if anyReportHasGates(allEntries)
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
            e = allEntries{idx};
            if isfield(e, 'isDashboard') && e.isDashboard
                app.ReportsTextArea.Visible       = 'off';
                app.ReportsDashboardPanel.Visible = 'on';
                renderDashboardPanel(app.ReportsDashboardPanel, ...
                    collectReportStructs(allEntries), ...
                    struct( ...
                        'onRefresh',        @() updateReportsTabImpl(app), ...
                        'onExport',         @() exportDashboardPNG(app, allEntries), ...
                        'onFailedRowClick', @(name) jumpToFileEntry(app, allEntries, name)));
            else
                app.ReportsDashboardPanel.Visible = 'off';
                app.ReportsTextArea.Visible       = 'on';
                if isfield(e, 'text')
                    app.ReportsTextArea.Value = e.text;
                end
            end
        end

        function reRenderReportsOnResize(app)
        % Repaint the Session Quality Dashboard after a window resize so its
        % absolute-positioned children (heatmap, table, histograms) reflow to
        % the new panel size. No-op unless the dashboard is the visible pane;
        % renderDashboardPanel itself clears and re-lays-out from parent size.
            if isempty(app.ReportsDashboardPanel) || ~isvalid(app.ReportsDashboardPanel)
                return
            end
            if ~strcmp(app.ReportsDashboardPanel.Visible, 'on'); return; end
            allEntries = [app.allPipelineReports, app.loadedReports];
            if anyReportHasGates(allEntries)
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
            renderDashboardPanel(fig, collectReportStructs(allEntries));
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

        function LoadReportsButtonPushed(app, ~)
        % Browse for a folder containing *_report_*.mat files and load them.
            folder = uigetdir(getpref('nestapp','lastDataFolder',''), ...
                'Select Folder with Pipeline Reports');
            if isequal(folder, 0); return; end

            matFiles = dir(fullfile(folder, '*_report_*.mat'));
            if isempty(matFiles)
                uialert(app.UIFigure, 'No *_report_*.mat files found in that folder.', ...
                    'No Reports Found');
                return
            end

            loaded = 0;
            for k = 1:numel(matFiles)
                fpath = fullfile(folder, matFiles(k).name);
                try
                    S = load(fpath, 'pipelineReport');
                    if ~isfield(S, 'pipelineReport'); continue; end
                    [txt, ~] = exportReport(S.pipelineReport, tempdir());
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

        function RefreshReportsButtonPushed(app, ~)
        % Reload reports from disk for any loadedReports entries, then refresh tab.
            if isempty(app.loadedReports)
                updateReportsTabImpl(app);
                return
            end
            % Re-derive the folder from the first loaded report
            firstPath = app.loadedReports{1}.report.inputFile;
            folder = fileparts(firstPath);
            if ~isfolder(folder)
                updateReportsTabImpl(app);
                return
            end
            % Clear disk-loaded reports and reload
            app.loadedReports = {};
            matFiles = dir(fullfile(folder, '*_report_*.mat'));
            for k = 1:numel(matFiles)
                fpath = fullfile(folder, matFiles(k).name);
                try
                    S = load(fpath, 'pipelineReport');
                    if ~isfield(S, 'pipelineReport'); continue; end
                    [txt, ~] = exportReport(S.pipelineReport, tempdir());
                    entry.text   = txt;
                    entry.report = S.pipelineReport;
                    app.loadedReports{end+1} = entry;
                catch
                end
            end
            updateReportsTabImpl(app);
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
        % Export the currently selected file's report + checkpoint PNGs
        % as a single PDF. Synthetic Summary / Dashboard entries skip.
            idx = app.ReportsListBox.Value;
            if isempty(idx)
                uialert(app.UIFigure, 'No report selected.', 'Export PDF');
                return
            end
            allEntries = [app.allPipelineReports, app.loadedReports];
            if anyReportHasGates(allEntries)
                allEntries{end+1} = struct('isDashboard', true, ...
                    'text', '', 'report', struct());
            end
            if ~isnumeric(idx) || idx < 1 || idx > numel(allEntries)
                uialert(app.UIFigure, 'No report selected.', 'Export PDF');
                return
            end
            e = allEntries{idx};
            if (isfield(e, 'isSummary') && e.isSummary) ...
                    || (isfield(e, 'isDashboard') && e.isDashboard)
                uialert(app.UIFigure, ...
                    'Pick a single file report (not the Summary or Dashboard entry) before exporting PDF.', ...
                    'Export PDF');
                return
            end
            [~, baseName] = fileparts(e.report.inputFile);
            [fname, fpath] = uiputfile('*.pdf', 'Export Report as PDF', ...
                [baseName, '_report.pdf']);
            if isequal(fname, 0); return; end
            try
                exportFileReportPDF(e.report, fullfile(fpath, fname));
                app.ReportsStatusLabel.Text = sprintf('PDF saved: %s', fname);
            catch err
                uialert(app.UIFigure, ...
                    sprintf('PDF export failed: %s', err.message), ...
                    'Export PDF', 'Icon', 'error');
            end
        end

        function CopyMethodsButtonPushed(app, ~)
        % Copy a methods paragraph to the clipboard. A single file report copies
        % that file's prose; the Session Summary / Dashboard copies the
        % cross-file aggregate (mean +/- SD). Shares the same builders the text
        % reports use (methodsParagraph / methodsParagraphAggregate).
            idx = app.ReportsListBox.Value;
            if isempty(idx); return; end
            allEntries = [app.allPipelineReports, app.loadedReports];
            if ~isnumeric(idx) || idx < 1 || idx > numel(allEntries); return; end
            e = allEntries{idx};

            if isfield(e, 'report') && isstruct(e.report) && ~isempty(e.report)
                methodsText = methodsParagraph(e.report);
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
        % Grey out UITable rows whose Value is a placeholder or literal '[]'.
        % Placeholders start with '(' by convention (e.g. '(all channels)').
        % Works with both cell array data (new format) and table (legacy).
            removeStyle(app.UITable);
            T = app.UITable.Data;
            if isempty(T); return; end
            grey = uistyle('FontColor', [0.6 0.6 0.6], 'FontAngle', 'italic');
            if iscell(T)
                nRows = size(T, 1);
                for row = 1:nRows
                    v = T{row, 2};
                    if isscalar(v) && (isstring(v) || ischar(v))
                        sv = string(v);
                        if (strlength(sv) > 0 && startsWith(sv, '(')) || strcmp(sv, '[]')
                            addStyle(app.UITable, grey, 'cell', [row, 2]);
                        end
                    end
                end
            elseif istable(T)
                for row = 1:height(T)
                    v = T.val{row};
                    if isscalar(v) && (isstring(v) || ischar(v))
                        sv = string(v);
                        if (strlength(sv) > 0 && startsWith(sv, '(')) || strcmp(sv, '[]')
                            addStyle(app.UITable, grey, 'cell', [row, 2]);
                        end
                    end
                end
            end
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

        function LoadSelecEEGdata(app)
            % Ensure EEGLAB functions are on the path.  eeglab('nogui') is
            % normally called by runPipelineCore, but the Visualizing tab can be
            % used independently, so initialise on demand if needed.
            if ~exist('pop_loadset', 'file')
                eeglab('nogui');
            end
            for nfile = 1:numel(app.SelectedFilesforTEP)
                EEGaux = pop_loadset('filename', app.SelectedFilesforTEP{nfile}, 'filepath', app.PathofSelectedFilesforTEP);
                app.EEGofAllSelectedFiles{nfile} = EEGaux;
                app.EEGtime = EEGaux.times;
            end
            app.EEG_SelectedTEPFiles_Loaded = true;
        end

        function LoadLabels(app)
            all_labels = cell(1,numel(app.SelectedFilesforTEP));
            if ~app.EEG_SelectedTEPFiles_Loaded
                LoadSelecEEGdata(app)
            end
            for nn=1:numel(app.SelectedFilesforTEP)
                EEG = app.EEGofAllSelectedFiles{nn};
                all_labels{nn} = {EEG.chanlocs.labels};
            end
            % Common labels
            app.Common_Labels.Items = app.elecList;
            for i = 1:length(all_labels)
                app.Common_Labels.Items = intersect(app.Common_Labels.Items, all_labels{i});
            end
            % All lables across selected files
            total_labels = app.elecList;
            for i = 1:length(all_labels)
                total_labels = union(total_labels, all_labels{i});
            end
            % Uncommon labels = total - common
            uncommon_labels.Items = intersect(app.elecList,setdiff(total_labels, app.Common_Labels.Items));
            for nn = 1:length(uncommon_labels.Items)
                propName = [upper(uncommon_labels.Items{nn}), 'Button'];
                if isprop(app, propName)
                    app.(propName).Enable = 'off';
                    app.(propName).Value  = 0;
                end
            end
            
        end
        
        function findTEPelecs(app)
            mm = 0; % Selected TEP elecs Counter
            app.ROIelecsLabels = []; % Stores the ROI elecs
            if app.DontfindcommonelectrodesCheckBox.Value
                for nn = 1:length(app.elecList)
                    if app.([upper(app.elecList{nn}),'Button']).Value
                        mm = mm+1;
                        app.ROIelecsLabels{mm} = app.elecList{nn};
                    end
                end
            else
                for nn = 1:length(app.Common_Labels.Items)
                    if app.([upper(app.Common_Labels.Items{nn}),'Button']).Value
                        mm = mm+1;
                        app.ROIelecsLabels{mm} = app.Common_Labels.Items{nn};
                    end
                end
            end
        end
        
        function plotTEP(app)
            if ~app.EEG_SelectedTEPFiles_Loaded
                LoadSelecEEGdata(app)
            end
            nFiles  = numel(app.EEGofAllSelectedFiles);
            nTimes  = numel(app.EEGtime);
            TEP_ROI = zeros(nFiles, nTimes);
            for nfile = 1:nFiles
                EEGaux = app.EEGofAllSelectedFiles{1, nfile};
                ROIind = find(ismember({EEGaux.chanlocs.labels}, app.ROIelecsLabels));
                TEP_ROI(nfile,:) = mean(mean(EEGaux.data(ROIind,:,:), 3, 'omitmissing'), 1, 'omitmissing');
            end

            SMOOTH_WIN_PTS = 5;       % 5-point moving average (~5 ms at 1 kHz)
            app.TEP2Export = TEP_ROI;
            grandMean = mean(TEP_ROI, 1, 'omitmissing');
            TEP_ROISD = std(TEP_ROI, 1, 1) / sqrt(nFiles);
            co    = app.UIAxes.ColorOrder;
            meanx = smoothdata(grandMean,  'movmean', SMOOTH_WIN_PTS);
            sdx   = smoothdata(TEP_ROISD, 'movmean', SMOOTH_WIN_PTS);
            xf = [app.EEGtime(1) app.EEGtime  app.EEGtime(end) app.EEGtime(end:-1:1)];
            yf = [meanx(1)-sdx(1)/2 meanx+sdx/2 meanx(end)-sdx(end)/2 meanx(end:-1:1)-sdx(end:-1:1)/2];

            % Legend label: base filename of first selected file
            if iscell(app.SelectedFilesforTEP) && ~isempty(app.SelectedFilesforTEP)
                [~, dispName, ~] = fileparts(app.SelectedFilesforTEP{1});
            else
                dispName = 'TEP';
            end

            if app.NewFigureButton.Value
                cla(app.UIAxes, 'reset');
                Colr = co(1, :);
                hold(app.UIAxes, 'on');
                fill(app.UIAxes, xf, yf, Colr, 'FaceAlpha', 0.5, 'LineStyle', 'none', 'HandleVisibility', 'off');
                plot(app.UIAxes, app.EEGtime, meanx, 'Color', Colr, 'LineWidth', 2, 'DisplayName', dispName);
                hold(app.UIAxes, 'off');
                xlim(app.UIAxes, app.DefaulTEPxLim);
            elseif app.AddtocurrentFigureButton.Value
                % Only count main TEP lines (HandleVisibility='on') to determine next color.
                mainLines = findobj(app.UIAxes, 'Type', 'Line', 'HandleVisibility', 'on');
                if isempty(mainLines)
                    usedColors = zeros(0, 3);
                else
                    usedColors = reshape([mainLines.Color], 3, [])';
                end
                Colr = co(1, :);
                for k = 1:size(co, 1)
                    candidate = co(k, :);
                    if ~any(all(abs(usedColors - candidate) < 1e-6, 2))
                        Colr = candidate;
                        break;
                    end
                end
                prevYLim = ylim(app.UIAxes);
                hold(app.UIAxes, 'on');
                fill(app.UIAxes, xf, yf, Colr, 'FaceAlpha', 0.5, 'LineStyle', 'none', 'HandleVisibility', 'off');
                plot(app.UIAxes, app.EEGtime, meanx, 'Color', Colr, 'LineWidth', 2, 'DisplayName', dispName);
                xlim(app.UIAxes, app.DefaulTEPxLim);
                % Expand y-axis to accommodate new data; never shrink existing range
                newYLim = ylim(app.UIAxes);
                ylim(app.UIAxes, [min(prevYLim(1), newYLim(1)), max(prevYLim(2), newYLim(2))]);
            end

            legend(app.UIAxes, 'show', 'Location', 'best');

            % Component detection runs on the SMOOTHED grand mean - the same
            % waveform that is plotted (meanx) and the same smoothing the batch
            % CSV path applies in batchTEPExtract. Feeding TESA the raw grandMean
            % let baseline noise wiggles register as spurious local extrema
            % (e.g. a "negative peak" on a flat baseline or a shoulder on a
            % monotonic rise), which defeats tesa_peakanalysis's own peak guard.
            try
                app.tepPeaks = tepPeakFinder(meanx, app.EEGtime, app.tepComponentDefs);
                if app.ShowComponentsButton.Value
                    overlayTEPComponents(app);
                end
                populateTEPComponentTable(app);
            catch ME
                if strcmp(ME.identifier, 'tepPeakFinder:noTESA')
                    uialert(app.UIFigure, ...
                        ['TESA not found. Add TESA to the MATLAB path to enable ' ...
                         'component detection (Show Components / Extract Peaks).'], ...
                        'TESA Required');
                    app.ShowComponentsButton.Value = false;
                else
                    rethrow(ME);
                end
            end
        end
        
        function EEG_topoplot(app)
            cla(app.UIAxes2)
            TOPOPLOT_INTRAD = 0.55;   % EEGLAB default interpolation radius
            SMOOTH_METHOD   = 'movmean';
            SMOOTH_WIN_PTS  = 5;      % 5-point moving average (~5 ms at 1 kHz)
            if ~app.EEG_SelectedTEPFiles_Loaded
                LoadSelecEEGdata(app)
            end
            LoadLabels(app);
            BIGEEG = zeros(numel(app.Common_Labels.Items), length(app.EEGtime),numel(app.EEGofAllSelectedFiles));
            for nfile = 1:numel(app.EEGofAllSelectedFiles)
                EEGaux = app.EEGofAllSelectedFiles{1,nfile};
                ChansLocs = EEGaux.chanlocs;
                commonElectrodsInd = ismember({ChansLocs.labels},app.Common_Labels.Items);
                BIGEEG(:,:,nfile) = mean(EEGaux.data(commonElectrodsInd,:,:),3,"omitmissing");
                
            end
            ChansLocs(~commonElectrodsInd) = [];
            yp = smoothdata(mean(BIGEEG,3,"omitmissing")',SMOOTH_METHOD,SMOOTH_WIN_PTS)'; % Smooth the EEGdata along subjects
            timepoint = app.TopoplottimeSpinner.Value;
            Topo_ind = [round(timepoint-app.WindowsizefortimeaveragedTopoplotEditField.Value/2),...
                round(timepoint+app.WindowsizefortimeaveragedTopoplotEditField.Value/2)];
            Topo_ind = [find(app.EEGtime==Topo_ind(1)), find(app.EEGtime==Topo_ind(2))];
            Topo = mean(yp(:,Topo_ind),2,"omitmissing");
            
            oldFig = gcf;

            % Create invisible figure to trick topoplot
            invisibleFig = figure('Visible', 'off');
            copyobj(app.UIAxes2, invisibleFig);  % clone axes
            newAx = findobj(invisibleFig, 'Type', 'Axes');

            % Plot into cloned axes
            axes(newAx);  % set as current
            topoplot(mean(Topo,2),ChansLocs,'electrodes','off',...
                'numcontour',5,'intsquare','on','style','map','conv', 'on', 'intrad',TOPOPLOT_INTRAD);axis auto
            colormap(app.UIAxes2,'hsv')
            % Copy contents back to app UIAxes
            cla(app.UIAxes2);
            copyobj(allchild(newAx), app.UIAxes2);

            % Clean up
            close(invisibleFig);  % Close hidden fig
            figure(oldFig);axis(app.UIAxes2,'auto')
            close
        end

        
        function tf = isFileSelected(app)
            tf = ~isempty(app.SelectedFilesforTEP);
        end

        function overlayTEPComponents(app)
        % Draw dashed vertical lines and text labels for each detected TEP component.
        % Assumes app.tepPeaks is already populated by tepPeakFinder.
            if isempty(app.tepPeaks)
                return
            end
            ax = app.UIAxes;
            yLims = ylim(ax);
            % Place labels near the top of the axes (80% height)
            labelY = yLims(1) + 0.80 * (yLims(2) - yLims(1));
            hold(ax, 'on');
            for i = 1:numel(app.tepPeaks)
                pk = app.tepPeaks(i);
                if ~pk.found
                    continue
                end
                xline(ax, pk.latencyMs, '--', 'Color', [0.4 0.4 0.4], ...
                    'LineWidth', 1, 'HandleVisibility', 'off');
                text(ax, pk.latencyMs, labelY, ...
                    sprintf('%s\n%.0f ms\n%.1f uV', pk.name, pk.latencyMs, pk.amplitudeUV), ...
                    'FontSize', 7, 'HorizontalAlignment', 'center', ...
                    'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'top');
            end
            hold(ax, 'off');
        end

        function populateTEPComponentTable(app)
        % Fill TEPComponentTable from app.tepPeaks.
        % Shows '-' for components not found.
            if isempty(app.tepPeaks)
                app.TEPComponentTable.Data = {};
                return
            end
            nComp = numel(app.tepPeaks);
            tableData = cell(nComp, 3);
            for i = 1:nComp
                pk = app.tepPeaks(i);
                tableData{i, 1} = pk.name;
                if pk.found
                    tableData{i, 2} = pk.latencyMs;
                    tableData{i, 3} = pk.amplitudeUV;
                else
                    tableData{i, 2} = '-';
                    tableData{i, 3} = '-';
                end
            end
            app.TEPComponentTable.Data = tableData;
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

            % Visualizing tab
            app.PLOTTEPButton.Tooltip         = 'Plot TMS-evoked potential waveforms for the selected files and electrodes';
            app.ShowComponentsButton.Tooltip          = 'Detect and overlay TEP component peaks on the TEP plot';
            app.EditComponentWindowsButton.Tooltip    = 'Customise the search windows used for each TEP component';
            app.TOPOPLOTButton.Tooltip        = 'Plot a scalp topographic map at the specified time point';
            app.ExportTEPFigureButton.Tooltip = 'Export the current TEP plot as PNG or PDF';
            app.ReLoadAvailableElectrodesButton.Tooltip = ...
                'Reload the electrode list from the currently selected files';
            app.SelectAllCheckBox.Tooltip   = 'Select all available files for TEP plotting';
            app.UseCurrentlyCleanedDataCheckBox.Tooltip = ...
                'Use the most recently processed output instead of selecting files manually';
            app.DontfindcommonelectrodesCheckBox.Tooltip = ...
                ['When checked: show all selected electrodes regardless of whether they ' ...
                'appear in every file. When unchecked: restrict to electrodes present across all selected files.'];
        end

        % Code that executes after component creation
        function startupFcn(app)
            clc
            steps = stepRegistry();
            app.StepsListBox.Items = {steps.name};

            app.info = {steps.info};

            app.spec = repmat(struct('name','','params',struct()), 0, 1);
            app.SelectedListBox.Items(:) = [];
            app.SelectedListBox.ItemsData(:) = [];
            app.UITable.Data = [];
            app.ItemNum = 1;
            app.originalSize      = app.UIFigure.Position(3:4);
            app.tepComponentDefs  = defaultTEPComponentDefs();
            applyTooltips(app);
            loadPrefs(app);
            buildRecentFilesMenu(app);
            buildRecentPipelinesMenu(app);
            updateStatusBar(app);
            clc
        end

        % Clicked callback: StepsListBox
        function StepsListBoxClicked(app, ~)
            % Detect double-click via inter-click interval (< 500 ms).
            % ListBoxInteraction has no NumClicks property in R2025b.
            t = datetime('now');
            if seconds(t - app.lastStepClick) < 0.5
                appendStep(app, app.StepsListBox.Value);
            end
            app.lastStepClick = t;
        end

        % Value changed function: StepsListBox
        function StepsListBoxValueChanged(app, ~)
            value = app.StepsListBox.Value;
            ind = find(ismember(app.StepsListBox.Items,value));
            app.InfoTextArea.Value = string(app.info{ind});
            app.selectedItem = [];
        end

        % Button pushed function: AddButton
        function AddButtonPushed(app, ~)
            stepName = app.StepsListBox.Value;
            appendStep(app, stepName);
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
        function SelectDataButtonPushed(app, ~)
            try
                startFolder = getpref('nestapp', 'lastDataFolder', '');
                [app.file,app.path] = uigetfile( ...
                    {'*.set;*.vhdr;*.cdt;*.cnt',...
                    'Data Files (*.set,*.vhdr,*.cdt,*.cnt)'; ...
                    '*.set','Set Files (*.set)'; ...
                    '*.vhdr','VHDR Files (*.vhdr)'; ...
                    '*.cdt','CDT Files (*.cdt)'; ...
                    '*.cnt','CNT Files (*.cnt)'; ...
                    '*.*',  'All Files (*.*)'}, ...
                    'Select File(s)', startFolder, 'multiSelect','on');
                if iscell(app.file)
                    app.NSelecFiles = numel(app.file);
                else
                    app.NSelecFiles = 1;
                    app.file = {app.file};
                end

                app.SelectedFilesListBox.Items = app.file;
                setpref('nestapp', 'lastDataFolder', app.path);
                pushRecent(app, 'recentFiles', app.path);
                buildRecentFilesMenu(app);
                updateStatusBar(app);

            catch
                warning('Please select at least one file!')
                app.SelectedFilesListBox.Items = {};
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
            end
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

            % Silently initialise EEGLAB if its plugins aren't on the path yet.
            global PLUGINLIST %#ok<GVMIS>
            if isempty(PLUGINLIST)
                try
                    evalc('eeglab nogui');
                catch ME
                    uialert(app.UIFigure, ...
                        ['EEGLAB could not be initialised: ' ME.message newline ...
                         'Verify the EEGLAB path in Preferences.'], ...
                        'EEGLAB Init Failed', 'Icon', 'error');
                    return
                end
            end

            filePaths = cellfun(@(f) fullfile(app.path, f), app.file, 'UniformOutput', false);

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

            opts.uiFigure     = app.UIFigure;
            opts.pipelineName = app.pipelineName;
            opts.statusBar    = app.StatusBar;
            opts.parallel     = app.ParallelCheckBox.Value;
            opts.chanLocFile  = app.preSelectedChanFile;

            try
                [allReports, allSummaries] = runPipelineCore(app.spec, filePaths, opts);
            catch err
                if strcmp(err.identifier, 'nestapp:cancelled')
                    return
                end
                uialert(app.UIFigure, err.message, 'Pipeline Error', 'Icon', 'error');
                return
            end

            if numel(allReports) > 1
                summEntry.text      = summarizeReports(allReports);
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
            if app.UseCurrentlyCleanedDataCheckBox.Value
                UseCurrentlyCleanedDataCheckBoxValueChanged(app)
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
                app.TextArea.Value = strjoin(val, ', ');
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
            drawnow limitrate  % throttle: skip redraws that arrive faster than screen refresh
            newSize = app.UIFigure.Position(3:4);
            minW = 650; minH = 420;
            if newSize(1) < minW || newSize(2) < minH
                newSize(1) = max(newSize(1), minW);
                newSize(2) = max(newSize(2), minH);
                app.UIFigure.Position(3:4) = newSize;
            end
            sX = newSize(1) / app.originalSize(1);
            sY = newSize(2) / app.originalSize(2);
            rescaleComponents(app, sX, sY);
            reRenderReportsOnResize(app);  % reflow the Quality Dashboard if it's showing
        end

        % Cell edit callback: UITable
        function UITableCellEdit(app, event)
            stepIdx = selectedStepIndex(app);
            if isempty(stepIdx) || stepIdx > numel(app.spec); return; end

            reg    = stepRegistry();
            step   = app.spec(stepIdx);
            regIdx = find(strcmp({reg.name}, step.name), 1);
            if isempty(regIdx); return; end

            app.spec = applyParamEdit(app.spec, stepIdx, event.Indices(1), event.NewData, reg(regIdx));
            app.pipelineDirty = true;
            updateStatusBar(app);
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

        % Button pushed function: PLOTTEPButton
        function PLOTTEPButtonPushed(app, ~)
            if ~isFileSelected(app)
                warning('Please select at least a file to plot the TEP!');
            else
                LoadLabels(app);
                findTEPelecs(app);
                plotTEP(app)
                app.TEPWindowSlider.Limits = [app.EEGtime(1) app.EEGtime(end)];
                app.TEPCreated = true;
                app.TEPWindowSlider.Value = app.DefaulTEPxLim;
                app.ExportTEPDataButton.Enable    = 'on';
                app.TEPvarNameEditFieldLabel.Enable = 'on';
                app.TEPvarNameEditField.Enable    = 'on';
            end
        end

        % Value changed function: UseCurrentlyCleanedDataCheckBox
        function UseCurrentlyCleanedDataCheckBoxValueChanged(app, ~)
            value = app.UseCurrentlyCleanedDataCheckBox.Value;
            if value
                if ~isempty(app.path) && ~isempty(app.cleanedName)
                    app.FilesListBox.Items = {};
                    for nn=1:app.NSelecFiles
                        dots = find(ismember(app.file{nn},'.'));
                        fname=[app.file{nn}(1:dots(end)-1),'_',app.cleanedName,'.set'];
                        app.TEPfiles{nn} = fname;
                        app.PathofSelectedFilesforTEP = app.path;
                        app.FolderEditField_2.Value = app.PathofSelectedFilesforTEP;

                    end
                    app.FilesListBox.Items =  app.TEPfiles';
                    app.TEPCreated = false;  % file selection changed - existing plot is stale

                elseif ~isempty(app.FilesListBox.Items)
                    warning('No files have been cleaned recently. Try selecting files!')
                end
                app.TOPOPLOTButton.Enable = 'on';
                app.ExportTEPFigureButton.Enable = "on";
                app.PLOTTEPButton.Enable = "on";
            else
                warning('Try selecting files!')
            end
        end

        % Button pushed function: SelectDataButton_2
        function SelectDataButton_2Pushed(app, ~)
            try
                [app.TEPfiles,app.PathofSelectedFilesforTEP] = uigetfile( ...
                    {'*.set',...
                    'Data Files (*.set)'; ...
                    '*.set','Set Files (*.set)'} , ...
                    'Select File(s)','multiSelect','on');
                if iscell(app.TEPfiles)
                    app.NumberOfSelecFilesforTEP = numel(app.TEPfiles);
                else
                    app.NSelecFiles = 1;
                    app.TEPfiles = {app.TEPfiles};
                end

                % Invalidate EEG cache - new files mean stale loaded data must be discarded.
                app.EEG_SelectedTEPFiles_Loaded = false;
                app.EEGofAllSelectedFiles = {};

                % app.FileEditField_2.Value = app.TEPfiles{1};
                app.FolderEditField_2.Value = app.PathofSelectedFilesforTEP;
                app.FilesListBox.Items =  app.TEPfiles';
                app.TOPOPLOTButton.Enable = 'on';
                app.ExportTEPFigureButton.Enable = "on";
                app.PLOTTEPButton.Enable = "on";
                app.PlotEEGdataButton.Enable = 'on';
                app.EEGDatasetDropDown.Enable = "on";
            catch
                warning('Please select at least one file!')
                if isempty(app.FilesListBox.Items)
                    app.FolderEditField_2.Value = '';
                    app.TOPOPLOTButton.Enable = 'off';
                    app.ExportTEPFigureButton.Enable = "off";
                    app.PLOTTEPButton.Enable = "off";
                    app.PlotEEGdataButton.Enable = 'off';
                    app.EEGDatasetDropDown.Enable = "off";
                end
            end
        end

        % Value changed function: FilesListBox
        function FilesListBoxValueChanged(app, event)
            app.EEG_SelectedTEPFiles_Loaded = false;
            app.EEGofAllSelectedFiles = []; % Every time the new file is checked clear the all loaded EEG data
            if ~isempty(event.Value) || event.Value ~= 0
                % fname = event.Value;
                % ind = event.ValueIndex;
                app.SelectedFilesforTEP = event.Value;
                app.SelectAllCheckBox.Value = 0;
            end
            app.EEGDatasetDropDown.Items = app.SelectedFilesforTEP;
        end

        % Button pushed function: TOPOPLOTButton
        function TOPOPLOTButtonPushed(app, ~)
            if isFileSelected(app)
                EEG_topoplot(app)
            end
        end

        % Value changed function: SelectAllCheckBox
        function SelectAllCheckBoxValueChanged(app, ~)
            value = app.SelectAllCheckBox.Value;
            if value
                app.SelectedFilesforTEP = app.TEPfiles;
                app.FilesListBox.ValueIndex = 1:max(size(app.TEPfiles));
                % ind = app.FilesListBox.ValueIndex;
            else
                app.SelectedFilesforTEP = [];
                app.FilesListBox.ValueIndex = [];
            end
        end

        % Value changed function: DontfindcommonelectrodesCheckBox
        function DontfindcommonelectrodesCheckBoxValueChanged(app, ~)
            value = app.DontfindcommonelectrodesCheckBox.Value;
            if ~value
                app.ReLoadAvailableElectrodesButton.Enable = 1;
            else
                app.ReLoadAvailableElectrodesButton.Enable = 0;
            end
        end

        % Button pushed function: ReLoadAvailableElectrodesButton
        function ReLoadAvailableElectrodesButtonPushed(app, ~)
            if isFileSelected(app)
                LoadLabels(app);
            end

        end

        % Button pushed function: ExportTEPFigureButton
        function ExportTEPFigureButtonPushed(app, ~)
            if ~app.TEPCreated
                uialert(app.UIFigure, 'Please plot a TEP first.', 'No figure');
                return
            end
            [fname, fpath] = uiputfile( ...
                {'*.png','PNG image';'*.pdf','PDF file';'*.fig','MATLAB figure'}, ...
                'Export TEP Figure', 'tep_figure');
            if isequal(fname, 0)
                return
            end
            exportgraphics(app.UIAxes, fullfile(fpath, fname), 'Resolution', 300);
            [~, nm, ~] = fileparts(fname);
            savefig(ancestor(app.UIAxes, 'figure'), fullfile(fpath, [nm '.fig']));
        end

        % Value changing function: TEPWindowSlider
        function TEPWindowSliderValueChanging(app, event)
            changingValue = event.Value;
            app.UIAxes.XLim = changingValue;
        end

        % Value changed function: TopoplottimeSpinner
        function TopoplottimeSpinnerValueChanged(app, event)
            value = event.Value;
            app.Slider.Value = value;
        end

        % Value changed function: EEGDatasetDropDown
        function EEGDatasetDropDownValueChanged(~, ~)
        end

        % Button pushed function: PlotEEGdataButton
        function PlotEEGdataButtonPushed(app, ~)
            subInd = strcmpi(app.SelectedFilesforTEP, app.EEGDatasetDropDown.Value);
            if isFileSelected(app)
                if ~app.EEG_SelectedTEPFiles_Loaded
                    LoadSelecEEGdata(app)
                end
                pop_eegplot(app.EEGofAllSelectedFiles{subInd},1,1,1)
            end
        end

        % Button pushed function: ExportTEPDataButton
        function ExportTEPDataButtonPushed(app, ~)
            assignin('base', app.TEPvarNameEditField.Value, app.TEP2Export)
        end

        % Value changed function: TEPvarNameEditField
        function TEPvarNameEditFieldValueChanged(~, ~)
        end

        % Value changed function: ShowComponentsButton
        function ShowComponentsButtonValueChanged(app, ~)
            if app.TEPCreated
                if app.ShowComponentsButton.Value
                    overlayTEPComponents(app);
                else
                    % Replot without overlays - cla then replot
                    plotTEP(app);
                end
            end
        end

        % Button pushed function: EditComponentWindowsButton
        function EditComponentWindowsButtonPushed(app, ~)
        % EDITCOMPONENTWINDOWSBUTTONPUSHED  Open a modal dialog for editing TEP component windows.
            defs = app.tepComponentDefs;
            nComp = numel(defs);

            fig = uifigure('Name', 'TEP Component Windows', ...
                'Position', [200 200 530 265], 'WindowStyle', 'modal');

            % Build editable table
            tableData = cell(nComp, 5);
            for i = 1:nComp
                tableData{i, 1} = defs(i).name;
                tableData{i, 2} = defs(i).polarity;
                tableData{i, 3} = defs(i).nomLatency;
                tableData{i, 4} = defs(i).winStart;
                tableData{i, 5} = defs(i).winEnd;
            end

            tbl = uitable(fig, ...
                'Position',       [10 55 510 195], ...
                'Data',           tableData, ...
                'ColumnName',     {'Component', 'Polarity', 'Nom. Latency (ms)', 'Win Start (ms)', 'Win End (ms)'}, ...
                'ColumnEditable', [false, false, true, true, true], ...
                'ColumnWidth',    {80, 65, 130, 110, 100}, ...
                'RowName',        {});

            uibutton(fig, 'Text', 'Reset Defaults', ...
                'Position', [10 12 130 30], ...
                'ButtonPushedFcn', @(~,~) resetDefaults(tbl));

            uibutton(fig, 'Text', 'Cancel', ...
                'Position', [300 12 100 30], ...
                'ButtonPushedFcn', @(~,~) close(fig));

            uibutton(fig, 'Text', 'Apply', ...
                'Position', [410 12 110 30], ...
                'ButtonPushedFcn', @(~,~) applyAndClose(tbl, fig));

            uiwait(fig);

            function resetDefaults(t)
                % app is accessible from the enclosing method scope
                defaults = defaultTEPComponentDefs();
                d = cell(numel(defaults), 5);
                for k = 1:numel(defaults)
                    d{k, 1} = defaults(k).name;
                    d{k, 2} = defaults(k).polarity;
                    d{k, 3} = defaults(k).nomLatency;
                    d{k, 4} = defaults(k).winStart;
                    d{k, 5} = defaults(k).winEnd;
                end
                t.Data = d;
            end

            function applyAndClose(t, f)
                % app is accessible from the enclosing method scope
                d = t.Data;
                for k = 1:size(d, 1)
                    app.tepComponentDefs(k).nomLatency = d{k, 3};
                    app.tepComponentDefs(k).winStart   = d{k, 4};
                    app.tepComponentDefs(k).winEnd     = d{k, 5};
                end
                % Re-detect and replot if TEP is already shown
                if app.TEPCreated
                    PLOTTEPButtonPushed(app, []);
                end
                close(f);
            end
        end

        % -- Analysis Tab callbacks ----------------------------------------

        % Button pushed function: ExtractPeaksCSVButton
        function ExtractPeaksCSVButtonPushed(app, ~)
        % Extract peaks across all selected files and save as CSV.
            findTEPelecs(app);   % refresh ROI from current electrode button state
            if isempty(app.ROIelecsLabels)
                uialert(app.UIFigure, ...
                    'No ROI electrodes selected. Choose electrodes on the Visualizing tab.', ...
                    'Extract Peaks');
                return
            end
            if isempty(app.SelectedFilesforTEP)
                uialert(app.UIFigure, ...
                    'No files selected. Select .set files on the Visualizing tab.', ...
                    'Extract Peaks');
                return
            end
            if isempty(which('tesa_peakanalysis'))
                uialert(app.UIFigure, ...
                    'TESA toolbox not found on path. Cannot run peak extraction.', ...
                    'Extract Peaks');
                return
            end

            [fname, fpath] = uiputfile('*.csv', 'Save TEP Peaks CSV', 'tep_peaks.csv');
            if isequal(fname, 0); return; end
            csvPath = fullfile(fpath, fname);

            filePaths = cellfun(@(f) fullfile(app.PathofSelectedFilesforTEP, f), ...
                app.SelectedFilesforTEP, 'UniformOutput', false);

            d = uiprogressdlg(app.UIFigure, ...
                'Title',          'Extracting TEP Peaks', ...
                'Message',        'Starting...', ...
                'Cancelable',     'off', ...
                'ShowPercentage', 'on');

            try
                [results, warnings] = batchTEPExtract(filePaths, app.ROIelecsLabels, ...
                    'compDefs',    app.tepComponentDefs, ...
                    'csvPath',     csvPath, ...
                    'progressFcn', @(i,n) updateExtractionProgress(d, i, n, filePaths));
            catch ME
                if isvalid(d); close(d); end
                uialert(app.UIFigure, ME.message, 'Extraction Error');
                return
            end
            if isvalid(d); close(d); end

            nRows = height(results);
            if isempty(warnings)
                app.AnalysisStatusLabel.Text = sprintf('Extracted %d rows -> %s', nRows, fname);
            else
                app.AnalysisStatusLabel.Text = sprintf( ...
                    'Extracted %d rows -> %s  (%d warning(s))', nRows, fname, numel(warnings));
                uialert(app.UIFigure, strjoin(warnings, newline), 'Extraction Warnings');
            end

            function updateExtractionProgress(dlg, iFile, nFiles, fps)
                [~, nm] = fileparts(fps{iFile});
                dlg.Value   = (iFile - 1) / nFiles;
                dlg.Message = sprintf('File %d / %d  -  %s', iFile, nFiles, nm);
                drawnow limitrate
            end
        end

        % Selection changed function: TabGroup
        function TabGroupSelectionChanged(app, event)
        % Refresh the Analysis tab selection summary whenever it becomes active.
            if event.NewValue == app.AnalysisTab
                updateAnalysisSelectionSummary(app);
            end
        end

        function updateAnalysisSelectionSummary(app)
        % Update the read-only summary label on the Analysis tab.
            findTEPelecs(app);   % refresh ROI from current electrode button state
            nFiles = numel(app.SelectedFilesforTEP);
            nROI   = numel(app.ROIelecsLabels);
            if nFiles == 0 && nROI == 0
                app.AnalysisSelectionLabel.Text = ...
                    'Select files and ROI electrodes on the Visualizing tab.';
                return
            end
            fileStr = sprintf('%d file(s) selected', nFiles);
            if nROI == 0
                roiStr = 'No ROI electrodes selected';
            elseif nROI <= 6
                roiStr = sprintf('ROI: %s', strjoin(app.ROIelecsLabels, ', '));
            else
                roiStr = sprintf('ROI: %s ... (%d electrodes total)', ...
                    strjoin(app.ROIelecsLabels(1:6), ', '), nROI);
            end
            app.AnalysisSelectionLabel.Text = sprintf('%s          %s', fileStr, roiStr);
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

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
