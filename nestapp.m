% WARNING: Do not open nestapp_designer.mlapp and save — App Designer will
% regenerate this file and overwrite startupFcn and other hand-edited methods.
% All edits must be made directly to nestapp.m.
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
        ExportReportsCSVButton          matlab.ui.control.Button
        CopyMethodsButton               matlab.ui.control.Button
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
        
    end
    properties (Access = public)
        % Tab Cleaning
        selectedItem % Selected Table Item Values
        info % Command Information and description
        path % File Path
        file % File Name
        steps2run % Stores the steps we select to run
        nstep; % Step Number in steps2run
        initialVars; % Stores initial vars generated by eeglab
        ChangedVal       % Display table per pipeline step: var=friendly label, val=value string
        DefaultsVal      % Default display tables (registry-indexed, parallel to StepsListBox)
        defaultParamKeys % Raw EEGLAB keys (registry-indexed, same row order as DefaultsVal)
        stepParamKeys    % Raw EEGLAB keys (pipeline-indexed, same row order as ChangedVal)
        originalSize     % [w h] of UIFigure at creation — used by UIFigureSizeChanged
        convert = 0;
        needchanloc = 1; % The EEG data must contain EEG channel Location
        NSelecFiles % Number of selcted Files for EEG preprocessing
        clickedItem = [];
        doubleClicked = 0;
        cleanedName % Name used to rename the save cleaned EEG data
        TEPfiles % File list for calculating TEPs
        
        % Tab Visualizing
        PathofSelectedFilesforTEP
        NumberOfSelecFilesforTEP
        SelectedFilesforTEP % Selected files to plot the TEP
        Common_Labels % Commong electrod name among files
        ROIelecsLabels % Selected electrodes as Region of Interest
        TEPCreated = 0; % If the TEP plot bottomn has been pressed once
        EEG_SelectedTEPFiles_Loaded = 0;
        EEGofAllSelectedFiles = [];
        DefaulTEPxLim = [-50 300]; % Default xLim for time in TEP
        EEGtime
        TEP2Export
        MenuRecentFiles     % Handle to 'Recent Files' submenu — rebuilt on open
        MenuRecentPipelines % Handle to 'Recent Pipelines' submenu — rebuilt on open
        StatusBar           % uilabel pinned to bottom of UIFigure — visible on both tabs
        pipelineDirty   = false    % true when pipeline has unsaved changes
        pipelineName    = ''       % filename of last saved/loaded pipeline
        lastStepClick   = NaT     % datetime of last StepsListBox click (double-click detection)
        tepPeaks        = struct([]) % struct array from tepPeakFinder; cached after each PLOT TEP
        allPipelineReports = {}    % cell array of report entry structs from current session
        loadedReports      = {}    % cell array of report entry structs loaded from disk
    end

    methods (Access = private)
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
                assignin('base', 'files', app.file);
                assignin('base', 'paths', app.path);
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
            % Temporarily set the file so LoadPipelineButtonPushed can read it.
            % We call the loader directly rather than through the button to avoid
            % opening the file dialog again.
            try
                pipeline = load(pPath, '-mat');
                app.SelectedListBox.Items     = pipeline.PLItems;
                app.SelectedListBox.ItemsData = pipeline.PLItemsData;
                app.ChangedVal                = pipeline.VarIns;
                app.UITable.Data              = [];
                if isfield(pipeline, 'ParamKeys')
                    app.stepParamKeys = pipeline.ParamKeys;
                else
                    % Trigger the old-format upgrade path
                    LoadPipelineButtonPushed(app, []);
                    return
                end
                setpref('nestapp', 'lastPipelineFolder', pFolder);
                pushRecent(app, 'recentPipelines', pPath);
                buildRecentPipelinesMenu(app);
            catch err
                uialert(app.UIFigure, err.message, 'Load Error', 'Icon', 'error');
            end
        end

        % MenuSelected callback wrappers (thin shims so uimenu can call private methods)
        function openPreferencesMenu(app, ~)
            openPreferences(app);
        end

        function showAboutMenu(app, ~)
            showAbout(app);
        end

        function loadPrefs(~)
        % LOADPREFS  Read persistent preferences and apply to app state.
        %   Called from startupFcn. Uses MATLAB getpref with 'nestapp' group.
        %   The app handle is accepted but not used — prefs apply globally
        %   (addpath) rather than writing to removed UI components.
            eeglabPath = getpref('nestapp', 'eeglabPath', '');
            if ~isempty(eeglabPath) && isfolder(eeglabPath)
                addpath(eeglabPath);
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
                'Position', [200 200 420 415], ...
                'WindowStyle', 'modal', 'Resize', 'off');

            % --- EEGLAB section ---
            uilabel(dlg, 'Text', 'EEGLAB', 'FontWeight', 'bold', ...
                'Position', [15 375 200 20]);
            uilabel(dlg, 'Text', 'Path:', ...
                'Position', [15 350 35 22], 'HorizontalAlignment', 'right');
            fEeglab = uieditfield(dlg, 'text', ...
                'Position', [55 350 275 22], 'Editable', 'on', ...
                'Value', getpref('nestapp','eeglabPath',''));
            uibutton(dlg, 'Text', 'Browse...', 'Position', [335 350 70 22], ...
                'ButtonPushedFcn', @(~,~) browseEeglab());

            % --- Default Locations section ---
            uilabel(dlg, 'Text', 'Default Locations', 'FontWeight', 'bold', ...
                'Position', [15 320 200 20]);
            uilabel(dlg, 'Text', 'Data folder:', ...
                'Position', [15 295 65 22], 'HorizontalAlignment', 'right');
            fData = uieditfield(dlg, 'text', ...
                'Position', [85 295 245 22], 'Editable', 'on', ...
                'Value', getpref('nestapp','lastDataFolder',''));
            uibutton(dlg, 'Text', 'Browse...', 'Position', [335 295 70 22], ...
                'ButtonPushedFcn', @(~,~) browseFolder(fData));
            uilabel(dlg, 'Text', 'Pipeline folder:', ...
                'Position', [15 267 80 22], 'HorizontalAlignment', 'right');
            fPipeline = uieditfield(dlg, 'text', ...
                'Position', [100 267 230 22], 'Editable', 'on', ...
                'Value', getpref('nestapp','lastPipelineFolder',''));
            uibutton(dlg, 'Text', 'Browse...', 'Position', [335 267 70 22], ...
                'ButtonPushedFcn', @(~,~) browseFolder(fPipeline));
            uilabel(dlg, 'Text', 'Reports folder:', ...
                'Position', [15 239 80 22], 'HorizontalAlignment', 'right');
            fReports = uieditfield(dlg, 'text', ...
                'Position', [100 239 230 22], 'Editable', 'on', ...
                'Value', getpref('nestapp','reportFolder',''));
            uibutton(dlg, 'Text', 'Browse...', 'Position', [335 239 70 22], ...
                'ButtonPushedFcn', @(~,~) browseFolder(fReports));

            % --- Behaviour section ---
            uilabel(dlg, 'Text', 'Behaviour', 'FontWeight', 'bold', ...
                'Position', [15 208 200 20]);
            cbReport = uicheckbox(dlg, 'Text', 'Switch to Reports tab after each run', ...
                'Position', [15 184 300 22], ...
                'Value', getpref('nestapp','showReport',true));
            cbConfirm = uicheckbox(dlg, 'Text', 'Confirm before clearing pipeline', ...
                'Position', [15 160 300 22], ...
                'Value', getpref('nestapp','confirmClear',true));
            cbOverwrite = uicheckbox(dlg, 'Text', 'Overwrite existing report files (no timestamp)', ...
                'Position', [15 136 320 22], ...
                'Value', getpref('nestapp','overwriteReports',false));

            % --- Experimental Features section ---
            uilabel(dlg, 'Text', 'Experimental Features', 'FontWeight', 'bold', ...
                'Position', [15 104 200 20]);
            cbQuality = uicheckbox(dlg, ...
                'Text', 'Compute TEP quality metrics after each run', ...
                'Position', [15 80 300 22], ...
                'Value', getpref('nestapp','computeQuality',false));
            uilabel(dlg, ...
                'Text', ['Results are approximate and under active development. ' ...
                         'See documentation for interpretation.'], ...
                'FontColor', [0.5 0.5 0.5], 'FontSize', 10, ...
                'WordWrap', 'on', 'Position', [15 52 390 26]);

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
                % EEGLAB path
                ep = strtrim(fEeglab.Value);
                if ~isempty(ep)
                    if ~isfolder(ep)
                        uialert(dlg, 'EEGLAB path does not exist.', 'Invalid Path');
                        return
                    end
                    addpath(ep);
                end
                setpref('nestapp', 'eeglabPath',          ep);
                setpref('nestapp', 'lastDataFolder',      strtrim(fData.Value));
                setpref('nestapp', 'lastPipelineFolder',  strtrim(fPipeline.Value));
                setpref('nestapp', 'reportFolder',        strtrim(fReports.Value));
                setpref('nestapp', 'showReport',          cbReport.Value);
                setpref('nestapp', 'confirmClear',        cbConfirm.Value);
                setpref('nestapp', 'overwriteReports',    cbOverwrite.Value);
                setpref('nestapp', 'computeQuality',      cbQuality.Value);
                close(dlg);
            end
        end

        function updateReportsTabImpl(app)
        % UPDATEREPORTSTABIMPL  Refresh the Reports tab listbox from session and loaded reports.
        %   Combines app.allPipelineReports (from current run) with app.loadedReports
        %   (loaded from disk). Updates listbox labels and status text.
            allEntries = [app.allPipelineReports, app.loadedReports];
            n = numel(allEntries);
            if n == 0
                app.ReportsListBox.Items = {};
                app.ReportsListBox.ItemsData = {};
                app.ReportsStatusLabel.Text = 'No reports loaded.';
                app.ReportsTextArea.Value = '';
                app.ExportReportsCSVButton.Enable = 'off';
                app.CopyMethodsButton.Enable = 'off';
                return
            end

            labels = cell(1, n);
            for i = 1:n
                e = allEntries{i};
                [~, baseName] = fileparts(e.report.inputFile);
                try
                    dateLabel = datestr(e.report.processedAt, 'yyyy-mm-dd HH:MM'); %#ok<DATST>
                catch
                    dateLabel = '?';
                end
                labels{i} = sprintf('%s (%s)', baseName, dateLabel);
            end

            % Preserve selection index across refresh if still valid
            prevIdx = app.ReportsListBox.Value;
            app.ReportsListBox.Items = labels;
            app.ReportsListBox.ItemsData = num2cell(1:n);

            if isnumeric(prevIdx) && ~isempty(prevIdx) && prevIdx >= 1 && prevIdx <= n
                app.ReportsListBox.Value = prevIdx;
                app.ReportsTextArea.Value = allEntries{prevIdx}.text;
            else
                app.ReportsListBox.Value = n;
                app.ReportsTextArea.Value = allEntries{n}.text;
            end

            nSess   = numel(app.allPipelineReports);
            nLoaded = numel(app.loadedReports);
            parts   = {};
            if nSess   > 0; parts{end+1} = sprintf('%d from session', nSess);   end
            if nLoaded > 0; parts{end+1} = sprintf('%d from disk', nLoaded); end
            app.ReportsStatusLabel.Text = strjoin(parts, ', ');
            app.ExportReportsCSVButton.Enable = 'on';
            app.CopyMethodsButton.Enable = 'on';
        end

        function ReportsListBoxValueChanged(app, ~)
        % Callback — show the report text for the newly selected entry.
            idx = app.ReportsListBox.Value;
            if isempty(idx); return; end
            allEntries = [app.allPipelineReports, app.loadedReports];
            if isnumeric(idx) && idx >= 1 && idx <= numel(allEntries)
                app.ReportsTextArea.Value = allEntries{idx}.text;
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
                'Retention,Artifact reduction,Background,' ...
                'Reproducibility,AEP likeness\n']);

            for i = 1:numel(allEntries)
                r = allEntries{i}.report;
                [~, baseName] = fileparts(r.inputFile);
                try
                    dStr = datestr(r.processedAt, 'yyyy-mm-dd HH:MM:SS'); %#ok<DATST>
                catch
                    dStr = '?';
                end
                safeVal = @(axField) safeAxisVal(r.teps, axField);

                fprintf(fid, '%s,%s,%d,%d,%d,%d,%d,%s,%s,%s,%s,%s\n', ...
                    baseName, dStr, ...
                    r.channels.original, r.channels.final, ...
                    r.trials.original, r.trials.final, ...
                    r.ica.nRejected, ...
                    safeVal('retention'), safeVal('artifactReduction'), ...
                    safeVal('bgRestoration'), safeVal('reproducibility'), ...
                    safeVal('aepLikeness'));
            end
            fclose(fid);
            app.ReportsStatusLabel.Text = sprintf('CSV saved: %s', fname);

            function s = safeAxisVal(teps, axField)
                if isstruct(teps) && isfield(teps, axField) && ~isnan(teps.(axField).value)
                    s = sprintf('%.3f', teps.(axField).value);
                else
                    s = '';
                end
            end
        end

        function CopyMethodsButtonPushed(app, ~)
        % Build a brief methods paragraph from the selected report and copy to clipboard.
            idx = app.ReportsListBox.Value;
            if isempty(idx); return; end
            allEntries = [app.allPipelineReports, app.loadedReports];
            if ~isnumeric(idx) || idx < 1 || idx > numel(allEntries); return; end
            r = allEntries{idx}.report;

            parts = {};
            if r.channels.original > 0
                nRej  = r.channels.nRejected;
                nIntp = r.channels.nInterpolated;
                if nRej > 0 && nIntp > 0
                    parts{end+1} = sprintf('%d of %d channels were retained (%d rejected, %d interpolated)', ...
                        r.channels.final, r.channels.original, nRej, nIntp);
                elseif nRej > 0
                    parts{end+1} = sprintf('%d of %d channels were retained (%d rejected)', ...
                        r.channels.final, r.channels.original, nRej);
                elseif nIntp > 0
                    parts{end+1} = sprintf('%d channels were retained (%d interpolated)', ...
                        r.channels.final, nIntp);
                else
                    parts{end+1} = sprintf('%d channels were retained', r.channels.final);
                end
            end
            if r.trials.original > 0
                parts{end+1} = sprintf('%d of %d epochs were retained (%d rejected)', ...
                    r.trials.final, r.trials.original, r.trials.rejected);
            end
            if r.ica.nComponents > 0
                parts{end+1} = sprintf('%d ICA components were identified and %d removed', ...
                    r.ica.nComponents, r.ica.nRejected);
            end

            if isempty(parts)
                methodsText = 'TMS-EEG data were preprocessed using nestapp.';
            else
                methodsText = sprintf('TMS-EEG data were preprocessed using nestapp. %s.', ...
                    strjoin(parts, '; '));
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
                'nestapp — TMS-EEG Processing\n\n' ...
                'EEGLAB:  %s\n' ...
                'MATLAB:  %s\n\n' ...
                'Please cite:\n' ...
                'Rogasch et al. (2017) NeuroImage — TESA toolbox\n' ...
                'Delorme & Makeig (2004) J Neurosci Methods — EEGLAB'], ...
                eeglabVer, version);
            uialert(app.UIFigure, msg, 'About nestapp', 'Icon', 'info');
        end

        function rescaleComponents(app, sX, sY)
        % RESCALECOMPONENTS  Rescale all UI components proportionally.
        %   sX = newWidth/origW, sY = newHeight/origH.
        %   Position helper p() scales [x y w h] by [sX sY sX sY].
        %   Font helper fs() uses the smaller scale factor to prevent distortion.
        %
        %   The figure coordinate space is 549px (529px TabGroup + 20px status bar).
        %   The uimenu bar renders outside the coordinate space — MATLAB shifts the
        %   window position when the menu is created, keeping the coordinate height
        %   unchanged. No MENU_BAR_H correction needed here.
            STATUS_H = 20;

            sf = min(sX, sY);
            p  = @(o) round(o .* [sX, sY, sX, sY]);
            fs = @(o) max(8, round(o * sf));

            % Status bar: fixed height, pinned to bottom, full width
            app.StatusBar.Position = round([0, 0, 867*sX, STATUS_H]);
            % TabGroup fills all coordinate space above the status bar.
            figH = round(sY * app.originalSize(2));
            tabH = figH - STATUS_H;
            app.TabGroup.Position = [1, STATUS_H, round(867*sX), tabH];

            %% Cleaning Tab
            app.StepsListBox.Position             = p([10 173 207 294]);
            app.StepsListBox.FontSize             = fs(11);
            app.CommandDescriptionLabel.Position  = p([12 152 31 22]);
            app.CommandDescriptionLabel.FontSize  = fs(14);
            app.InfoTextArea.Position             = p([10 10 207 143]);
            app.StepsListBoxLabel.Position        = p([89 475 49 22]);
            app.StepsListBoxLabel.FontSize        = fs(16);
            app.SelectedListBox.Position          = p([230 104 215 360]);
            app.SelectedListBox.FontSize          = fs(11);
            app.MoveUpButton.Position             = p([306 56 66 36]);
            app.MoveDownButton.Position           = p([305 12 66 36]);
            app.AddButton.Position                = p([233 56 66 36]);
            app.RemoveButton.Position             = p([232 12 66 36]);
            app.SelectedListBoxLabel.Position     = p([278 472 119 25]);
            app.SelectedListBoxLabel.FontSize     = fs(16);
            app.UITable.Position                  = p([450 104 188 363]);
            app.DefaultValueButton.Position       = p([485 15 110 23]);
            app.TextArea.Position                 = p([450 46 188 56]);
            app.SelectedListBoxLabel_2.Position   = p([492 472 102 25]);
            app.SelectedListBoxLabel_2.FontSize   = fs(16);

            % Select Data panel + children (children coords are panel-relative)
            app.SelectDatatoPerformAnalysisPanel.Position = p([649 135 208 206]);
            app.SelectedFilesListBox.Position     = p([5 30 195 145]);
            app.SelectDataButton.Position         = p([5 5 195 23]);

            app.RunAnalysisButton.Position        = p([657 15 201 60]);
            app.RunAnalysisButton.FontSize        = fs(18);
            app.Image.Position                    = p([653 453 203 44]);
            app.NESTAPPLabel.Position             = p([785 448 71 22]);
            app.NESTAPPLabel.FontSize             = fs(14);

            app.ReStartStepsButton.Position       = p([658 91 201 36]);
            app.ReStartStepsButton.FontSize       = fs(18);

            %% Visualizing Tab
            % Three zones: left (x:0-340 electrode map), center (x:340-648 TEP/topo),
            % right (x:651-867 file selection). Bottom strip (y:0-165) holds controls.
            % TEP window slider lives above UIAxes; topoplot controls sit right of UIAxes2.

            % Left action column — bottom strip (x:5-145)
            app.PLOTTEPButton.Position            = p([5 130 140 30]);
            app.ShowComponentsButton.Position     = p([5 104 140 23]);
            app.ExportTEPFigureButton.Position    = p([5 76 140 23]);
            app.PlotEEGdataButton.Position        = p([5 50 108 23]);

            % Center-left controls — bottom strip (x:152-340)
            app.PlottingModeButtonGroup.Position  = p([152 88 150 67]);
            app.NewFigureButton.Position          = p([11 21 83 22]);
            app.AddtocurrentFigureButton.Position = p([11 -1 135 22]);
            app.EEGDatasetDropDownLabel.Position  = p([152 58 75 22]);
            app.EEGDatasetDropDown.Position       = p([230 58 100 22]);

            % Center column — TEP window slider above the TEP plot
            app.TEPWindowSliderLabel.Position     = p([669 420 85 16]);
            app.TEPWindowSlider.Position          = p([669 388 183 3]);
            app.UIAxes.Position                   = p([340 319 308 186]);
            app.TEPComponentTable.Position        = p([340 215 308 104]);
            app.TEPComponentTable.ColumnWidth     = num2cell(round([80, 100, 128] * sX));

            % Center column — topoplot with controls to its right
            app.UIAxes2.Position                  = p([340 65 200 148]);
            app.TOPOPLOTButton.Position           = p([543 167 104 44]);
            app.TopoplottimeSpinnerLabel.Position = p([543 145 103 20]);
            app.TopoplottimeSpinner.Position      = p([543 120 103 22]);
            app.WindowsizeforTopoplotLabel.Position = p([543 90 103 26]);
            app.WindowsizefortimeaveragedTopoplotEditField.Position = p([543 65 70 22]);
            app.Slider.Position                   = p([340 48 200 3]);

            % Head image (electrode map) — unchanged
            app.Image2.Position                   = p([-1 165 350 336]);

            % Right column — top: export/output; bottom: data selection (prototype layout)
            app.ExportTEPDataButton.Position      = p([669 468 183 23]);
            app.TEPvarNameEditFieldLabel.Position = p([669 438 80 22]);
            app.TEPvarNameEditField.Position      = p([754 438 98 22]);
            app.UseCurrentlyCleanedDataCheckBox.Position = p([671 325 180 22]);
            app.SelectDatatoVisulaizeTEPsPanel.Position  = p([651 230 208 90]);
            app.FolderEditField_2Label.Position   = p([1 41 40 22]);
            app.FolderEditField_2.Position        = p([49 41 145 22]);
            app.SelectDataButton_2.Position       = p([13 10 183 23]);
            app.FilesListBoxLabel.Position        = p([740 220 30 22]);
            app.FilesListBox.Position             = p([669 71 183 149]);
            app.SelectAllCheckBox.Position        = p([670 46 71 22]);
            app.DontfindcommonelectrodesCheckBox.Position = p([670 28 180 22]);
            app.ReLoadAvailableElectrodesButton.Position  = p([686 7 153 23]);

            %% Electrode buttons (Visualizing Tab — 64 buttons)
            app.AF3Button.Position   = p([108 410 25 23]);
            app.FP1Button.Position   = p([131 433 25 23]);
            app.FPZButton.Position   = p([161 439 25 23]);
            app.FP2Button.Position   = p([191 433 25 23]);
            app.AF4Button.Position   = p([215 409 25 23]);
            app.F8Button.Position    = p([266 390 25 23]);
            app.F6Button.Position    = p([240 385 25 23]);
            app.F4Button.Position    = p([214 380 25 23]);
            app.F2Button.Position    = p([187 385 25 23]);
            app.F5Button.Position    = p([80 385 25 23]);
            app.F3Button.Position    = p([107 380 25 23]);
            app.FZButton.Position    = p([161 386 25 23]);
            app.FC2Button.Position   = p([192 348 25 23]);
            app.FC4Button.Position   = p([222 348 25 23]);
            app.FC6Button.Position   = p([252 350 25 23]);
            app.F1Button.Position    = p([134 385 25 23]);
            app.C4Button.Position    = p([229 316 25 23]);
            app.C6Button.Position    = p([262 316 25 23]);
            app.FT8Button.Position   = p([285 354 25 23]);
            app.F7Button.Position    = p([55 391 25 23]);
            app.FC1Button.Position   = p([130 348 25 23]);
            app.FCZButton.Position   = p([161 349 25 23]);
            app.FC3Button.Position   = p([99 348 25 23]);
            app.C1Button.Position    = p([128 316 25 23]);
            app.CZButton.Position    = p([161 316 25 23]);
            app.C2Button.Position    = p([195 316 25 23]);
            app.CP3Button.Position   = p([97 283 25 23]);
            app.CP1Button.Position   = p([128 284 25 23]);
            app.CP2Button.Position   = p([192 284 25 23]);
            app.T8Button.Position    = p([293 316 25 23]);
            app.FT7Button.Position   = p([36 354 25 23]);
            app.FC5Button.Position   = p([68 350 25 23]);
            app.C5Button.Position    = p([59 316 25 23]);
            app.C3Button.Position    = p([93 316 25 23]);
            app.T7Button.Position    = p([26 316 25 23]);
            app.TP7Button.Position   = p([34 277 25 23]);
            app.CP5Button.Position   = p([62 279 25 23]);
            app.CPZButton.Position   = p([161 286 25 23]);
            app.CP4Button.Position   = p([224 283 25 23]);
            app.CP6Button.Position   = p([258 279 25 23]);
            app.TP8Button.Position   = p([288 277 25 23]);
            app.P8Button.Position    = p([275 237 25 23]);
            app.P3Button.Position    = p([105 251 25 23]);
            app.P1Button.Position    = p([132 251 25 23]);
            app.P2Button.Position    = p([188 251 25 23]);
            app.P7Button.Position    = p([46 237 25 23]);
            app.P5Button.Position    = p([76 246 25 23]);
            app.PZButton.Position    = p([161 250 25 23]);
            app.P4Button.Position    = p([216 251 25 23]);
            app.P6Button.Position    = p([245 246 25 23]);
            app.O1Button.Position    = p([128 179 25 23]);
            app.PO3Button.Position   = p([106 216 25 23]);
            app.POZButton.Position   = p([161 211 25 23]);
            app.PO4Button.Position   = p([217 217 25 23]);
            app.PO7Button.Position   = p([52 204 25 23]);
            app.PO5Button.Position   = p([78 215 25 23]);
            app.PO2Button.Position   = p([189 212 25 23]);
            app.PO8Button.Position   = p([270 202 25 23]);
            app.CB1Button.Position   = p([98 189 25 23]);
            app.OZButton.Position    = p([160 177 25 23]);
            app.O2Button.Position    = p([192 179 25 23]);
            app.CB2Button.Position   = p([224 189 25 23]);
            app.TP10Button.Position  = p([302 253 25 23]);
            app.TP9Button.Position   = p([20 253 25 23]);
            app.AFZButton.Position   = p([161 412 25 23]);
            app.AF7Button.Position   = p([79 419 25 23]);
            app.AF8Button.Position   = p([245 419 25 23]);
            app.PO1Button.Position   = p([133 212 25 23]);
            app.PO6Button.Position   = p([243 215 25 23]);

            %% Reports Tab
            app.ReportsListBoxLabel.Position    = p([5 472 205 22]);
            app.ReportsListBoxLabel.FontSize     = fs(16);
            app.ReportsListBox.Position          = p([5 73 205 393]);
            app.ReportsListBox.FontSize          = fs(11);
            app.LoadReportsButton.Position       = p([5 45 100 25]);
            app.RefreshReportsButton.Position    = p([110 45 100 25]);
            app.ReportsFolderLabel.Position      = p([5 25 205 18]);
            app.ReportsStatusLabel.Position      = p([5 5 205 18]);
            app.ExportReportsCSVButton.Position  = p([580 470 130 24]);
            app.CopyMethodsButton.Position       = p([715 470 147 24]);
            app.ReportsTextArea.Position         = p([220 10 637 457]);
            app.ReportsTextArea.FontSize         = fs(10);
        end

        function styleParamTable(app)
        % Grey out UITable rows whose Value is a placeholder or literal '[]'.
        % Placeholders start with '(' by convention (e.g. '(all channels)').
        % '[]' from old-format saved pipelines is treated the same way.
        % Call this after any assignment of UITable.Data to a parameter table.
            removeStyle(app.UITable);
            T = app.UITable.Data;
            if isempty(T) || ~istable(T); return; end
            grey = uistyle('FontColor', [0.6 0.6 0.6], 'FontAngle', 'italic');
            for row = 1:height(T)
                v = T.val{row};
                % Only scalar strings can be placeholders; skip numerics and arrays.
                if isscalar(v) && (isstring(v) || ischar(v))
                    sv = string(v);
                    isPlaceholder = strlength(sv) > 0 && startsWith(sv, '(');
                    isEmptyBracket = strcmp(sv, '[]');
                    if isPlaceholder || isEmptyBracket
                        addStyle(app.UITable, grey, 'cell', [row, 2]);
                    end
                end
            end
        end

        function LoadSelecEEGdata(app)
            % Ensure EEGLAB functions are on the path.  eeglab('nogui') is
            % normally called by runPipeline, but the Visualizing tab can be
            % used independently, so initialise on demand if needed.
            if ~exist('pop_loadset', 'file')
                eeglab('nogui');
            end
            for nfile = 1:numel(app.SelectedFilesforTEP)
                EEGaux = pop_loadset(app.SelectedFilesforTEP(nfile),app.PathofSelectedFilesforTEP);
                app.EEGofAllSelectedFiles{nfile} = EEGaux;
                app.EEGtime = EEGaux.times;
            end
            app.EEG_SelectedTEPFiles_Loaded = 1;
        end

        function app = LoadLabels(app)
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
                app.([upper(uncommon_labels.Items{nn}),'Button']).Enable = 'off';
                app.([upper(uncommon_labels.Items{nn}),'Button']).Value = 0;
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

            swin   = 5;
            app.TEP2Export = TEP_ROI;
            grandMean = mean(TEP_ROI, 1, 'omitmissing');
            TEP_ROISD = std(TEP_ROI, 1, 1) / sqrt(nFiles);
            Colr  = rand(1,3);
            meanx = smoothdata(grandMean,   'movmean', swin);
            sdx   = smoothdata(TEP_ROISD,  'movmean', swin);
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
                hold(app.UIAxes, 'on');
                fill(app.UIAxes, xf, yf, Colr(1,:), 'FaceAlpha', 0.5, 'LineStyle', 'none', 'HandleVisibility', 'off');
                plot(app.UIAxes, app.EEGtime, meanx, 'Color', Colr(1,:), 'LineWidth', 2, 'DisplayName', dispName);
                hold(app.UIAxes, 'off');
                xlim(app.UIAxes, app.DefaulTEPxLim);
            elseif app.AddtocurrentFigureButton.Value
                prevYLim = ylim(app.UIAxes);
                hold(app.UIAxes, 'on');
                fill(app.UIAxes, xf, yf, Colr(1,:), 'FaceAlpha', 0.5, 'LineStyle', 'none', 'HandleVisibility', 'off');
                plot(app.UIAxes, app.EEGtime, meanx, 'Color', Colr(1,:), 'LineWidth', 2, 'DisplayName', dispName);
                xlim(app.UIAxes, app.DefaulTEPxLim);
                % Expand y-axis to accommodate new data; never shrink existing range
                newYLim = ylim(app.UIAxes);
                ylim(app.UIAxes, [min(prevYLim(1), newYLim(1)), max(prevYLim(2), newYLim(2))]);
            end

            legend(app.UIAxes, 'show', 'Location', 'best');

            % Component detection on grand mean (runs regardless of toggle to cache peaks)
            app.tepPeaks = tepPeakFinder(grandMean, app.EEGtime);
            if app.ShowComponentsButton.Value
                overlayTEPComponents(app);
            end
            populateTEPComponentTable(app);
        end
        
        function EEG_topoplot(app)
            cla(app.UIAxes2)
            BIGEEG = [];
            IntRad = 0.55;
            avgType = 'movmean';
            swin = 5;
            if ~app.EEG_SelectedTEPFiles_Loaded
                LoadSelecEEGdata(app)
            end
            app = LoadLabels(app);
            BIGEEG = zeros(numel(app.Common_Labels.Items), length(app.EEGtime),numel(app.EEGofAllSelectedFiles));
            for nfile = 1:numel(app.EEGofAllSelectedFiles)
                EEGaux = app.EEGofAllSelectedFiles{1,nfile};
                ChansLocs = EEGaux.chanlocs;
                commonElectrodsInd = ismember({ChansLocs.labels},app.Common_Labels.Items);
                BIGEEG(:,:,nfile) = mean(EEGaux.data(commonElectrodsInd,:,:),3,"omitmissing");
                
            end
            ChansLocs(~commonElectrodsInd) = [];
            yp = smoothdata(mean(BIGEEG,3,"omitmissing")',avgType,swin)'; % Smooth the EEGdata along subjects
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
                'numcontour',5,'intsquare','on','style','map','conv', 'on', 'intrad',IntRad);axis auto
            colormap(app.UIAxes2,'hsv')
            % Copy contents back to app UIAxes
            cla(app.UIAxes2);
            copyobj(allchild(newAx), app.UIAxes2);

            % Clean up
            close(invisibleFig);  % Close hidden fig
            figure(oldFig);axis(app.UIAxes2,'auto')
            close
        end

        
        function selected = CheckifanyFileSelected(app)
            selected = 0;
            if isempty(app.SelectedFilesforTEP)
                warning('Please select at least a file to plot')
            else
                selected = 1;
            end


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
                    sprintf('%s\n%.0f ms\n%.1f µV', pk.name, pk.latencyMs, pk.amplitudeUV), ...
                    'FontSize', 7, 'HorizontalAlignment', 'center', ...
                    'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'top');
            end
            hold(ax, 'off');
        end

        function populateTEPComponentTable(app)
        % Fill TEPComponentTable from app.tepPeaks.
        % Shows '—' for components not found.
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
                    tableData{i, 2} = '—';
                    tableData{i, 3} = '—';
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
            app.ShowComponentsButton.Tooltip  = 'Detect and overlay N15/P30/N45/P60/N100/P180 peaks on the TEP plot';
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

            for i = 1:numel(steps)
                fields = fieldnames(steps(i).defaults);
                var = cell(numel(fields), 1);
                val = cell(numel(fields), 1);

                % Build key→placeholder lookup for this step
                paramKeys   = {steps(i).params.key};
                paramPH     = {steps(i).params.placeholder};

                for j = 1:numel(fields)
                    % Show friendlyName (unit) when available; raw key otherwise.
                    fn = fields{j};
                    pIdx = find(strcmp(paramKeys, fn), 1);
                    if ~isempty(pIdx)
                        p = steps(i).params(pIdx);
                        if isempty(p.unit)
                            var{j} = string(p.friendlyName);
                        else
                            var{j} = string([p.friendlyName, ' (', p.unit, ')']);
                        end
                    else
                        var{j} = string(fn);
                    end
                    V = steps(i).defaults.(fields{j});
                    if ischar(V)
                        val{j} = string(V);
                    elseif isempty(V)
                        % Show a human-readable placeholder instead of '[]'.
                        % Look up param-specific text; fall back to generic.
                        pIdx = find(strcmp(paramKeys, fields{j}), 1);
                        if ~isempty(pIdx) && ~isempty(paramPH{pIdx})
                            val{j} = string(paramPH{pIdx});
                        else
                            val{j} = "(not set)";
                        end
                    elseif iscell(V)
                        val{j} = [string(V)];
                    else
                        % Numeric scalar or array. Store as a mat2str string so
                        % UITable can render it — arrays like [1 6] otherwise
                        % display as '[]'. The preprocessing loop in runPipeline
                        % converts these strings back to numeric via str2num.
                        val{j} = string(mat2str(V));
                    end
                end
                app.info{i} = steps(i).info;
                app.DefaultsVal{i} = table(var, val);
                app.defaultParamKeys{i} = fields;  % raw EEGLAB keys, same row order
            end

            app.SelectedListBox.Items(:) = [];
            app.SelectedListBox.ItemsData(:) = [];
            app.UITable.Data = [];
            app.ItemNum = 1;
            app.needchanloc = 1;
            app.originalSize = app.UIFigure.Position(3:4);
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
                AddButtonPushed(app, []);
            end
            app.lastStepClick = t;
        end

        % Value changed function: StepsListBox
        function StepsListBoxValueChanged(app, event)
            value = app.StepsListBox.Value;
            ind = find(ismember(app.StepsListBox.Items,value));
            app.InfoTextArea.Value = string(app.info{ind});
            app.selectedItem = [];
        end

        % Button pushed function: AddButton
        function AddButtonPushed(app, event)
            if isempty(app.SelectedListBox.Items) || isempty(app.SelectedListBox.Items{1})
                ind = find(ismember(app.StepsListBox.Items,app.StepsListBox.Value));
                app.SelectedListBox.Items{1} = app.StepsListBox.Value;
                app.SelectedListBox.ItemsData{1} = ['item',num2str(1)];
                app.ChangedVal{1}       = app.DefaultsVal{ind};
                app.stepParamKeys{1}    = app.defaultParamKeys{ind};
            % elseif isempty(app.SelectedListBox.Items{1})
            %     ind = find(ismember(app.StepsListBox.Items,app.StepsListBox.Value));
            %     app.SelectedListBox.Items{1} = app.StepsListBox.Value;
            %     app.SelectedListBox.ItemsData{1} = ['item',num2str(1)];
            %     app.ChangedVal{1}=app.DefaultsVal{ind};
            else
                ind = find(ismember(app.StepsListBox.Items,app.StepsListBox.Value));
                app.ItemNum = numel(app.SelectedListBox.Items);
                app.SelectedListBox.Items{app.ItemNum+1} = app.StepsListBox.Value;
                app.SelectedListBox.ItemsData{app.ItemNum+1} = ['item',num2str(app.ItemNum+1)];
                app.ChangedVal{app.ItemNum+1}    = app.DefaultsVal{ind};
                app.stepParamKeys{app.ItemNum+1} = app.defaultParamKeys{ind};
            end
            app.pipelineDirty = true;
            updateStatusBar(app);
        end

        % Button pushed function: MoveUpButton
        function MoveUpButtonPushed(app, event)
            value = app.SelectedListBox.Value;
            ind1 = find(ismember(app.SelectedListBox.ItemsData,value));
            if ind1 ~= 1
                A = app.SelectedListBox.Items{ind1-1};
                B = app.SelectedListBox.Items{ind1};
                app.SelectedListBox.Items{ind1-1} = B;
                app.SelectedListBox.Items{ind1} = A;

                app.SelectedListBox.ItemsData{ind1-1} = strcat('Item',num2str(ind1-1));
                app.SelectedListBox.ItemsData{ind1} = strcat('Item',num2str(ind1));
                C = app.ChangedVal{ind1-1};
                D = app.ChangedVal{ind1};
                app.ChangedVal{ind1-1} = D;
                app.ChangedVal{ind1}   = C;
                Ck = app.stepParamKeys{ind1-1};
                Dk = app.stepParamKeys{ind1};
                app.stepParamKeys{ind1-1} = Dk;
                app.stepParamKeys{ind1}   = Ck;
            end
            app.SelectedListBox.Value=app.SelectedListBox.ItemsData{ind1-1};
            app.pipelineDirty = true;
            updateStatusBar(app);
        end

        % Button pushed function: SavePipelineButton
        function SavePipelineButtonPushed(app, event) %#ok<INUSD>
            PLItems = app.SelectedListBox.Items;     %#ok<NASGU>
            PLItemsData = app.SelectedListBox.ItemsData; %#ok<NASGU>
            VarIns = app.ChangedVal;                 %#ok<NASGU>
            ParamKeys = app.stepParamKeys;           %#ok<NASGU>
            startFolder = getpref('nestapp', 'lastPipelineFolder', '');
            uisave({'PLItems','PLItemsData','VarIns','ParamKeys'}, ...
                fullfile(startFolder, '*.mat'));
            % uisave does not return the chosen path; dirty flag and name update
            % deferred until Save Pipeline moves to a proper uiputfile dialog.
        end

        % Button pushed function: RemoveButton
        function RemoveButtonPushed(app, event)
            ind = find(ismember(app.SelectedListBox.ItemsData,app.SelectedListBox.Value));
            app.SelectedListBox.Items(ind) = [];
            app.SelectedListBox.ItemsData(ind) = [];
            app.ChangedVal(ind)    = [];
            app.stepParamKeys(ind) = [];
            for i = ind : length(app.SelectedListBox.ItemsData)
                app.SelectedListBox.ItemsData{i} = ['Item',num2str(i)];
            end
            app.pipelineDirty = true;
            updateStatusBar(app);
        end

        % Button pushed function: MoveDownButton
        function MoveDownButtonPushed(app, event)
            value = app.SelectedListBox.Value;
            ind1 = find(ismember(app.SelectedListBox.ItemsData,value));
            if ind1 ~= numel(app.SelectedListBox.Items)
                A = app.SelectedListBox.Items{ind1+1};
                B = app.SelectedListBox.Items{ind1};
                app.SelectedListBox.Items{ind1+1} = B;
                app.SelectedListBox.Items{ind1} = A;

                app.SelectedListBox.ItemsData{ind1+1} = strcat('Item',num2str(ind1+1));
                app.SelectedListBox.ItemsData{ind1} = strcat('Item',num2str(ind1));
                C = app.ChangedVal{ind1+1};
                D = app.ChangedVal{ind1};
                app.ChangedVal{ind1+1} = D;
                app.ChangedVal{ind1}   = C;
                Ck = app.stepParamKeys{ind1+1};
                Dk = app.stepParamKeys{ind1};
                app.stepParamKeys{ind1+1} = Dk;
                app.stepParamKeys{ind1}   = Ck;
            end
            app.SelectedListBox.Value=app.SelectedListBox.ItemsData{ind1+1};
            app.pipelineDirty = true;
            updateStatusBar(app);
        end

        % Button pushed function: LoadPipelineButton
        function LoadPipelineButtonPushed(app, event)
            startFolder = getpref('nestapp', 'lastPipelineFolder', '');
            [pName,pPath] = uigetfile('*.mat', 'Load Pipeline', startFolder);
            pipeline = load([pPath,pName],'-mat');
            app.SelectedListBox.Items    = pipeline.PLItems;
            app.SelectedListBox.ItemsData= pipeline.PLItemsData;
            app.ChangedVal               = pipeline.VarIns;
            app.UITable.Data             = [];

            % Restore or rebuild stepParamKeys.
            % New-format files have ParamKeys saved directly.
            % Old-format files have raw EEGLAB keys in VarIns{k}.var — upgrade them
            % to friendly names while extracting the raw keys into stepParamKeys.
            if isfield(pipeline, 'ParamKeys')
                app.stepParamKeys = pipeline.ParamKeys;
            else
                reg = stepRegistry();
                app.stepParamKeys = cell(numel(pipeline.PLItems), 1);
                for k = 1:numel(pipeline.PLItems)
                    stepName = pipeline.PLItems{k};
                    regIdx   = find(strcmp({reg.name}, stepName), 1);
                    T        = app.ChangedVal{k};
                    nRows    = height(T);
                    rawKeys  = cell(nRows, 1);
                    for row = 1:nRows
                        rawKey = char(T.var{row});
                        rawKeys{row} = rawKey;
                        if ~isempty(regIdx)
                            pIdx = find(strcmp({reg(regIdx).params.key}, rawKey), 1);
                            if ~isempty(pIdx)
                                p = reg(regIdx).params(pIdx);
                                % Upgrade label: raw key → friendly name (unit)
                                if isempty(p.unit)
                                    T.var{row} = string(p.friendlyName);
                                else
                                    T.var{row} = string([p.friendlyName, ' (', p.unit, ')']);
                                end
                                % Upgrade value: '[]' → descriptive placeholder
                                if isscalar(T.val{row}) && strcmp(string(T.val{row}), '[]')
                                    if ~isempty(p.placeholder)
                                        T.val{row} = string(p.placeholder);
                                    else
                                        T.val{row} = "(not set)";
                                    end
                                end
                            end
                        end
                    end
                    app.stepParamKeys{k} = rawKeys;
                    app.ChangedVal{k}    = T;
                end
            end
            if ~isequal(pPath, 0)
                setpref('nestapp', 'lastPipelineFolder', pPath);
                pushRecent(app, 'recentPipelines', fullfile(pPath, pName));
                buildRecentPipelinesMenu(app);
                [~, nm, ~] = fileparts(pName);
                app.pipelineName  = nm;
                app.pipelineDirty = false;
                updateStatusBar(app);
            end
        end

        % Button pushed function: SelectDataButton
        function SelectDataButtonPushed(app, event)
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

                assignin('base','files',app.file)
                assignin('base','paths',app.path)

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
            app.SelectedListBox.Items(:)  = [];
            app.SelectedListBox.ItemsData(:) = [];
            app.UITable.Data    = [];
            app.ChangedVal      = {};
            app.stepParamKeys   = {};
            app.ItemNum         = 0;
            app.nstep           = 1;
            app.needchanloc     = 1;
        end

        % Menu selected function: Load Template...
        function LoadTemplateMenuSelected(app, ~)
        % LOADTEMPLATEMENUSELECTED  Show a template picker and load the chosen template.
        %   Clears the current pipeline (without confirmation) and populates
        %   the step list from the template, applying any parameter overrides.
            templates = pipelineTemplates();
            if isempty(templates); return; end

            % Build a small modal picker
            names = {templates.name};
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
                t = templates(idx);

                % Clear pipeline without confirmation
                app.SelectedListBox.Items(:)  = [];
                app.SelectedListBox.ItemsData(:) = [];
                app.UITable.Data    = [];
                app.ChangedVal      = {};
                app.stepParamKeys   = {};
                app.ItemNum         = 0;
                app.nstep           = 1;
                app.needchanloc     = 1;
                clc

                % Add each step and apply overrides
                for si = 1:numel(t.steps)
                    stepLabel = t.steps{si};
                    regIdx = find(strcmp(app.StepsListBox.Items, stepLabel), 1);
                    if isempty(regIdx)
                        warning('nestapp:template', ...
                            'Template step "%s" not found in registry — skipped.', stepLabel);
                        continue
                    end
                    % Add the step using the same logic as AddButtonPushed
                    n = numel(app.SelectedListBox.Items);
                    if n == 0 || (n == 1 && isempty(app.SelectedListBox.Items{1}))
                        pos = 1;
                    else
                        pos = n + 1;
                    end
                    app.SelectedListBox.Items{pos}    = stepLabel;
                    app.SelectedListBox.ItemsData{pos} = ['item' num2str(pos)];
                    app.ChangedVal{pos}    = app.DefaultsVal{regIdx};
                    app.stepParamKeys{pos} = app.defaultParamKeys{regIdx};

                    % Apply parameter overrides for this step
                    if si <= numel(t.overrides) && ~isempty(fieldnames(t.overrides{si}))
                        ov = t.overrides{si};
                        ovFields = fieldnames(ov);
                        keys = app.stepParamKeys{pos};
                        T    = app.ChangedVal{pos};
                        for fi = 1:numel(ovFields)
                            rawKey = ovFields{fi};
                            row = find(strcmp(keys, rawKey), 1);
                            if isempty(row); continue; end
                            v = ov.(rawKey);
                            if isnumeric(v)
                                T.val{row} = v;
                            else
                                T.val{row} = string(v);
                            end
                        end
                        app.ChangedVal{pos} = T;
                    end
                end

                app.pipelineName  = t.name;
                app.pipelineDirty = true;
                updateStatusBar(app);
                % Show the first step's params if any steps were added
                if ~isempty(app.SelectedListBox.Items) && ~isempty(app.SelectedListBox.Items{1})
                    app.SelectedListBox.Value = app.SelectedListBox.ItemsData{1};
                    app.UITable.Data = app.ChangedVal{1};
                    styleParamTable(app);
                end
            end
        end

        % Button pushed function: RunAnalysisButton
        function RunAnalysisButtonPushed(app, event)
            app.RunAnalysisButton.Text = {'Run';'Analysis'};
            app.needchanloc = 1;
            try
                runPipeline(app);
            catch err
                uialert(app.UIFigure, err.message, 'Pipeline Error', 'Icon', 'error');
                return
            end
            if app.UseCurrentlyCleanedDataCheckBox.Value
                UseCurrentlyCleanedDataCheckBoxValueChanged(app)
            end
        end

        % Value changed function: TextArea
        function TextAreaValueChanged(app, event)
            value = app.TextArea.Value;
            value(strcmp(value,''))=[];
            if app.convert
                app.UITable.Data{app.selectedItem(1),app.selectedItem(2)} = {str2double(value(:))};
            else
                if numel(value)>1
                    app.UITable.Data{app.selectedItem(1),app.selectedItem(2)} = {value};
                else
                    app.UITable.Data{app.selectedItem(1),app.selectedItem(2)} = value;
                end
            end
            val = app.SelectedListBox.Value;
            indNum2 = find(ismember(app.SelectedListBox.ItemsData,val));
            app.ChangedVal{indNum2} = app.UITable.Data;

        end

        % Cell selection callback: UITable
        function UITableCellSelection(app, event)
            if isempty(event.Indices); return; end
            app.selectedItem = event.Indices;
            if app.selectedItem(1) > height(app.UITable.Data); return; end
            x = app.UITable.Data{app.selectedItem(1),app.selectedItem(2)};
            y = x{:};
            app.convert = isnumeric(y);
            if app.convert
                app.TextArea.Value = string(y);
            else
                app.TextArea.Value = y;
            end
        end

        % Button pushed function: DefaultValueButton
        function DefaultValueButtonPushed(app, event)
            ind2 = find(ismember(app.SelectedListBox.ItemsData,app.SelectedListBox.Value));
            ind1 = find(ismember(app.StepsListBox.Items,app.SelectedListBox.Items{ind2}));
            app.ChangedVal{ind2}    = app.DefaultsVal{ind1};
            app.stepParamKeys{ind2} = app.defaultParamKeys{ind1};
            app.UITable.Data        = app.DefaultsVal{ind1};
            styleParamTable(app);
            app.pipelineDirty = true;
            updateStatusBar(app);
        end

        % Size changed function: UIFigure
        function UIFigureSizeChanged(app, ~)
            if isempty(app.originalSize); return; end
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
        end

        % Cell edit callback: UITable
        function UITableCellEdit(app, event)
            indices = event.Indices; %#ok<NASGU>
            newData = event.NewData; %#ok<NASGU>

            value = app.SelectedListBox.Value;
            indNum2 = find(ismember(app.SelectedListBox.ItemsData,value));
            app.ChangedVal{indNum2} = app.UITable.Data;
            app.pipelineDirty = true;
            updateStatusBar(app);
        end

        % Value changed function: SelectedListBox
        function SelectedListBoxValueChanged(app, event)
            value = app.SelectedListBox.Value;
            indNum2 = find(ismember(app.SelectedListBox.ItemsData,value));
            ItemName = app.SelectedListBox.Items{indNum2};
            indNum1 = find(ismember(app.StepsListBox.Items,ItemName));
            app.UITable.Data = app.DefaultsVal{indNum1};

            if ~isequal(table2struct(app.ChangedVal{indNum2}),table2struct(app.DefaultsVal{indNum1}))
                app.UITable.Data = app.ChangedVal{indNum2};
            end

            styleParamTable(app);
            app.TextArea.Value='';

        end

        % Button pushed function: PLOTTEPButton
        function PLOTTEPButtonPushed(app, event)
            if ~CheckifanyFileSelected(app)
                warning('Please select at least a file to plot the TEP!');
            else
                app = LoadLabels(app);
                findTEPelecs(app);
                plotTEP(app)
                app.Slider.Limits = [app.EEGtime(1) app.EEGtime(end)];
                app.TEPWindowSlider.Limits = [app.EEGtime(1) app.EEGtime(end)];
                app.TEPCreated = 1;
                app.TEPWindowSlider.Value = app.DefaulTEPxLim;
                
                app.ExportTEPDataButton.Enable = 'on';
                app.TEPvarNameEditField.Enable = 'on';
            end
            
            % app.Slider.Value = 
        end

        % Value changed function: UseCurrentlyCleanedDataCheckBox
        function UseCurrentlyCleanedDataCheckBoxValueChanged(app, event)
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

        % Value changed function: FolderEditField_2
        function FolderEditField_2ValueChanged(app, event)
            if ~isempty(app.cleanedName) && ~isempty(app.path)
                app.FolderEditField_2.Value = app.path;
            else
                app.FolderEditField_2.Value = '';
            end

        end

        % Value changing function: Slider
        function SliderValueChanging(app, event)
            changingValue = event.Value;
            app.TopoplottimeSpinner.Value = round(changingValue);

        end

        % Button pushed function: SelectDataButton_2
        function SelectDataButton_2Pushed(app, event)
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
            app.EEG_SelectedTEPFiles_Loaded = 0;
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
        function TOPOPLOTButtonPushed(app, event)
            if CheckifanyFileSelected(app)
                EEG_topoplot(app)
            end
        end

        % Value changed function: SelectAllCheckBox
        function SelectAllCheckBoxValueChanged(app, event)
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
        function DontfindcommonelectrodesCheckBoxValueChanged(app, event)
            value = app.DontfindcommonelectrodesCheckBox.Value;
            if ~value
                app.ReLoadAvailableElectrodesButton.Enable = 1;
            else
                app.ReLoadAvailableElectrodesButton.Enable = 0;
            end
        end

        % Button pushed function: ReLoadAvailableElectrodesButton
        function ReLoadAvailableElectrodesButtonPushed(app, event)
            if CheckifanyFileSelected(app)
                app = LoadLabels(app);
            end

        end

        % Button pushed function: ExportTEPFigureButton
        function ExportTEPFigureButtonPushed(app, event)
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
        function EEGDatasetDropDownValueChanged(app, event)
            value = app.EEGDatasetDropDown.Value;
            
        end

        % Button pushed function: PlotEEGdataButton
        function PlotEEGdataButtonPushed(app, event)
            subInd = strcmpi(app.SelectedFilesforTEP, app.EEGDatasetDropDown.Value);
            if CheckifanyFileSelected(app)
                if ~app.EEG_SelectedTEPFiles_Loaded
                    LoadSelecEEGdata(app)
                    pop_eegplot(app.EEGofAllSelectedFiles{subInd},1,1,1)
                else
                    pop_eegplot(app.EEGofAllSelectedFiles{subInd},1,1,1)
                end
            end
        end

        % Button pushed function: ExportTEPDataButton
        function ExportTEPDataButtonPushed(app, event)
            assignin('base', app.TEPvarNameEditField.Value, app.TEP2Export)
        end

        % Value changed function: TEPvarNameEditField
        function TEPvarNameEditFieldValueChanged(app, event)
            value = app.TEPvarNameEditField.Value;

        end

        % Value changed function: ShowComponentsButton
        function ShowComponentsButtonValueChanged(app, ~)
            if app.TEPCreated
                if app.ShowComponentsButton.Value
                    overlayTEPComponents(app);
                else
                    % Replot without overlays — cla then replot
                    plotTEP(app);
                end
            end
        end

    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            % Height = 529px tab area + 20px status bar = 549px.
            % uimenu renders outside the coordinate space (MATLAB shifts the window
            % upward when the menu is created; coordinate height stays unchanged).
            app.UIFigure.Position = [100 100 867 549];
            app.UIFigure.Name = 'nestapp — TMS-EEG Processing';
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @UIFigureSizeChanged, true);

            % Create menu bar
            mFile = uimenu(app.UIFigure, 'Text', 'File');
            uimenu(mFile, 'Text', 'Open Data...', 'Accelerator', 'O', ...
                'MenuSelectedFcn', createCallbackFcn(app, @SelectDataButtonPushed, true));
            app.MenuRecentFiles = uimenu(mFile, 'Text', 'Recent Files');
            uimenu(mFile, 'Text', 'Load Pipeline...', 'Accelerator', 'L', 'Separator', 'on', ...
                'MenuSelectedFcn', createCallbackFcn(app, @LoadPipelineButtonPushed, true));
            uimenu(mFile, 'Text', 'Save Pipeline', 'Accelerator', 'S', ...
                'MenuSelectedFcn', createCallbackFcn(app, @SavePipelineButtonPushed, true));
            app.MenuRecentPipelines = uimenu(mFile, 'Text', 'Recent Pipelines');
            uimenu(mFile, 'Text', 'Load Template...', 'Separator', 'on', ...
                'MenuSelectedFcn', createCallbackFcn(app, @LoadTemplateMenuSelected, true));
            uimenu(mFile, 'Text', 'Exit', 'Separator', 'on', ...
                'MenuSelectedFcn', @(~,~) delete(app));

            mSettings = uimenu(app.UIFigure, 'Text', 'Settings');
            uimenu(mSettings, 'Text', 'Preferences...', ...
                'MenuSelectedFcn', createCallbackFcn(app, @openPreferencesMenu, true));

            mHelp = uimenu(app.UIFigure, 'Text', 'Help');
            uimenu(mHelp, 'Text', 'About nestapp', ...
                'MenuSelectedFcn', createCallbackFcn(app, @showAboutMenu, true));

            % Create status bar — pinned to bottom of UIFigure, visible on both tabs
            app.StatusBar = uilabel(app.UIFigure);
            app.StatusBar.Position = [0 0 867 20];
            app.StatusBar.BackgroundColor = [0.90 0.90 0.90];
            app.StatusBar.FontSize = 10;
            app.StatusBar.Text = '  Ready';
            app.StatusBar.HorizontalAlignment = 'left';

            % Create TabGroup — starts at y=20 to leave room for status bar
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.AutoResizeChildren = 'off';
            app.TabGroup.Position = [1 20 867 529];

            % Create CleaningTab
            app.CleaningTab = uitab(app.TabGroup);
            app.CleaningTab.AutoResizeChildren = 'off';
            app.CleaningTab.Title = 'Cleaning';

            % Create StepsListBox
            app.StepsListBox = uilistbox(app.CleaningTab);
            app.StepsListBox.Items = {'Load Data', 'Load Channel Location', 'Save New Set', 'Choose Data Set', 'Remove un-needed Channels', 'Remove Baseline', 'Remove Bad Channels', 'Clean Artifacts', 'Automatic Continuous Rejection', 'Automatic Cleaning Data', 'De-Trend Epoch', 'TESA De-Trend', 'Re-Sample', 'Re-Reference', 'Frequency Filter (CleanLine)', 'Frequency Filter (TESA)', 'Frequency Filter', 'Run ICA', 'Run TESA ICA', 'Label ICA Components', 'Flag ICA Components for Rejection', 'Remove Flagged ICA Components', 'Remove ICA Components (TESA)', 'Epoching', 'Remove Bad Epoch', 'Find TMS Pulses (TESA)', 'Remove TMS Artifacts (TESA)', 'Fix TMS Pulse (TESA)', 'Interpolate Channels', 'Interpolate Missing Data (TESA)', 'Find Artifacts EDM (TESA)', 'SSP SIR', 'Median Filter 1D', 'Extract TEP (TESA)', 'Find TEP Peaks (TESA)', 'TEP Peak Output', 'Remove Recording Noise (SOUND)', 'Visualize EEG Data', 'Manual Command'};
            app.StepsListBox.ValueChangedFcn = createCallbackFcn(app, @StepsListBoxValueChanged, true);
            app.StepsListBox.FontSize = 11;
            app.StepsListBox.ClickedFcn = createCallbackFcn(app, @StepsListBoxClicked, true);
            app.StepsListBox.Position = [10 173 207 294];
            app.StepsListBox.Value = 'Load Data';

            % Create CommandDescriptionLabel
            app.CommandDescriptionLabel = uilabel(app.CleaningTab);
            app.CommandDescriptionLabel.FontSize = 14;
            app.CommandDescriptionLabel.FontWeight = 'bold';
            app.CommandDescriptionLabel.Position = [12 152 31 22];
            app.CommandDescriptionLabel.Text = 'Info';

            % Create InfoTextArea
            app.InfoTextArea = uitextarea(app.CleaningTab);
            app.InfoTextArea.Editable = 'off';
            app.InfoTextArea.Position = [10 10 207 143];

            % Create StepsListBoxLabel
            app.StepsListBoxLabel = uilabel(app.CleaningTab);
            app.StepsListBoxLabel.FontSize = 16;
            app.StepsListBoxLabel.FontWeight = 'bold';
            app.StepsListBoxLabel.Position = [89 475 49 22];
            app.StepsListBoxLabel.Text = 'Steps';

            % Create SelectedListBox
            app.SelectedListBox = uilistbox(app.CleaningTab);
            app.SelectedListBox.Items = {''};
            app.SelectedListBox.ValueChangedFcn = createCallbackFcn(app, @SelectedListBoxValueChanged, true);
            app.SelectedListBox.FontSize = 11;
            app.SelectedListBox.Position = [230 104 215 360];
            app.SelectedListBox.Value = '';

            % Create MoveUpButton
            app.MoveUpButton = uibutton(app.CleaningTab, 'push');
            app.MoveUpButton.ButtonPushedFcn = createCallbackFcn(app, @MoveUpButtonPushed, true);
            app.MoveUpButton.BackgroundColor = [0.8 0.8 0.8];
            app.MoveUpButton.Position = [306 56 66 36];
            app.MoveUpButton.Text = {'Move'; 'Up'};

            % Create MoveDownButton
            app.MoveDownButton = uibutton(app.CleaningTab, 'push');
            app.MoveDownButton.ButtonPushedFcn = createCallbackFcn(app, @MoveDownButtonPushed, true);
            app.MoveDownButton.BackgroundColor = [0.8 0.8 0.8];
            app.MoveDownButton.Position = [305 12 66 36];
            app.MoveDownButton.Text = {'Move'; 'Down'};

            % Create AddButton
            app.AddButton = uibutton(app.CleaningTab, 'push');
            app.AddButton.ButtonPushedFcn = createCallbackFcn(app, @AddButtonPushed, true);
            app.AddButton.BackgroundColor = [0.8 0.8 0.8];
            app.AddButton.Position = [233 56 66 36];
            app.AddButton.Text = 'Add';

            % Create RemoveButton
            app.RemoveButton = uibutton(app.CleaningTab, 'push');
            app.RemoveButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveButtonPushed, true);
            app.RemoveButton.BackgroundColor = [0.8 0.8 0.8];
            app.RemoveButton.Position = [232 12 66 36];
            app.RemoveButton.Text = 'Remove';

            % Create SelectedListBoxLabel
            app.SelectedListBoxLabel = uilabel(app.CleaningTab);
            app.SelectedListBoxLabel.FontSize = 16;
            app.SelectedListBoxLabel.FontWeight = 'bold';
            app.SelectedListBoxLabel.Position = [278 472 119 25];
            app.SelectedListBoxLabel.Text = 'Selected Steps';

            % Create UITable
            app.UITable = uitable(app.CleaningTab);
            app.UITable.ColumnName = {'Properties'; 'Value'};
            app.UITable.RowName = {};
            app.UITable.CellEditCallback = createCallbackFcn(app, @UITableCellEdit, true);
            app.UITable.CellSelectionCallback = createCallbackFcn(app, @UITableCellSelection, true);
            app.UITable.Position = [450 104 188 363];

            % Create DefaultValueButton
            app.DefaultValueButton = uibutton(app.CleaningTab, 'push');
            app.DefaultValueButton.ButtonPushedFcn = createCallbackFcn(app, @DefaultValueButtonPushed, true);
            app.DefaultValueButton.BackgroundColor = [0.8 0.8 0.8];
            app.DefaultValueButton.Position = [485 15 110 23];
            app.DefaultValueButton.Text = 'Default Value';

            % Create TextArea
            app.TextArea = uitextarea(app.CleaningTab);
            app.TextArea.ValueChangedFcn = createCallbackFcn(app, @TextAreaValueChanged, true);
            app.TextArea.Position = [450 46 188 56];

            % Create SelectedListBoxLabel_2
            app.SelectedListBoxLabel_2 = uilabel(app.CleaningTab);
            app.SelectedListBoxLabel_2.FontSize = 16;
            app.SelectedListBoxLabel_2.FontWeight = 'bold';
            app.SelectedListBoxLabel_2.Position = [492 472 102 25];
            app.SelectedListBoxLabel_2.Text = 'Parameter(s)';

            % Create SelectDatatoPerformAnalysisPanel
            % Panel expanded to show file listbox (was 116px tall, now 206px)
            app.SelectDatatoPerformAnalysisPanel = uipanel(app.CleaningTab);
            app.SelectDatatoPerformAnalysisPanel.AutoResizeChildren = 'off';
            app.SelectDatatoPerformAnalysisPanel.BorderType = 'none';
            app.SelectDatatoPerformAnalysisPanel.Title = 'Select Data to Perform Analysis';
            app.SelectDatatoPerformAnalysisPanel.Position = [649 135 208 206];

            % Create SelectedFilesListBox — shows all queued files
            app.SelectedFilesListBox = uilistbox(app.SelectDatatoPerformAnalysisPanel);
            app.SelectedFilesListBox.Items = {};
            app.SelectedFilesListBox.Position = [5 30 195 145];
            app.SelectedFilesListBox.FontSize = 10;

            % Create SelectDataButton (full-width, at bottom of panel)
            app.SelectDataButton = uibutton(app.SelectDatatoPerformAnalysisPanel, 'push');
            app.SelectDataButton.ButtonPushedFcn = createCallbackFcn(app, @SelectDataButtonPushed, true);
            app.SelectDataButton.Position = [5 5 195 23];
            app.SelectDataButton.Text = 'Select Data Files...';

            % Create RunAnalysisButton
            app.RunAnalysisButton = uibutton(app.CleaningTab, 'push');
            app.RunAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @RunAnalysisButtonPushed, true);
            app.RunAnalysisButton.BackgroundColor = [0.20 0.55 0.20];
            app.RunAnalysisButton.FontColor = [1 1 1];
            app.RunAnalysisButton.FontSize = 18;
            app.RunAnalysisButton.FontWeight = 'bold';
            app.RunAnalysisButton.Position = [657 15 201 60];
            app.RunAnalysisButton.Text = 'Run Analysis';

            % Create Image
            app.Image = uiimage(app.CleaningTab);
            app.Image.Position = [653 453 203 44];
            app.Image.ImageSource = fullfile(pathToMLAPP, 'LogoNest.jpg');

            % Create NESTAPPLabel
            app.NESTAPPLabel = uilabel(app.CleaningTab);
            app.NESTAPPLabel.FontSize = 14;
            app.NESTAPPLabel.FontWeight = 'bold';
            app.NESTAPPLabel.FontAngle = 'italic';
            app.NESTAPPLabel.Position = [785 448 71 22];
            app.NESTAPPLabel.Text = 'NESTAPP';

            % Create ReStartStepsButton
            app.ReStartStepsButton = uibutton(app.CleaningTab, 'push');
            app.ReStartStepsButton.ButtonPushedFcn = createCallbackFcn(app, @ReStartStepsButtonPushed, true);
            app.ReStartStepsButton.BackgroundColor = [0.651 0.651 0.651];
            app.ReStartStepsButton.FontSize = 18;
            app.ReStartStepsButton.FontWeight = 'bold';
            app.ReStartStepsButton.Position = [658 91 201 36];
            app.ReStartStepsButton.Text = 'ReStart Steps';

            % Create VisualizingTab
            app.VisualizingTab = uitab(app.TabGroup);
            app.VisualizingTab.AutoResizeChildren = 'off';
            app.VisualizingTab.Title = 'Visualizing';

            % Create UIAxes
            app.UIAxes = uiaxes(app.VisualizingTab);
            title(app.UIAxes, 'TMS Evoked Potential')
            xlabel(app.UIAxes, 'Time')
            ylabel(app.UIAxes, 'TEP')
            app.UIAxes.TickDir = 'both';
            app.UIAxes.Position = [340 319 308 186];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.VisualizingTab);
            app.UIAxes2.TickLabelInterpreter = 'none';
            app.UIAxes2.XAxisLocation = 'origin';
            app.UIAxes2.XTick = [];
            app.UIAxes2.YAxisLocation = 'origin';
            app.UIAxes2.YTick = [];
            app.UIAxes2.ZTick = [];
            app.UIAxes2.Position = [340 65 200 148];

            % Create TEPComponentTable
            app.TEPComponentTable = uitable(app.VisualizingTab);
            app.TEPComponentTable.ColumnName = {'Component', 'Latency (ms)', 'Amplitude (µV)'};
            app.TEPComponentTable.ColumnWidth = {80, 100, 128};
            app.TEPComponentTable.RowName = {};
            app.TEPComponentTable.Enable = 'on';
            app.TEPComponentTable.Position = [340 215 308 104];

            % Create ShowComponentsButton
            app.ShowComponentsButton = uibutton(app.VisualizingTab, 'state');
            app.ShowComponentsButton.ValueChangedFcn = createCallbackFcn(app, @ShowComponentsButtonValueChanged, true);
            app.ShowComponentsButton.Text = 'Show Components';
            app.ShowComponentsButton.Position = [5 104 140 23];

            % Create PLOTTEPButton
            app.PLOTTEPButton = uibutton(app.VisualizingTab, 'push');
            app.PLOTTEPButton.ButtonPushedFcn = createCallbackFcn(app, @PLOTTEPButtonPushed, true);
            app.PLOTTEPButton.Enable = 'off';
            app.PLOTTEPButton.Position = [5 130 140 30];
            app.PLOTTEPButton.Text = 'PLOT TEP';

            % Create SelectDatatoVisulaizeTEPsPanel
            app.SelectDatatoVisulaizeTEPsPanel = uipanel(app.VisualizingTab);
            app.SelectDatatoVisulaizeTEPsPanel.AutoResizeChildren = 'off';
            app.SelectDatatoVisulaizeTEPsPanel.BorderType = 'none';
            app.SelectDatatoVisulaizeTEPsPanel.Title = 'Select Data to Visualize TEPs';
            app.SelectDatatoVisulaizeTEPsPanel.Position = [651 230 208 90];

            % Create FolderEditField_2Label
            app.FolderEditField_2Label = uilabel(app.SelectDatatoVisulaizeTEPsPanel);
            app.FolderEditField_2Label.HorizontalAlignment = 'right';
            app.FolderEditField_2Label.Position = [1 41 40 22];
            app.FolderEditField_2Label.Text = 'Folder';

            % Create FolderEditField_2
            app.FolderEditField_2 = uieditfield(app.SelectDatatoVisulaizeTEPsPanel, 'text');
            app.FolderEditField_2.ValueChangedFcn = createCallbackFcn(app, @FolderEditField_2ValueChanged, true);
            app.FolderEditField_2.Editable = 'off';
            app.FolderEditField_2.Position = [49 41 145 22];

            % Create SelectDataButton_2
            app.SelectDataButton_2 = uibutton(app.SelectDatatoVisulaizeTEPsPanel, 'push');
            app.SelectDataButton_2.ButtonPushedFcn = createCallbackFcn(app, @SelectDataButton_2Pushed, true);
            app.SelectDataButton_2.Position = [13 10 183 23];
            app.SelectDataButton_2.Text = 'Select Data';

            % Create UseCurrentlyCleanedDataCheckBox
            app.UseCurrentlyCleanedDataCheckBox = uicheckbox(app.VisualizingTab);
            app.UseCurrentlyCleanedDataCheckBox.ValueChangedFcn = createCallbackFcn(app, @UseCurrentlyCleanedDataCheckBoxValueChanged, true);
            app.UseCurrentlyCleanedDataCheckBox.Text = 'Use Currently Cleaned Data';
            app.UseCurrentlyCleanedDataCheckBox.FontWeight = 'bold';
            app.UseCurrentlyCleanedDataCheckBox.Position = [671 325 180 22];

            % Create FilesListBoxLabel
            app.FilesListBoxLabel = uilabel(app.VisualizingTab);
            app.FilesListBoxLabel.HorizontalAlignment = 'right';
            app.FilesListBoxLabel.Position = [740 220 30 22];
            app.FilesListBoxLabel.Text = 'Files';

            % Create FilesListBox
            app.FilesListBox = uilistbox(app.VisualizingTab);
            app.FilesListBox.Items = {};
            app.FilesListBox.Multiselect = 'on';
            app.FilesListBox.ValueChangedFcn = createCallbackFcn(app, @FilesListBoxValueChanged, true);
            app.FilesListBox.Position = [669 71 183 149];
            app.FilesListBox.Value = {};

            % Create Image2
            app.Image2 = uiimage(app.VisualizingTab);
            app.Image2.Position = [-1 165 350 336];
            app.Image2.ImageSource = fullfile(pathToMLAPP, 'Head.png');

            % Create Slider
            app.Slider = uislider(app.VisualizingTab);
            app.Slider.ValueChangingFcn = createCallbackFcn(app, @SliderValueChanging, true);
            app.Slider.Position = [340 48 200 3];
            app.Slider.Value = 60;

            % Create WindowsizeforTopoplotLabel
            app.WindowsizeforTopoplotLabel = uilabel(app.VisualizingTab);
            app.WindowsizeforTopoplotLabel.HorizontalAlignment = 'center';
            app.WindowsizeforTopoplotLabel.Position = [543 90 103 26];
            app.WindowsizeforTopoplotLabel.Text = {'Window size for'; ' time averaged'; ' Topoplot'};

            % Create WindowsizefortimeaveragedTopoplotEditField
            app.WindowsizefortimeaveragedTopoplotEditField = uieditfield(app.VisualizingTab, 'numeric');
            app.WindowsizefortimeaveragedTopoplotEditField.ValueDisplayFormat = '%.0f';
            app.WindowsizefortimeaveragedTopoplotEditField.Position = [543 65 70 22];

            % Create TOPOPLOTButton
            app.TOPOPLOTButton = uibutton(app.VisualizingTab, 'push');
            app.TOPOPLOTButton.ButtonPushedFcn = createCallbackFcn(app, @TOPOPLOTButtonPushed, true);
            app.TOPOPLOTButton.Enable = 'off';
            app.TOPOPLOTButton.Position = [543 167 104 44];
            app.TOPOPLOTButton.Text = 'TOPOPLOT';

            % Create ExportTEPFigureButton
            app.ExportTEPFigureButton = uibutton(app.VisualizingTab, 'push');
            app.ExportTEPFigureButton.ButtonPushedFcn = createCallbackFcn(app, @ExportTEPFigureButtonPushed, true);
            app.ExportTEPFigureButton.Enable = 'off';
            app.ExportTEPFigureButton.Position = [5 76 140 23];
            app.ExportTEPFigureButton.Text = 'Export TEP Figure';

            % Create PlottingModeButtonGroup
            app.PlottingModeButtonGroup = uibuttongroup(app.VisualizingTab);
            app.PlottingModeButtonGroup.AutoResizeChildren = 'off';
            app.PlottingModeButtonGroup.BorderType = 'none';
            app.PlottingModeButtonGroup.Title = 'Plotting Mode';
            app.PlottingModeButtonGroup.Position = [152 88 150 67];

            % Create NewFigureButton
            app.NewFigureButton = uiradiobutton(app.PlottingModeButtonGroup);
            app.NewFigureButton.Text = 'New Figure';
            app.NewFigureButton.Position = [11 21 83 22];
            app.NewFigureButton.Value = true;

            % Create AddtocurrentFigureButton
            app.AddtocurrentFigureButton = uiradiobutton(app.PlottingModeButtonGroup);
            app.AddtocurrentFigureButton.Text = 'Add to current Figure';
            app.AddtocurrentFigureButton.Position = [11 -1 135 22];

            % Create AF3Button
            app.AF3Button = uibutton(app.VisualizingTab, 'state');
            app.AF3Button.IconAlignment = 'center';
            app.AF3Button.HorizontalAlignment = 'left';
            app.AF3Button.Text = 'AF3';
            app.AF3Button.FontSize = 8;
            app.AF3Button.FontWeight = 'bold';
            app.AF3Button.Position = [108 410 25 23];
            app.AF3Button.Value = true;

            % Create FP1Button
            app.FP1Button = uibutton(app.VisualizingTab, 'state');
            app.FP1Button.IconAlignment = 'center';
            app.FP1Button.HorizontalAlignment = 'left';
            app.FP1Button.Text = 'FP1';
            app.FP1Button.FontSize = 8;
            app.FP1Button.FontWeight = 'bold';
            app.FP1Button.Position = [131 433 25 23];

            % Create FPZButton
            app.FPZButton = uibutton(app.VisualizingTab, 'state');
            app.FPZButton.IconAlignment = 'center';
            app.FPZButton.HorizontalAlignment = 'left';
            app.FPZButton.Text = 'FPZ';
            app.FPZButton.FontSize = 8;
            app.FPZButton.FontWeight = 'bold';
            app.FPZButton.Position = [161 439 25 23];

            % Create FP2Button
            app.FP2Button = uibutton(app.VisualizingTab, 'state');
            app.FP2Button.IconAlignment = 'center';
            app.FP2Button.HorizontalAlignment = 'left';
            app.FP2Button.Text = 'FP2';
            app.FP2Button.FontSize = 8;
            app.FP2Button.FontWeight = 'bold';
            app.FP2Button.Position = [191 433 25 23];

            % Create AF4Button
            app.AF4Button = uibutton(app.VisualizingTab, 'state');
            app.AF4Button.IconAlignment = 'center';
            app.AF4Button.HorizontalAlignment = 'left';
            app.AF4Button.Text = 'AF4';
            app.AF4Button.FontSize = 8;
            app.AF4Button.FontWeight = 'bold';
            app.AF4Button.Position = [215 409 25 23];

            % Create F8Button
            app.F8Button = uibutton(app.VisualizingTab, 'state');
            app.F8Button.IconAlignment = 'center';
            app.F8Button.HorizontalAlignment = 'left';
            app.F8Button.Text = 'F8';
            app.F8Button.FontSize = 8;
            app.F8Button.FontWeight = 'bold';
            app.F8Button.Position = [266 390 25 23];

            % Create F6Button
            app.F6Button = uibutton(app.VisualizingTab, 'state');
            app.F6Button.IconAlignment = 'center';
            app.F6Button.HorizontalAlignment = 'left';
            app.F6Button.Text = 'F6';
            app.F6Button.FontSize = 8;
            app.F6Button.FontWeight = 'bold';
            app.F6Button.Position = [240 385 25 23];

            % Create F4Button
            app.F4Button = uibutton(app.VisualizingTab, 'state');
            app.F4Button.IconAlignment = 'center';
            app.F4Button.HorizontalAlignment = 'left';
            app.F4Button.Text = 'F4';
            app.F4Button.FontSize = 8;
            app.F4Button.FontWeight = 'bold';
            app.F4Button.Position = [214 380 25 23];

            % Create F2Button
            app.F2Button = uibutton(app.VisualizingTab, 'state');
            app.F2Button.IconAlignment = 'center';
            app.F2Button.HorizontalAlignment = 'left';
            app.F2Button.Text = 'F2';
            app.F2Button.FontSize = 8;
            app.F2Button.FontWeight = 'bold';
            app.F2Button.Position = [187 385 25 23];

            % Create F5Button
            app.F5Button = uibutton(app.VisualizingTab, 'state');
            app.F5Button.IconAlignment = 'center';
            app.F5Button.HorizontalAlignment = 'left';
            app.F5Button.Text = 'F5';
            app.F5Button.FontSize = 8;
            app.F5Button.FontWeight = 'bold';
            app.F5Button.Position = [80 385 25 23];

            % Create F3Button
            app.F3Button = uibutton(app.VisualizingTab, 'state');
            app.F3Button.IconAlignment = 'center';
            app.F3Button.HorizontalAlignment = 'left';
            app.F3Button.Text = 'F3';
            app.F3Button.FontSize = 8;
            app.F3Button.FontWeight = 'bold';
            app.F3Button.Position = [107 380 25 23];
            app.F3Button.Value = true;

            % Create FZButton
            app.FZButton = uibutton(app.VisualizingTab, 'state');
            app.FZButton.IconAlignment = 'center';
            app.FZButton.HorizontalAlignment = 'left';
            app.FZButton.Text = 'FZ';
            app.FZButton.FontSize = 8;
            app.FZButton.FontWeight = 'bold';
            app.FZButton.Position = [161 386 25 23];

            % Create FC2Button
            app.FC2Button = uibutton(app.VisualizingTab, 'state');
            app.FC2Button.IconAlignment = 'center';
            app.FC2Button.HorizontalAlignment = 'left';
            app.FC2Button.Text = 'FC2';
            app.FC2Button.FontSize = 8;
            app.FC2Button.FontWeight = 'bold';
            app.FC2Button.Position = [192 348 25 23];

            % Create FC4Button
            app.FC4Button = uibutton(app.VisualizingTab, 'state');
            app.FC4Button.IconAlignment = 'center';
            app.FC4Button.HorizontalAlignment = 'left';
            app.FC4Button.Text = 'FC4';
            app.FC4Button.FontSize = 8;
            app.FC4Button.FontWeight = 'bold';
            app.FC4Button.Position = [222 348 25 23];

            % Create FC6Button
            app.FC6Button = uibutton(app.VisualizingTab, 'state');
            app.FC6Button.IconAlignment = 'center';
            app.FC6Button.HorizontalAlignment = 'left';
            app.FC6Button.Text = 'FC6';
            app.FC6Button.FontSize = 8;
            app.FC6Button.FontWeight = 'bold';
            app.FC6Button.Position = [252 350 25 23];

            % Create F1Button
            app.F1Button = uibutton(app.VisualizingTab, 'state');
            app.F1Button.IconAlignment = 'center';
            app.F1Button.HorizontalAlignment = 'left';
            app.F1Button.Text = 'F1';
            app.F1Button.FontSize = 8;
            app.F1Button.FontWeight = 'bold';
            app.F1Button.Position = [134 385 25 23];
            app.F1Button.Value = true;

            % Create C4Button
            app.C4Button = uibutton(app.VisualizingTab, 'state');
            app.C4Button.IconAlignment = 'center';
            app.C4Button.HorizontalAlignment = 'left';
            app.C4Button.Text = 'C4';
            app.C4Button.FontSize = 8;
            app.C4Button.FontWeight = 'bold';
            app.C4Button.Position = [229 316 25 23];

            % Create C6Button
            app.C6Button = uibutton(app.VisualizingTab, 'state');
            app.C6Button.IconAlignment = 'center';
            app.C6Button.HorizontalAlignment = 'left';
            app.C6Button.Text = 'C6';
            app.C6Button.FontSize = 8;
            app.C6Button.FontWeight = 'bold';
            app.C6Button.Position = [262 316 25 23];

            % Create FT8Button
            app.FT8Button = uibutton(app.VisualizingTab, 'state');
            app.FT8Button.IconAlignment = 'center';
            app.FT8Button.HorizontalAlignment = 'left';
            app.FT8Button.Text = 'FT8';
            app.FT8Button.FontSize = 8;
            app.FT8Button.FontWeight = 'bold';
            app.FT8Button.Position = [285 354 25 23];

            % Create F7Button
            app.F7Button = uibutton(app.VisualizingTab, 'state');
            app.F7Button.IconAlignment = 'center';
            app.F7Button.HorizontalAlignment = 'left';
            app.F7Button.Text = 'F7';
            app.F7Button.FontSize = 8;
            app.F7Button.FontWeight = 'bold';
            app.F7Button.Position = [55 391 25 23];

            % Create FC1Button
            app.FC1Button = uibutton(app.VisualizingTab, 'state');
            app.FC1Button.IconAlignment = 'center';
            app.FC1Button.HorizontalAlignment = 'left';
            app.FC1Button.Text = 'FC1';
            app.FC1Button.FontSize = 8;
            app.FC1Button.FontWeight = 'bold';
            app.FC1Button.Position = [130 348 25 23];
            app.FC1Button.Value = true;

            % Create FCZButton
            app.FCZButton = uibutton(app.VisualizingTab, 'state');
            app.FCZButton.IconAlignment = 'center';
            app.FCZButton.HorizontalAlignment = 'left';
            app.FCZButton.Text = 'FCZ';
            app.FCZButton.FontSize = 8;
            app.FCZButton.FontWeight = 'bold';
            app.FCZButton.Position = [161 349 25 23];

            % Create FC3Button
            app.FC3Button = uibutton(app.VisualizingTab, 'state');
            app.FC3Button.IconAlignment = 'center';
            app.FC3Button.HorizontalAlignment = 'left';
            app.FC3Button.Text = 'FC3';
            app.FC3Button.FontSize = 8;
            app.FC3Button.FontWeight = 'bold';
            app.FC3Button.Position = [99 348 25 23];
            app.FC3Button.Value = true;

            % Create C1Button
            app.C1Button = uibutton(app.VisualizingTab, 'state');
            app.C1Button.IconAlignment = 'center';
            app.C1Button.HorizontalAlignment = 'left';
            app.C1Button.Text = 'C1';
            app.C1Button.FontSize = 8;
            app.C1Button.FontWeight = 'bold';
            app.C1Button.Position = [128 316 25 23];

            % Create CZButton
            app.CZButton = uibutton(app.VisualizingTab, 'state');
            app.CZButton.IconAlignment = 'center';
            app.CZButton.HorizontalAlignment = 'left';
            app.CZButton.Text = 'CZ';
            app.CZButton.FontSize = 8;
            app.CZButton.FontWeight = 'bold';
            app.CZButton.Position = [161 316 25 23];

            % Create C2Button
            app.C2Button = uibutton(app.VisualizingTab, 'state');
            app.C2Button.IconAlignment = 'center';
            app.C2Button.HorizontalAlignment = 'left';
            app.C2Button.Text = 'C2';
            app.C2Button.FontSize = 8;
            app.C2Button.FontWeight = 'bold';
            app.C2Button.Position = [195 316 25 23];

            % Create CP3Button
            app.CP3Button = uibutton(app.VisualizingTab, 'state');
            app.CP3Button.IconAlignment = 'center';
            app.CP3Button.HorizontalAlignment = 'left';
            app.CP3Button.Text = 'CP3';
            app.CP3Button.FontSize = 8;
            app.CP3Button.FontWeight = 'bold';
            app.CP3Button.Position = [97 283 25 23];

            % Create CP1Button
            app.CP1Button = uibutton(app.VisualizingTab, 'state');
            app.CP1Button.IconAlignment = 'center';
            app.CP1Button.HorizontalAlignment = 'left';
            app.CP1Button.Text = 'CP1';
            app.CP1Button.FontSize = 8;
            app.CP1Button.FontWeight = 'bold';
            app.CP1Button.Position = [128 284 25 23];

            % Create CP2Button
            app.CP2Button = uibutton(app.VisualizingTab, 'state');
            app.CP2Button.IconAlignment = 'center';
            app.CP2Button.HorizontalAlignment = 'left';
            app.CP2Button.Text = 'CP2';
            app.CP2Button.FontSize = 8;
            app.CP2Button.FontWeight = 'bold';
            app.CP2Button.Position = [192 284 25 23];

            % Create T8Button
            app.T8Button = uibutton(app.VisualizingTab, 'state');
            app.T8Button.IconAlignment = 'center';
            app.T8Button.HorizontalAlignment = 'left';
            app.T8Button.Text = 'T8';
            app.T8Button.FontSize = 8;
            app.T8Button.FontWeight = 'bold';
            app.T8Button.Position = [293 316 25 23];

            % Create FT7Button
            app.FT7Button = uibutton(app.VisualizingTab, 'state');
            app.FT7Button.IconAlignment = 'center';
            app.FT7Button.HorizontalAlignment = 'left';
            app.FT7Button.Text = 'FT7';
            app.FT7Button.FontSize = 8;
            app.FT7Button.FontWeight = 'bold';
            app.FT7Button.Position = [36 354 25 23];

            % Create FC5Button
            app.FC5Button = uibutton(app.VisualizingTab, 'state');
            app.FC5Button.IconAlignment = 'center';
            app.FC5Button.HorizontalAlignment = 'left';
            app.FC5Button.Text = 'FC5';
            app.FC5Button.FontSize = 8;
            app.FC5Button.FontWeight = 'bold';
            app.FC5Button.Position = [68 350 25 23];

            % Create C5Button
            app.C5Button = uibutton(app.VisualizingTab, 'state');
            app.C5Button.IconAlignment = 'center';
            app.C5Button.HorizontalAlignment = 'left';
            app.C5Button.Text = 'C5';
            app.C5Button.FontSize = 8;
            app.C5Button.FontWeight = 'bold';
            app.C5Button.Position = [59 316 25 23];

            % Create C3Button
            app.C3Button = uibutton(app.VisualizingTab, 'state');
            app.C3Button.IconAlignment = 'center';
            app.C3Button.HorizontalAlignment = 'left';
            app.C3Button.Text = 'C3';
            app.C3Button.FontSize = 8;
            app.C3Button.FontWeight = 'bold';
            app.C3Button.Position = [93 316 25 23];

            % Create T7Button
            app.T7Button = uibutton(app.VisualizingTab, 'state');
            app.T7Button.IconAlignment = 'center';
            app.T7Button.HorizontalAlignment = 'left';
            app.T7Button.Text = 'T7';
            app.T7Button.FontSize = 8;
            app.T7Button.FontWeight = 'bold';
            app.T7Button.Position = [26 316 25 23];

            % Create TP7Button
            app.TP7Button = uibutton(app.VisualizingTab, 'state');
            app.TP7Button.IconAlignment = 'center';
            app.TP7Button.HorizontalAlignment = 'left';
            app.TP7Button.Text = 'TP7';
            app.TP7Button.FontSize = 8;
            app.TP7Button.FontWeight = 'bold';
            app.TP7Button.Position = [34 277 25 23];

            % Create CP5Button
            app.CP5Button = uibutton(app.VisualizingTab, 'state');
            app.CP5Button.IconAlignment = 'center';
            app.CP5Button.HorizontalAlignment = 'left';
            app.CP5Button.Text = 'CP5';
            app.CP5Button.FontSize = 8;
            app.CP5Button.FontWeight = 'bold';
            app.CP5Button.Position = [62 279 25 23];

            % Create CPZButton
            app.CPZButton = uibutton(app.VisualizingTab, 'state');
            app.CPZButton.IconAlignment = 'center';
            app.CPZButton.HorizontalAlignment = 'left';
            app.CPZButton.Text = 'CPZ';
            app.CPZButton.FontSize = 8;
            app.CPZButton.FontWeight = 'bold';
            app.CPZButton.Position = [161 286 25 23];

            % Create CP4Button
            app.CP4Button = uibutton(app.VisualizingTab, 'state');
            app.CP4Button.IconAlignment = 'center';
            app.CP4Button.HorizontalAlignment = 'left';
            app.CP4Button.Text = 'CP4';
            app.CP4Button.FontSize = 8;
            app.CP4Button.FontWeight = 'bold';
            app.CP4Button.Position = [224 283 25 23];

            % Create CP6Button
            app.CP6Button = uibutton(app.VisualizingTab, 'state');
            app.CP6Button.IconAlignment = 'center';
            app.CP6Button.HorizontalAlignment = 'left';
            app.CP6Button.Text = 'CP6';
            app.CP6Button.FontSize = 8;
            app.CP6Button.FontWeight = 'bold';
            app.CP6Button.Position = [258 279 25 23];

            % Create TP8Button
            app.TP8Button = uibutton(app.VisualizingTab, 'state');
            app.TP8Button.IconAlignment = 'center';
            app.TP8Button.HorizontalAlignment = 'left';
            app.TP8Button.Text = 'TP8';
            app.TP8Button.FontSize = 8;
            app.TP8Button.FontWeight = 'bold';
            app.TP8Button.Position = [288 277 25 23];

            % Create P8Button
            app.P8Button = uibutton(app.VisualizingTab, 'state');
            app.P8Button.IconAlignment = 'center';
            app.P8Button.HorizontalAlignment = 'left';
            app.P8Button.Text = 'P8';
            app.P8Button.FontSize = 8;
            app.P8Button.FontWeight = 'bold';
            app.P8Button.Position = [275 237 25 23];

            % Create P3Button
            app.P3Button = uibutton(app.VisualizingTab, 'state');
            app.P3Button.IconAlignment = 'center';
            app.P3Button.HorizontalAlignment = 'left';
            app.P3Button.Text = 'P3';
            app.P3Button.FontSize = 8;
            app.P3Button.FontWeight = 'bold';
            app.P3Button.Position = [105 251 25 23];

            % Create P1Button
            app.P1Button = uibutton(app.VisualizingTab, 'state');
            app.P1Button.IconAlignment = 'center';
            app.P1Button.HorizontalAlignment = 'left';
            app.P1Button.Text = 'P1';
            app.P1Button.FontSize = 8;
            app.P1Button.FontWeight = 'bold';
            app.P1Button.Position = [132 251 25 23];

            % Create P2Button
            app.P2Button = uibutton(app.VisualizingTab, 'state');
            app.P2Button.IconAlignment = 'center';
            app.P2Button.HorizontalAlignment = 'left';
            app.P2Button.Text = 'P2';
            app.P2Button.FontSize = 8;
            app.P2Button.FontWeight = 'bold';
            app.P2Button.Position = [188 251 25 23];

            % Create P7Button
            app.P7Button = uibutton(app.VisualizingTab, 'state');
            app.P7Button.IconAlignment = 'center';
            app.P7Button.HorizontalAlignment = 'left';
            app.P7Button.Text = 'P7';
            app.P7Button.FontSize = 8;
            app.P7Button.FontWeight = 'bold';
            app.P7Button.Position = [46 237 25 23];

            % Create P5Button
            app.P5Button = uibutton(app.VisualizingTab, 'state');
            app.P5Button.IconAlignment = 'center';
            app.P5Button.HorizontalAlignment = 'left';
            app.P5Button.Text = 'P5';
            app.P5Button.FontSize = 8;
            app.P5Button.FontWeight = 'bold';
            app.P5Button.Position = [76 246 25 23];

            % Create PZButton
            app.PZButton = uibutton(app.VisualizingTab, 'state');
            app.PZButton.IconAlignment = 'center';
            app.PZButton.HorizontalAlignment = 'left';
            app.PZButton.Text = 'PZ';
            app.PZButton.FontSize = 8;
            app.PZButton.FontWeight = 'bold';
            app.PZButton.Position = [161 250 25 23];

            % Create P4Button
            app.P4Button = uibutton(app.VisualizingTab, 'state');
            app.P4Button.IconAlignment = 'center';
            app.P4Button.HorizontalAlignment = 'left';
            app.P4Button.Text = 'P4';
            app.P4Button.FontSize = 8;
            app.P4Button.FontWeight = 'bold';
            app.P4Button.Position = [216 251 25 23];

            % Create P6Button
            app.P6Button = uibutton(app.VisualizingTab, 'state');
            app.P6Button.IconAlignment = 'center';
            app.P6Button.HorizontalAlignment = 'left';
            app.P6Button.Text = 'P6';
            app.P6Button.FontSize = 8;
            app.P6Button.FontWeight = 'bold';
            app.P6Button.Position = [245 246 25 23];

            % Create O1Button
            app.O1Button = uibutton(app.VisualizingTab, 'state');
            app.O1Button.IconAlignment = 'center';
            app.O1Button.HorizontalAlignment = 'left';
            app.O1Button.Text = 'O1';
            app.O1Button.FontSize = 8;
            app.O1Button.FontWeight = 'bold';
            app.O1Button.Position = [128 179 25 23];

            % Create PO3Button
            app.PO3Button = uibutton(app.VisualizingTab, 'state');
            app.PO3Button.IconAlignment = 'center';
            app.PO3Button.HorizontalAlignment = 'left';
            app.PO3Button.Text = 'PO3';
            app.PO3Button.FontSize = 8;
            app.PO3Button.FontWeight = 'bold';
            app.PO3Button.Position = [106 216 25 23];

            % Create POZButton
            app.POZButton = uibutton(app.VisualizingTab, 'state');
            app.POZButton.IconAlignment = 'center';
            app.POZButton.HorizontalAlignment = 'left';
            app.POZButton.Text = 'POZ';
            app.POZButton.FontSize = 8;
            app.POZButton.FontWeight = 'bold';
            app.POZButton.Position = [161 211 25 23];

            % Create PO4Button
            app.PO4Button = uibutton(app.VisualizingTab, 'state');
            app.PO4Button.IconAlignment = 'center';
            app.PO4Button.HorizontalAlignment = 'left';
            app.PO4Button.Text = 'PO4';
            app.PO4Button.FontSize = 8;
            app.PO4Button.FontWeight = 'bold';
            app.PO4Button.Position = [217 217 25 23];

            % Create PO7Button
            app.PO7Button = uibutton(app.VisualizingTab, 'state');
            app.PO7Button.IconAlignment = 'center';
            app.PO7Button.HorizontalAlignment = 'left';
            app.PO7Button.Text = 'PO7';
            app.PO7Button.FontSize = 8;
            app.PO7Button.FontWeight = 'bold';
            app.PO7Button.Position = [52 204 25 23];

            % Create PO5Button
            app.PO5Button = uibutton(app.VisualizingTab, 'state');
            app.PO5Button.IconAlignment = 'center';
            app.PO5Button.HorizontalAlignment = 'left';
            app.PO5Button.Text = 'PO5';
            app.PO5Button.FontSize = 8;
            app.PO5Button.FontWeight = 'bold';
            app.PO5Button.Position = [78 215 25 23];

            % Create PO2Button
            app.PO2Button = uibutton(app.VisualizingTab, 'state');
            app.PO2Button.IconAlignment = 'center';
            app.PO2Button.HorizontalAlignment = 'left';
            app.PO2Button.Text = 'PO2';
            app.PO2Button.FontSize = 8;
            app.PO2Button.FontWeight = 'bold';
            app.PO2Button.Position = [189 212 25 23];

            % Create PO8Button
            app.PO8Button = uibutton(app.VisualizingTab, 'state');
            app.PO8Button.IconAlignment = 'center';
            app.PO8Button.HorizontalAlignment = 'left';
            app.PO8Button.Text = 'PO8';
            app.PO8Button.FontSize = 8;
            app.PO8Button.FontWeight = 'bold';
            app.PO8Button.Position = [270 202 25 23];

            % Create CB1Button
            app.CB1Button = uibutton(app.VisualizingTab, 'state');
            app.CB1Button.IconAlignment = 'center';
            app.CB1Button.HorizontalAlignment = 'left';
            app.CB1Button.Text = 'CB1';
            app.CB1Button.FontSize = 8;
            app.CB1Button.FontWeight = 'bold';
            app.CB1Button.Position = [98 189 25 23];

            % Create OZButton
            app.OZButton = uibutton(app.VisualizingTab, 'state');
            app.OZButton.IconAlignment = 'center';
            app.OZButton.HorizontalAlignment = 'left';
            app.OZButton.Text = 'OZ';
            app.OZButton.FontSize = 8;
            app.OZButton.FontWeight = 'bold';
            app.OZButton.Position = [160 177 25 23];

            % Create O2Button
            app.O2Button = uibutton(app.VisualizingTab, 'state');
            app.O2Button.IconAlignment = 'center';
            app.O2Button.HorizontalAlignment = 'left';
            app.O2Button.Text = 'O2';
            app.O2Button.FontSize = 8;
            app.O2Button.FontWeight = 'bold';
            app.O2Button.Position = [192 179 25 23];

            % Create CB2Button
            app.CB2Button = uibutton(app.VisualizingTab, 'state');
            app.CB2Button.IconAlignment = 'center';
            app.CB2Button.HorizontalAlignment = 'left';
            app.CB2Button.Text = 'CB2';
            app.CB2Button.FontSize = 8;
            app.CB2Button.FontWeight = 'bold';
            app.CB2Button.Position = [224 189 25 23];

            % Create TP10Button
            app.TP10Button = uibutton(app.VisualizingTab, 'state');
            app.TP10Button.IconAlignment = 'center';
            app.TP10Button.HorizontalAlignment = 'left';
            app.TP10Button.Text = 'TP10';
            app.TP10Button.FontSize = 7;
            app.TP10Button.FontWeight = 'bold';
            app.TP10Button.Position = [302 253 25 23];

            % Create TP9Button
            app.TP9Button = uibutton(app.VisualizingTab, 'state');
            app.TP9Button.IconAlignment = 'center';
            app.TP9Button.HorizontalAlignment = 'left';
            app.TP9Button.Text = 'TP9';
            app.TP9Button.FontSize = 7;
            app.TP9Button.FontWeight = 'bold';
            app.TP9Button.Position = [20 253 25 23];

            % Create AFZButton
            app.AFZButton = uibutton(app.VisualizingTab, 'state');
            app.AFZButton.IconAlignment = 'center';
            app.AFZButton.HorizontalAlignment = 'left';
            app.AFZButton.Text = 'AFZ';
            app.AFZButton.FontSize = 8;
            app.AFZButton.FontWeight = 'bold';
            app.AFZButton.Position = [161 412 25 23];

            % Create AF7Button
            app.AF7Button = uibutton(app.VisualizingTab, 'state');
            app.AF7Button.IconAlignment = 'center';
            app.AF7Button.HorizontalAlignment = 'left';
            app.AF7Button.Text = 'AF7';
            app.AF7Button.FontSize = 8;
            app.AF7Button.FontWeight = 'bold';
            app.AF7Button.Position = [79 419 25 23];

            % Create AF8Button
            app.AF8Button = uibutton(app.VisualizingTab, 'state');
            app.AF8Button.IconAlignment = 'center';
            app.AF8Button.HorizontalAlignment = 'left';
            app.AF8Button.Text = 'AF8';
            app.AF8Button.FontSize = 8;
            app.AF8Button.FontWeight = 'bold';
            app.AF8Button.Position = [245 419 25 23];

            % Create SelectAllCheckBox
            app.SelectAllCheckBox = uicheckbox(app.VisualizingTab);
            app.SelectAllCheckBox.ValueChangedFcn = createCallbackFcn(app, @SelectAllCheckBoxValueChanged, true);
            app.SelectAllCheckBox.Text = 'Select All';
            app.SelectAllCheckBox.Position = [670 46 71 22];

            % Create DontfindcommonelectrodesCheckBox
            app.DontfindcommonelectrodesCheckBox = uicheckbox(app.VisualizingTab);
            app.DontfindcommonelectrodesCheckBox.ValueChangedFcn = createCallbackFcn(app, @DontfindcommonelectrodesCheckBoxValueChanged, true);
            app.DontfindcommonelectrodesCheckBox.Text = 'Don''t find common electrodes';
            app.DontfindcommonelectrodesCheckBox.Position = [670 28 180 22];
            app.DontfindcommonelectrodesCheckBox.Value = true;

            % Create PO1Button
            app.PO1Button = uibutton(app.VisualizingTab, 'state');
            app.PO1Button.IconAlignment = 'center';
            app.PO1Button.HorizontalAlignment = 'left';
            app.PO1Button.Text = 'PO1';
            app.PO1Button.FontSize = 8;
            app.PO1Button.FontWeight = 'bold';
            app.PO1Button.Position = [133 212 25 23];

            % Create PO6Button
            app.PO6Button = uibutton(app.VisualizingTab, 'state');
            app.PO6Button.IconAlignment = 'center';
            app.PO6Button.HorizontalAlignment = 'left';
            app.PO6Button.Text = 'PO6';
            app.PO6Button.FontSize = 8;
            app.PO6Button.FontWeight = 'bold';
            app.PO6Button.Position = [243 215 25 23];

            % Create ReLoadAvailableElectrodesButton
            app.ReLoadAvailableElectrodesButton = uibutton(app.VisualizingTab, 'push');
            app.ReLoadAvailableElectrodesButton.ButtonPushedFcn = createCallbackFcn(app, @ReLoadAvailableElectrodesButtonPushed, true);
            app.ReLoadAvailableElectrodesButton.BackgroundColor = [0 0.4471 0.7412];
            app.ReLoadAvailableElectrodesButton.FontSize = 10;
            app.ReLoadAvailableElectrodesButton.FontWeight = 'bold';
            app.ReLoadAvailableElectrodesButton.FontColor = [0.9294 0.6941 0.1255];
            app.ReLoadAvailableElectrodesButton.Enable = 'off';
            app.ReLoadAvailableElectrodesButton.Position = [669 7 183 23];
            app.ReLoadAvailableElectrodesButton.Text = 'Re/Load Available Electrodes';

            % Create TEPWindowSliderLabel
            app.TEPWindowSliderLabel = uilabel(app.VisualizingTab);
            app.TEPWindowSliderLabel.HorizontalAlignment = 'right';
            app.TEPWindowSliderLabel.Position = [669 420 85 16];
            app.TEPWindowSliderLabel.Text = 'TEP Window';

            % Create TEPWindowSlider
            app.TEPWindowSlider = uislider(app.VisualizingTab, 'range');
            app.TEPWindowSlider.Limits = [-100 300];
            app.TEPWindowSlider.ValueChangingFcn = createCallbackFcn(app, @TEPWindowSliderValueChanging, true);
            app.TEPWindowSlider.Position = [669 388 183 3];
            app.TEPWindowSlider.Value = [-50 300];

            % Create TopoplottimeSpinnerLabel
            app.TopoplottimeSpinnerLabel = uilabel(app.VisualizingTab);
            app.TopoplottimeSpinnerLabel.HorizontalAlignment = 'right';
            app.TopoplottimeSpinnerLabel.Position = [543 145 103 20];
            app.TopoplottimeSpinnerLabel.Text = 'Topoplot time';

            % Create TopoplottimeSpinner
            app.TopoplottimeSpinner = uispinner(app.VisualizingTab);
            app.TopoplottimeSpinner.RoundFractionalValues = 'on';
            app.TopoplottimeSpinner.ValueDisplayFormat = '%.0f';
            app.TopoplottimeSpinner.ValueChangedFcn = createCallbackFcn(app, @TopoplottimeSpinnerValueChanged, true);
            app.TopoplottimeSpinner.Position = [543 120 103 22];
            app.TopoplottimeSpinner.Value = 60;

            % Create EEGDatasetDropDownLabel
            app.EEGDatasetDropDownLabel = uilabel(app.VisualizingTab);
            app.EEGDatasetDropDownLabel.HorizontalAlignment = 'right';
            app.EEGDatasetDropDownLabel.Position = [152 58 75 22];
            app.EEGDatasetDropDownLabel.Text = 'EEG Dataset';

            % Create EEGDatasetDropDown
            app.EEGDatasetDropDown = uidropdown(app.VisualizingTab);
            app.EEGDatasetDropDown.Items = {'Select a file'};
            app.EEGDatasetDropDown.ValueChangedFcn = createCallbackFcn(app, @EEGDatasetDropDownValueChanged, true);
            app.EEGDatasetDropDown.Enable = 'off';
            app.EEGDatasetDropDown.Position = [230 58 100 22];
            app.EEGDatasetDropDown.Value = 'Select a file';

            % Create PlotEEGdataButton
            app.PlotEEGdataButton = uibutton(app.VisualizingTab, 'push');
            app.PlotEEGdataButton.ButtonPushedFcn = createCallbackFcn(app, @PlotEEGdataButtonPushed, true);
            app.PlotEEGdataButton.Enable = 'off';
            app.PlotEEGdataButton.Position = [5 50 108 23];
            app.PlotEEGdataButton.Text = 'Plot EEG data';

            % Create ExportTEPDataButton — right column, below files listbox
            app.ExportTEPDataButton = uibutton(app.VisualizingTab, 'push');
            app.ExportTEPDataButton.ButtonPushedFcn = createCallbackFcn(app, @ExportTEPDataButtonPushed, true);
            app.ExportTEPDataButton.Enable = 'off';
            app.ExportTEPDataButton.Position = [669 468 183 23];
            app.ExportTEPDataButton.Text = 'Export TEP Data';

            % Create TEPvarNameEditFieldLabel — right column, below Export TEP Data
            app.TEPvarNameEditFieldLabel = uilabel(app.VisualizingTab);
            app.TEPvarNameEditFieldLabel.HorizontalAlignment = 'right';
            app.TEPvarNameEditFieldLabel.Enable = 'off';
            app.TEPvarNameEditFieldLabel.Position = [669 438 80 22];
            app.TEPvarNameEditFieldLabel.Text = 'TEP var Name';

            % Create TEPvarNameEditField
            app.TEPvarNameEditField = uieditfield(app.VisualizingTab, 'text');
            app.TEPvarNameEditField.ValueChangedFcn = createCallbackFcn(app, @TEPvarNameEditFieldValueChanged, true);
            app.TEPvarNameEditField.Enable = 'off';
            app.TEPvarNameEditField.Position = [754 438 98 22];

            % Create ReportsTab
            app.ReportsTab = uitab(app.TabGroup);
            app.ReportsTab.AutoResizeChildren = 'off';
            app.ReportsTab.Title = 'Reports';

            % Reports tab — left column: session list
            app.ReportsListBoxLabel = uilabel(app.ReportsTab);
            app.ReportsListBoxLabel.FontSize = 16;
            app.ReportsListBoxLabel.FontWeight = 'bold';
            app.ReportsListBoxLabel.Position = [5 472 205 22];
            app.ReportsListBoxLabel.Text = 'Session Reports';

            app.ReportsListBox = uilistbox(app.ReportsTab);
            app.ReportsListBox.Items = {};
            app.ReportsListBox.Position = [5 73 205 393];
            app.ReportsListBox.ValueChangedFcn = createCallbackFcn(app, @ReportsListBoxValueChanged, true);

            app.LoadReportsButton = uibutton(app.ReportsTab, 'push');
            app.LoadReportsButton.ButtonPushedFcn = createCallbackFcn(app, @LoadReportsButtonPushed, true);
            app.LoadReportsButton.Position = [5 45 100 25];
            app.LoadReportsButton.Text = 'Load from Folder';
            app.LoadReportsButton.Tooltip = 'Load pipeline reports from a folder on disk';

            app.RefreshReportsButton = uibutton(app.ReportsTab, 'push');
            app.RefreshReportsButton.ButtonPushedFcn = createCallbackFcn(app, @RefreshReportsButtonPushed, true);
            app.RefreshReportsButton.Position = [110 45 100 25];
            app.RefreshReportsButton.Text = 'Refresh';
            app.RefreshReportsButton.Tooltip = 'Reload reports from the current folder';

            app.ReportsFolderLabel = uilabel(app.ReportsTab);
            app.ReportsFolderLabel.FontSize = 9;
            app.ReportsFolderLabel.FontColor = [0.5 0.5 0.5];
            app.ReportsFolderLabel.Position = [5 25 205 18];
            app.ReportsFolderLabel.Text = '';

            app.ReportsStatusLabel = uilabel(app.ReportsTab);
            app.ReportsStatusLabel.FontSize = 9;
            app.ReportsStatusLabel.FontColor = [0.5 0.5 0.5];
            app.ReportsStatusLabel.Position = [5 5 205 18];
            app.ReportsStatusLabel.Text = 'No reports loaded.';

            % Reports tab — right column: report text + actions
            app.ExportReportsCSVButton = uibutton(app.ReportsTab, 'push');
            app.ExportReportsCSVButton.ButtonPushedFcn = createCallbackFcn(app, @ExportReportsCSVButtonPushed, true);
            app.ExportReportsCSVButton.Position = [580 470 130 24];
            app.ExportReportsCSVButton.Text = 'Export CSV';
            app.ExportReportsCSVButton.Tooltip = 'Export a CSV summary of all loaded reports';
            app.ExportReportsCSVButton.Enable = 'off';

            app.CopyMethodsButton = uibutton(app.ReportsTab, 'push');
            app.CopyMethodsButton.ButtonPushedFcn = createCallbackFcn(app, @CopyMethodsButtonPushed, true);
            app.CopyMethodsButton.Position = [715 470 147 24];
            app.CopyMethodsButton.Text = 'Copy Methods Text';
            app.CopyMethodsButton.Tooltip = 'Copy a methods paragraph for the selected report to the clipboard';
            app.CopyMethodsButton.Enable = 'off';

            app.ReportsTextArea = uitextarea(app.ReportsTab);
            app.ReportsTextArea.Editable = 'off';
            app.ReportsTextArea.FontName = 'Courier New';
            app.ReportsTextArea.FontSize = 10;
            app.ReportsTextArea.Position = [220 10 637 457];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % Public methods callable from external functions (e.g. runPipeline)
    methods (Access = public)

        function updateReportsTab(app)
        % UPDATEREPORTSTAB  Public entry point — refreshes the Reports tab.
        %   Delegates to the private implementation. Exposed as public so
        %   runPipeline.m can call it after each processing run.
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