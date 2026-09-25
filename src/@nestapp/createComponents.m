
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function createComponents(app)
% CREATECOMPONENTS  Create all UI components for nestapp.
%   Sets every property on app.UIFigure and all child controls.
%   Called from the nestapp constructor via App Designer's createComponents hook.
%   This file lives in src/@nestapp/ - a proper class method with access to
%   protected methods like createCallbackFcn.
%
% WARNING: Do not open nestapp_designer.mlapp and save - App Designer will
% regenerate nestapp.m and may overwrite the createComponents call path.
% All layout edits belong in this file.

            % Get the file path for locating images
            imageDir = fullfile(nestappRoot(), 'src');

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            % Height = 529px tab area + 20px status bar = 549px.
            % uimenu renders outside the coordinate space (MATLAB shifts the window
            % upward when the menu is created; coordinate height stays unchanged).
            app.UIFigure.Position = [100 100 867 549];
            app.UIFigure.Name = 'nestapp - TMS-EEG Processing';
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.SizeChangedFcn    = createCallbackFcn(app, @UIFigureSizeChanged, true);
            app.UIFigure.CloseRequestFcn   = createCallbackFcn(app, @UIFigureCloseRequest, true);
            % Drives the dwell-delayed step legend (see stepsTreeLegend).
            app.UIFigure.WindowButtonMotionFcn = createCallbackFcn(app, @UIFigureMouseMoved, true);

            % Create menu bar
            mFile = uimenu(app.UIFigure, 'Text', 'File');
            uimenu(mFile, 'Text', 'Open Data...', 'Accelerator', 'O', ...
                'MenuSelectedFcn', createCallbackFcn(app, @SelectDataButtonPushed, true));
            app.MenuRecentFiles = uimenu(mFile, 'Text', 'Recent Files');
            uimenu(mFile, 'Text', 'Load Pipeline...', 'Accelerator', 'L', 'Separator', 'on', ...
                'MenuSelectedFcn', createCallbackFcn(app, @LoadPipelineButtonPushed, true));
            uimenu(mFile, 'Text', 'Save Pipeline', 'Accelerator', 'S', ...
                'MenuSelectedFcn', createCallbackFcn(app, @SavePipelineButtonPushed, true));
            uimenu(mFile, 'Text', 'Copy Pipeline Description', ...
                'MenuSelectedFcn', createCallbackFcn(app, @copyPipelineDescriptionMenu, true));
            app.MenuRecentPipelines = uimenu(mFile, 'Text', 'Recent Pipelines');
            uimenu(mFile, 'Text', 'Load Analysis...', 'Separator', 'on', ...
                'MenuSelectedFcn', createCallbackFcn(app, @LoadAnalysisMenuSelected, true));
            uimenu(mFile, 'Text', 'Load Template...', 'Separator', 'on', ...
                'MenuSelectedFcn', createCallbackFcn(app, @LoadTemplateMenuSelected, true));
            uimenu(mFile, 'Text', 'Exit', 'Separator', 'on', ...
                'MenuSelectedFcn', createCallbackFcn(app, @UIFigureCloseRequest, true));

            mSettings = uimenu(app.UIFigure, 'Text', 'Settings');
            uimenu(mSettings, 'Text', 'Preferences...', ...
                'MenuSelectedFcn', createCallbackFcn(app, @openPreferencesMenu, true));

            mTools = uimenu(app.UIFigure, 'Text', 'Tools');
            uimenu(mTools, 'Text', 'Browse EEG...', ...
                'MenuSelectedFcn', createCallbackFcn(app, @browseRawEegMenu, true));

            mHelp = uimenu(app.UIFigure, 'Text', 'Help');
            uimenu(mHelp, 'Text', 'About nestapp', ...
                'MenuSelectedFcn', createCallbackFcn(app, @showAboutMenu, true));
            uimenu(mHelp, 'Text', 'Copy Diagnostics to Clipboard', 'Separator', 'on', ...
                'MenuSelectedFcn', createCallbackFcn(app, @copyDiagnosticsMenu, true));
            uimenu(mHelp, 'Text', 'Collect Support Bundle...', ...
                'MenuSelectedFcn', createCallbackFcn(app, @collectSupportBundleMenu, true));
            uimenu(mHelp, 'Text', 'Check My Install', ...
                'MenuSelectedFcn', createCallbackFcn(app, @selfTestMenu, true));
            uimenu(mHelp, 'Text', 'Install AARATEP Helpers...', 'Separator', 'on', ...
                'MenuSelectedFcn', createCallbackFcn(app, @installAaratepMenu, true));

            % Create status bar - pinned to bottom of UIFigure, visible on both tabs
            app.StatusBar = uilabel(app.UIFigure);
            app.StatusBar.Position = [0 0 867 20];
            app.StatusBar.BackgroundColor = [0.90 0.90 0.90];
            app.StatusBar.FontSize = 10;
            app.StatusBar.Text = '  Ready';
            app.StatusBar.HorizontalAlignment = 'left';

            % Create TabGroup - starts at y=20 to leave room for status bar
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.AutoResizeChildren = 'off';
            app.TabGroup.Position = [1 20 867 529];
            app.TabGroup.SelectionChangedFcn = createCallbackFcn(app, @TabGroupSelectionChanged, true);

            % Create CleaningTab
            app.CleaningTab = uitab(app.TabGroup);
            app.CleaningTab.AutoResizeChildren = 'off';
            app.CleaningTab.Title = 'Cleaning';

            % Create StepsTree - a stage-grouped tree of the available steps
            % (see stepTaxonomy / populateStepsTree). Nodes are filled in
            % startupFcn via populateStepsTree so the empty tree here and the
            % populated one share one source. NodeData on each leaf is the exact
            % registry step name, which the Add/selection callbacks read.
            app.StepsTree = uitree(app.CleaningTab);
            app.StepsTree.FontSize = 11;
            app.StepsTree.Position = [10 173 207 294];
            app.StepsTree.SelectionChangedFcn = createCallbackFcn(app, @StepsTreeSelectionChanged, true);
            % Double-click a step to add it. DoubleClickedFcn is a recent uitree
            % addition; guard so older releases just fall back to the Add button.
            if isprop(app.StepsTree, 'DoubleClickedFcn')
                app.StepsTree.DoubleClickedFcn = createCallbackFcn(app, @StepsTreeDoubleClicked, true);
            end

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
            app.MoveUpButton.Position = [340 56 105 36];
            app.MoveUpButton.Text = 'Move Up';

            % Create MoveDownButton
            app.MoveDownButton = uibutton(app.CleaningTab, 'push');
            app.MoveDownButton.ButtonPushedFcn = createCallbackFcn(app, @MoveDownButtonPushed, true);
            app.MoveDownButton.BackgroundColor = [0.8 0.8 0.8];
            app.MoveDownButton.Position = [340 12 105 36];
            app.MoveDownButton.Text = 'Move Down';

            % Create AddButton
            app.AddButton = uibutton(app.CleaningTab, 'push');
            app.AddButton.ButtonPushedFcn = createCallbackFcn(app, @AddButtonPushed, true);
            app.AddButton.BackgroundColor = [0.8 0.8 0.8];
            app.AddButton.Position = [230 56 105 36];
            app.AddButton.Text = 'Add';

            % Create RemoveButton
            app.RemoveButton = uibutton(app.CleaningTab, 'push');
            app.RemoveButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveButtonPushed, true);
            app.RemoveButton.BackgroundColor = [0.8 0.8 0.8];
            app.RemoveButton.Position = [230 12 105 36];
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
            app.SelectDatatoPerformAnalysisPanel.Position = [649 237 208 206];

            % Create SelectedFilesListBox - shows all queued files
            app.SelectedFilesListBox = uilistbox(app.SelectDatatoPerformAnalysisPanel);
            app.SelectedFilesListBox.Items = {};
            app.SelectedFilesListBox.Position = [5 30 195 145];
            app.SelectedFilesListBox.FontSize = 10;

            % Create SelectDataButton (full-width "Browse..." at bottom of panel)
            app.SelectDataButton = uibutton(app.SelectDatatoPerformAnalysisPanel, 'push');
            app.SelectDataButton.ButtonPushedFcn = createCallbackFcn(app, @SelectDataButtonPushed, true);
            app.SelectDataButton.Position = [5 5 195 23];
            app.SelectDataButton.Text = 'Browse...';
            app.SelectDataButton.Tooltip = {'Browse a folder tree with checkboxes and a path filter to pick data files across many subject folders'};

            % Create RunAnalysisButton
            app.RunAnalysisButton = uibutton(app.CleaningTab, 'push');
            app.RunAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @RunAnalysisButtonPushed, true);
            app.RunAnalysisButton.BackgroundColor = [0.20 0.55 0.20];
            app.RunAnalysisButton.FontColor = [1 1 1];
            app.RunAnalysisButton.FontSize = 18;
            app.RunAnalysisButton.FontWeight = 'bold';
            app.RunAnalysisButton.Position = [657 117 201 60];
            app.RunAnalysisButton.Text = 'Run Analysis';

            % Create Image
            app.Image = uiimage(app.CleaningTab);
            app.Image.Position = [653 453 203 44];
            app.Image.ImageSource = fullfile(imageDir, 'LogoNest.jpg');

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
            app.ReStartStepsButton.Position = [658 193 201 36];
            app.ReStartStepsButton.Text = 'ReStart Steps';

            % Create ParallelCheckBox
            app.ParallelCheckBox = uicheckbox(app.CleaningTab);
            app.ParallelCheckBox.Text = 'Parallel Processing';
            app.ParallelCheckBox.Position = [657 85 201 24];
            app.ParallelCheckBox.Value = false;
            if license('test', 'Distrib_Computing_Toolbox')
                app.ParallelCheckBox.Enable = 'on';
            else
                app.ParallelCheckBox.Enable = 'off';
                app.ParallelCheckBox.Tooltip = 'Requires Parallel Computing Toolbox';
            end

            % Create ReportsTab
            app.ReportsTab = uitab(app.TabGroup);
            app.ReportsTab.AutoResizeChildren = 'off';
            app.ReportsTab.Title = 'Reports';

            % Reports tab - left column: session list
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

            app.ClearReportsButton = uibutton(app.ReportsTab, 'push');
            app.ClearReportsButton.ButtonPushedFcn = createCallbackFcn(app, @ClearReportsButtonPushed, true);
            app.ClearReportsButton.Position = [110 45 100 25];
            app.ClearReportsButton.Text = 'Clear List';
            app.ClearReportsButton.Tooltip = ['Empty the report list - both this session''s runs and ' ...
                'anything loaded from disk. Nothing on disk is deleted; Load from Folder brings ' ...
                'saved reports back'];

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

            % Reports tab - right column: what to show for the selected file.
            % Text and QC images describe the same file, so they share the
            % pane rather than competing for room in it.
            app.ReportsViewGroup = uibuttongroup(app.ReportsTab);
            app.ReportsViewGroup.BorderType = 'none';
            app.ReportsViewGroup.Position = [220 470 156 24];
            app.ReportsViewGroup.SelectionChangedFcn = createCallbackFcn(app, @ReportsViewChanged, true);

            app.ReportsTextViewButton = uitogglebutton(app.ReportsViewGroup);
            app.ReportsTextViewButton.Position = [0 0 60 24];
            app.ReportsTextViewButton.Text = 'Text';
            app.ReportsTextViewButton.Value = true;

            app.ReportsImageViewButton = uitogglebutton(app.ReportsViewGroup);
            app.ReportsImageViewButton.Position = [60 0 96 24];
            app.ReportsImageViewButton.Text = 'QC images';
            app.ReportsImageViewButton.Tooltip = ['The quality figures written during the run: ' ...
                'channel x trial attributes, ICA components, butterfly and PSD'];

            app.OpenReportSetButton = uibutton(app.ReportsTab, 'push');
            app.OpenReportSetButton.ButtonPushedFcn = createCallbackFcn(app, @OpenReportSetButtonPushed, true);
            app.OpenReportSetButton.Position = [382 470 68 24];
            app.OpenReportSetButton.Text = 'Open...';
            app.OpenReportSetButton.Tooltip = ['Open the cleaned recording this report describes ' ...
                'in the EEG scrolling viewer'];
            app.OpenReportSetButton.Enable = 'off';

            app.ExportReportsCSVButton = uibutton(app.ReportsTab, 'push');
            app.ExportReportsCSVButton.ButtonPushedFcn = createCallbackFcn(app, @ExportReportsCSVButtonPushed, true);
            app.ExportReportsCSVButton.Position = [455 470 145 24];
            app.ExportReportsCSVButton.Text = 'Export Metrics CSV';
            app.ExportReportsCSVButton.Tooltip = ['Write a CSV file: one row per file (channels and ' ...
                'trials retained, ICA components removed, quality verdict) for every report listed ' ...
                'here, including ones loaded from disk - so it can span several batch runs'];
            app.ExportReportsCSVButton.Enable = 'off';

            app.ExportPDFButton = uibutton(app.ReportsTab, 'push');
            app.ExportPDFButton.ButtonPushedFcn = createCallbackFcn(app, @ExportPDFButtonPushed, true);
            app.ExportPDFButton.Position = [605 470 110 24];
            app.ExportPDFButton.Text = 'Export PDF...';
            app.ExportPDFButton.Tooltip = ['Report text plus QC checkpoint images as a single PDF, ' ...
                'for the selected file or for every listed report. Not needed when ' ...
                'Auto-export PDF is on in Settings'];
            app.ExportPDFButton.Enable = 'off';

            app.CopyMethodsButton = uibutton(app.ReportsTab, 'push');
            app.CopyMethodsButton.ButtonPushedFcn = createCallbackFcn(app, @CopyMethodsButtonPushed, true);
            app.CopyMethodsButton.Position = [720 470 142 24];
            app.CopyMethodsButton.Text = 'Copy Methods Text';
            app.CopyMethodsButton.Tooltip = ['Copies the full parameterized methods paragraph for the ' ...
                'selected file - longer than the one-sentence note shown in the report - or the ' ...
                'cross-file aggregate when the Session Summary is selected'];
            app.CopyMethodsButton.Enable = 'off';

            app.ReportsTextArea = uitextarea(app.ReportsTab);
            app.ReportsTextArea.Editable = 'off';
            app.ReportsTextArea.FontName = 'Courier New';
            app.ReportsTextArea.FontSize = 10;
            app.ReportsTextArea.Position = [220 10 637 457];

            % Quality Dashboard panel - same rectangle as the text area,
            % hidden by default. Visible when the user picks the
            % synthetic "Session Quality Dashboard" entry in the listbox.
            app.ReportsDashboardPanel = uipanel(app.ReportsTab);
            app.ReportsDashboardPanel.Position = [220 10 637 457];
            app.ReportsDashboardPanel.BorderType = 'none';
            app.ReportsDashboardPanel.AutoResizeChildren = 'off';
            app.ReportsDashboardPanel.Visible = 'off';

            % QC images panel - same rectangle again, hidden by default.
            % Visible when the view switch is on QC images for a file entry.
            app.ReportsImagePanel = uipanel(app.ReportsTab);
            app.ReportsImagePanel.Position = [220 10 637 457];
            app.ReportsImagePanel.BorderType = 'none';
            app.ReportsImagePanel.AutoResizeChildren = 'off';
            app.ReportsImagePanel.Visible = 'off';

            % Floating hover tip for the Steps tree. Parented to the figure
            % rather than a tab so it is never clipped, and created last so it
            % paints over the TabGroup. See stepsTreeLegend for why this is not
            % the native Tooltip.
            app.StepsTipPanel = uipanel(app.UIFigure);
            app.StepsTipPanel.BorderType = 'line';
            app.StepsTipPanel.BackgroundColor = [1 1 0.88];
            app.StepsTipPanel.AutoResizeChildren = 'off';
            % Tall enough for the wrapped legend at this width; the text is
            % fixed, so a fixed box that fits it is simpler than measuring.
            app.StepsTipPanel.Position = [0 0 300 84];
            app.StepsTipPanel.Visible = 'off';

            app.StepsTipLabel = uilabel(app.StepsTipPanel);
            app.StepsTipLabel.Position = [8 5 284 74];
            app.StepsTipLabel.VerticalAlignment = 'top';
            app.StepsTipLabel.WordWrap = 'on';
            app.StepsTipLabel.FontSize = 11;


            %% ---- Explore tab -------------------------------------------
            % The workspace that replaces Visualizing + Analysis: groups on the
            % left, one plot from the registry in the middle, four ways out
            % along the bottom. Built last for now so the existing two tabs are
            % untouched while this is proven; they come out, and this moves into
            % their place, once it carries their work.
            app.ExploreTab = uitab(app.TabGroup);
            app.ExploreTab.AutoResizeChildren = 'off';
            app.ExploreTab.Title = 'Explore';

            % The rail is wide enough for FOUR columns. At 197 it was not: the
            % results view showed Win/Mean/Peak ms/Peak uV behind a horizontal
            % scrollbar, and the define view had no room for the polarity that
            % decides which way a peak is read. The width comes out of the
            % canvas, which has it to spare - the right edge is unchanged.
            RAIL_X = 8;
            RAIL_W = 250;
            MAIN_X = 268;
            MAIN_W = 592;

            HALF_W = floor((RAIL_W - 6) / 2);   % two buttons across the rail
            THIRD_W = floor((RAIL_W - 12) / 3); % three buttons across it

            % -- groups --------------------------------------------------
            app.ExploreGroupsLabel = uilabel(app.ExploreTab, ...
                'Text', 'GROUPS', 'FontWeight', 'bold', 'FontSize', 10, ...
                'Position', [RAIL_X 472 RAIL_W 18]);

            app.ExploreGroupsListBox = uilistbox(app.ExploreTab);
            app.ExploreGroupsListBox.Items = {};
            app.ExploreGroupsListBox.ItemsData = {};
            app.ExploreGroupsListBox.Position = [RAIL_X 404 RAIL_W 66];
            app.ExploreGroupsListBox.Tooltip = {'Each group is a set of recordings compared as one condition. n counts SUBJECTS - see Files... for which file belongs to whom.'};

            app.ExploreAddGroupButton = uibutton(app.ExploreTab, 'push');
            app.ExploreAddGroupButton.ButtonPushedFcn = createCallbackFcn(app, @ExploreAddGroupButtonPushed, true);
            app.ExploreAddGroupButton.Text = 'Add group...';
            app.ExploreAddGroupButton.Position = [RAIL_X 374 HALF_W 26];

            app.ExploreRemoveGroupButton = uibutton(app.ExploreTab, 'push');
            app.ExploreRemoveGroupButton.ButtonPushedFcn = createCallbackFcn(app, @ExploreRemoveGroupButtonPushed, true);
            app.ExploreRemoveGroupButton.Enable = 'off';
            app.ExploreRemoveGroupButton.Text = 'Remove';
            app.ExploreRemoveGroupButton.Position = [RAIL_X + HALF_W + 6 374 HALF_W 26];

            app.ExploreFilesButton = uibutton(app.ExploreTab, 'push');
            app.ExploreFilesButton.ButtonPushedFcn = createCallbackFcn(app, @ExploreFilesButtonPushed, true);
            app.ExploreFilesButton.Enable = 'off';
            app.ExploreFilesButton.Text = 'Files, subjects, groups...';
            app.ExploreFilesButton.Position = [RAIL_X 344 RAIL_W 26];
            app.ExploreFilesButton.Tooltip = {'See and correct which file belongs to which subject and group. This is where n comes from.'};

            % -- design --------------------------------------------------
            % An explicit control, not an inference. It used to be derived from
            % guessed subject ids, so a naming coincidence could switch to
            % paired and narrow every interval without saying so.
            app.ExploreDesignLabel = uilabel(app.ExploreTab, ...
                'Text', 'DESIGN', 'FontWeight', 'bold', 'FontSize', 10, ...
                'Position', [RAIL_X 320 RAIL_W 18]);

            app.ExploreDesignGroup = uibuttongroup(app.ExploreTab);
            app.ExploreDesignGroup.AutoResizeChildren = 'off';
            app.ExploreDesignGroup.BorderType = 'none';
            app.ExploreDesignGroup.SelectionChangedFcn = createCallbackFcn(app, @ExploreDesignChanged, true);
            app.ExploreDesignGroup.Position = [RAIL_X 294 RAIL_W 24];

            app.ExploreUnpairedButton = uiradiobutton(app.ExploreDesignGroup);
            app.ExploreUnpairedButton.Text = 'unpaired';
            app.ExploreUnpairedButton.Position = [0 1 78 22];
            app.ExploreUnpairedButton.Value = true;

            app.ExplorePairedButton = uiradiobutton(app.ExploreDesignGroup);
            app.ExplorePairedButton.Text = 'paired';
            app.ExplorePairedButton.Position = [96 1 78 22];
            app.ExplorePairedButton.Enable = 'off';

            app.ExploreDesignNoteLabel = uilabel(app.ExploreTab, ...
                'Position', [RAIL_X 274 RAIL_W 18], 'FontSize', 11, ...
                'FontColor', [0.35 0.38 0.43]);

            % -- region of interest --------------------------------------
            app.ExploreRoiLabel = uilabel(app.ExploreTab, ...
                'Text', 'REGION OF INTEREST', 'FontWeight', 'bold', ...
                'FontSize', 10, 'Position', [RAIL_X 250 RAIL_W 18]);

            app.ExploreRoiDropDown = uidropdown(app.ExploreTab);
            app.ExploreRoiDropDown.Items = {};
            app.ExploreRoiDropDown.ValueChangedFcn = createCallbackFcn(app, @ExploreRoiDropDownValueChanged, true);
            app.ExploreRoiDropDown.Position = [RAIL_X 224 RAIL_W 24];

            app.ExploreRoiEditButton = uibutton(app.ExploreTab, 'push');
            app.ExploreRoiEditButton.ButtonPushedFcn = createCallbackFcn(app, @ExploreRoiEditButtonPushed, true);
            app.ExploreRoiEditButton.Text = 'Edit electrodes...';
            app.ExploreRoiEditButton.Position = [RAIL_X 196 RAIL_W 24];

            app.ExploreRoiSummaryLabel = uilabel(app.ExploreTab, ...
                'Position', [RAIL_X 174 RAIL_W 20], 'FontSize', 11, ...
                'FontColor', [0.35 0.38 0.43]);

            % -- windows of interest -------------------------------------
            app.ExploreWindowsLabel = uilabel(app.ExploreTab, ...
                'Text', 'WINDOWS', 'FontWeight', 'bold', ...
                'FontSize', 10, 'Position', [RAIL_X 150 70 18]);

            % Define or measure, in the same space. The Analysis tab showed the
            % window bounds AND their measures in one table; the rail is too
            % narrow for six columns, so it switches instead of dropping three.
            app.ExploreWindowsModeDropDown = uidropdown(app.ExploreTab);
            app.ExploreWindowsModeDropDown.Items = {'define', 'results'};
            app.ExploreWindowsModeDropDown.ValueChangedFcn = createCallbackFcn(app, @ExploreWindowsModeChanged, true);
            app.ExploreWindowsModeDropDown.Position = [RAIL_X + RAIL_W - 108 148 108 22];

            app.ExploreWindowsTable = uitable(app.ExploreTab);
            app.ExploreWindowsTable.ColumnName = {'Name'; 'T1'; 'T2'; 'Peak'};
            app.ExploreWindowsTable.ColumnWidth = exploreWindowColWidths(app);
            app.ExploreWindowsTable.ColumnEditable = [true true true true];
            app.ExploreWindowsTable.ColumnFormat = {'char', 'numeric', 'numeric', ...
                                                    {'auto', 'pos', 'neg'}};
            app.ExploreWindowsTable.CellEditCallback = createCallbackFcn(app, @ExploreWindowsTableCellEdit, true);
            app.ExploreWindowsTable.RowName = {};
            app.ExploreWindowsTable.Position = [RAIL_X 42 RAIL_W 104];

            app.ExploreWindowsAddButton = uibutton(app.ExploreTab, 'push');
            app.ExploreWindowsAddButton.ButtonPushedFcn = createCallbackFcn(app, @ExploreWindowsAddButtonPushed, true);
            app.ExploreWindowsAddButton.Text = 'Add';
            app.ExploreWindowsAddButton.Position = [RAIL_X 14 THIRD_W 24];

            app.ExploreWindowsRemoveButton = uibutton(app.ExploreTab, 'push');
            app.ExploreWindowsRemoveButton.ButtonPushedFcn = createCallbackFcn(app, @ExploreWindowsRemoveButtonPushed, true);
            app.ExploreWindowsRemoveButton.Text = 'Remove';
            app.ExploreWindowsRemoveButton.Position = [RAIL_X + THIRD_W + 6 14 THIRD_W 24];

            app.ExploreWindowsResetButton = uibutton(app.ExploreTab, 'push');
            app.ExploreWindowsResetButton.ButtonPushedFcn = createCallbackFcn(app, @ExploreWindowsResetButtonPushed, true);
            app.ExploreWindowsResetButton.Text = 'Reset';
            app.ExploreWindowsResetButton.Position = [RAIL_X + 2*THIRD_W + 12 14 THIRD_W 24];

            % -- plot picker ---------------------------------------------
            app.ExplorePlotLabel = uilabel(app.ExploreTab, ...
                'Text', 'Plot', 'FontWeight', 'bold', ...
                'Position', [MAIN_X 470 32 22]);

            app.ExplorePlotDropDown = uidropdown(app.ExploreTab);
            app.ExplorePlotDropDown.Items = {};
            app.ExplorePlotDropDown.ValueChangedFcn = createCallbackFcn(app, @ExplorePlotDropDownValueChanged, true);
            app.ExplorePlotDropDown.Position = [MAIN_X + 36 470 300 24];

            % Plot-specific settings live here rather than on the rail. The
            % rail is for things that describe the DATASET - groups, ROI,
            % windows - and it is already crowded; a setting that means
            % something only while one plot is chosen does not belong beside
            % them, and putting it there would grow the rail with every plot
            % added.
            app.ExplorePlotOptionsButton = uibutton(app.ExploreTab, ...
                'Text', 'Options...', 'Position', [MAIN_X + 344 470 84 24], ...
                'ButtonPushedFcn', createCallbackFcn(app, @ExplorePlotOptionsButtonPushed, true));

            app.ExplorePlotInfoLabel = uilabel(app.ExploreTab, ...
                'Position', [MAIN_X + 436 466 MAIN_W - 436 28], ...
                'FontSize', 11, 'WordWrap', 'on', ...
                'VerticalAlignment', 'center', 'FontColor', [0.55 0.33 0.10]);

            % -- canvas --------------------------------------------------
            % A panel rather than a fixed axes: a topography draws one map per
            % group and the bars one panel per window, so the number of axes is
            % a property of the chosen plot. renderExplorePlot fills this.
            app.ExploreCanvas = uipanel(app.ExploreTab);
            app.ExploreCanvas.BorderType = 'none';
            app.ExploreCanvas.AutoResizeChildren = 'off';
            app.ExploreCanvas.Position = [MAIN_X 66 MAIN_W 396];

            app.ExploreEmptyLabel = uilabel(app.ExploreCanvas, ...
                'Position', [20 180 MAIN_W - 40 40], ...
                'HorizontalAlignment', 'center', 'FontSize', 13, ...
                'FontColor', [0.45 0.48 0.53], 'WordWrap', 'on', ...
                'Text', 'Add a group to begin. A group is a set of recordings compared as one condition - pre and post, or one cohort against another.');

            % -- the four ways out ---------------------------------------
            EXIT_W = 155;
            EXIT_GAP = 8;
            app.ExploreFigureButton = uibutton(app.ExploreTab, 'push');
            app.ExploreFigureButton.ButtonPushedFcn = createCallbackFcn(app, @ExploreFigureButtonPushed, true);
            app.ExploreFigureButton.Text = 'Figure...';
            app.ExploreFigureButton.Enable = 'off';
            app.ExploreFigureButton.Position = [MAIN_X 30 EXIT_W 26];
            app.ExploreFigureButton.Tooltip = {'Open this plot in a standard MATLAB figure, for editing and saving at publication resolution'};

            app.ExploreCsvButton = uibutton(app.ExploreTab, 'push');
            app.ExploreCsvButton.ButtonPushedFcn = createCallbackFcn(app, @ExploreCsvButtonPushed, true);
            app.ExploreCsvButton.Text = 'Measures -> CSV...';
            app.ExploreCsvButton.Enable = 'off';
            app.ExploreCsvButton.Position = [MAIN_X + (EXIT_W + EXIT_GAP) 30 EXIT_W 26];
            app.ExploreCsvButton.Tooltip = {'One row per group x subject x window - the small tabular form, for R, JASP or Prism'};

            app.ExploreResultsButton = uibutton(app.ExploreTab, 'push');
            app.ExploreResultsButton.ButtonPushedFcn = createCallbackFcn(app, @ExploreResultsButtonPushed, true);
            app.ExploreResultsButton.Text = 'Results -> MATLAB...';
            app.ExploreResultsButton.Enable = 'off';
            app.ExploreResultsButton.Position = [MAIN_X + 2*(EXIT_W + EXIT_GAP) 30 EXIT_W 26];
            app.ExploreResultsButton.Tooltip = {'The whole result as a struct - curves at sampling rate, intervals, provenance - saved or sent to the workspace'};

            app.ExploreStatusLabel = uilabel(app.ExploreTab, ...
                'Position', [MAIN_X 6 MAIN_W 20], 'FontSize', 11, ...
                'FontColor', [0.35 0.38 0.43], 'Text', 'Ready.');

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
