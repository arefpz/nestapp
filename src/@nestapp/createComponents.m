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
            pathToMLAPP = fileparts(fileparts(mfilename('fullpath')));

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
                'MenuSelectedFcn', createCallbackFcn(app, @UIFigureCloseRequest, true));

            mSettings = uimenu(app.UIFigure, 'Text', 'Settings');
            uimenu(mSettings, 'Text', 'Preferences...', ...
                'MenuSelectedFcn', createCallbackFcn(app, @openPreferencesMenu, true));

            mTools = uimenu(app.UIFigure, 'Text', 'Tools');
            uimenu(mTools, 'Text', 'Browse Raw EEG...', ...
                'MenuSelectedFcn', createCallbackFcn(app, @PlotEEGdataButtonPushed, true));

            mHelp = uimenu(app.UIFigure, 'Text', 'Help');
            uimenu(mHelp, 'Text', 'About nestapp', ...
                'MenuSelectedFcn', createCallbackFcn(app, @showAboutMenu, true));

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

            % Create StepsListBox - items derived from stepRegistry, not hardcoded
            reg_init = stepRegistry();
            app.StepsListBox = uilistbox(app.CleaningTab);
            app.StepsListBox.Items = {reg_init.name};
            app.StepsListBox.ValueChangedFcn = createCallbackFcn(app, @StepsListBoxValueChanged, true);
            app.StepsListBox.FontSize = 11;
            app.StepsListBox.ClickedFcn = createCallbackFcn(app, @StepsListBoxClicked, true);
            app.StepsListBox.Position = [10 173 207 294];
            app.StepsListBox.Value = reg_init(1).name;

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
            app.SelectDatatoPerformAnalysisPanel.Position = [649 237 208 206];

            % Create SelectedFilesListBox - shows all queued files
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
            app.RunAnalysisButton.Position = [657 117 201 60];
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
            app.UIAxes.Position = [340 230 308 270];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.VisualizingTab);
            app.UIAxes2.TickLabelInterpreter = 'none';
            app.UIAxes2.XAxisLocation = 'origin';
            app.UIAxes2.XTick = [];
            app.UIAxes2.YAxisLocation = 'origin';
            app.UIAxes2.YTick = [];
            app.UIAxes2.ZTick = [];
            app.UIAxes2.Position = [340 7 308 179];

            % Create ShowComponentsButton
            app.ShowComponentsButton = uibutton(app.VisualizingTab, 'state');
            app.ShowComponentsButton.ValueChangedFcn = createCallbackFcn(app, @ShowComponentsButtonValueChanged, true);
            app.ShowComponentsButton.Text = 'Show Components';
            app.ShowComponentsButton.Position = [5 61 140 23];

            % Create PLOTTEPButton
            app.PLOTTEPButton = uibutton(app.VisualizingTab, 'push');
            app.PLOTTEPButton.ButtonPushedFcn = createCallbackFcn(app, @PLOTTEPButtonPushed, true);
            app.PLOTTEPButton.Enable = 'off';
            app.PLOTTEPButton.Position = [5 88 140 30];
            app.PLOTTEPButton.Text = 'PLOT TEP';

            % Create SelectDatatoVisulaizeTEPsPanel
            app.SelectDatatoVisulaizeTEPsPanel = uipanel(app.VisualizingTab);
            app.SelectDatatoVisulaizeTEPsPanel.AutoResizeChildren = 'off';
            app.SelectDatatoVisulaizeTEPsPanel.BorderType = 'none';
            app.SelectDatatoVisulaizeTEPsPanel.Title = 'Select Data to Visualize TEPs';
            app.SelectDatatoVisulaizeTEPsPanel.Position = [651 406 208 90];

            % Create FolderEditField_2Label
            app.FolderEditField_2Label = uilabel(app.SelectDatatoVisulaizeTEPsPanel);
            app.FolderEditField_2Label.HorizontalAlignment = 'right';
            app.FolderEditField_2Label.Position = [1 41 40 22];
            app.FolderEditField_2Label.Text = 'Folder';

            % Create FolderEditField_2
            app.FolderEditField_2 = uieditfield(app.SelectDatatoVisulaizeTEPsPanel, 'text');
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
            app.UseCurrentlyCleanedDataCheckBox.Visible = 'off';

            % Create FilesListBoxLabel
            app.FilesListBoxLabel = uilabel(app.VisualizingTab);
            app.FilesListBoxLabel.HorizontalAlignment = 'right';
            app.FilesListBoxLabel.Position = [740 382 30 22];
            app.FilesListBoxLabel.Text = 'Files';

            % Create FilesListBox
            app.FilesListBox = uilistbox(app.VisualizingTab);
            app.FilesListBox.Items = {};
            app.FilesListBox.Multiselect = 'on';
            app.FilesListBox.ValueChangedFcn = createCallbackFcn(app, @FilesListBoxValueChanged, true);
            app.FilesListBox.Position = [669 71 183 306];
            app.FilesListBox.Value = {};

            % Create Image2
            app.Image2 = uiimage(app.VisualizingTab);
            app.Image2.Position = [-1 165 350 336];
            app.Image2.ImageSource = fullfile(pathToMLAPP, 'Head.png');

            % Create WindowsizeforTopoplotLabel
            app.WindowsizeforTopoplotLabel = uilabel(app.VisualizingTab);
            app.WindowsizeforTopoplotLabel.HorizontalAlignment = 'center';
            app.WindowsizeforTopoplotLabel.Position = [245 10 35 22];
            app.WindowsizeforTopoplotLabel.Text = 'Win';

            % Create WindowsizefortimeaveragedTopoplotEditField
            app.WindowsizefortimeaveragedTopoplotEditField = uieditfield(app.VisualizingTab, 'numeric');
            app.WindowsizefortimeaveragedTopoplotEditField.ValueDisplayFormat = '%.0f';
            app.WindowsizefortimeaveragedTopoplotEditField.Position = [282 10 52 22];

            % Create TOPOPLOTButton
            app.TOPOPLOTButton = uibutton(app.VisualizingTab, 'push');
            app.TOPOPLOTButton.ButtonPushedFcn = createCallbackFcn(app, @TOPOPLOTButtonPushed, true);
            app.TOPOPLOTButton.Enable = 'off';
            app.TOPOPLOTButton.Position = [5 7 140 23];
            app.TOPOPLOTButton.Text = 'TOPOPLOT';

            % Create ExportTEPFigureButton
            app.ExportTEPFigureButton = uibutton(app.VisualizingTab, 'push');
            app.ExportTEPFigureButton.ButtonPushedFcn = createCallbackFcn(app, @ExportTEPFigureButtonPushed, true);
            app.ExportTEPFigureButton.Enable = 'off';
            app.ExportTEPFigureButton.Position = [5 34 140 23];
            app.ExportTEPFigureButton.Text = 'Export TEP Figure';

            % Create PlottingModeButtonGroup
            app.PlottingModeButtonGroup = uibuttongroup(app.VisualizingTab);
            app.PlottingModeButtonGroup.AutoResizeChildren = 'off';
            app.PlottingModeButtonGroup.BorderType = 'none';
            app.PlottingModeButtonGroup.Title = 'Plotting Mode';
            app.PlottingModeButtonGroup.Position = [152 36 150 67];

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
            app.TEPWindowSliderLabel.Position = [380 204 130 16];
            app.TEPWindowSliderLabel.Text = 'TEP Window';

            % Create TEPWindowSlider
            app.TEPWindowSlider = uislider(app.VisualizingTab, 'range');
            app.TEPWindowSlider.Limits = [-100 300];
            app.TEPWindowSlider.ValueChangingFcn = createCallbackFcn(app, @TEPWindowSliderValueChanging, true);
            app.TEPWindowSlider.Position = [380 193 268 3];
            app.TEPWindowSlider.Value = [-50 300];

            % Create TopoplottimeSpinnerLabel
            app.TopoplottimeSpinnerLabel = uilabel(app.VisualizingTab);
            app.TopoplottimeSpinnerLabel.HorizontalAlignment = 'right';
            app.TopoplottimeSpinnerLabel.Position = [152 10 35 22];
            app.TopoplottimeSpinnerLabel.Text = 'Time';

            % Create TopoplottimeSpinner
            app.TopoplottimeSpinner = uispinner(app.VisualizingTab);
            app.TopoplottimeSpinner.RoundFractionalValues = 'on';
            app.TopoplottimeSpinner.ValueDisplayFormat = '%.0f';
            app.TopoplottimeSpinner.ValueChangedFcn = createCallbackFcn(app, @TopoplottimeSpinnerValueChanged, true);
            app.TopoplottimeSpinner.Position = [189 10 52 22];
            app.TopoplottimeSpinner.Value = 60;

            % Create AnalysisTab
            app.AnalysisTab = uitab(app.TabGroup);
            app.AnalysisTab.AutoResizeChildren = 'off';
            app.AnalysisTab.Title = 'Analysis';

            % Analysis tab - current selection summary panel (near top)
            app.AnalysisSelPanel = uipanel(app.AnalysisTab, 'Title', 'Current Selection', ...
                'AutoResizeChildren', 'off', ...
                'Position', [10 430 847 55]);
            app.AnalysisSelectionLabel = uilabel(app.AnalysisSelPanel, ...
                'Position', [10 5 820 32], ...
                'Text', 'Select files and ROI electrodes on the Visualizing tab.', ...
                'WordWrap', 'on', 'FontSize', 11);

            % Analysis tab - LEFT column: component windows
            app.AnalysisCompWindowsLabel = uilabel(app.AnalysisTab, 'Position', [10 407 300 18], ...
                'Text', 'COMPONENT WINDOWS', 'FontWeight', 'bold', 'FontSize', 10);

            % TEPComponentTable - taller to show all 6 components without scrolling
            app.TEPComponentTable = uitable(app.AnalysisTab);
            app.TEPComponentTable.ColumnName  = {'Component', 'Latency (ms)', 'Amplitude (uV)'};
            app.TEPComponentTable.ColumnWidth = {'auto', 'auto', 'auto'};
            app.TEPComponentTable.RowName     = {};
            app.TEPComponentTable.Enable      = 'on';
            app.TEPComponentTable.Position    = [10 225 360 178];

            app.EditComponentWindowsButton = uibutton(app.AnalysisTab, 'push');
            app.EditComponentWindowsButton.ButtonPushedFcn = createCallbackFcn(app, @EditComponentWindowsButtonPushed, true);
            app.EditComponentWindowsButton.Text     = 'Edit Component Windows...';
            app.EditComponentWindowsButton.Position = [10 196 220 25];

            % Analysis tab - RIGHT column: workspace export + batch extraction grouped
            app.AnalysisWorkspaceLabel = uilabel(app.AnalysisTab, 'Position', [450 407 380 18], ...
                'Text', 'WORKSPACE EXPORT', 'FontWeight', 'bold', 'FontSize', 10);

            app.ExportTEPDataButton = uibutton(app.AnalysisTab, 'push');
            app.ExportTEPDataButton.ButtonPushedFcn = createCallbackFcn(app, @ExportTEPDataButtonPushed, true);
            app.ExportTEPDataButton.Enable   = 'off';
            app.ExportTEPDataButton.Text     = 'Export TEP to Workspace';
            app.ExportTEPDataButton.Position = [450 374 220 28];

            app.TEPvarNameEditFieldLabel = uilabel(app.AnalysisTab);
            app.TEPvarNameEditFieldLabel.HorizontalAlignment = 'right';
            app.TEPvarNameEditFieldLabel.Enable   = 'off';
            app.TEPvarNameEditFieldLabel.Position = [450 348 60 22];
            app.TEPvarNameEditFieldLabel.Text     = 'Variable:';

            app.TEPvarNameEditField = uieditfield(app.AnalysisTab, 'text');
            app.TEPvarNameEditField.ValueChangedFcn = createCallbackFcn(app, @TEPvarNameEditFieldValueChanged, true);
            app.TEPvarNameEditField.Enable   = 'off';
            app.TEPvarNameEditField.Position = [515 348 155 22];
            app.TEPvarNameEditField.Value    = 'TEPdata';

            app.AnalysisBatchLabel = uilabel(app.AnalysisTab, 'Position', [450 313 380 18], ...
                'Text', 'BATCH EXTRACTION', 'FontWeight', 'bold', 'FontSize', 10);
            app.AnalysisBatchDescLabel = uilabel(app.AnalysisTab, 'Position', [450 291 380 18], ...
                'Text', ['Extract peak latency and amplitude from each file. ' ...
                    'Results saved as CSV for import into R/SPSS/Excel.'], ...
                'WordWrap', 'on', 'FontSize', 9, 'FontColor', [0.4 0.4 0.4]);

            app.ExtractPeaksCSVButton = uibutton(app.AnalysisTab, 'push');
            app.ExtractPeaksCSVButton.ButtonPushedFcn = createCallbackFcn(app, @ExtractPeaksCSVButtonPushed, true);
            app.ExtractPeaksCSVButton.Text    = 'Extract Peaks  ->  CSV';
            app.ExtractPeaksCSVButton.Position = [450 254 220 32];
            app.ExtractPeaksCSVButton.Tooltip = ...
                'Run peak detection across all selected files and save results as a CSV table';

            app.AnalysisStatusLabel = uilabel(app.AnalysisTab);
            app.AnalysisStatusLabel.Position   = [10 15 847 22];
            app.AnalysisStatusLabel.Text       = 'Ready.';
            app.AnalysisStatusLabel.FontSize   = 10;
            app.AnalysisStatusLabel.FontColor  = [0.4 0.4 0.4];

            % Create EEGDatasetDropDownLabel (kept for PlotEEGdataButtonPushed callback)
            app.EEGDatasetDropDownLabel = uilabel(app.VisualizingTab);
            app.EEGDatasetDropDownLabel.HorizontalAlignment = 'right';
            app.EEGDatasetDropDownLabel.Position = [152 58 75 22];
            app.EEGDatasetDropDownLabel.Text = 'EEG Dataset';
            app.EEGDatasetDropDownLabel.Visible = 'off';

            % Create EEGDatasetDropDown (kept for PlotEEGdataButtonPushed; hidden)
            app.EEGDatasetDropDown = uidropdown(app.VisualizingTab);
            app.EEGDatasetDropDown.Items = {'Select a file'};
            app.EEGDatasetDropDown.ValueChangedFcn = createCallbackFcn(app, @EEGDatasetDropDownValueChanged, true);
            app.EEGDatasetDropDown.Enable = 'off';
            app.EEGDatasetDropDown.Position = [230 58 100 22];
            app.EEGDatasetDropDown.Value = 'Select a file';
            app.EEGDatasetDropDown.Visible = 'off';

            % Create PlotEEGdataButton (hidden; triggered via Tools menu)
            app.PlotEEGdataButton = uibutton(app.VisualizingTab, 'push');
            app.PlotEEGdataButton.ButtonPushedFcn = createCallbackFcn(app, @PlotEEGdataButtonPushed, true);
            app.PlotEEGdataButton.Enable  = 'off';
            app.PlotEEGdataButton.Visible = 'off';
            app.PlotEEGdataButton.Position = [5 26 108 23];
            app.PlotEEGdataButton.Text = 'Plot EEG data';

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

            % Reports tab - right column: report text + actions
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
