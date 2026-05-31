function rescaleComponents(app, sX, sY)
% RESCALECOMPONENTS  Rescale all nestapp UI components proportionally.
%   Called from UIFigureSizeChanged. This file lives in src/@nestapp/ - a
%   proper class method with access to protected App Designer methods.
            STATUS_H = 20;

            sf = min(sX, sY);
            p  = @(o) round(o .* [sX, sY, sX, sY]);
            % ph: scale x/y/w but keep original height - for fixed-height controls
            % (uicheckbox, uieditfield, uispinner, uirangedslider) that emit
            % "height cannot be changed" warnings when h is set to a scaled value.
            ph = @(o) [round(o(1)*sX), round(o(2)*sY), round(o(3)*sX), o(4)];
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
            app.SelectDatatoPerformAnalysisPanel.Position = p([649 237 208 206]);
            app.SelectedFilesListBox.Position     = p([5 30 195 145]);
            app.SelectDataButton.Position         = p([5 5 195 23]);

            app.RunAnalysisButton.Position        = p([657 117 201 60]);
            app.RunAnalysisButton.FontSize        = fs(18);
            app.Image.Position                    = p([653 453 203 44]);
            app.NESTAPPLabel.Position             = p([785 448 71 22]);
            app.NESTAPPLabel.FontSize             = fs(14);

            app.ReStartStepsButton.Position       = p([658 193 201 36]);
            app.ReStartStepsButton.FontSize       = fs(18);
            app.ParallelCheckBox.Position         = ph([657 85 201 24]);

            %% Visualizing Tab
            % Three zones: left (x:0-340 electrode map), center (x:340-648 TEP/topo),
            % right (x:651-867 file selection). Bottom strip (y:0-165) holds controls.
            % TEP window slider lives above UIAxes; topoplot controls sit right of UIAxes2.

            % Left action column - bottom-aligned with ReLoad button (y=7)
            % TOPOPLOT joins this group as the lowest button
            app.PLOTTEPButton.Position                = p([5 88 140 30]);
            app.ShowComponentsButton.Position         = p([5 61 140 23]);
            app.ExportTEPFigureButton.Position        = p([5 34 140 23]);
            app.TOPOPLOTButton.Position               = p([5 7 140 23]);

            % Center-left controls - bottom strip (x:152-340)
            % PlottingModeButtonGroup sits above the single-row topoplot controls
            app.PlottingModeButtonGroup.Position  = p([152 36 150 67]);
            app.NewFigureButton.Position          = p([11 21 83 22]);
            app.AddtocurrentFigureButton.Position = p([11 -1 135 22]);
            % Topoplot time and window on one line - 3-digit fields
            app.TopoplottimeSpinnerLabel.Position = p([152 10 35 22]);
            app.TopoplottimeSpinner.Position      = ph([189 10 52 22]);
            app.WindowsizeforTopoplotLabel.Position = p([245 10 35 22]);
            app.WindowsizefortimeaveragedTopoplotEditField.Position = ph([282 10 52 22]);

            % Center column - TEP window slider above the TEP plot
            % TEP plot (60%) and topoplot (40%) of 448px available - slider in gap
            app.UIAxes.Position                   = p([340 230 308 270]);
            app.TEPWindowSliderLabel.Position     = p([380 204 130 16]);
            app.TEPWindowSlider.Position          = ph([380 193 268 3]);
            app.UIAxes2.Position                  = p([340 7 308 179]);

            % Head image (electrode map) - unchanged
            app.Image2.Position                   = p([-1 165 350 336]);

            % Right column - data selection
            app.SelectDatatoVisulaizeTEPsPanel.Position  = p([651 406 208 90]);
            app.FolderEditField_2Label.Position   = p([1 41 40 22]);
            app.FolderEditField_2.Position        = ph([49 41 145 22]);
            app.SelectDataButton_2.Position       = p([13 10 183 23]);
            app.FilesListBoxLabel.Position        = p([740 382 30 22]);
            app.FilesListBox.Position             = p([669 71 183 306]);
            app.SelectAllCheckBox.Position        = ph([670 46 71 22]);
            app.DontfindcommonelectrodesCheckBox.Position = ph([670 28 180 22]);
            app.ReLoadAvailableElectrodesButton.Position  = p([686 7 153 23]);

            %% Electrode buttons (Visualizing Tab - 64 buttons)
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

            %% Analysis Tab
            app.AnalysisSelPanel.Position          = p([10 430 847 55]);
            app.AnalysisSelectionLabel.Position    = p([10 5 820 32]);
            app.AnalysisCompWindowsLabel.Position  = p([10 407 300 18]);
            app.TEPComponentTable.Position         = p([10 225 360 178]);
            app.TEPComponentTable.ColumnWidth      = {'auto', 'auto', 'auto'};
            app.EditComponentWindowsButton.Position = p([10 196 220 25]);
            app.AnalysisWorkspaceLabel.Position    = p([450 407 380 18]);
            app.ExportTEPDataButton.Position       = p([450 374 220 28]);
            app.TEPvarNameEditFieldLabel.Position  = p([450 348 60 22]);
            app.TEPvarNameEditField.Position       = ph([515 348 155 22]);
            app.AnalysisBatchLabel.Position        = p([450 313 380 18]);
            app.AnalysisBatchDescLabel.Position    = p([450 291 380 18]);
            app.ExtractPeaksCSVButton.Position     = p([450 254 220 32]);
            app.AnalysisStatusLabel.Position       = p([10 15 847 22]);

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
