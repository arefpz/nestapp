% WARNING: Do not open nestapp_designer.mlapp and save — App Designer will
% regenerate this file and overwrite startupFcn and other hand-edited methods.
% All edits must be made directly to nestapp.m.
classdef nestapp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        CleaningTab                     matlab.ui.container.Tab
        RunningstepEditField            matlab.ui.control.EditField
        RunningstepEditFieldLabel       matlab.ui.control.Label
        ProcessingfileEditField         matlab.ui.control.EditField
        ProcessingfileEditFieldLabel    matlab.ui.control.Label
        outofEditField                  matlab.ui.control.EditField
        outofEditFieldLabel             matlab.ui.control.Label
        ReStartStepsButton              matlab.ui.control.Button
        EEGLABpathifalreadyisnotinpathPanel  matlab.ui.container.Panel
        SelectEEGLABFolderButton        matlab.ui.control.Button
        PathEditField                   matlab.ui.control.EditField
        PathEditFieldLabel              matlab.ui.control.Label
        NESTAPPLabel                    matlab.ui.control.Label
        Image                           matlab.ui.control.Image
        RunAnalysisButton               matlab.ui.control.Button
        SelectDatatoPerformAnalysisPanel  matlab.ui.container.Panel
        SelectDataButton                matlab.ui.control.Button
        FileEditField                   matlab.ui.control.EditField
        FileEditFieldLabel              matlab.ui.control.Label
        FolderEditField                 matlab.ui.control.EditField
        FolderEditFieldLabel            matlab.ui.control.Label
        SelectedListBoxLabel_2          matlab.ui.control.Label
        TextArea                        matlab.ui.control.TextArea
        DefaultValueButton              matlab.ui.control.Button
        UITable                         matlab.ui.control.Table
        SelectedListBoxLabel            matlab.ui.control.Label
        SavePipelineButton              matlab.ui.control.Button
        LoadPipelineButton              matlab.ui.control.Button
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
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes                          matlab.ui.control.UIAxes
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
        ChangedVal % Non-Default Value for each eeglab function
        DefaultsVal % Storing default values for each eeglab function
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
        DefaulTEPxLim = [-100 300]; % Defaul xLim for time in TEP
        EEGtime
        TEP2Export
    end

    methods (Access = private)
        function styleParamTable(app)
        % Grey out UITable rows whose Value is a display placeholder.
        % Placeholders start with '(' by convention, e.g. '(all channels)'.
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
                    if strlength(sv) > 0 && startsWith(sv, '(')
                        addStyle(app.UITable, grey, 'cell', [row, 2]);
                    end
                end
            end
        end

        function LoadSelecEEGdata(app)
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
        
        function plotTEP(app) %%%%%% LEGEND IS NOT GOOD
            if ~app.EEG_SelectedTEPFiles_Loaded
                LoadSelecEEGdata(app)
            end
            for nfile = 1:numel(app.SelectedFilesforTEP)
                EEGaux = app.EEGofAllSelectedFiles{1,nfile};
                ROIind = find(ismember({EEGaux.chanlocs.labels},app.ROIelecsLabels));
                TEP_ROI(nfile,:)=mean(mean(EEGaux.data(ROIind,:,:),3,'omitmissing'),1,"omitmissing");
            end
            swin = 5;
            app.TEP2Export = TEP_ROI;
            TEP_ROISD = std(TEP_ROI,1,1)/sqrt(size(TEP_ROI,1));
            Colr = rand(1,3);
            meanx = smoothdata(mean(TEP_ROI,1,"omitmissing"),'movmean',swin); sdx = smoothdata(TEP_ROISD,'movmean',swin);
            xf=[app.EEGtime(1) app.EEGtime  app.EEGtime(end) app.EEGtime(end:-1:1)];
            yf=[meanx(1)-sdx(1)/2 meanx+sdx/2 meanx(end)-sdx(end)/2 meanx(end:-1:1)-sdx(end:-1:1)/2];
            

            if app.NewFigureButton.Value
                cla(app.UIAxes, 'reset');
                hold(app.UIAxes, 'on');
                fill(app.UIAxes,xf,yf,Colr(1,:),'FaceAlpha',0.5,'LineStyle','none')
                plot(app.UIAxes,app.EEGtime,meanx,'Color',Colr(1,:),'LineWidth',2);
                hold(app.UIAxes, 'off');
                xlim(app.UIAxes,app.DefaulTEPxLim);
            elseif app.AddtocurrentFigureButton.Value
                hold(app.UIAxes, 'on');
                fill(app.UIAxes,xf,yf,Colr(1,:),'FaceAlpha',0.5,'LineStyle','none')
                plot(app.UIAxes,app.EEGtime,meanx,'Color',Colr(1,:),'LineWidth',2);
                xlim(app.UIAxes,app.DefaulTEPxLim);
            end
                
            
            
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
    end


    % Callbacks that handle component events
    methods (Access = private)

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
                    var{j} = string(fields{j});
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
            end

            app.SelectedListBox.Items(:) = [];
            app.SelectedListBox.ItemsData(:) = [];
            app.UITable.Data = [];
            app.ItemNum = 1;
            app.needchanloc = 1;
            clc
        end

        % Clicked callback: StepsListBox
        function StepsListBoxClicked(app, event)
            item = event.InteractionInformation.Item;

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
                app.ChangedVal{1}=app.DefaultsVal{ind};
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
                app.ChangedVal{app.ItemNum+1}=app.DefaultsVal{ind};
            end
            % app.updateAppLayout
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
                app.ChangedVal{ind1} = C;
            end
            app.SelectedListBox.Value=app.SelectedListBox.ItemsData{ind1-1};
        end

        % Button pushed function: SavePipelineButton
        function SavePipelineButtonPushed(app, event)
            PLItems = app.SelectedListBox.Items;
            PLItemsData = app.SelectedListBox.ItemsData;
            VarIns = app.ChangedVal;
            uisave({'PLItems','PLItemsData','VarIns'},'*.mat')
        end

        % Button pushed function: RemoveButton
        function RemoveButtonPushed(app, event)
            ind = find(ismember(app.SelectedListBox.ItemsData,app.SelectedListBox.Value));
            app.SelectedListBox.Items(ind) = [];
            app.SelectedListBox.ItemsData(ind) =[];
            app.ChangedVal(ind)=[];
            for i = ind : length(app.SelectedListBox.ItemsData)
                app.SelectedListBox.ItemsData{i} = ['Item',num2str(i)];
            end
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
                app.ChangedVal{ind1} = C;
            end
            app.SelectedListBox.Value=app.SelectedListBox.ItemsData{ind1+1};
        end

        % Button pushed function: LoadPipelineButton
        function LoadPipelineButtonPushed(app, event)
            [pName,pPath] = uigetfile('*.mat');
            pipeline = load([pPath,pName],'-mat','PLItems','PLItemsData','VarIns');
            app.SelectedListBox.Items = pipeline.PLItems;
            app.SelectedListBox.ItemsData = pipeline.PLItemsData;
            app.ChangedVal = pipeline.VarIns;
            app.UITable.Data = [];
        end

        % Button pushed function: SelectEEGLABFolderButton
        function SelectEEGLABFolderButtonPushed(app, event)
            try
                eeglabpath = uigetdir('*.*','Select EEGLAB Folder');
                app.PathEditField.Value=eeglabpath;
                addpath(eeglabpath)
            catch
                warning('Please select at least one file!')
                app.PathEditField.Value= '';

            end
        end

        % Button pushed function: SelectDataButton
        function SelectDataButtonPushed(app, event)
            try
                [app.file,app.path] = uigetfile( ...
                    {'*.set;*.vhdr;*.cdt;*.cnt',...
                    'Data Files (*.set,*.vhdr,*.cdt,*.cnt)'; ...
                    '*.set','Set Files (*.set)'; ...
                    '*.vhdr','VHDR Files (*.vhdr)'; ...
                    '*.cdt','CDT Files (*.cdt)'; ...
                    '*.cnt','CNT Files (*.cnt)'; ...
                    '*.*',  'All Files (*.*)'}, ...
                    'Select File(s)','multiSelect','on');
                if iscell(app.file)
                    app.NSelecFiles = numel(app.file);
                else
                    app.NSelecFiles = 1;
                    app.file = {app.file};
                end

                assignin('base','files',app.file)
                assignin('base','paths',app.path)

                app.FileEditField.Value = app.file{1};
                app.FolderEditField.Value = app.path;
                app.outofEditField.Value = num2str(max(size(app.file)));

            catch
                warning('Please select at least one file!')
                app.FileEditField.Value = '';
                app.FolderEditField.Value = '';
            end
        end

        % Button pushed function: ReStartStepsButton
        function ReStartStepsButtonPushed(app, event)
            clc
            app.SelectedListBox.Items(:)=[];
            app.SelectedListBox.ItemsData(:)=[];
            app.UITable.Data = [];
            app.ItemNum = 0;
            app.nstep = 1;
            app.needchanloc = 1;
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
            app.ChangedVal{ind2} = app.DefaultsVal{ind1};
            app.UITable.Data = app.DefaultsVal{ind1};
            styleParamTable(app);
        end

        % Value changed function: outofEditField
        function outofEditFieldValueChanged(app, event)
            value = app.FolderEditField.Value;

        end

        % Value changed function: RunningstepEditField
        function RunningstepEditFieldValueChanged(app, event)
            app.RunningstepEditField.Value=Step;

        end

        % Cell edit callback: UITable
        function UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;

            value = app.SelectedListBox.Value;
            indNum2 = find(ismember(app.SelectedListBox.ItemsData,value));
            app.ChangedVal{indNum2} = app.UITable.Data;

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
                app.FolderEditField_2.Value = app.FolderEditField.Value;
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
                warning('Please plot TEP first!')
            else
                % Create a new figure and axes
                fig = figure;
                ax = axes('Parent', fig);

                % Copy all children (plots, lines, etc.) from app.UIAxes to the new axes
                copyobj(allchild(app.UIAxes), ax);

                % Copy axis labels, title, and limits
                ax.XLabel.String = app.UIAxes.XLabel.String;
                ax.YLabel.String = app.UIAxes.YLabel.String;
                ax.Title.String  = app.UIAxes.Title.String;

                % Optional: Match axis limits
                ax.XLim = app.UIAxes.XLim;
                ax.YLim = app.UIAxes.YLim;

                % Optional: Copy legend if exists
                lgd = findobj(app.UIAxes.Parent, 'Type', 'Legend');
                if ~isempty(lgd)
                    legend(ax, '');
                end
            end
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
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 867 529];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.AutoResizeChildren = 'off';
            app.TabGroup.Position = [1 1 867 529];

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

            % Create LoadPipelineButton
            app.LoadPipelineButton = uibutton(app.CleaningTab, 'push');
            app.LoadPipelineButton.ButtonPushedFcn = createCallbackFcn(app, @LoadPipelineButtonPushed, true);
            app.LoadPipelineButton.BackgroundColor = [0.8 0.8 0.8];
            app.LoadPipelineButton.Position = [377 12 66 36];
            app.LoadPipelineButton.Text = {'Load'; 'Pipeline'};

            % Create SavePipelineButton
            app.SavePipelineButton = uibutton(app.CleaningTab, 'push');
            app.SavePipelineButton.ButtonPushedFcn = createCallbackFcn(app, @SavePipelineButtonPushed, true);
            app.SavePipelineButton.BackgroundColor = [0.8 0.8 0.8];
            app.SavePipelineButton.Position = [378 56 66 36];
            app.SavePipelineButton.Text = {'Save'; 'Pipeline'};

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
            app.SelectDatatoPerformAnalysisPanel = uipanel(app.CleaningTab);
            app.SelectDatatoPerformAnalysisPanel.AutoResizeChildren = 'off';
            app.SelectDatatoPerformAnalysisPanel.BorderType = 'none';
            app.SelectDatatoPerformAnalysisPanel.Title = 'Select Data to Perform Analysis';
            app.SelectDatatoPerformAnalysisPanel.Position = [649 225 208 116];

            % Create FolderEditFieldLabel
            app.FolderEditFieldLabel = uilabel(app.SelectDatatoPerformAnalysisPanel);
            app.FolderEditFieldLabel.HorizontalAlignment = 'right';
            app.FolderEditFieldLabel.Position = [5 64 40 22];
            app.FolderEditFieldLabel.Text = 'Folder';

            % Create FolderEditField
            app.FolderEditField = uieditfield(app.SelectDatatoPerformAnalysisPanel, 'text');
            app.FolderEditField.Editable = 'off';
            app.FolderEditField.Position = [53 64 145 22];

            % Create FileEditFieldLabel
            app.FileEditFieldLabel = uilabel(app.SelectDatatoPerformAnalysisPanel);
            app.FileEditFieldLabel.HorizontalAlignment = 'right';
            app.FileEditFieldLabel.Position = [6 39 25 22];
            app.FileEditFieldLabel.Text = 'File';

            % Create FileEditField
            app.FileEditField = uieditfield(app.SelectDatatoPerformAnalysisPanel, 'text');
            app.FileEditField.Editable = 'off';
            app.FileEditField.Position = [53 38 145 22];

            % Create SelectDataButton
            app.SelectDataButton = uibutton(app.SelectDatatoPerformAnalysisPanel, 'push');
            app.SelectDataButton.ButtonPushedFcn = createCallbackFcn(app, @SelectDataButtonPushed, true);
            app.SelectDataButton.Position = [15 5 183 23];
            app.SelectDataButton.Text = 'Select Data';

            % Create RunAnalysisButton
            app.RunAnalysisButton = uibutton(app.CleaningTab, 'push');
            app.RunAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @RunAnalysisButtonPushed, true);
            app.RunAnalysisButton.BackgroundColor = [0.651 0.651 0.651];
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

            % Create EEGLABpathifalreadyisnotinpathPanel
            app.EEGLABpathifalreadyisnotinpathPanel = uipanel(app.CleaningTab);
            app.EEGLABpathifalreadyisnotinpathPanel.AutoResizeChildren = 'off';
            app.EEGLABpathifalreadyisnotinpathPanel.BorderType = 'none';
            app.EEGLABpathifalreadyisnotinpathPanel.Title = 'EEGLAB path (if already is not in path)';
            app.EEGLABpathifalreadyisnotinpathPanel.Position = [651 342 208 90];

            % Create PathEditFieldLabel
            app.PathEditFieldLabel = uilabel(app.EEGLABpathifalreadyisnotinpathPanel);
            app.PathEditFieldLabel.HorizontalAlignment = 'right';
            app.PathEditFieldLabel.Position = [13 40 30 22];
            app.PathEditFieldLabel.Text = 'Path';

            % Create PathEditField
            app.PathEditField = uieditfield(app.EEGLABpathifalreadyisnotinpathPanel, 'text');
            app.PathEditField.Editable = 'off';
            app.PathEditField.Position = [50 40 145 22];

            % Create SelectEEGLABFolderButton
            app.SelectEEGLABFolderButton = uibutton(app.EEGLABpathifalreadyisnotinpathPanel, 'push');
            app.SelectEEGLABFolderButton.ButtonPushedFcn = createCallbackFcn(app, @SelectEEGLABFolderButtonPushed, true);
            app.SelectEEGLABFolderButton.Position = [15 10 183 23];
            app.SelectEEGLABFolderButton.Text = 'Select EEGLAB Folder';

            % Create ReStartStepsButton
            app.ReStartStepsButton = uibutton(app.CleaningTab, 'push');
            app.ReStartStepsButton.ButtonPushedFcn = createCallbackFcn(app, @ReStartStepsButtonPushed, true);
            app.ReStartStepsButton.BackgroundColor = [0.651 0.651 0.651];
            app.ReStartStepsButton.FontSize = 18;
            app.ReStartStepsButton.FontWeight = 'bold';
            app.ReStartStepsButton.Position = [658 91 201 36];
            app.ReStartStepsButton.Text = 'ReStart Steps';

            % Create outofEditFieldLabel
            app.outofEditFieldLabel = uilabel(app.CleaningTab);
            app.outofEditFieldLabel.HorizontalAlignment = 'center';
            app.outofEditFieldLabel.FontSize = 11;
            app.outofEditFieldLabel.Position = [772 191 36 27];
            app.outofEditFieldLabel.Text = 'out of';

            % Create outofEditField
            app.outofEditField = uieditfield(app.CleaningTab, 'text');
            app.outofEditField.ValueChangedFcn = createCallbackFcn(app, @outofEditFieldValueChanged, true);
            app.outofEditField.FontSize = 11;
            app.outofEditField.Position = [808 193 32 22];

            % Create ProcessingfileEditFieldLabel
            app.ProcessingfileEditFieldLabel = uilabel(app.CleaningTab);
            app.ProcessingfileEditFieldLabel.HorizontalAlignment = 'center';
            app.ProcessingfileEditFieldLabel.FontSize = 11;
            app.ProcessingfileEditFieldLabel.Position = [661 191 80 27];
            app.ProcessingfileEditFieldLabel.Text = 'Processing file ';

            % Create ProcessingfileEditField
            app.ProcessingfileEditField = uieditfield(app.CleaningTab, 'text');
            app.ProcessingfileEditField.Editable = 'off';
            app.ProcessingfileEditField.Position = [739 193 32 22];

            % Create RunningstepEditFieldLabel
            app.RunningstepEditFieldLabel = uilabel(app.CleaningTab);
            app.RunningstepEditFieldLabel.HorizontalAlignment = 'right';
            app.RunningstepEditFieldLabel.Position = [660 163 79 22];
            app.RunningstepEditFieldLabel.Text = 'Running step:';

            % Create RunningstepEditField
            app.RunningstepEditField = uieditfield(app.CleaningTab, 'text');
            app.RunningstepEditField.ValueChangedFcn = createCallbackFcn(app, @RunningstepEditFieldValueChanged, true);
            app.RunningstepEditField.Position = [664 140 181 22];

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
            app.UIAxes.Position = [344 299 300 200];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.VisualizingTab);
            app.UIAxes2.TickLabelInterpreter = 'none';
            app.UIAxes2.XAxisLocation = 'origin';
            app.UIAxes2.XTick = [];
            app.UIAxes2.YAxisLocation = 'origin';
            app.UIAxes2.YTick = [];
            app.UIAxes2.ZTick = [];
            app.UIAxes2.Position = [359 87 176 168];

            % Create PLOTTEPButton
            app.PLOTTEPButton = uibutton(app.VisualizingTab, 'push');
            app.PLOTTEPButton.ButtonPushedFcn = createCallbackFcn(app, @PLOTTEPButtonPushed, true);
            app.PLOTTEPButton.Enable = 'off';
            app.PLOTTEPButton.Position = [21 113 108 23];
            app.PLOTTEPButton.Text = 'PLOT TEP';

            % Create SelectDatatoVisulaizeTEPsPanel
            app.SelectDatatoVisulaizeTEPsPanel = uipanel(app.VisualizingTab);
            app.SelectDatatoVisulaizeTEPsPanel.AutoResizeChildren = 'off';
            app.SelectDatatoVisulaizeTEPsPanel.BorderType = 'none';
            app.SelectDatatoVisulaizeTEPsPanel.Title = 'Select Data to Visulaize TEPs';
            app.SelectDatatoVisulaizeTEPsPanel.Position = [651 342 208 90];

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
            app.UseCurrentlyCleanedDataCheckBox.Position = [671 455 180 22];

            % Create FilesListBoxLabel
            app.FilesListBoxLabel = uilabel(app.VisualizingTab);
            app.FilesListBoxLabel.HorizontalAlignment = 'right';
            app.FilesListBoxLabel.Position = [740 325 30 22];
            app.FilesListBoxLabel.Text = 'Files';

            % Create FilesListBox
            app.FilesListBox = uilistbox(app.VisualizingTab);
            app.FilesListBox.Items = {};
            app.FilesListBox.Multiselect = 'on';
            app.FilesListBox.ValueChangedFcn = createCallbackFcn(app, @FilesListBoxValueChanged, true);
            app.FilesListBox.Position = [669 71 183 259];
            app.FilesListBox.Value = {};

            % Create Image2
            app.Image2 = uiimage(app.VisualizingTab);
            app.Image2.Position = [-1 165 350 336];
            app.Image2.ImageSource = fullfile(pathToMLAPP, 'Head.png');

            % Create Slider
            app.Slider = uislider(app.VisualizingTab);
            app.Slider.ValueChangingFcn = createCallbackFcn(app, @SliderValueChanging, true);
            app.Slider.Position = [371 38 254 3];
            app.Slider.Value = 60;

            % Create WindowsizeforTopoplotLabel
            app.WindowsizeforTopoplotLabel = uilabel(app.VisualizingTab);
            app.WindowsizeforTopoplotLabel.HorizontalAlignment = 'center';
            app.WindowsizeforTopoplotLabel.Position = [353 44 90 44];
            app.WindowsizeforTopoplotLabel.Text = {'Window size for'; ' time averaged'; ' Topoplot'};

            % Create WindowsizefortimeaveragedTopoplotEditField
            app.WindowsizefortimeaveragedTopoplotEditField = uieditfield(app.VisualizingTab, 'numeric');
            app.WindowsizefortimeaveragedTopoplotEditField.ValueDisplayFormat = '%.0f';
            app.WindowsizefortimeaveragedTopoplotEditField.Position = [453 58 41 22];

            % Create TOPOPLOTButton
            app.TOPOPLOTButton = uibutton(app.VisualizingTab, 'push');
            app.TOPOPLOTButton.ButtonPushedFcn = createCallbackFcn(app, @TOPOPLOTButtonPushed, true);
            app.TOPOPLOTButton.Enable = 'off';
            app.TOPOPLOTButton.Position = [541 150 100 42];
            app.TOPOPLOTButton.Text = 'TOPOPLOT';

            % Create ExportTEPFigureButton
            app.ExportTEPFigureButton = uibutton(app.VisualizingTab, 'push');
            app.ExportTEPFigureButton.ButtonPushedFcn = createCallbackFcn(app, @ExportTEPFigureButtonPushed, true);
            app.ExportTEPFigureButton.Enable = 'off';
            app.ExportTEPFigureButton.Position = [21 84 108 23];
            app.ExportTEPFigureButton.Text = 'Export TEP Figure';

            % Create PlottingModeButtonGroup
            app.PlottingModeButtonGroup = uibuttongroup(app.VisualizingTab);
            app.PlottingModeButtonGroup.AutoResizeChildren = 'off';
            app.PlottingModeButtonGroup.BorderType = 'none';
            app.PlottingModeButtonGroup.Title = 'Plotting Mode';
            app.PlottingModeButtonGroup.Position = [174 88 149 67];

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
            app.ReLoadAvailableElectrodesButton.Position = [686 7 153 23];
            app.ReLoadAvailableElectrodesButton.Text = 'Re/Load Available Electrodes';

            % Create TEPWindowSliderLabel
            app.TEPWindowSliderLabel = uilabel(app.VisualizingTab);
            app.TEPWindowSliderLabel.HorizontalAlignment = 'right';
            app.TEPWindowSliderLabel.Position = [350 266 74 22];
            app.TEPWindowSliderLabel.Text = 'TEP Window';

            % Create TEPWindowSlider
            app.TEPWindowSlider = uislider(app.VisualizingTab, 'range');
            app.TEPWindowSlider.Limits = [-100 300];
            app.TEPWindowSlider.ValueChangingFcn = createCallbackFcn(app, @TEPWindowSliderValueChanging, true);
            app.TEPWindowSlider.Position = [441 285 193 3];
            app.TEPWindowSlider.Value = [-50 250];

            % Create TopoplottimeSpinnerLabel
            app.TopoplottimeSpinnerLabel = uilabel(app.VisualizingTab);
            app.TopoplottimeSpinnerLabel.HorizontalAlignment = 'right';
            app.TopoplottimeSpinnerLabel.Position = [496 57 76 22];
            app.TopoplottimeSpinnerLabel.Text = 'Topoplot time';

            % Create TopoplottimeSpinner
            app.TopoplottimeSpinner = uispinner(app.VisualizingTab);
            app.TopoplottimeSpinner.RoundFractionalValues = 'on';
            app.TopoplottimeSpinner.ValueDisplayFormat = '%.0f';
            app.TopoplottimeSpinner.ValueChangedFcn = createCallbackFcn(app, @TopoplottimeSpinnerValueChanged, true);
            app.TopoplottimeSpinner.Position = [574 57 66 22];
            app.TopoplottimeSpinner.Value = 60;

            % Create EEGDatasetDropDownLabel
            app.EEGDatasetDropDownLabel = uilabel(app.VisualizingTab);
            app.EEGDatasetDropDownLabel.HorizontalAlignment = 'right';
            app.EEGDatasetDropDownLabel.Position = [132 21 75 22];
            app.EEGDatasetDropDownLabel.Text = 'EEG Dataset';

            % Create EEGDatasetDropDown
            app.EEGDatasetDropDown = uidropdown(app.VisualizingTab);
            app.EEGDatasetDropDown.Items = {'Select a file'};
            app.EEGDatasetDropDown.ValueChangedFcn = createCallbackFcn(app, @EEGDatasetDropDownValueChanged, true);
            app.EEGDatasetDropDown.Enable = 'off';
            app.EEGDatasetDropDown.Position = [214 21 110 22];
            app.EEGDatasetDropDown.Value = 'Select a file';

            % Create PlotEEGdataButton
            app.PlotEEGdataButton = uibutton(app.VisualizingTab, 'push');
            app.PlotEEGdataButton.ButtonPushedFcn = createCallbackFcn(app, @PlotEEGdataButtonPushed, true);
            app.PlotEEGdataButton.Enable = 'off';
            app.PlotEEGdataButton.Position = [21 21 108 23];
            app.PlotEEGdataButton.Text = 'Plot EEG data';

            % Create ExportTEPDataButton
            app.ExportTEPDataButton = uibutton(app.VisualizingTab, 'push');
            app.ExportTEPDataButton.ButtonPushedFcn = createCallbackFcn(app, @ExportTEPDataButtonPushed, true);
            app.ExportTEPDataButton.Enable = 'off';
            app.ExportTEPDataButton.Position = [21 56 108 23];
            app.ExportTEPDataButton.Text = 'Export TEP Data';

            % Create TEPvarNameEditFieldLabel
            app.TEPvarNameEditFieldLabel = uilabel(app.VisualizingTab);
            app.TEPvarNameEditFieldLabel.HorizontalAlignment = 'right';
            app.TEPvarNameEditFieldLabel.Enable = 'off';
            app.TEPvarNameEditFieldLabel.Position = [128 56 83 22];
            app.TEPvarNameEditFieldLabel.Text = 'TEP var Name';

            % Create TEPvarNameEditField
            app.TEPvarNameEditField = uieditfield(app.VisualizingTab, 'text');
            app.TEPvarNameEditField.ValueChangedFcn = createCallbackFcn(app, @TEPvarNameEditFieldValueChanged, true);
            app.TEPvarNameEditField.Enable = 'off';
            app.TEPvarNameEditField.Position = [212 56 111 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
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