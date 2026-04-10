function runPipeline(app)
% RUNPIPELINE  Execute the selected cleaning pipeline on all chosen files.
%   runPipeline(app) runs each step in app.SelectedListBox against every
%   file in app.file. Called from nestapp.RunAnalysisButtonPushed.
%   app is a handle object; all UI state is read and written through it.
%
% Copyright (C) 2023  Aref Pariz, University of Ottawa & The Royal
% Institute for Mental Health, Ottawa, Ontario, Canada.
% apariz@uottawa.ca
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
%
% STEPS TO PERFORM ON EEG DATA
%
% This file contains the commands and functions already available in eeglab
% package. This function is called by the app "nestapp".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION
clc
global EEG ALLEEG CURRENTSET ALLCOM
app.steps2run=cell(1,2 * numel(app.SelectedListBox.Items));
if isempty(app.file)
    error('--------------Please select at least one data!----------------')
else
    for i = 1:numel(app.SelectedListBox.Items)
        app.steps2run{2*(i-1)+1} = app.SelectedListBox.Items(i);
        app.steps2run{2*i} = app.ChangedVal(i);
    end

    for i = 2:2:size(app.steps2run,2)
        x = table2cell(app.steps2run{i}{:});
        inputvals = cell(1,2*size(x,1));
        for j = 1:size(x,1)
            inputvals{2*j-1} = x{j,1};
            v = x{j,2};
            if ischar(v) || isstring(v)
                sv = string(v);
                if isscalar(sv)
                    % Scalar string: check for placeholder or numeric restoration.
                    if strlength(sv) > 0 && startsWith(sv, '(')
                        % Display placeholder — treat as empty for EEGLAB.
                        v = '[]';
                    else
                        % Restore mat2str strings back to numeric
                        % (e.g. '[1 6]' -> [1 6], '250' -> 250).
                        num = str2num(char(sv)); %#ok<ST2NM>
                        if ~isempty(num)
                            v = num;
                        end
                    end
                end
                % Non-scalar string arrays (e.g. ["TP9" "TP10"]) pass through unchanged.
            end
            inputvals{2*j} = v;
        end
        app.steps2run{i} = [];
        app.steps2run{i} = inputvals;
    end
    % Below will export the Steps and their name to MATLAB workspace
    assignin('base','steps2run',app.steps2run);
    assignin('base','stepsName',app.SelectedListBox.Items);

    app.nstep = 1; % The Steps starting point
    chName = []; % No File Selected for Channel Location
    dstep = 2; % DO NOT CHANGE THIS.
end

%% MAIN

nFiles = app.NSelecFiles;
nSteps = numel(app.steps2run) / 2;
allSummaries = {};

dlg = uiprogressdlg(app.UIFigure, ...
    'Title',           'Running Pipeline', ...
    'Message',         'Initialising...', ...
    'Cancelable',      'on', ...
    'ShowPercentage',  'on');

for nfile = 1:nFiles
    app.ProcessingfileEditField.Value = num2str(nfile);
    % To avoid any unforseen error, for each data new eeglab window will be
    % used.
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');
    app.initialVars = who; % Variables available at the begining of the analysis. Will be used to save them.
    pathName=app.path; % Path to data folder
    fileName = app.file{nfile}; % Data name(s) to be analyzed.

    disp(['!--------FILE ',fileName,' IS BEING PROCESSED--------!'])

    % Per-file step log — accumulated across the step loop and written to disk.
    stepLog = struct('step',{},'duration_s',{},'chanBefore',{},'chanAfter',{}, ...
                     'epochBefore',{},'epochAfter',{},'error',{});

    % Per-file pipeline report (channels, trials, ICA, TEP metrics).
    fileReport = initPipelineReport(fullfile(pathName, fileName));
    % ICA stats captured before pop_subcomp/pop_tesa_compselect (indices invalid after removal).
    pendingICAStats = struct();

    % In below loop, all assigned steps will be evaluated.
    for Step=app.nstep:dstep:numel(app.steps2run)
        stepName = app.steps2run{Step}{:};
        stepIdx  = (Step - 1) / 2 + 1;

        % Update progress dialog and check for user cancellation.
        dlg.Value   = ((nfile - 1) * nSteps + stepIdx - 1) / (nFiles * nSteps);
        dlg.Message = sprintf('File %d / %d  \x2014  %s', nfile, nFiles, stepName);
        if dlg.CancelRequested
            writeSessionLog(pathName, fileName, stepLog);
            close(dlg);
            return
        end

        app.RunningstepEditField.Value = stepName; pause(0.1);
        varin = app.steps2run{Step+1};
        disp(strcat('step ',num2str(stepIdx), ': "',stepName,'" is running!'));

        % Capture EEG state before the step runs.
        % EEG is [] (double) until Load Data runs — guard against that.
        if isstruct(EEG) && ~isempty(EEG)
            nChanBefore  = EEG.nbchan;
            nEpochBefore = size(EEG.data, 3);
        else
            nChanBefore  = 0;
            nEpochBefore = 0;
        end
        t0 = tic;

        try
            switch app.steps2run{Step}{:}
                case 'Load Channel Location'
                    %% Loading Channel Locations

                    vars = convertContainedStringsToChars(varin);
                    ind1 = find(strcmpi(vars,'eachFilediffPath'));
                    eachFilediffPath = vars{ind1+1};
                    ind2 = find(strcmpi(vars,'needchanloc'));
                    needchanloc = vars{ind2+1};

                    if strcmp(eachFilediffPath,'yes')
                        needchanloc='yes';
                        chName = [];
                    end
                    if strcmp(needchanloc,'yes')
                        if isempty(chName)
                            pathEEGLAB = which('eeglab');
                            if isunix
                                pathEEGLAB = replace(pathEEGLAB,'\','/');
                                pathEEGLAB = replace(pathEEGLAB,'eeglab.m','');
                                D = dir([pathEEGLAB,'plugins/dipfit*']);
                                lookforchnlocs = [D.folder,'/',D.name,'/standard_BEM/elec/standard_1005.elc'];
                            elseif ispc
                                pathEEGLAB = replace(pathEEGLAB,'/','\');
                                pathEEGLAB = replace(pathEEGLAB,'eeglab.m','');
                                D = dir([pathEEGLAB,'plugins\dipfit*']);
                                lookforchnlocs = [D.folder,'\',D.name,'\standard_BEM\elec\standard_1005.elc'];
                            end
                            [chName,chPath] = uigetfile('*.*','Select a file');
                        end
                        EEG = pop_chanedit(EEG, 'lookup', lookforchnlocs, ...
                            'load', {[chPath,chName],'filetype','autodetect'});
                        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
                    end

                case 'Load Data'
                    %% Loading the Files

                    mode = varin{1,2};
                    if   strcmpi(fileName(end-2:end),'set')
                        EEG = pop_loadset( [pathName fileName]);
                        fileFormat = 'set';
                    elseif strcmpi(fileName(end-2:end),'cnt')
                        fileFormat = 'cnt';
                        EEG = pop_loadcnt([pathName fileName] , 'dataformat', 'int32' );
                    elseif strcmpi(fileName(end-2:end),'cdt')
                        fileFormat = 'cdt';
                        EEG = loadcurry([pathName fileName], 'CurryLocations', 'False');
                    elseif strcmpi(fileName(end-3:end),'vhdr')
                        fileFormat = 'vhdr';
                        EEG  = pop_loadbv(pathName , fileName );
                    end
                    EEG.filename=fileName;
                    % Set the mode as 'on' to redraw loaded file to eeglab
                    if strcmpi(mode,'on')
                        postVars = who;
                        assignin('base','postVars',postVars)
                        for i = 1:numel(postVars)
                            if find(ismember(postVars{i}, app.initialVars))
                                assignin('base',postVars{i},eval(postVars{i}))
                            end
                        end
                        % eeglab redraw
                    end
                    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

                case 'Save New Set'
                    %% newSet

                    EEG = eeg_checkset( EEG );
                    vars = convertContainedStringsToChars(varin);
                    inds = find(strcmpi(vars,'includeFileName'));
                    IFN = vars{inds + 1};
                    vars([inds, inds+1]) =[];
                    fname='';
                    if strcmp(IFN,'yes')
                        fname = [app.path,fileName];
                        fname = fname(1:end-numel(fileFormat)-1);
                        fname = replace(fname,' ','_');
                        fname = replace(fname,'-','_');
                        fname = strcat(fname,'_');
                    end
                    ind1 = find(strcmp(vars,'savenew'));
                    if ~strcmp(vars{ind1+1},'[]') %|| ~isempty(vars{ind+1})
                        vars{ind1+1} = [fname,vars{ind1+1}];
                    end
                    ind2 = find(strcmp(vars,'saveold'));
                    if ~strcmp(vars{ind2+1},'[]') %|| ~isempty(vars{ind+1})
                        vars{ind2+1} = [fname,vars{ind2+1}];
                    end
                    inds = find(strcmp(vars,'[]'));
                    vars([inds,inds-1])=[];
                    
                    ind3 = find(strcmp(vars,'retrieve'));
                    if isnan(vars{ind3+1})
                        vars([ind3, ind3+1])=[];
                    end
                    EEG = eeg_checkset(EEG);
                    EEG.pipeline=app.steps2run;
                    % CURRENTSET = CURRENTSET + 1;
                    % assignin('base','EEG',EEG)
                    % assignin('base','ALLEEG',ALLEEG)
                    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,vars{:});
                    
                    % eeglab redraw
                    
                case 'Manual Command'
                    %% Manual Command
                    vars = convertContainedStringsToChars(varin);
                    ind=find(strcmpi(vars,'command'));
                    if size(vars{ind+1},1)>1
                        for nstep=1:numel(vars{ind+1})
                            eval(vars{ind+1}{nstep});
                        end
                    else
                        eval(vars{ind+1});
                    end
                case 'Choose Data Set'
                    vars = convertContainedStringsToChars(varin);
                    ind = find(strcmp(vars,'dataSetInd'));
                    setIndex = vars{ind+1};
                    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',setIndex);
                case 'Visualize EEG Data'
                    %% Visulaising the EEG data
                    if sum(isnan(varin{1,2}))
                        pop_eegplot(EEG);
                    else
                        pop_eegplot(EEG, varin{1,2}(1),varin{1,2}(2),varin{1,2}(3));
                    end
                    uiconfirm(app.UIFigure,'Press OK when done viewing the EEG plot.','Visualize EEG','Options',{'OK'},'DefaultOption',1);
                case 'Remove un-needed Channels'
                    %% Remove un-needed channels

                    %select here the channel you want to remove if they have a bad impedence
                    vars = convertContainedStringsToChars(varin);
                    inds = find(strcmp(vars,'[]'));
                    vars([inds,inds-1])=[];
                    assignin('base',"vars",vars)
                    EEG = pop_select( EEG,vars{:});
                    EEG = eeg_checkset( EEG );
                
                case 'Remove Bad Channels'
                    %%  Remove bad channels

                    vars = convertContainedStringsToChars(varin);
                    EEGelecNames = {EEG.chanlocs(1:end).labels};
                    ind1 = find(strcmpi(vars,'impelec'));
                    AuximportantElects = vars{ind1+1};
                    importantElects=matches(EEGelecNames, AuximportantElects,"IgnoreCase",true);
                    vars([ind1, ind1+1]) = [];

                    ind2 = find(strcmpi(vars,'elec'));
                    if strcmp(vars{1,ind2+1},'[]')
                        vars{1,ind2+1}= 1:EEG.nbchan;
                    elseif iscell(vars{1,ind2+1})
                        vars{1,ind2+1}=find(ismember(EEGelecNames, vars{1,ind2+1}));
                    else
                        vars{1,ind2+1} = 1:EEG.nbchan;
                    end

                    ind3 = find(strcmpi(vars,'freqrange'));
                    if sum(isnan(vars{ind3+1})) || strcmpi(vars{ind3},'[]')
                        vars([ind3, ind3+1]) = [];
                    end

                    if sum(importantElects)
                        vars{1,ind2+1} = find(~importantElects);
                        EEG = pop_rejchan(EEG, vars{:});
                    else
                        EEG = pop_rejchan(EEG, vars{:});
                    end
                case 'Automatic Continuous Rejection'
                    vars = convertContainedStringsToChars(varin);
                    ind=find(strcmpi(vars,'elecrange'));
                    if max(vars{ind+1})>EEG.nbchan
                        elecrange = 1:EEG.nbchan;
                    else
                        elecrange = vars{ind+1}(1):vars{ind+1}(end);
                    end
                    vars([ind,ind+1])=[];
                    EEG = pop_rejcont(EEG,'elecrange',elecrange,vars{:});
                case 'Clean Artifacts'
                    %% Clean Artifacts
                    vars = convertContainedStringsToChars(varin);
                    ind1=find(strcmpi(vars,'Channels'));
                    ind2=find(strcmpi(vars,'Channels_ignore'));
                    if sum(strcmpi(vars{ind1+1},'[]')) && sum(strcmpi(vars{ind2+1},'[]'))
                        chans={EEG.chanlocs.labels};
                        vars{ind1+1}=chans;
                        vars{ind2+1}=[];
                    elseif sum(strcmpi(vars{ind1+1},'[]')) && sum(~strcmpi(vars{ind2+1},'[]'))
                        if size(vars{ind2+1},1)>size(vars{ind2},2)
                            vars{ind2+1}=vars{ind2+1}';
                        end
                        vars([ind1 ind1+1])=[];
                    elseif sum(~strcmpi(vars{ind1+1},'[]')) && sum(strcmpi(vars{ind2+1},'[]'))
                        if size(vars{ind1+1},1)>size(vars{ind1},2)
                            vars{ind1+1}=vars{ind1+1}';
                        end
                        vars([ind2 ind2+1])=[];
                    end
                    EEG = clean_artifacts(EEG,vars{:});
                    EEG = eeg_checkset( EEG );

               case 'Automatic Cleaning Data'
                    %% Automatic Cleaning Raw Data

                    vars = convertContainedStringsToChars(varin);
                    ind = find(strcmp(vars,'Highpass'));
                    if ~strcmpi(vars{ind+1},'off')
                        highpass = [str2double(vars{ind+1})];
                        if size(highpass,2)<size(highpass,1)
                            highpass = highpass';
                        end
                        vars{ind+1} = highpass;
                    end
                    ind = find(strcmp(vars,'[]'));
                    vars([ind, ind-1])=[];
                    EEG = pop_clean_rawdata(EEG, vars{:});
                    EEG = eeg_checkset( EEG );
                case 'Remove Baseline'
                    %% Remove Baseline Offset

                    % The DC offset may cause problem in analyzing and cleaning data. So, it's
                    % better to remove it before other steps. However if we want to study the
                    % ERP, it is recomended to do not remove it, since it may affect the
                    % results.
                    vars = convertContainedStringsToChars(varin);
                    ind = find(strcmp(vars,'[]'), 1);
                    if ~isempty(ind)
                        timerange=vars{1,2}; % Default [] -> all
                        if ~strcmp(timerange,'[]')
                            timerange = str2double(timerange);
                            EEG = pop_rmbase( EEG, timerange);
                        end
                    else
                        timerange=vars{1,2}; % Default [] -> all
                        pointrange=vars{1,4}; % Default [] -> all
                        chanlist=vars{1,6}; % Default [] -> all
                        pointrange = str2double(pointrange);
                        EEG = pop_rmbase( EEG, timerange, pointrange, chanlist);
                    end
                    EEG = eeg_checkset( EEG );

                case 'De-Trend Epoch'
                    %%  DeTrending Epochs

                    % Built in dtrend function
                    % The default value is to use npoly = 1 but based on the stacked data, some
                    % channels show wronf trend and I used npoly 2.
                    vars = convertContainedStringsToChars(varin);
                    for elecc=1:size(EEG.data,1)
                        for epoo=1:size(EEG.data,3)
                            EEG.data(elecc,:,epoo)=detrend(EEG.data(elecc,:,epoo),vars{1,2});
                        end
                    end
                    EEG = eeg_checkset( EEG );
                case 'TESA De-Trend'
                    vars = convertContainedStringsToChars(varin);
                    ind1 = find(strcmpi(vars,'detrend'));
                    Tdetrend = vars{ind1+1};
                    ind2 = find(strcmpi(vars,'timeWin'));
                    TtimeWin = vars{ind2+1};
                    pop_tesa_detrend(EEG, Tdetrend, TtimeWin)
                case 'Re-Sample'
                    %% Re-Sampleing
                    vars = convertContainedStringsToChars(varin);

                    EEG = pop_resample(EEG,vars{2:2:end});
                    EEG = eeg_checkset( EEG );

                case 'Re-Reference'
                    %% Re-referencing

                    vars = convertContainedStringsToChars(varin);
                    ind = find(strcmp(vars,'ref'));
                    ref = vars{ind+1};
                    if ~ismember(ref,{EEG.chanlocs.labels}) & ~strcmp(ref,'[]')
                        answer = inputdlg('The reference channel is not in the data. Enter a new reference channel label:','Re-Reference',[1 50],{''});
                        if isempty(answer) || isempty(answer{1})
                            error('Re-Reference cancelled: no reference channel provided.');
                        end
                        ref = answer{1};
                    end
                    if strcmp(ref,'[]')
                        ref =eval(ref);
                    end
                    vars([ind,ind+1])=[];

                    inds = find(strcmp(vars,'[]'));
                    vars([inds,inds-1])=[];
                    vars([ind,ind+1])=[];
                    EEG = pop_reref(EEG,ref, vars{:}); % Re-reference to average reference
                    EEG = eeg_checkset( EEG );

                case 'Frequency Filter (CleanLine)'
                    %%  Frequency Filtering

                    % Below line will remove the notch frequency using cleanline extension.
                    % Notch frequency, 60 Hz and its first harmonic, 120 Hz is being filtered
                    % from all channels, 1:63
                    vars = convertContainedStringsToChars(varin);
                    assignin('base','vars',vars);
                    ind = find(strcmp(vars,'chanlist'));

                    if vars{ind+1}(2) > EEG.nbchan
                        vars{ind+1} = 1:EEG.nbchan-1;
                    else
                        vars{ind+1} = 1:vars{1,ind+1};
                    end

                    EEG = pop_cleanline(EEG, vars{:});
                    % EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:63] ,'computepower',1,'linefreqs',60,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
                    EEG = eeg_checkset( EEG );

                case 'Frequency Filter (TESA)'
                    %% Butterworth filter TESA

                    vars = convertContainedStringsToChars(varin);
                    ind1=find(strcmp(vars,'high'));
                    high = vars{ind1+1};
                    ind2=find(strcmp(vars,'low'));
                    low = vars{ind2+1};
                    ind3=find(strcmp(vars,'ord'));
                    ord = vars{ind3+1};
                    ind4=find(strcmp(vars,'type'));
                    type = vars{ind4+1};
                    vars([ind1,ind1+1,ind2,ind2+2,ind3,ind3+1])=[]; %#ok<NASGU> intentional cleanup; named values already extracted above

                    EEG = pop_tesa_filtbutter( EEG, high, low, ord, type ); % Zero-phase, 4th-order band pass butterworth filter between 1-100 Hz.
                    EEG = eeg_checkset( EEG );
                case 'Frequency Filter'
                    vars = convertContainedStringsToChars(varin);
                    ind = find(strcmp(vars,'filtorder'));
                    if mod(vars{ind+1},2)~=0 || isstring(vars{ind}) 
                        error('The Filtorder should be an even number!')
                    elseif vars{ind+1}==0
                        vars([ind, ind+1])=[];
                    end
                    EEG = pop_eegfiltnew(EEG, vars{:});
                    
                
                

                case 'Remove Bad Epoch'
                    %% Remove bad Epoch
                    % EEG = pop_jointprob(EEG,1,1:EEG.nbchan ,5,5,0,0);
                    vars = convertContainedStringsToChars(varin);
                    inds = find(strcmp(vars,'[]'));
                    vars([inds,inds-1])=[];
                    Ind1 = find(strcmpi(vars,'startprob'));
                    if ~isempty(Ind1)
                        vars{Ind1+1} = str2double(vars{Ind1+1});
                    end
                    Ind2 = find(strcmpi(vars,'maxrej'));
                    if ~isempty(Ind2)
                        vars{Ind2+1} = str2double(vars{Ind2+1});
                    end
                    
                    [EEG, rejepochs] = pop_autorej(EEG, vars{:});
                    EEG.rejEpochs=rejepochs;
                    EEG = eeg_checkset( EEG );

                case 'Run ICA'
                    %%  Runing ICA
                    EEG.data=double(EEG.data);
                    vars = convertContainedStringsToChars(varin);
                    EEG = pop_runica(EEG,vars{:});
                    EEG = eeg_checkset( EEG );

                case 'Label ICA Components'
                    %% Label ICA Components
                    vars = convertContainedStringsToChars(varin);
                    version = vars{2};
                    EEG = pop_iclabel(EEG, version);
                    EEG = eeg_checkset( EEG );

                case 'Flag ICA Components for Rejection'
                    %% Flag ICA Components for Rejection
                    vars = convertContainedStringsToChars(varin);
                    threshold = zeros(7,2);
                    for nflag = 2:2:14
                        threshold(nflag/2,:)=vars{nflag};
                    end
                    EEG = pop_icflag(EEG, threshold);
                    EEG = eeg_checkset( EEG );

                case 'Remove Flagged ICA Components'
                    %% Remove flagged components
                    vars = convertContainedStringsToChars(varin);
                    if strcmp(vars{2},'[]')
                        var_comp=[];
                    end
                    if ~exist('ICA_Rejected_Comp','var')
                        Rej = EEG.reject.gcompreject;
                        ICA_Rejected_Comp{1}=reshape(Rej,1,numel(Rej));
                    else
                        Rej = EEG.reject.gcompreject;
                        Rej = reshape(Rej, 1, numel(Rej));
                        ICA_Rejected_Comp{end+1}=Rej; %#ok<AGROW> count of ICA rounds unknown at call time
                    end

                    % Capture rejection mask and per-component variance BEFORE pop_subcomp.
                    % After removal icaweights loses the rejected rows, making the full
                    % logical mask invalid for indexing component activations.
                    rejMask = logical(reshape(EEG.reject.gcompreject, 1, []));
                    pendingICAStats = struct('rejMask', rejMask);
                    if ~isempty(EEG.icaweights) && size(EEG.icaweights,1) == numel(rejMask)
                        act2D = computeICAActivation(EEG);
                        data2D = reshape(EEG.data(EEG.icachansind,:,:), numel(EEG.icachansind), []);
                        totalVar = sum(var(data2D, 0, 2));
                        if totalVar > 0
                            pendingICAStats.compVarPct = (var(act2D, 0, 2) / totalVar * 100)';
                        end
                    end
                    if isfield(EEG,'etc') && isfield(EEG.etc,'ic_classification') && ...
                            isfield(EEG.etc.ic_classification,'ICLabel')
                        pendingICAStats.iclabelProbs = EEG.etc.ic_classification.ICLabel;
                    end

                    EEG = pop_subcomp( EEG, var_comp, vars{4}, vars{6});
                    EEG.ICA_Rejected_Comp=ICA_Rejected_Comp;
                    EEG = eeg_checkset( EEG );

                case 'Interpolate Channels'
                    %% Interpolate channels
                    vars = convertContainedStringsToChars(varin);
                    ind1 = find(strcmp(vars,'method'));
                    method = vars{ind1+1};
                    ind2 = find(strcmp(vars,'trange'));
                    if strcmp(vars{ind2+1},'[]')
                        trange = [];
                    end
                    
                    EEG = pop_interp(EEG, EEG.chaninfo.removedchans(1:size(EEG.chaninfo.removedchans,2)), method, trange);
                    if ~exist('interpElecs','var')
                        interpElecs = EEG.chaninfo.removedchans;
                    else
                        interpElecs = [interpElecs;EEG.chaninfo.removedchans]; %#ok<AGROW> accumulates across interpolation rounds
                    end
                    EEG.interpElecs = interpElecs;
                    EEG.setname = [EEG.setname '_interp'];
                    EEG.filename=[EEG.setname '.set'];
                    EEG.datfile=[EEG.setname '.fdt'];
                    EEG = eeg_checkset( EEG );

                    % make sure to modify the initial chanlocs file to have the final
                    % electrodes that you want to include in your study

                case 'Find TMS Pulses (TESA)'
                    %% TESA Finding TMS puls locations
                    % Refractory period (i.e. the time it takes for the TMS pulse to recover)
                    % The rate of change for detecting the TMS artifact (in uV/ms). If too many non-TMS pulse artifacts are being incorrectly labeled, increase this number.
                    % Label for single TMS pulses
                    vars = convertContainedStringsToChars(varin);
                    ind = find(strcmp(vars,'elec'));
                    elec = vars{ind+1};
                    vars([ind,ind+1])=[];
                    if iscell(elec)
                        elec=elec{:};
                    end
                    EEG = pop_tesa_findpulse( EEG, elec, vars{:});
                    EEG = eeg_checkset( EEG );

                case 'Fix TMS Pulse (TESA)'
                    %% TESA Fixing TMS pulse latencies
                    % This function finds TMS pulses by detecting the large TMS artifacts present in already epoched data.
                    % This script is designed for instances when the recorded events do not correspond with when the TMS pulse was given.
                    % The script works by extracting a single channel and finding the time points in which the first derivatives exceed a certain threshold (defined by 'rate')
                    vars = convertContainedStringsToChars(varin);
                    ind1 = find(strcmp(vars,'elec'));
                    elec = vars{1,ind1+1};
                    ind2 = find(strcmp(vars,'epoch_len'));
                    epoch_len = vars{1,ind2+1};
                    ind3 = find(strcmp(vars,'type'));
                    type = vars{1,ind3+1};
                    vars([ind1,ind1+1,ind2,ind2+1,ind3,ind3+1])=[];

                    EEG = tesa_fixevent( EEG, elec, epoch_len, type, vars{:} );
                    EEG = eeg_checkset( EEG );

                case 'Remove TMS Artifacts (TESA)'
                    %% TESA Remove TMS pulse artifact
                    vars = convertContainedStringsToChars(varin);
                    cutTimesTMS = vars{1,find(strcmp(vars,'cutTimesTMS'))+1}; % Time period in ms to be replaced by 0
                    replaceTimes = vars{1,find(strcmp(vars,'replaceTimes'))+1}; %[-500, -100]; % if not empty, Period used to be replaced TMS pulse
                    cutEvent = vars{1,find(strcmp(vars,'cutEvent'))+1}; % If not empty, Replace with 0s around event 'TMS'
                    if ~iscell(cutEvent)
                        cutEvent = {cutEvent};
                    end
                    if strcmp(replaceTimes,'[]')
                        replaceTimes = eval(replaceTimes);
                    end
                    if size(cutTimesTMS,1)>size(cutTimesTMS,2)
                        cutTimesTMS = cutTimesTMS';
                    end
                    EEG = pop_tesa_removedata(EEG, cutTimesTMS, replaceTimes, cutEvent);
                    EEG = eeg_checkset( EEG );

                case 'Epoching'
                    %% Epoching
                    vars = convertContainedStringsToChars(varin);
                    ind1 = find(strcmp(vars,'types'));
                    type = vars{1,ind1+1};
                    ind2 = find(strcmp(vars,'timelim'));
                    timelim = vars{1,ind2+1};
                    vars([ind1,ind2])=[];

                    % Using eeglab function to epoch the signal
                    EEG = pop_epoch( EEG, type, timelim, vars{:});
                    EEG = eeg_checkset( EEG );

                case 'Interpolate Missing Data (TESA)'
                    %% TESA Interpolate missing data LINEAR
                    % replaces missing data with linear interpolation.
                    % Linear function is fitted on data point before and after missing data.
                    vars = convertContainedStringsToChars(varin);
                    interpolation = vars{1,find(strcmp(vars,'interpolation'))+1};
                    interpWin = vars{1,find(strcmp(vars,'interpWin'))+1};
                    if size(interpWin,1) > size(interpWin,2)
                        interpWin = interpWin';
                    end
                    EEG = pop_tesa_interpdata( EEG, interpolation,interpWin);
                    % EEG = eeg_checkset( EEG );

                case 'Run TESA ICA'
                    %% TESA fastICA
                    % Uses the gauss contrast function and turns on the stabilized
                    % FastICA version to aid with convergence. g -> 'tanh' or 'gauss' or 'pow3' or 'skew'
                    vars = convertContainedStringsToChars(varin);
                    EEG = pop_tesa_fastica( EEG, vars{:} );
                    EEG = eeg_checkset( EEG );

                case 'Remove ICA Components (TESA)'
                    %% TESA Is Removing Components
                    % Turn off electrode noise detection, change threshold for
                    % blinks to 3, change electrodes used to AF3 and AF4 and turn
                    % on the feedback of blink threhsolds for individual components
                    % in the command window.
                    vars = convertContainedStringsToChars(varin);
                    ind1 = find(strcmpi(vars,'plotTimeX'));
                    TP = vars{ind1+1};
                    if TP(1) ~= EEG.times(1) && TP(2) ~= EEG.times(end)
                        vars{ind1+1} = [EEG.times(1) EEG.times(end)];
                    end
                    for nInd=2:2:numel(vars)
                        if size(vars{nInd},1)>size(vars{nInd},2)
                            vars{nInd}=vars{nInd}';
                        end
                    end
                    % Capture component count and variance BEFORE pop_tesa_compselect.
                    % TESA's compselect doesn't expose a rejection mask beforehand,
                    % so we can only attribute count change, not per-component stats.
                    nCompBeforeTESA = size(EEG.icaweights, 1);
                    pendingICAStats = struct('nCompBefore', nCompBeforeTESA);
                    if ~isempty(EEG.icaweights)
                        act2D_tesa = computeICAActivation(EEG);
                        data2D_tesa = reshape(EEG.data(EEG.icachansind,:,:), numel(EEG.icachansind), []);
                        totalVar_tesa = sum(var(data2D_tesa, 0, 2));
                        if totalVar_tesa > 0
                            pendingICAStats.compVarPct = (var(act2D_tesa, 0, 2) / totalVar_tesa * 100)';
                        end
                    end

                    EEG = pop_tesa_compselect( EEG,vars{:});
                    EEG = eeg_checkset( EEG );

                case 'Find Artifacts EDM (TESA)'
                    %% TESA Enhanced deflation method (EDM)
                    % This function finds artifactual components automatically by
                    % using the enhanced deflation method (EDM)

                    % change the threshold for artefact detection to 10, change the window for comparions
                    % to 11-50 ms and return threshold values for each component
                    % in the command window.
                    % nc = 30; % Number of component to be find. [] for all
                    % sr = 1000; % Sampling frequency in Hz. ([] if Sf is in EEG structure)
                    % chanl = []; % Channel locations. [] is cahnlocation is already in EEG structure
                    % cmp =10; % describes the number of components to perform selection on (e.g. first 10 components). Leave empty for all components.
                    % 'tmsMuscleThresh' the threshold for detecting components representing TMS-evoked muscle activity.
                    % 'tmsMuscleWin' [start,end] Vector describing the target window for TMS-evoked muscle activity (in ms).
                    % 'tmsMuscleFeedback' 'on' or 'off' turning on/off feedback of TMS-evoked muscle threshold value for each component in the command window.

                    vars = convertContainedStringsToChars(varin);
                    ind1 = find(strcmpi(vars,'chanl'));
                    chanlocations = vars{ind1+1};
                    vars([ind1, ind1+1]) = [];

                    ind2 = find(strcmpi(vars,'nc'));
                    nc = vars{ind2+1};
                    vars([ind2, ind2+1]) = [];

                    ind3 = find(strcmpi(vars,'sf'));
                    sf = vars{ind3+1};
                    vars([ind3, ind3+1]) = [];
                    if sf ~= EEG.srate
                        sf = EEG.srate;
                    end

                    EEG = pop_tesa_edm( EEG,chanlocations, nc, sf, vars{:});
                    EEG = eeg_checkset( EEG );

                case 'SSP SIR'
                    %% SSP SIR
                    % Both artefacts, TMS-evoked potentials (TEPs), and
                    % peripherally evoked potentials (PEPs) can be suppressed by
                    % SSP–SIR implemented in TESA.

                    %  NOTE: This command is better to run if the data has not been
                    %  highpass filltered below 200 Hz

                    % Suppresses control data by removing the first n
                    % principal components of controlResponse.
                    vars = convertContainedStringsToChars(varin);
                    EEG = pop_tesa_sspsir(EEG, vars{:});
                    EEG = eeg_checkset( EEG );

                case 'Remove Recording Noise (SOUND)'
                    %% Remove recording noise
                    % The SOUND algorithm automatically detect and remove recording noise.
                    vars = convertContainedStringsToChars(varin);
                    EEG = pop_tesa_sound(EEG, vars{:} ); %Run SOUND using customised input values
                    EEG = eeg_checkset( EEG );

                case 'Median Filter 1D'
                    %% 1-dimensional median filter of nth-order to remove artifacts such as spikes and muscle artifacts
                    % Vector with time range for applying median filter in ms.
                    % Note that t1 must be 0 or negative and t2 positive
                    vars = convertContainedStringsToChars(varin);
                    ind1 = find(strcmp(vars,'timeWin'));
                    timeWin=vars{1,ind1+1};
                    ind2 = find(strcmp(vars,'mdorder'));
                    mdorder=vars{1,ind2+1};
                    ind3 = find(strcmp(vars,'event_type'));
                    event_type=vars{1,ind3+1};
                    EEG = tesa_filtmedian( EEG, timeWin, mdorder, event_type );
                    EEG = eeg_checkset( EEG );

                case 'Remove Bad Trials'
                    %% Remove bad Trials
                    vars = convertContainedStringsToChars(varin);
                    indLoc = find(strcmpi(vars,'localThresh'));
                    indGlb = find(strcmpi(vars,'globalThresh'));
                    localThresh  = str2double(vars{indLoc+1});
                    globalThresh = str2double(vars{indGlb+1});
                    EEG = pop_jointprob(EEG, 1, 1:size(EEG.data,1), localThresh, globalThresh, 0, 0);
                    pop_rejmenu(EEG, 1);
                    uiconfirm(app.UIFigure,'Highlight bad trials in the rejection menu, then press OK to continue.','Remove Bad Trials','Options',{'OK'},'DefaultOption',1);
                    EEG.BadTr = unique([find(EEG.reject.rejjp==1) find(EEG.reject.rejmanual==1)]);
                    EEG = pop_rejepoch( EEG, EEG.BadTr ,0);
                    EEG = eeg_checkset( EEG );

                case 'Extract TEP (TESA)'
                    %% Extract TEP
                    vars = convertContainedStringsToChars(varin);
                    ind1 = find(strcmp(vars,'type'));
                    type=vars{1,ind1+1};
                    vars([ind1, ind1+1]) = [];
                    ind2 = find(strcmpi(vars,'pairCorrect'));
                    if ~strcmp(vars{ind2+1},'on')
                        ind3 = find(strcmpi(vars,'ISI'));
                        vars([ind3, ind3 +1]) = [];
                    end
                    inds = find(strcmp(vars,'[]'));
                    vars([inds,inds-1])=[];
                    EEG = pop_tesa_tepextract( EEG, type, vars );

                case 'Find TEP Peaks (TESA)'
                    %% Find TEP Peaks
                    vars = convertContainedStringsToChars(varin);
                    ind1 = find(strcmpi(vars,'input'));
                    input = vars{ind1+1};
                    vars([ind1, ind1+1]) = [];

                    ind2 = find(strcmpi(vars,'direction'));
                    direction = vars{ind2+1};
                    vars([ind2, ind2+1]) = [];

                    ind3 = find(strcmpi(vars,'peak'));
                    peak = vars{ind3+1};
                    vars([ind3, ind3+1]) = [];

                    ind4 = find(strcmpi(vars,'peakWin'));
                    peakWin = vars{ind4+1};
                    vars([ind4, ind4+1]) = [];

                    inds = find(strcmp(vars,'[]'));
                    vars([inds,inds-1])=[];

                    EEG = pop_tesa_peakanalysis( EEG, input, ...
                        direction, peak, peakWin, ...
                        vars(:) );
                case 'TEP Peak Output'
                    vars = convertContainedStringsToChars(varin);
                    inds = find(strcmp(vars,'[]'));
                    vars([inds,inds-1]) = [];
                    tepoutput = pop_tesa_peakoutput( EEG, vars{:} );
                    assignin('base','tep_output',tepoutput)

            end

            % Capture post-step channel/epoch counts (must precede report update).
            if isstruct(EEG) && ~isempty(EEG)
                nChanAfter  = EEG.nbchan;
                nEpochAfter = size(EEG.data, 3);
            else
                nChanAfter  = nChanBefore;
                nEpochAfter = nEpochBefore;
            end

            % Update fileReport with post-step metrics.
            if isstruct(EEG) && ~isempty(EEG)
                % Channel counts
                if strcmp(stepName, 'Load Data')
                    fileReport.channels.original = EEG.nbchan;
                end
                fileReport.channels.final = EEG.nbchan;
                if nChanAfter < nChanBefore
                    if any(strcmp(stepName, {'Remove Bad Channels','Remove un-needed Channels', ...
                            'Automatic Cleaning Data','Clean Artifacts', ...
                            'Automatic Continuous Rejection'}))
                        fileReport.channels.nRejected = fileReport.channels.nRejected + ...
                            (nChanBefore - nChanAfter);
                    end
                end
                if any(strcmp(stepName, {'Interpolate Channels','Interpolate Missing Data (TESA)'})) ...
                        && nChanAfter > nChanBefore
                    fileReport.channels.nInterpolated = fileReport.channels.nInterpolated + ...
                        (nChanAfter - nChanBefore);
                end

                % Trial/epoch counts
                if strcmp(stepName, 'Epoching') && fileReport.trials.original == 0
                    fileReport.trials.original = size(EEG.data, 3);
                end
                if size(EEG.data, 3) > 1
                    fileReport.trials.final = nEpochAfter;
                    if nEpochBefore > 1 && nEpochAfter < nEpochBefore
                        fileReport.trials.rejected = fileReport.trials.rejected + ...
                            (nEpochBefore - nEpochAfter);
                    end
                end

                % ICA — component count after decomposition
                if any(strcmp(stepName, {'Run ICA','Run TESA ICA'})) && ~isempty(EEG.icaweights)
                    fileReport.ica.nComponents = size(EEG.icaweights, 1);
                end

                % ICA — standard removal: use rejMask captured before pop_subcomp
                if strcmp(stepName, 'Remove Flagged ICA Components') && ...
                        isfield(pendingICAStats, 'rejMask')
                    rMask = pendingICAStats.rejMask;
                    nRej  = sum(rMask);
                    fileReport.ica.nRejected = fileReport.ica.nRejected + nRej;
                    fileReport.ica.nKept     = fileReport.ica.nComponents - fileReport.ica.nRejected;
                    if isfield(pendingICAStats, 'compVarPct') && ...
                            numel(pendingICAStats.compVarPct) == numel(rMask)
                        rejPct = pendingICAStats.compVarPct(rMask);
                        if isnan(fileReport.ica.varRemoved)
                            fileReport.ica.varRemoved = sum(rejPct);
                        else
                            fileReport.ica.varRemoved = fileReport.ica.varRemoved + sum(rejPct);
                        end
                        fileReport.ica.varMin = min([fileReport.ica.varMin, rejPct]);
                        fileReport.ica.varMax = max([fileReport.ica.varMax, rejPct]);
                    end
                    if isfield(pendingICAStats, 'iclabelProbs')
                        probs = pendingICAStats.iclabelProbs;
                        [~, bestCat] = max(probs, [], 2);
                        for ci = 1:7
                            inCat = (bestCat == ci) & rMask(:);
                            fileReport.ica.categories.nRemoved(ci) = ...
                                fileReport.ica.categories.nRemoved(ci) + sum(inCat);
                            if isfield(pendingICAStats, 'compVarPct')
                                fileReport.ica.categories.varShare(ci) = ...
                                    fileReport.ica.categories.varShare(ci) + ...
                                    sum(pendingICAStats.compVarPct(inCat));
                            end
                        end
                    end
                    pendingICAStats = struct();
                end

                % ICA — TESA removal: use component count delta (no per-component breakdown)
                if strcmp(stepName, 'Remove ICA Components (TESA)') && ...
                        isfield(pendingICAStats, 'nCompBefore')
                    nRejTESA = pendingICAStats.nCompBefore - size(EEG.icaweights, 1);
                    if nRejTESA > 0
                        fileReport.ica.nRejected = fileReport.ica.nRejected + nRejTESA;
                        fileReport.ica.nKept     = fileReport.ica.nComponents - fileReport.ica.nRejected;
                    end
                    pendingICAStats = struct();
                end
            end

            % Append step record to fileReport.
            stepRec.name         = stepName;
            stepRec.chansBefore  = nChanBefore;
            stepRec.chansAfter   = nChanAfter;
            stepRec.trialsBefore = nEpochBefore;
            stepRec.trialsAfter  = nEpochAfter;
            stepRec.duration     = toc(t0);
            stepRec.timestamp    = datetime('now');
            fileReport.steps{end+1} = stepRec;

            % Record successful step in session log.
            stepLog(end+1) = struct( ...
                'step',        stepName, ...
                'duration_s',  toc(t0), ...
                'chanBefore',  nChanBefore, ...
                'chanAfter',   nChanAfter, ...
                'epochBefore', nEpochBefore, ...
                'epochAfter',  nEpochAfter, ...
                'error',       ''); %#ok<AGROW>

        catch err
            stepLog(end+1) = struct( ...
                'step',        stepName, ...
                'duration_s',  toc(t0), ...
                'chanBefore',  nChanBefore, ...
                'chanAfter',   nChanBefore, ...
                'epochBefore', nEpochBefore, ...
                'epochAfter',  nEpochBefore, ...
                'error',       err.message); %#ok<AGROW>

            disp(err.message)
            warning(strcat('An error acoured at file ',fileName,...
                ' at step ',num2str(stepIdx), ': ',stepName));
            toContinue = uiconfirm(app.UIFigure, ...
                sprintf('Error at step %d (%s):\n%s\n\nContinue to next step?', ...
                    stepIdx, stepName, err.message), ...
                'Step Failed','Options',{'Continue','Abort'},'DefaultOption','Continue','CancelOption','Abort');
            if strcmp(toContinue,'Abort')
                postVars = who;
                assignin('base','postVars',postVars)
                for i = 1:numel(postVars)
                    if find(ismember(postVars{i}, app.initialVars))
                        assignin('base',postVars{i},eval(postVars{i}))
                    end
                end
                assignin('base','lastVarin',varin);
                assignin('base','lastStepInd',Step);
                writeSessionLog(pathName, fileName, stepLog);
                if isvalid(dlg); close(dlg); end
                return
            end
        end
    end

    writeSessionLog(pathName, fileName, stepLog);
    [summaryText, ~] = exportReport(fileReport, pathName);
    allSummaries{end+1} = summaryText; %#ok<AGROW>
    eeglab redraw
    disp('-----------------Data processed!-----------------')
end

if isvalid(dlg); close(dlg); end
showReportDialog(app, allSummaries);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local helpers

function writeSessionLog(pathName, fileName, stepLog)
% WRITESESSIONLOG  Write a plain-text processing log alongside the data file.
%   Saved as <fileName_noext>_nestapp_log.txt in pathName.
    [~, baseName] = fileparts(fileName);
    logPath = fullfile(pathName, [baseName, '_nestapp_log.txt']);
    fid = fopen(logPath, 'w');
    if fid == -1; return; end  % silently skip if path is not writable

    fprintf(fid, '=== nestapp session log ===\n');
    fprintf(fid, 'File:      %s\n', fileName);
    fprintf(fid, 'Processed: %s\n', datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    fprintf(fid, 'MATLAB:    %s\n', version);
    fprintf(fid, '\n%-4s  %-35s  %7s  %9s  %9s  %s\n', ...
        '#', 'Step', 'Time(s)', 'Ch before', 'Ch after', 'Note');
    fprintf(fid, '%s\n', repmat('-', 1, 80));

    for k = 1:numel(stepLog)
        s = stepLog(k);
        note = s.error;
        if s.epochBefore ~= s.epochAfter && isempty(s.error)
            note = sprintf('epochs %d → %d', s.epochBefore, s.epochAfter);
        end
        fprintf(fid, '%-4d  %-35s  %7.2f  %9d  %9d  %s\n', ...
            k, s.step, s.duration_s, s.chanBefore, s.chanAfter, note);
    end

    fprintf(fid, '\nTotal steps: %d\n', numel(stepLog));
    errSteps = sum(~cellfun(@isempty, {stepLog.error}));
    if errSteps > 0
        fprintf(fid, 'Steps with errors: %d\n', errSteps);
    end
    fclose(fid);
end

function act2D = computeICAActivation(EEG)
% COMPUTEICAACTIVATION  Return 2-D ICA activations (nComp x nSamples).
%   Computes activations on demand if EEG.icaact is empty.
if ~isempty(EEG.icaact)
    act2D = reshape(EEG.icaact, size(EEG.icaact,1), []);
else
    data2D = reshape(EEG.data(EEG.icachansind,:,:), numel(EEG.icachansind), []);
    act2D  = (EEG.icaweights * EEG.icasphere) * data2D;
end
end

function showReportDialog(~, summaries)
% SHOWREPORTDIALOG  Display per-file pipeline report summaries in a modal window.
if isempty(summaries); return; end
separator = [newline repmat('=', 1, 60) newline];
combined  = strjoin(summaries, separator);
fig = uifigure('Name', 'Pipeline Report', 'Position', [150 100 660 540], ...
    'WindowStyle', 'modal');
uitextarea(fig, 'Position', [10 50 640 480], ...
    'Editable', 'off', 'FontName', 'Courier New', 'FontSize', 10, ...
    'Value', combined);
uibutton(fig, 'Text', 'Close', 'Position', [280 10 100 30], ...
    'ButtonPushedFcn', @(~,~) close(fig));
uiwait(fig);
end