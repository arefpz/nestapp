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
% STEPS TO PEFORM ON EEG DATA
%
% This file contains the commands and functions already available in eeglab
% package. This program is being called by the app "nestapp".

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
            inputvals{2*j} = x{j,2};
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

for nfile = 1:app.NSelecFiles
    app.ProcessingfileEditField.Value = num2str(nfile);
    % To avoid any unforseen error, for each data new eeglab window will be
    % used.
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    app.initialVars = who; % Variables available at the begining of the analysis. Will be used to save them.
    pathName=app.path; % Path to data folder
    fileName = app.file{nfile}; % Data name(s) to be analyzed.

    disp(['!--------FILE ',fileName,' IS BEING PROCESSED--------!'])

    % In below loop, all assigned steps will be evaluated.
    for Step=app.nstep:dstep:numel(app.steps2run)
        varin = app.steps2run{Step+1};
        disp(strcat('step ',num2str(fix(Step/2)+1), ': "',app.steps2run{Step},'" is running!'));
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
                    if strcmp(needchanloc,'yes') && isempty(chName)
                        pathEEGLAB = which('eeglab');
                        pathEEGLAB = replace(pathEEGLAB,'\','/');
                        dashes = find(pathEEGLAB=='/');
                        pathEEGLAB(dashes(end)+1:end)=[];
                        [chName,chPath] = uigetfile('*.*','Select a file');
                    end

                    EEG=pop_chanedit(EEG, 'lookup',...
                        [pathEEGLAB,'plugins/dipfit/standard_BEM/elec/standard_1005.elc'],...
                        'load',{[chPath,chName],'filetype','autodetect'});

                    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

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
                        eeglab redraw
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
                        replace(fname,' ','_');
                        replace(fname,'-','_');
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
                    input('Please press enter to continue!')
                case 'Remove high-std Channels'
                    %% Remove high std channels.

                    % Here based on the standard deviation of signal for each channel, we can
                    % remove the high std channels. However in future steps, we will use
                    % different function to remove the channels.
                    period = EEG.pnts/10+1:2*EEG.pnts/10+1;
                    SD=zeros(1,EEG.nbchan);
                    for ii=1:EEG.nbchan
                        SD(ii) = std(EEG.data(ii,period));
                    end
                    bad_SD_threshold=varin{1,2};
                    bad_channels_SD=find(SD>bad_SD_threshold);
                    remove_channels = cell(1,numel(bad_channels_SD));

                    for ii=1:numel(bad_channels_SD)
                        remove_channels{ii}=EEG.chanlocs(1,bad_channels_SD(ii)).labels;
                    end
                    disp(['Channels ',remove_channels,' are removed!'])
                    EEG = pop_select( EEG,'nochannel',remove_channels);
                    EEG = eeg_checkset( EEG );

                case 'Remove un-needed Channels'
                    %% Remove un-needed channels

                    %select here the channel you want to remove if they have a bad impedence
                    vars = convertContainedStringsToChars(varin);
                    inds = find(strcmp(vars,'[]'));
                    vars([inds,inds-1])=[];
                    assignin('base',"vars",vars)
                    EEG = pop_select( EEG,vars{:});
                    EEG = eeg_checkset( EEG );

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

                case 'Remove Baseline'
                    %% Remove Baseline Offset

                    % The DC offset may cause problem in analyzing and cleaning data. So, it's
                    % better to remove it before other steps. However if we want to study the
                    % ERP, it is recomended to do not remove it, since it may affect the
                    % results.
                    vars = convertContainedStringsToChars(varin);

                    timerange=vars{1,2}; % Default [] -> all
                    pointrange=vars{1,4}; % Default [] -> all
                    chanlist=vars{1,6}; % Default [] -> all
                    EEG = pop_rmbase( EEG, eval(timerange), eval(pointrange), eval(chanlist));
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
                        ref=input('The referenced channel has been removed, please enter new ref: ','s');
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
                    vars([ind1,ind1+1,ind2,ind2+2,ind3,ind3+1])=[];

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


                case 'Remove Bad Channels'
                    %%  Remove bad channels

                    vars = convertContainedStringsToChars(varin);
                    EEGelecNames = {EEG.chanlocs(1:end).labels};
                    ind1 = find(strcmpi(vars,'impelec'));
                    AuximportantElects = vars{ind1+1};
                    importantElects=ismember(EEGelecNames, AuximportantElects);
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

                case 'Remove Bad Epoch'
                    %% Remove bad Epoch
                    % EEG = pop_jointprob(EEG,1,1:EEG.nbchan ,5,5,0,0);
                    vars = convertContainedStringsToChars(varin);
                    inds = find(strcmp(vars,'[]'));
                    vars([inds,inds-1])=[];
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
                    EEG = pop_subcomp( EEG, var_comp, vars{4}, vars{6});
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

                case 'Interpolate Missing Data TESA'
                    %% TESA Interpolate missing data LINEAR
                    % replaces missing data with linear interpolation.
                    % Linear function is fitted on data point before and after missing data.
                    vars = convertContainedStringsToChars(varin);
                    interpolation = vars{1,find(strcmp(vars,'interpolation'))+1};
                    interpWin = vars{1,find(strcmp(vars,'interpWin'))+1};

                    EEG = pop_tesa_interpdata( EEG, interpolation,interpWin);
                    EEG = eeg_checkset( EEG );

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
                    ind1 = find(strcmp{vars,'timeWin'});
                    timeWin=vars{1,ind1+1};
                    ind2 = find(strcmp{vars,'mdorder'});
                    mdorder=vars{1,ind2+1};
                    ind3 = find(strcmp{vars,'event_type'});
                    event_type=vars{1,ind3+1};
                    EEG = tesa_filtmedian( EEG, timeWin, mdorder, event_type );
                    EEG = eeg_checkset( EEG );

                case 'Remove Bad Trials'
                    %% Remove bad Trials
                    EEG = pop_jointprob(EEG,1,1:size(EEG.data,1) ,5,5,0,0);
                    pop_rejmenu(EEG,1);
                    pause_script = input('Highlight bad trials, update marks and then press enter');
                    EEG.BadTr = unique([find(EEG.reject.rejjp==1) find(EEG.reject.rejmanual==1)]);
                    EEG = pop_rejepoch( EEG, EEG.BadTr ,0);
                    EEG = eeg_checkset( EEG );

                case 'Extract TEP (TESA)'
                    %% Extract TEP
                    vars = convertContainedStringsToChars(varin);
                    ind1 = find(strcmp{vars,'type'});
                    type=vars{1,ind1+1};
                    vars([ind1, ind1+1]) = [];
                    ind2 = find(strcmpi(vars,'pairCorrect'));
                    if ~strcmp(vars{ind2+1},'on')
                        ind3 = faind(strcmpi(vars,'ISI'));
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

                    tepoutput = pop_tesa_peakoutput( EEG, vars(:) );
                    assignin('base','tep_output',tepoutput)

            end
        catch err
            disp(err.message)
            warning(strcat('An error acoured at file ',fileName,...
                ' at step ',num2str(round(Step/2)), ': ',app.steps2run{Step}{:}));
            toContinue = input('Do you want to continue? (y/n)? ','s');
            if strcmpi(toContinue,'n')
                EEG.nestappSteps = app.steps2run;
                postVars = who;
                assignin('base','postVars',postVars)
                for i = 1:numel(postVars)
                    if find(ismember(postVars{i}, app.initialVars))
                        assignin('base',postVars{i},eval(postVars{i}))
                    end
                end
                assignin('base','lastVarin',varin);
                assignin('base','lastStepInd',Step);
                eeglab redraw
                break
            end
        end
    end
    
    % assignin('base','EEG',EEG)
    % assignin('base','ALLEEG',ALLEEG)
    % eeglab redraw
    pause(2)
    disp('-----------------Data processed!-----------------')
end