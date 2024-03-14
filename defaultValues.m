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
% This app is designed to use EEGLAB, TESA, and other available packages/
% plugins, to apply to EEG data. No new analysis is being introduced other 
% than the mentioned functions. 
% 
% DEFAULT VALUES
% 
% This file contains the default values used in the original functions. You
% may change them for your frequently used values.
% 
% NOTE 1: Naming is the same as the steps name. ( Spaces, " ", and 
% parenthesis "(", ")" are replaced with underline "_").
% 
% NOTE 2: You may find some fields that are not introduced within the 
% original functions.
% 
% NOTE 3: DO NOT ADD extra fields for each structure unless it has been 
% already introduced within original functions. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data Default Values

Load_Data.redraw = 'on';

%% Load Channel Location

Load_Channel_Location.needchanloc = 'no';
Load_Channel_Location.eachFilediffPath = 'no';

%% Save New Set

Save_New_Set.setname = 'New Set';
Save_New_Set.comments = '[]';
Save_New_Set.overwrite = 'off';
Save_New_Set.saveold = '[]';
Save_New_Set.savenew = '[]';
Save_New_Set.retrieve = NaN;
Save_New_Set.gui = 'off';
Save_New_Set.includeFileName = 'yes';

%% De Trend Default Values

De_Trend_Epoch.npoly = 1;

%% TESA De-Trend

TESA_De_Trend.detrend = 'linear';
TESA_De_Trend.timeWin = [11, 500];

%% Re Referencing Default Values
Re_Reference.ref = 'Cz';
Re_Reference.exclude = [];
Re_Reference.interpchan = 'off';
Re_Reference.keepref = 'off';
Re_Reference.refloc = [];

%% Re Sampling Default Values

Re_Sample.freq = 1000;
Re_Sample.fc = 0.9;
Re_Sample.df = 0.2;

%% Remove Baseline Default Values

Remove_Baseline.timerange = [];
Remove_Baseline.pointrange = [];
Remove_Baseline.chanlist = [];

%% Remove high std Channels Default Values

Remove_high_std_Channels.bad_SD_threshold = 500;

%% Remove un-needed Channels
Remove_un_needed_Channels.time=[];
Remove_un_needed_Channels.rmtime=[];
Remove_un_needed_Channels.point=[];
Remove_un_needed_Channels.rmpoint=[];
Remove_un_needed_Channels.trial=[];
Remove_un_needed_Channels.rmtrial=[];
Remove_un_needed_Channels.sorttrial='on';
Remove_un_needed_Channels.nochannel={'TP9','TP10'};
Remove_un_needed_Channels.channel = [];
Remove_un_needed_Channels.chantype=[];
Remove_un_needed_Channels.rmchantype=[];

%% Frequency Filter EEGLAB Default Values
Frequency_Filter.locutoff=1;
Frequency_Filter.hicutoff=80;
Frequency_Filter.filtorder=0;
Frequency_Filter.revfilt=0;
Frequency_Filter.usefft=0;
Frequency_Filter.plotfreqz=0;
Frequency_Filter.minphase=0;
Frequency_Filter.usefftfilt=0;


%% Frequency Filter TESA Default Values

Frequency_Filter_TESA.high = 1;
Frequency_Filter_TESA.low = 80;
Frequency_Filter_TESA.ord = 4;
Frequency_Filter_TESA.type = 'bandpass';

%% Frequency Filter CleanLine Default Values

Frequency_Filter_CleanLine.bandwidth =  2;
Frequency_Filter_CleanLine.chanlist =  [1 64];
Frequency_Filter_CleanLine.computepower =  1;
Frequency_Filter_CleanLine.linefreqs =  [60 120];
Frequency_Filter_CleanLine.newversion =  0;
Frequency_Filter_CleanLine.normSpectrum =  0;
Frequency_Filter_CleanLine.p =  0.0100;
Frequency_Filter_CleanLine.pad =  2;
Frequency_Filter_CleanLine.plotfigures =  0;
Frequency_Filter_CleanLine.scanforlines =  1;
Frequency_Filter_CleanLine.sigtype =  'Channels';
Frequency_Filter_CleanLine.taperbandwidth =  2;
Frequency_Filter_CleanLine.tau =  100;
Frequency_Filter_CleanLine.verb =  1;
Frequency_Filter_CleanLine.winsize =  4;
Frequency_Filter_CleanLine.winstep =  1;

%% Visualize EEG Data Default Values

Visualize_EEG_Data.state = [1, 1, 1];

%% Visualize ICA Component Default Values

% Select_ICA_Comp.compnum = [1:63];

%% Visualize ICA Component Default Values

Visualize_ICA_Comp_TESA.figSize = 'small';
Visualize_ICA_Comp_TESA.freqScale = 'db';
Visualize_ICA_Comp_TESA.plotFreqX = [2, 45];
Visualize_ICA_Comp_TESA.plotTimeX = [-600, 600];
Visualize_ICA_Comp_TESA.saveWeights = 'on';

%%  Runing ICA Default Values

Run_ICA.icatype = 'fastica';
Run_ICA.approach = 'symm';
Run_ICA.g = 'tanh';
Run_ICA.stabilization = 'on';

%% Runing TESA ICA
Run_TESA_ICA.approach = 'symm';
Run_TESA_ICA.g =  'gauss';
Run_TESA_ICA.stabilization  = 'on';

%% Plot Component Spectra and maps
% Plot_Component_Spectra_and_maps.dataflag = 1;
% Plot_Component_Spectra_and_maps.Epoch_time_range = [0, 10000];
Plot_Component_Spectra_and_Maps.freq = 10;
Plot_Component_Spectra_and_Maps.plotchan = 0;
Plot_Component_Spectra_and_Maps.percent = 20;
Plot_Component_Spectra_and_Maps.icacomps = 1:32;
Plot_Component_Spectra_and_Maps.nicamaps = 5;
Plot_Component_Spectra_and_Maps.freqrange = [2 25];
Plot_Component_Spectra_and_Maps.electrodes = 'off';

%% Find TMS Pulses (TESA)

Find_TMS_Pulses_TESA.elec = 'Cz';
Find_TMS_Pulses_TESA.refract = 3;
Find_TMS_Pulses_TESA.rate = 10000;
Find_TMS_Pulses_TESA.tmslabel = 'TMS';
Find_TMS_Pulses_TESA.plots = 'off';
Find_TMS_Pulses_TESA.paired = 'no';
Find_TMS_Pulses_TESA.repetitive = 'no';
Find_TMS_Pulses_TESA.ISI = [];
Find_TMS_Pulses_TESA.pairlabel = {'pp'};
Find_TMS_Pulses_TESA.ITI = [];
Find_TMS_Pulses_TESA.pulseNum = [];

%% Interpolate Channels Default Values

Interpolate_Channels.method = 'spherical';
Interpolate_Channels.trange = [];


%% Remove ICA components TESA
Remove_ICA_Components_TESA.compCheck = 'off';
Remove_ICA_Components_TESA.comps = 15;
Remove_ICA_Components_TESA.tmsMuscle = 'on';
Remove_ICA_Components_TESA.tmsMuscleThresh = 8;
Remove_ICA_Components_TESA.tmsMuscleWin = [11 30];
Remove_ICA_Components_TESA.tmsMuscleFeedback = 'off';
Remove_ICA_Components_TESA.blink = 'off';
Remove_ICA_Components_TESA.blinkThresh = 2.5;
Remove_ICA_Components_TESA.blinkElecs = {'Fp1','Fp2'};
Remove_ICA_Components_TESA.blinkFeedback = 'off';
Remove_ICA_Components_TESA.move = 'off';
Remove_ICA_Components_TESA.moveThresh = 2;
Remove_ICA_Components_TESA.moveElecs = {'F7','F8'};
Remove_ICA_Components_TESA.moveFeedback = 'off';
Remove_ICA_Components_TESA.muscle = 'off';
Remove_ICA_Components_TESA.muscleThresh = 0.6;
Remove_ICA_Components_TESA.muscleFreqIn = [30 100];
Remove_ICA_Components_TESA.muscleFeedback = 'off';
Remove_ICA_Components_TESA.elecNoise = 'off';
Remove_ICA_Components_TESA.elecNoiseThresh = 4;
Remove_ICA_Components_TESA.elecNoiseFeedback = 'off';

%% Find Artifacts EDM (TESA)
Find_Artifacts_EDM_TESA.chanl = [];
Find_Artifacts_EDM_TESA.nc = 30;
Find_Artifacts_EDM_TESA.sf = 1000;
Find_Artifacts_EDM_TESA.comps = 10;
Find_Artifacts_EDM_TESA.tmsMuscleThresh = 10;
Find_Artifacts_EDM_TESA.tmsMuscleWin = [11,50];
Find_Artifacts_EDM_TESA.tmsMuscleFeedback = 'on';

%% SSP SIR
SSP_SIR.artScale =  'automatic';
SSP_SIR.timeRange =  [5,50];
SSP_SIR.PC = {'data', 90};

%% Median Filter
Median_Filter_1D.timeWin = [-2,2];
Median_Filter_1D.mdorder = 5; % % Apply a X order median filter to remove small spike artifacts
Median_Filter_1D.event_type = 'spike'; % Or 'TMS'
 
%% Interpolate Missing Data TESA
Interpolate_Missing_Data_TESA.interpolation = 'linear';
Interpolate_Missing_Data_TESA.interpWin = [20, 20];

%% Epoching
Epoching.types = {'TMS'};
Epoching.timelim = [-2, 2];
Epoching.epochinfo = 'yes';

%% Remove TMS Artifacts (TESA)
Remove_TMS_Artifacts_TESA.cutTimesTMS = [-2, 10];
Remove_TMS_Artifacts_TESA.replaceTimes = [];
Remove_TMS_Artifacts_TESA.cutEvent = {'TMS'};


%% Automatic Cleaning Data
Automatic_Cleaning_Data.FlatlineCriterion=5;
Automatic_Cleaning_Data.ChannelCriterion=0.8;
Automatic_Cleaning_Data.LineNoiseCriterion=4;
Automatic_Cleaning_Data.Highpass='off';
Automatic_Cleaning_Data.BurstCriterion=20;
Automatic_Cleaning_Data.WindowCriterion=0.25;
Automatic_Cleaning_Data.BurstRejection='on';
Automatic_Cleaning_Data.Distance='Euclidian';
Automatic_Cleaning_Data.channels_ignore={''};
Automatic_Cleaning_Data.WindowCriterionTolerances=[-Inf 7];

%% Remove bad channels 
Remove_Bad_Channels.impelec = {''};
Remove_Bad_Channels.elec=[];
Remove_Bad_Channels.threshold=5;
Remove_Bad_Channels.norm='on';
Remove_Bad_Channels.measure='kurt';
Remove_Bad_Channels.freqrange=[1 6];

%% Remove Bad Epoch
Remove_Bad_Epoch.threshold = 1000; % in uV
Remove_Bad_Epoch.electrodes = [];
Remove_Bad_Epoch.icacomps = [];
Remove_Bad_Epoch.startprob = [];
Remove_Bad_Epoch.maxrej = [];
Remove_Bad_Epoch.nogui='on'; 
Remove_Bad_Epoch.eegplot='off';

%% Fix TMS Pulse (TESA)
Fix_TMS_Pulse_TESA.elec='Cz';
Fix_TMS_Pulse_TESA.epoch_len=[-1,1];
Fix_TMS_Pulse_TESA.type = 'TMS';
Fix_TMS_Pulse_TESA.refract= 4;
Fix_TMS_Pulse_TESA.rate=1e4;
Fix_TMS_Pulse_TESA.paired= 'no';
Fix_TMS_Pulse_TESA.ISI=[];

%% Remove Recording Noise
Remove_Recording_Noise_SOUND.lambdaValue=0.2;
Remove_Recording_Noise_SOUND.iter=10;

%% Extract TEP (TESA)
Extract_TEP_TESA.type = 'ROI';
Extract_TEP_TESA.elecs = 'all';
Extract_TEP_TESA.tepName = [];
Extract_TEP_TESA.pairCorrect ='on';
Extract_TEP_TESA.ISI = [];

%% Find_TEP_Peaks_TESA
Find_TEP_Peaks_TESA.input = 'ROI';
Find_TEP_Peaks_TESA.direction = 'positive';
Find_TEP_Peaks_TESA.peak = [30 50;70 90;180 220];
Find_TEP_Peaks_TESA.peakWin = 0;
Find_TEP_Peaks_TESA.method = 'largest';
Find_TEP_Peaks_TESA.samples = 5;
Find_TEP_Peaks_TESA.tepName = 'TEP';

%% Plot TEP TESA
Plot_TEP_TESA.tepType = 'ROI';
Plot_TEP_TESA.tepName = 'R1';
Plot_TEP_TESA.xlim = [-100 500];
Plot_TEP_TESA.ylim = [];
Plot_TEP_TESA.CI = 'off';
Plot_TEP_TESA.plotPeak = 'on';

%% TEP peak output
TEP_Peak_Output.tepName = 'R1';
TEP_Peak_Output.calcType = 'amplitude';
TEP_Peak_Output.winType = 'individual';
TEP_Peak_Output.averageWin = 5;
TEP_Peak_Output.fixedPeak = [];
TEP_Peak_Output.tablePlot = 'on';

%% Flag ICA Components
Flag_ICA_Components_for_Rejection.Brain=[NaN,NaN];
Flag_ICA_Components_for_Rejection.Muscle=[0.9,1];
Flag_ICA_Components_for_Rejection.Eye=[0.9,1];
Flag_ICA_Components_for_Rejection.Heart=[NaN,NaN];
Flag_ICA_Components_for_Rejection.LineNoise=[NaN,NaN];
Flag_ICA_Components_for_Rejection.ChannelNoise=[NaN,NaN];
Flag_ICA_Components_for_Rejection.Other=[NaN,NaN];

%% Label ICA Components
Label_ICA_Components.version='default';

%% Remove Flagged ICA Components
Remove_Flagged_ICA_Components.components=[];
Remove_Flagged_ICA_Components.plotag=0;
Remove_Flagged_ICA_Components.keepcomp=0;

%% Manual Commands
Manual_Command.command='';

%% Edit Data Set Info
% Edit_Data_Set_info.setname="Sets' Name" ;
% Edit_Data_Set_info.data="" ;
% Edit_Data_Set_info.dataformat="" ;
% Edit_Data_Set_info.subject="" ;
% Edit_Data_Set_info.condition="" ;
% Edit_Data_Set_info.group="" ;
% Edit_Data_Set_info.run="";
% Edit_Data_Set_info.session="";
% Edit_Data_Set_info.chanlocs="";
% Edit_Data_Set_info.nbchan="";
% Edit_Data_Set_info.xmin="";
% Edit_Data_Set_info.pnts="";
% Edit_Data_Set_info.srate="";
% Edit_Data_Set_info.ref="";
% Edit_Data_Set_info.icaweight="";
% Edit_Data_Set_info.icasphere="";
% Edit_Data_Set_info.comments="";

%% Choose Data Set

Choose_Data_Set.dataSetInd = '';

%% Automatic Continuous Rejection
Automatic_Continuous_Rejection.elecrange = 1:64;
Automatic_Continuous_Rejection.epochlength = 0.5;
Automatic_Continuous_Rejection.overlap = 0.25;
Automatic_Continuous_Rejection.freqlimit = [35 128];
Automatic_Continuous_Rejection.mode = 'max';
Automatic_Continuous_Rejection.correct = 'remove';
Automatic_Continuous_Rejection.threshold = 10;
Automatic_Continuous_Rejection.contiguous = 4;
Automatic_Continuous_Rejection.addlength = 0.25;
Automatic_Continuous_Rejection.eegplot = 'off';
Automatic_Continuous_Rejection.onlyreturnselection = 'off';
Automatic_Continuous_Rejection.verbose = 'off';
Automatic_Continuous_Rejection.taper = 'hamming';
