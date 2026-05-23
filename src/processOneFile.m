function [fileReport, stepLog] = processOneFile(spec, fullPath, opts)
% PROCESSONEFILE  Execute a typed pipeline spec against a single EEG file.
%   [fileReport, stepLog] = PROCESSONEFILE(spec, fullPath, opts)
%
%   spec      - struct array of PipelineStep (name + params struct)
%   fullPath  - absolute path to one EEG data file
%   opts      - struct with fields:
%     .pipelineName   - pipeline name string for EEG.history provenance
%     .chanLocFile    - pre-selected channel location file path (or '')
%     .progressFcn    - @(si, stepName) called before each step (serial mode)
%     .progressQueue  - parallel.pool.DataQueue (parallel mode) or []
%     .onStepError    - @(si, stepName, err) -> 'Continue'|'Abort', or []
%     .onPickChanFile - @() -> charPath for per-file channel picking (serial only)
%     .fileIndex      - integer, used in progressQueue messages
%     .uiFigure       - UIFigure handle (needed by interactive steps in serial mode)
%
%   fileReport - PipelineReport struct (see initPipelineReport)
%   stepLog    - struct array with per-step timing and channel counts
%
%   On any unrecoverable error, throws 'nestapp:cancelled'. The caller is
%   responsible for closing progress dialogs and collecting partial results.
%
%   See also: runPipelineCore, initPipelineReport, exportReport

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'pipelineName'),    opts.pipelineName    = ''; end
if ~isfield(opts, 'chanLocFile'),     opts.chanLocFile     = ''; end
if ~isfield(opts, 'progressFcn'),    opts.progressFcn    = []; end
if ~isfield(opts, 'progressQueue'),  opts.progressQueue  = []; end
if ~isfield(opts, 'onStepError'),    opts.onStepError    = []; end
if ~isfield(opts, 'onPickChanFile'), opts.onPickChanFile = []; end
if ~isfield(opts, 'fileIndex'),      opts.fileIndex      = 0; end
if ~isfield(opts, 'uiFigure'),       opts.uiFigure       = []; end
if ~isfield(opts, 'logQueue'),       opts.logQueue       = []; end
if ~isfield(opts, 'nWorkers'),      opts.nWorkers       = 1;  end
if ~isfield(opts, 'batchCtx'),          opts.batchCtx          = []; end
if ~isfield(opts, 'autoQualityReport'), opts.autoQualityReport = false; end
if ~isfield(opts, 'qcAttribute'),       opts.qcAttribute       = 'minmax_no_tms'; end
if ~isfield(opts, 'qcTmsWindow'),       opts.qcTmsWindow       = [0 25];          end
if ~isfield(opts, 'qcTmsAutoDetect'),   opts.qcTmsAutoDetect   = true;            end
if ~isfield(opts, 'skipOnQualityFail'), opts.skipOnQualityFail = false;           end
if ~isfield(opts, 'autoExportPDF'),     opts.autoExportPDF     = false;           end

% eeglab('nogui') is expensive (plugin scan, path setup); run it once per
% worker then just reset globals for subsequent files on the same worker.
persistent eeglabWorkerReady lastThreadCount
global EEG ALLEEG CURRENTSET ALLCOM %#ok<GVMIS>

[pathDir, fileBase, fileExt] = fileparts(fullPath);
pathName = [pathDir, filesep];
fileName = [fileBase, fileExt];

wLabel   = sprintf('FILE-%d %s', opts.fileIndex, fileBase);
fileTic  = tic;

% Limit each worker to its fair share of BLAS threads to prevent
% over-subscription when N workers each default to all cores.
% Pool reuse can change nWorkers between pipeline runs, so re-check each call
% but skip the BLAS call when the value is already correct.
if ~isempty(opts.progressQueue) && opts.nWorkers > 1
    desired = max(1, floor(feature('numcores') / opts.nWorkers));
    if isempty(lastThreadCount) || lastThreadCount ~= desired
        maxNumCompThreads(desired);
        lastThreadCount = desired;
    end
end

if isempty(eeglabWorkerReady)
    sendWorkerLog(opts.logQueue, wLabel, 'eeglab(''nogui'') - first file on this worker, initializing...');
    t0eeg = tic;
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');
    eeglabWorkerReady = true;
    sendWorkerLog(opts.logQueue, wLabel, 'eeglab(''nogui'') done (%.2fs)', toc(t0eeg));
else
    sendWorkerLog(opts.logQueue, wLabel, 'EEGLAB already initialized on this worker - resetting globals only');
    EEG = []; ALLEEG = []; CURRENTSET = 0; ALLCOM = {};
end

ICA_Rejected_Comp = {};
interpElecs       = {};
pendingICAStats   = struct();
histLenBefore     = 0;

nSteps = numel(spec);

sendWorkerLog(opts.logQueue, wLabel, 'START  %s  (%d steps)', fileName, nSteps);

% Stagger the first worker wave to avoid simultaneous disk reads at batch start.
% Files beyond the first wave are already staggered naturally by worker completion times.
if ~isempty(opts.progressQueue) && opts.fileIndex > 1 && opts.fileIndex <= opts.nWorkers
    pause(0.25 * (opts.fileIndex - 1));
end

stepLog = struct('step',{},'duration_s',{},'chanBefore',{},'chanAfter',{}, ...
                 'epochBefore',{},'epochAfter',{},'error',{});
fileReport = initPipelineReport(fullPath);

for si = 1:nSteps
    step     = spec(si);
    stepName = step.name;
    varin    = paramsToVarin(step.params);

    % Progress notification - may throw nestapp:cancelled if user cancelled.
    if ~isempty(opts.progressFcn)
        opts.progressFcn(si, stepName);
    elseif ~isempty(opts.progressQueue)
        send(opts.progressQueue, struct( ...
            'fi', opts.fileIndex, 'si', si, ...
            'nSteps', nSteps, 'stepName', stepName));
    end

    sendWorkerLog(opts.logQueue, wLabel, 'Step %d/%d START  "%s"', si, nSteps, stepName);

    if isstruct(EEG) && ~isempty(EEG)
        nChanBefore  = EEG.nbchan;
        nEpochBefore = size(EEG.data, 3);
    else
        nChanBefore  = 0;
        nEpochBefore = 0;
    end
    t0 = tic;

    try
        switch stepName
            case 'Load Channel Location'
                vars = convertContainedStringsToChars(varin);
                ind1 = find(strcmpi(vars,'eachFilediffPath'));
                eachFilediffPath = vars{ind1+1};
                ind2 = find(strcmpi(vars,'needchanloc'));
                needchanloc = vars{ind2+1};

                if strcmp(eachFilediffPath,'yes')
                    needchanloc = 'yes';
                end
                if strcmp(needchanloc,'yes')
                    if strcmp(eachFilediffPath,'yes') && ~isempty(opts.onPickChanFile)
                        chanLocFile = opts.onPickChanFile();
                    elseif ~isempty(opts.chanLocFile)
                        chanLocFile = opts.chanLocFile;
                    else
                        error('nestapp:noChanFile', ...
                            'Load Channel Location: no channel file available.');
                    end
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
                    [chPath, chBase, chExt] = fileparts(chanLocFile);
                    chPath = [chPath, filesep];
                    EEG = pop_chanedit(EEG, 'lookup', lookforchnlocs, ...
                        'load', {[chPath, chBase, chExt], 'filetype', 'autodetect'});
                    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
                end

            case 'Load Data'
                if   strcmpi(fileName(end-2:end),'set')
                    EEG = pop_loadset( [pathName fileName]);
                elseif strcmpi(fileName(end-2:end),'cnt')
                    EEG = pop_loadcnt([pathName fileName] , 'dataformat', 'int32' );
                elseif strcmpi(fileName(end-2:end),'cdt')
                    EEG = loadcurry([pathName fileName], 'CurryLocations', 'False');
                elseif strcmpi(fileName(end-3:end),'vhdr')
                    EEG = pop_loadbv(pathName , fileName );
                else
                    error('nestapp:unknownFormat', ...
                        'Load Data: unrecognized file extension in "%s". Supported: .set, .cnt, .cdt, .vhdr', ...
                        fileName);
                end
                EEG.filename = fileName;
                [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
                if isfield(EEG, 'history')
                    histLenBefore = numel(EEG.history);
                end

            case 'Save New Set'
                EEG = eeg_checkset( EEG );
                vars = convertContainedStringsToChars(varin);
                inds = find(strcmpi(vars,'includeFileName'));
                IFN = vars{inds + 1};
                vars([inds, inds+1]) = [];
                fname = '';
                if strcmp(IFN,'yes')
                    [~, fbase, ~] = fileparts(fullfile(pathName, fileName));
                    fbase = replace(fbase, ' ', '_');
                    fbase = replace(fbase, '-', '_');
                    % .set destination now lives under the batch folder
                    % (data/ for typeBased, <stem>/ for perInput).
                    if ~isempty(opts.batchCtx)
                        targetDir = outputPaths(opts.batchCtx, 'data', fbase);
                    else
                        targetDir = pathDir;
                    end
                    fname = fullfile(targetDir, [fbase, '_']);
                end
                ind1 = find(strcmp(vars,'savenew'));
                sv1 = vars{ind1+1};
                if ischar(sv1) && ~isempty(sv1) && ~strcmp(sv1,'[]')
                    vars{ind1+1} = [fname, sv1];
                end
                ind2 = find(strcmp(vars,'saveold'));
                sv2 = vars{ind2+1};
                if ischar(sv2) && ~isempty(sv2) && ~strcmp(sv2,'[]')
                    vars{ind2+1} = [fname, sv2];
                end
                vars = stripEmptyVarin(vars);
                EEG = eeg_checkset(EEG);
                [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,vars{:});

            case 'Manual Command'
                cmd = step.params.command;
                if ischar(cmd) && ~isrow(cmd)
                    cmd = strjoin(cellstr(cmd), newline);
                elseif iscell(cmd)
                    cmd = strjoin(cmd(:)', newline);
                end
                eval(cmd);

            case 'Choose Data Set'
                vars = convertContainedStringsToChars(varin);
                ind = find(strcmp(vars,'dataSetInd'));
                setIndex = vars{ind+1};
                [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',setIndex);

            case 'Visualize EEG Data'
                vars  = convertContainedStringsToChars(varin);
                ind   = find(strcmpi(vars, 'state'));
                state = vars{ind+1};
                if sum(isnan(state))
                    pop_eegplot(EEG);
                else
                    pop_eegplot(EEG, state(1), state(2), state(3));
                end
                uiconfirm(opts.uiFigure,'Press OK when done viewing the EEG plot.','Visualize EEG','Options',{'OK'},'DefaultOption',1);

            case 'Remove un-needed Channels'
                vars = convertContainedStringsToChars(varin);
                vars = stripEmptyVarin(vars);
                EEG = pop_select( EEG,vars{:});
                EEG = eeg_checkset( EEG );

            case 'Remove Bad Channels'
                vars = convertContainedStringsToChars(varin);
                EEGelecNames = {EEG.chanlocs(1:end).labels};
                ind1 = find(strcmpi(vars,'impelec'));
                AuximportantElects = vars{ind1+1};
                importantElects = matches(EEGelecNames, AuximportantElects,"IgnoreCase",true);
                vars([ind1, ind1+1]) = [];

                ind2 = find(strcmpi(vars,'elec'));
                if strcmp(vars{1,ind2+1},'[]')
                    vars{1,ind2+1} = 1:EEG.nbchan;
                elseif iscell(vars{1,ind2+1})
                    vars{1,ind2+1} = find(ismember(EEGelecNames, vars{1,ind2+1}));
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
                ind = find(strcmpi(vars,'elecrange'));
                if max(vars{ind+1})>EEG.nbchan
                    elecrange = 1:EEG.nbchan;
                else
                    elecrange = vars{ind+1}(1):vars{ind+1}(end);
                end
                vars([ind,ind+1]) = [];
                EEG = pop_rejcont(EEG,'elecrange',elecrange,vars{:});

            case 'Clean Artifacts'
                vars = convertContainedStringsToChars(varin);
                ind1 = find(strcmpi(vars,'Channels'));
                ind2 = find(strcmpi(vars,'Channels_ignore'));
                if sum(strcmpi(vars{ind1+1},'[]')) && sum(strcmpi(vars{ind2+1},'[]'))
                    chans = {EEG.chanlocs.labels};
                    vars{ind1+1} = chans;
                    vars{ind2+1} = [];
                elseif sum(strcmpi(vars{ind1+1},'[]')) && sum(~strcmpi(vars{ind2+1},'[]'))
                    if size(vars{ind2+1},1)>size(vars{ind2},2)
                        vars{ind2+1} = vars{ind2+1}';
                    end
                    vars([ind1 ind1+1]) = [];
                elseif sum(~strcmpi(vars{ind1+1},'[]')) && sum(strcmpi(vars{ind2+1},'[]'))
                    if size(vars{ind1+1},1)>size(vars{ind1},2)
                        vars{ind1+1} = vars{ind1+1}';
                    end
                    vars([ind2 ind2+1]) = [];
                end
                EEG = clean_artifacts(EEG,vars{:});
                EEG = eeg_checkset( EEG );

            case 'Automatic Cleaning Data'
                vars = convertContainedStringsToChars(varin);
                ind = find(strcmp(vars,'Highpass'));
                if ~strcmpi(vars{ind+1},'off')
                    highpass = vars{ind+1};
                    if ischar(highpass) || isstring(highpass)
                        highpass = str2double(highpass);
                    end
                    if size(highpass,2)<size(highpass,1)
                        highpass = highpass';
                    end
                    vars{ind+1} = highpass;
                end
                vars = stripEmptyVarin(vars);
                EEG = pop_clean_rawdata(EEG, vars{:});
                EEG = eeg_checkset( EEG );

            case 'Remove Baseline'
                vars = convertContainedStringsToChars(varin);
                ind  = find(strcmpi(vars, 'timerange'));
                timerange = vars{ind+1};
                if ischar(timerange) && strcmp(timerange, '[]')
                    timerange = [];
                elseif ischar(timerange) || isstring(timerange)
                    timerange = str2num(char(timerange)); %#ok<ST2NM>
                end
                if isnumeric(timerange) && numel(timerange) == 2
                    timerange(1) = max(timerange(1), EEG.times(1));
                    timerange(2) = min(timerange(2), EEG.times(end));
                end
                EEG = pop_rmbase(EEG, timerange);
                EEG = eeg_checkset(EEG);

            case 'De-Trend Epoch'
                vars = convertContainedStringsToChars(varin);
                [nCh, nT, nEp] = size(EEG.data);
                d2 = reshape(permute(EEG.data, [2 1 3]), nT, nCh*nEp);
                ind = find(strcmpi(vars, 'npoly'));
                d2 = detrend(d2, vars{ind+1});
                EEG.data = permute(reshape(d2, nT, nCh, nEp), [2 1 3]);
                EEG = eeg_checkset( EEG );

            case 'TESA De-Trend'
                vars = convertContainedStringsToChars(varin);
                ind1 = find(strcmpi(vars,'detrend'));
                Tdetrend = vars{ind1+1};
                ind2 = find(strcmpi(vars,'timeWin'));
                TtimeWin = vars{ind2+1};
                pop_tesa_detrend(EEG, Tdetrend, TtimeWin)

            case 'Re-Sample'
                vars = convertContainedStringsToChars(varin);
                EEG = pop_resample(EEG,vars{2:2:end});
                EEG = eeg_checkset( EEG );

            case 'Re-Reference'
                vars = convertContainedStringsToChars(varin);
                ind = find(strcmp(vars,'ref'));
                ref = vars{ind+1};
                if ~ismember(ref,{EEG.chanlocs.labels}) && ~strcmp(ref,'[]')
                    answer = inputdlg('The reference channel is not in the data. Enter a new reference channel label:','Re-Reference',[1 50],{''});
                    if isempty(answer) || isempty(answer{1})
                        error('Re-Reference cancelled: no reference channel provided.');
                    end
                    ref = answer{1};
                end
                if strcmp(ref,'[]')
                    ref = eval(ref);
                end
                vars([ind,ind+1]) = [];
                vars = stripEmptyVarin(vars);
                EEG = pop_reref(EEG,ref, vars{:});
                EEG = eeg_checkset( EEG );

            case 'Frequency Filter (CleanLine)'
                if ~exist('hlp_varargin2struct', 'file')
                    cleanlineRoot = fileparts(which('pop_cleanline'));
                    if ~isempty(cleanlineRoot)
                        addpath(genpath(cleanlineRoot));
                    end
                end
                vars = convertContainedStringsToChars(varin);
                ind = find(strcmp(vars,'chanlist'));
                if vars{ind+1}(2) > EEG.nbchan
                    vars{ind+1} = 1:EEG.nbchan-1;
                else
                    vars{ind+1} = 1:vars{1,ind+1};
                end
                EEG = pop_cleanline(EEG, vars{:});
                EEG = eeg_checkset( EEG );

            case 'Frequency Filter (TESA)'
                vars = convertContainedStringsToChars(varin);
                ind1 = find(strcmp(vars,'high'));
                high = vars{ind1+1};
                ind2 = find(strcmp(vars,'low'));
                low = vars{ind2+1};
                ind3 = find(strcmp(vars,'ord'));
                ord = vars{ind3+1};
                ind4 = find(strcmp(vars,'type'));
                type = vars{ind4+1};
                vars([ind1,ind1+1,ind2,ind2+2,ind3,ind3+1]) = []; %#ok<NASGU>
                EEG = pop_tesa_filtbutter( EEG, high, low, ord, type );
                EEG = eeg_checkset( EEG );

            case 'Frequency Filter'
                vars = convertContainedStringsToChars(varin);
                ind = find(strcmp(vars,'filtorder'));
                if mod(vars{ind+1},2)~=0 || isstring(vars{ind})
                    error('The Filtorder should be an even number!')
                elseif vars{ind+1}==0
                    vars([ind, ind+1]) = [];
                end
                EEG = pop_eegfiltnew(EEG, vars{:});

            case 'Remove Bad Epoch'
                vars = convertContainedStringsToChars(varin);
                vars = stripEmptyVarin(vars);
                [EEG, rejepochs] = pop_autorej(EEG, vars{:});
                EEG.rejEpochs = rejepochs;
                EEG = eeg_checkset( EEG );

            case 'Run ICA'
                EEG.data = double(EEG.data);
                vars = convertContainedStringsToChars(varin);
                EEG = pop_runica(EEG,vars{:});
                EEG = eeg_checkset( EEG );

            case 'Label ICA Components'
                vars = convertContainedStringsToChars(varin);
                ind = find(strcmpi(vars, 'iclabelVersion'));
                iclabelVersion = vars{ind+1};
                EEG = pop_iclabel(EEG, iclabelVersion);
                EEG = eeg_checkset( EEG );

            case 'Flag ICA Components for Rejection'
                vars = convertContainedStringsToChars(varin);
                threshold = zeros(7,2);
                for nflag = 2:2:14
                    threshold(nflag/2,:) = vars{nflag};
                end
                EEG = pop_icflag(EEG, threshold);
                EEG = eeg_checkset( EEG );

            case 'Remove Flagged ICA Components'
                var_comp = step.params.components;
                plotag   = step.params.plotag;
                keepcomp = step.params.keepcomp;

                Rej = EEG.reject.gcompreject;
                Rej = reshape(Rej, 1, numel(Rej));
                ICA_Rejected_Comp{end+1} = Rej; %#ok<AGROW>

                gcr = EEG.reject.gcompreject;
                if isnumeric(gcr) || islogical(gcr)
                    rejMask = logical(reshape(gcr, 1, []));
                else
                    rejMask = false(1, max(size(EEG.icaweights, 1), 0));
                end
                pendingICAStats = struct('rejMask', rejMask);
                if ~isempty(EEG.icaweights) && size(EEG.icaweights,1) == numel(rejMask)
                    act2D = computeICAActivation(EEG);
                    % Cache activations so pop_subcomp's eeg_getica doesn't recompute them.
                    if isempty(EEG.icaact)
                        EEG.icaact = reshape(act2D, size(EEG.icaweights,1), ...
                            size(EEG.data,2), size(EEG.data,3));
                    end
                    data2D = reshape(EEG.data(EEG.icachansind,:,:), numel(EEG.icachansind), []);
                    totalVar = sum(var(data2D, 0, 2));
                    if totalVar > 0
                        pendingICAStats.compVarPct = double((var(act2D, 0, 2) / totalVar * 100)');
                    end
                end
                if isfield(EEG,'etc') && isfield(EEG.etc,'ic_classification') && ...
                        isfield(EEG.etc.ic_classification,'ICLabel') && ...
                        isfield(EEG.etc.ic_classification.ICLabel,'classifications')
                    pendingICAStats.iclabelProbs = ...
                        EEG.etc.ic_classification.ICLabel.classifications;
                end

                if ~(isnumeric(EEG.reject.gcompreject) || islogical(EEG.reject.gcompreject))
                    EEG.reject.gcompreject = zeros(1, size(EEG.icaweights, 1));
                end
                EEG = pop_subcomp( EEG, var_comp, plotag, keepcomp);
                EEG.ICA_Rejected_Comp = ICA_Rejected_Comp;
                EEG = eeg_checkset( EEG );

            case 'Interpolate Channels'
                method = step.params.method;
                trange = step.params.trange;
                EEG = pop_interp(EEG, EEG.chaninfo.removedchans(1:size(EEG.chaninfo.removedchans,2)), method, trange);
                interpElecs = [interpElecs; num2cell(EEG.chaninfo.removedchans)]; %#ok<AGROW>
                EEG.interpElecs = interpElecs;
                EEG.setname = [EEG.setname '_interp'];
                EEG.filename = [EEG.setname '.set'];
                EEG.datfile  = [EEG.setname '.fdt'];
                EEG = eeg_checkset( EEG );

            case 'Find TMS Pulses (TESA)'
                vars = convertContainedStringsToChars(varin);
                ind = find(strcmp(vars,'elec'));
                elec = vars{ind+1};
                vars([ind,ind+1]) = [];
                if iscell(elec)
                    elec = elec{:};
                end
                EEG = pop_tesa_findpulse( EEG, elec, vars{:});
                EEG = eeg_checkset( EEG );

            case 'Fix TMS Pulse (TESA)'
                vars = convertContainedStringsToChars(varin);
                ind1 = find(strcmp(vars,'elec'));
                elec = vars{1,ind1+1};
                ind2 = find(strcmp(vars,'epoch_len'));
                epoch_len = vars{1,ind2+1};
                ind3 = find(strcmp(vars,'type'));
                type = vars{1,ind3+1};
                vars([ind1,ind1+1,ind2,ind2+1,ind3,ind3+1]) = [];
                EEG = tesa_fixevent( EEG, elec, epoch_len, type, vars{:} );
                EEG = eeg_checkset( EEG );

            case 'Remove TMS Artifacts (TESA)'
                vars = convertContainedStringsToChars(varin);
                cutTimesTMS = vars{1,find(strcmp(vars,'cutTimesTMS'))+1};
                replaceTimes = vars{1,find(strcmp(vars,'replaceTimes'))+1};
                cutEvent = vars{1,find(strcmp(vars,'cutEvent'))+1};
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
                vars = convertContainedStringsToChars(varin);
                vars = stripEmptyVarin(vars);
                ind1 = find(strcmp(vars,'types'));
                type = vars{ind1+1};
                ind2 = find(strcmp(vars,'timelim'));
                timelim = vars{ind2+1};
                vars([ind1, ind1+1, ind2, ind2+1]) = [];
                EEG = pop_epoch( EEG, type, timelim, vars{:});
                EEG = eeg_checkset( EEG );

            case 'Interpolate Missing Data (TESA)'
                vars = convertContainedStringsToChars(varin);
                interpolation = vars{1,find(strcmp(vars,'interpolation'))+1};
                interpWin = vars{1,find(strcmp(vars,'interpWin'))+1};
                if size(interpWin,1) > size(interpWin,2)
                    interpWin = interpWin';
                end
                EEG = pop_tesa_interpdata( EEG, interpolation, interpWin);

            case 'Run TESA ICA'
                EEG.data = double(EEG.data);
                vars = convertContainedStringsToChars(varin);
                EEG = pop_tesa_fastica( EEG, vars{:} );
                EEG = eeg_checkset( EEG );

            case 'Remove ICA Components (TESA)'
                vars = convertContainedStringsToChars(varin);
                compsIdx = find(strcmpi(vars, 'comps'), 1);
                if ~isempty(compsIdx)
                    vars([compsIdx, compsIdx+1]) = [];
                end
                ind1 = find(strcmpi(vars,'plotTimeX'));
                TP = vars{ind1+1};
                if TP(1) ~= EEG.times(1) && TP(2) ~= EEG.times(end)
                    vars{ind1+1} = [EEG.times(1) EEG.times(end)];
                end
                for nInd=2:2:numel(vars)
                    if size(vars{nInd},1)>size(vars{nInd},2)
                        vars{nInd} = vars{nInd}';
                    end
                end
                EEG = pop_tesa_compselect( EEG,vars{:});
                EEG = eeg_checkset( EEG );

            case 'Find Artifacts EDM (TESA)'
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
                EEG = pop_tesa_edm( EEG, chanlocations, nc, sf, vars{:});
                EEG = eeg_checkset( EEG );

            case 'SSP SIR'
                vars = convertContainedStringsToChars(varin);
                EEG = pop_tesa_sspsir(EEG, vars{:});
                EEG = eeg_checkset( EEG );

            case 'Remove Recording Noise (SOUND)'
                vars = convertContainedStringsToChars(varin);
                EEG = pop_tesa_sound(EEG, vars{:} );
                EEG = eeg_checkset( EEG );

            case 'Median Filter 1D'
                vars = convertContainedStringsToChars(varin);
                ind1 = find(strcmp(vars,'timeWin'));
                timeWin = vars{1,ind1+1};
                ind2 = find(strcmp(vars,'mdorder'));
                mdorder = vars{1,ind2+1};
                ind3 = find(strcmp(vars,'event_type'));
                event_type = vars{1,ind3+1};
                EEG = tesa_filtmedian( EEG, timeWin, mdorder, event_type );
                EEG = eeg_checkset( EEG );

            case 'Remove Bad Trials'
                localThresh  = step.params.localThresh;
                globalThresh = step.params.globalThresh;
                EEG = pop_jointprob(EEG, 1, 1:size(EEG.data,1), localThresh, globalThresh, 0, 0);
                pop_rejmenu(EEG, 1);
                uiconfirm(opts.uiFigure,'Highlight bad trials in the rejection menu, then press OK to continue.','Remove Bad Trials','Options',{'OK'},'DefaultOption',1);
                EEG.BadTr = unique([find(EEG.reject.rejjp==1) find(EEG.reject.rejmanual==1)]);
                EEG = pop_rejepoch( EEG, EEG.BadTr ,0);
                EEG = eeg_checkset( EEG );

            case 'Extract TEP (TESA)'
                vars = convertContainedStringsToChars(varin);
                ind1 = find(strcmp(vars,'type'));
                type = vars{1,ind1+1};
                vars([ind1, ind1+1]) = [];
                ind2 = find(strcmpi(vars,'pairCorrect'));
                if ~strcmp(vars{ind2+1},'on')
                    ind3 = find(strcmpi(vars,'ISI'));
                    vars([ind3, ind3+1]) = [];
                end
                vars = stripEmptyVarin(vars);
                EEG = pop_tesa_tepextract( EEG, type, vars );

            case 'Find TEP Peaks (TESA)'
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
                vars = stripEmptyVarin(vars);
                EEG = pop_tesa_peakanalysis( EEG, input, direction, peak, peakWin, vars(:) );

            case 'TEP Peak Output'
                vars = convertContainedStringsToChars(varin);
                vars = stripEmptyVarin(vars);
                pop_tesa_peakoutput( EEG, vars{:} );

            case 'Quality Gate'
                gate = qualityGate(EEG, step.params);
                gate.stepIndex = si;
                gate.stepName  = stepName;
                if ~isfield(fileReport, 'quality') ...
                        || ~isfield(fileReport.quality, 'gates')
                    if ~isfield(fileReport, 'quality')
                        fileReport.quality = struct( ...
                            'figures', {{}}, 'gates', {{}}, 'worstVerdict', 'Pass');
                    else
                        fileReport.quality.gates = {};
                        fileReport.quality.worstVerdict = 'Pass';
                    end
                end
                fileReport.quality.gates{end+1} = gate;
                fileReport.quality.worstVerdict = worseVerdict( ...
                    fileReport.quality.worstVerdict, gate.verdict);
                sendWorkerLog(opts.logQueue, wLabel, ...
                    'Quality Gate "%s" -> %s', gate.label, gate.verdict);
                % Stream the verdict back to the progress dialog so a
                % long parallel run shows live Pass/Marg/Fail chips per
                % file instead of only after the run completes.
                if ~isempty(opts.progressQueue)
                    send(opts.progressQueue, struct( ...
                        'fi',          opts.fileIndex, ...
                        'si',          si, ...
                        'nSteps',      nSteps, ...
                        'stepName',    stepName, ...
                        'gateVerdict', gate.verdict, ...
                        'gateLabel',   gate.label));
                end
                if strcmp(gate.verdict, 'Fail') && opts.skipOnQualityFail
                    error('nestapp:qualityFail', ...
                        'Step %d (Quality Gate "%s") failed: %s', ...
                        si, gate.label, strjoin(gate.reasons, '; '));
                end

        end % switch

        %% Post-step metrics and report update
        if isstruct(EEG) && ~isempty(EEG)
            nChanAfter  = EEG.nbchan;
            nEpochAfter = size(EEG.data, 3);
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
            if any(strcmp(stepName, {'Run ICA','Run TESA ICA'})) && ~isempty(EEG.icaweights)
                fileReport.ica.nComponents = size(EEG.icaweights, 1);
            end
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
            if strcmp(stepName, 'Remove ICA Components (TESA)') && ...
                    isfield(EEG, 'icaCompClass') && isstruct(EEG.icaCompClass) && ...
                    ~isempty(fieldnames(EEG.icaCompClass))
                tesaKeys = fieldnames(EEG.icaCompClass);
                cl = EEG.icaCompClass.(tesaKeys{end});
                TESA_CATS  = tesaICACategories();
                TESA_CODES = tesaICAClassCodes();
                rejIdx   = cl.compClass > 1;
                nRejTESA = sum(rejIdx);
                rnd.roundNum    = numel(tesaKeys);
                rnd.nComponents = numel(cl.compClass);
                rnd.nRejected   = nRejTESA;
                rnd.varRemoved  = NaN;
                rnd.varMin      = NaN;
                rnd.varMax      = NaN;
                rnd.categories.names    = TESA_CATS;
                rnd.categories.nRemoved = zeros(1, numel(TESA_CATS));
                rnd.categories.varShare = zeros(1, numel(TESA_CATS));
                hasVars = isfield(cl, 'compVars') && numel(cl.compVars) >= numel(cl.compClass);
                if hasVars && nRejTESA > 0
                    rejPct         = double(cl.compVars(rejIdx));
                    rnd.varRemoved = sum(rejPct);
                    rnd.varMin     = min(rejPct);
                    rnd.varMax     = max(rejPct);
                end
                for ci = 1:numel(TESA_CODES)
                    inCat = (cl.compClass == TESA_CODES(ci));
                    rnd.categories.nRemoved(ci) = sum(inCat);
                    if hasVars
                        rnd.categories.varShare(ci) = sum(cl.compVars(inCat));
                    end
                end
                fileReport.ica.rounds{end+1} = rnd;
                fileReport.ica.nRejected = fileReport.ica.nRejected + nRejTESA;
                fileReport.ica.nKept     = fileReport.ica.nComponents - fileReport.ica.nRejected;
                if ~strcmp(fileReport.ica.categories.names{1}, 'TMS Muscle')
                    fileReport.ica.categories.names    = TESA_CATS;
                    fileReport.ica.categories.nRemoved = zeros(1, numel(TESA_CATS));
                    fileReport.ica.categories.varShare = zeros(1, numel(TESA_CATS));
                end
                fileReport.ica.categories.nRemoved = ...
                    fileReport.ica.categories.nRemoved + rnd.categories.nRemoved;
                if isscalar(fileReport.ica.rounds)
                    fileReport.ica.varRemoved = rnd.varRemoved;
                    fileReport.ica.varMin     = rnd.varMin;
                    fileReport.ica.varMax     = rnd.varMax;
                    fileReport.ica.categories.varShare = rnd.categories.varShare;
                end
            end
        else
            nChanAfter  = nChanBefore;
            nEpochAfter = nEpochBefore;
        end % if isstruct(EEG)

        stepRec.name         = stepName;
        stepRec.chansBefore  = nChanBefore;
        stepRec.chansAfter   = nChanAfter;
        stepRec.trialsBefore = nEpochBefore;
        stepRec.trialsAfter  = nEpochAfter;
        elapsed              = toc(t0);
        stepRec.duration     = elapsed;
        stepRec.timestamp    = datetime('now');
        fileReport.steps{end+1} = stepRec;

        sendWorkerLog(opts.logQueue, wLabel, 'Step %d/%d END    "%s"  (%.2fs)  ch:%d->%d  ep:%d->%d', ...
            si, nSteps, stepName, elapsed, nChanBefore, nChanAfter, nEpochBefore, nEpochAfter);

        stepLog(end+1) = struct( ...
            'step',        stepName, ...
            'duration_s',  elapsed, ...
            'chanBefore',  nChanBefore, ...
            'chanAfter',   nChanAfter, ...
            'epochBefore', nEpochBefore, ...
            'epochAfter',  nEpochAfter, ...
            'error',       ''); %#ok<AGROW>

        % Auto Quality Report: render a per-(file, gate) QC PNG after
        % every successful Quality Gate. Wrapped in try/catch so a
        % rendering bug never aborts the actual data run.
        if ~isempty(opts.batchCtx) && opts.autoQualityReport ...
                && strcmp(stepName, 'Quality Gate')
            pngName = sprintf('%02d_%s.png', si, sanitizeForPath(stepName));
            outPath = fullfile(outputPaths(opts.batchCtx, 'qc', fileBase), pngName);
            try
                if opts.qcTmsAutoDetect
                    tmsWin = inferTmsWindow(EEG, opts.qcTmsWindow);
                else
                    tmsWin = opts.qcTmsWindow;
                end
                qcOpts = struct( ...
                    'title',     fileBase, ...
                    'stepLabel', sprintf('Step %d / %s', si, stepName), ...
                    'attribute', opts.qcAttribute, ...
                    'tmsWindow', tmsWin, ...
                    'size',      [1600 1200]);
                qcOpts.panels = struct('attribMatrix',true,'icaGrid',true, ...
                                       'butterfly',true,'psd',true);
                renderQualityFigure(EEG, outPath, qcOpts);
                if ~isfield(fileReport, 'quality') ...
                        || ~isfield(fileReport.quality, 'figures')
                    fileReport.quality = struct('figures', {{}});
                end
                fileReport.quality.figures{end+1} = outPath;
                sendWorkerLog(opts.logQueue, wLabel, 'QC %s', outPath);
            catch qcErr
                sendWorkerLog(opts.logQueue, wLabel, ...
                    'QC render FAILED (continuing): %s', qcErr.message);
            end
        end

    catch err
        elapsed = toc(t0);
        sendWorkerLog(opts.logQueue, wLabel, 'Step %d/%d ERROR  "%s"  (%.2fs)  %s', ...
            si, nSteps, stepName, elapsed, err.message);

        stepLog(end+1) = struct( ...
            'step',        stepName, ...
            'duration_s',  elapsed, ...
            'chanBefore',  nChanBefore, ...
            'chanAfter',   nChanBefore, ...
            'epochBefore', nEpochBefore, ...
            'epochAfter',  nEpochBefore, ...
            'error',       err.message); %#ok<AGROW>

        % Quality Gate fails are user-defined skip points, not pipeline
        % errors. Don't prompt the user (skipping IS the chosen behavior)
        % and don't re-wrap the identifier; runPipelineCore.parseFailure
        % needs the original 'nestapp:qualityFail' message shape to tag
        % the file as 'skipped' rather than 'errored'.
        if strcmp(err.identifier, 'nestapp:qualityFail')
            if ~isempty(opts.progressQueue)
                send(opts.progressQueue, struct( ...
                    'fi', opts.fileIndex, 'si', 0, ...
                    'nSteps', nSteps, 'stepName', 'Skipped', 'failed', true));
            end
            rethrow(err);
        end

        warning('An error occurred at file %s at step %d: %s', fileName, si, stepName);

        shouldContinue = false;
        if ~isempty(opts.onStepError)
            choice = opts.onStepError(si, stepName, err);
            shouldContinue = strcmp(choice, 'Continue');
        end

        if ~shouldContinue
            if ~isempty(opts.onStepError)
                % Serial mode: user chose Abort on the step prompt - cancel whole run.
                error('nestapp:cancelled', 'Run aborted at step %d (%s): %s', ...
                    si, stepName, err.message);
            else
                % Parallel mode (no prompt): mark this file failed with
                % structured step info so runPipelineCore can summarize.
                % Release the dialog slot first - the success sentinel at the
                % bottom of this function is unreachable from here, and without
                % a release later files' progress messages get dropped.
                if ~isempty(opts.progressQueue)
                    send(opts.progressQueue, struct( ...
                        'fi', opts.fileIndex, 'si', 0, ...
                        'nSteps', nSteps, 'stepName', 'Failed', 'failed', true));
                end
                error('nestapp:stepFailed', 'Step %d (%s) failed: %s', ...
                    si, stepName, err.message);
            end
        end
    end
end % step loop

%% Write pipeline provenance to EEG.history
if isstruct(EEG) && isfield(EEG, 'history')
    EEG.history = [EEG.history, newline, buildHistoryEntry(spec, opts.pipelineName)];
    assignin('base', 'EEG', EEG);
    newHist  = EEG.history(histLenBefore + 1 : end);
    newLines = strtrim(strsplit(newHist, newline));
    newLines = newLines(~cellfun('isempty', newLines));
    ALLCOM = [newLines(:)', ALLCOM];
end

% Auto-export per-file PDF when the pref is on. Failure is logged
% but never aborts the pipeline.
if opts.autoExportPDF
    try
        if ~isempty(opts.batchCtx)
            pdfPath = exportFileReportPDF(fileReport, opts.batchCtx);
        else
            pdfPath = exportFileReportPDF(fileReport, pathName);
        end
        sendWorkerLog(opts.logQueue, wLabel, 'PDF report written: %s', pdfPath);
    catch pdfErr
        sendWorkerLog(opts.logQueue, wLabel, ...
            'PDF export FAILED (continuing): %s', pdfErr.message);
    end
end

sendWorkerLog(opts.logQueue, wLabel, 'DONE   %s  (total %.2fs)', fileName, toc(fileTic));

% Sentinel: si=0 signals the file is truly done (after all cleanup).
% updateParallelProgress uses this - not si==nSteps - to go green, so the
% bar stays blue until the worker has actually finished, not just started the
% last step.
if ~isempty(opts.progressQueue)
    sendWorkerLog(opts.logQueue, wLabel, 'Sending done sentinel to DataQueue');
    send(opts.progressQueue, struct( ...
        'fi', opts.fileIndex, 'si', 0, ...
        'nSteps', nSteps, 'stepName', 'Done'));
end
end

% -- local helpers ---------------------------------------------------------

function cats = tesaICACategories()
% TESA component category names, in the order TESA's icaCompClass.compClass codes map to.
cats = {'TMS Muscle','Blink','Eye Move','Muscle','Elec Noise','Sensory','Reject'};
end

function codes = tesaICAClassCodes()
% TESA compClass integer codes, matched positionally to tesaICACategories().
codes = [3, 4, 5, 6, 7, 8, 2];
end

function v = worseVerdict(a, b)
% Return the worst severity of two verdicts: Fail > Marginal > Pass > NotChecked.
order = {'NotChecked', 'Pass', 'Marginal', 'Fail', 'Pending'};
ia = find(strcmp(a, order));
ib = find(strcmp(b, order));
if isempty(ia), ia = 0; end
if isempty(ib), ib = 0; end
v = order{max(ia, ib)};
end

function sendWorkerLog(q, label, fmt, varargin)
% Route timestamped log lines through the DataQueue so they print in real-time
% on the client. PCT buffers plain fprintf until the future completes.
if isempty(q)
    nestLog(label, fmt, varargin{:});
    return
end
ts = char(datetime('now', 'Format', 'HH:mm:ss.SSS'));
send(q, struct('log', true, 'ts', ts, 'label', label, ...
               'text', sprintf(fmt, varargin{:})));
end
