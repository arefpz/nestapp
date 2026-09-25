
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
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
if ~isfield(opts, 'aaratepKeepIntermediates'), opts.aaratepKeepIntermediates = true; end
if ~isfield(opts, 'qcAttribute'),       opts.qcAttribute       = 'minmax_no_tms'; end
if ~isfield(opts, 'qcTmsWindow'),       opts.qcTmsWindow       = [0 25];          end
if ~isfield(opts, 'qcTmsAutoDetect'),   opts.qcTmsAutoDetect   = true;            end
if ~isfield(opts, 'skipOnQualityFail'), opts.skipOnQualityFail = false;           end
if ~isfield(opts, 'autoExportPDF'),     opts.autoExportPDF     = false;           end
if ~isfield(opts, 'saveErrorBundle'),   opts.saveErrorBundle   = false;           end

% eeglab('nogui') is expensive (plugin scan, path setup); run it once per
% worker then just reset globals for subsequent files on the same worker.
persistent eeglabWorkerReady lastThreadCount
global EEG ALLEEG CURRENTSET ALLCOM %#ok<GVMIS>
% Channel for harvesting a manual selection out of an eegplot callback.
% eegplot evaluates its 'command' string in the base workspace, so a global
% is the only way the marks reach this scope - see Remove Bad Trials.
global NESTAPP_TMPREJ %#ok<GVMIS>

[pathDir, fileBase, fileExt] = fileparts(fullPath);
pathName = [pathDir, filesep];
fileName = [fileBase, fileExt];

wLabel   = sprintf('FILE-%d %s', opts.fileIndex, fileBase);
fileTic  = tic;

% Give each worker its fair share of BLAS threads. Process workers START at
% maxNumCompThreads == 1 (the Processes profile's NumThreads default is 1), so
% this RAISES the allowance rather than capping it - there is no default
% over-subscription to prevent. floor(numcores/nWorkers) can never exceed the
% core budget in aggregate, and the raise is worth having: 1 -> 4 threads took a
% 3000x3000 matmul from ~1.2s to ~0.35s per worker (R2026a, 16 cores).
% NB feature('numcores') reports the machine's cores even on a worker, not the
% worker's own allowance, which is why the share is computed against it.
% Done here, per worker, rather than by setting the cluster's NumThreads: that
% property is saved profile state, and a pipeline run should not edit the user's
% MATLAB configuration as a side effect.
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
latestICASnapshot = struct([]);   % most recent pre-rejection ICA state for QC
histLenBefore     = 0;

nSteps = numel(spec);

sendWorkerLog(opts.logQueue, wLabel, 'START  %s  (%d steps)', fileName, nSteps);

% Stagger the first worker wave to avoid simultaneous disk reads at batch start.
% Files beyond the first wave are already staggered naturally by worker completion times.
if ~isempty(opts.progressQueue) && opts.fileIndex > 1 && opts.fileIndex <= opts.nWorkers
    pause(0.25 * (opts.fileIndex - 1));
end

% Files that may be removed once this file completes - see the AARATEP case.
pendingDeletes = {};
stepLog = struct('step',{},'duration_s',{},'chanBefore',{},'chanAfter',{}, ...
                 'epochBefore',{},'epochAfter',{},'error',{});
fileReport = initPipelineReport(fullPath);
fileReport.pipelineName = opts.pipelineName;  % provenance (citations come from steps)

for si = 1:nSteps
    step     = spec(si);
    stepName = canonicalStepName(step.name);   % migrate legacy step names before dispatch
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
        if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
            labelsBefore = {EEG.chanlocs.labels};
        else
            labelsBefore = {};
        end
    else
        nChanBefore  = 0;
        nEpochBefore = 0;
        labelsBefore = {};
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
                coordCheck = 'off';
                ci = find(strcmpi(vars,'coordCheck'),1);
                if ~isempty(ci), coordCheck = vars{ci+1}; end

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
                    % 'load' must come BEFORE 'lookup': pop_chanedit processes
                    % its args in order, and the 'load' case replaces the whole
                    % chanlocs struct (and resets chaninfo). With 'lookup' first
                    % its work was discarded, so the standard_1005 coordinates
                    % were never actually merged. Load the file, THEN look up.
                    EEG = pop_chanedit(EEG, ...
                        'load', {[chPath, chBase, chExt], 'filetype', 'autodetect'}, ...
                        'lookup', lookforchnlocs);
                    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
                end
                % Validate coordinates (after any lookup): invalid coords
                % (missing/NaN/0,0,0) silently break the downstream spatial
                % RANSAC bad-channel and interpolation steps.
                EEG = validateChannelCoords(EEG, coordCheck);

            case 'Load Data'
                EEG = loadEegFile([pathName fileName]);
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

                % Remember where it landed. pop_saveset stamps EEG.filename /
                % .filepath, which is authoritative because it reflects the
                % extension and any name normalisation EEGLAB applied; the
                % composed name is the fallback for a save that did not go
                % through it. Without this the report describes a file it
                % cannot point at.
                fileReport.outputFile = savedSetPath(EEG, vars);

            case 'Find TMS Pulses (AARATEP)'
                ensureAaratepOnPath();
                fp = varinToStruct(varin);
                % Upstream asserts on this itself, but from here the message
                % names the step and the fix.
                if size(EEG.data, 3) > 1
                    error('nestapp:aaratepFindPulsesNeedsContinuous', ...
                        ['Find TMS Pulses (AARATEP) detects pulses in continuous ' ...
                         'data. Run it before Epoching.']);
                end
                % Upstream builds its new events with
                % c_struct_createEmptyCopy(EEG.event), which asserts isstruct.
                % EEGLAB leaves EEG.event as [] when a recording carries no
                % events - which is precisely the case this detector exists
                % for - so seed the struct it expects. Cannot be fixed in the
                % vendored copy, so it is prepared here.
                if ~isstruct(EEG.event)
                    EEG.event = struct('type', {}, 'latency', {}, ...
                                       'duration', {}, 'urevent', {});
                end
                % minNumPulses is left at its default of 0 on purpose - see the
                % registry entry; any other value can reach a `keyboard`.
                EEG = c_TMSEEG_findTMSPulses(EEG, ...
                    'addEventsOfType',    fp.addEventsOfType, ...
                    'filterMethod',       fp.filterMethod, ...
                    'maxPulseRate',       fp.maxPulseRate, ...
                    'minPulseThreshold',  fp.minPulseThreshold, ...
                    'maxPulseThreshold',  fp.maxPulseThreshold, ...
                    'doRefineOnsetTimes', strcmpi(fp.doRefineOnsetTimes, 'on'));
                EEG = eeg_checkset(EEG, 'eventconsistency');

            case 'AARATEP Pipeline (whole)'
                ensureAaratepOnPath();
                o = varinToStruct(varin);

                % Upstream expects CONTINUOUS data with pulse events already
                % present - it epochs itself. Fail here with something the
                % user can act on rather than deep inside the orchestrator.
                if size(EEG.data, 3) > 1
                    error('nestapp:aaratepNeedsContinuous', ...
                        ['AARATEP Pipeline (whole) does its own epoching, so it ' ...
                         'needs continuous data. Remove the Epoching step - the ' ...
                         'pipeline epochs using its own epochTimespan.']);
                end
                if isempty(o.pulseEvents) || all(cellfun(@isempty, cellstr(o.pulseEvents)))
                    error('nestapp:aaratepNoPulseEvents', ...
                        ['AARATEP Pipeline (whole) needs the event type marking ' ...
                         'each TMS pulse. Run "Find TMS Pulses (TESA)" first and ' ...
                         'set Pulse event type(s) to match.']);
                end
                % Left blank, the output folder comes from the same output
                % root every other step writes under, in a per-file subfolder.
                % Per-file is not cosmetic: upstream MOVES an existing output
                % folder to <folder>_old# before writing, so one shared folder
                % would have each file displace the previous file's results.
                if isempty(o.outputDir)
                    o.outputDir = aaratepOutputDir(opts.batchCtx, fullPath);
                end
                if isempty(o.outputDir)
                    error('nestapp:aaratepNoOutputDir', ...
                        ['AARATEP Pipeline (whole) writes its results to a folder. ' ...
                         'Normally this comes from the output root, but this run ' ...
                         'has no batch context, so set Output folder explicitly.']);
                end

                % Passed through in upstream's own units and defaults. Note
                % maximizePlotsToMonitor is deliberately absent - its
                % validator (@isschar) is not a function, so setting it at all
                % raises an undefined-function error.
                EEG = c_TMSEEG_Preprocess_AARATEPPipeline(EEG, ...
                    'pulseEvents',                  cellstr(o.pulseEvents), ...
                    'outputDir',                    o.outputDir, ...
                    'epochTimespan',                o.epochTimespan, ...
                    'outputFilePrefix',             o.outputFilePrefix, ...
                    'artifactTimespan',             o.artifactTimespan, ...
                    'baselineTimespan',             o.baselineTimespan, ...
                    'filterPrePostExtrapolationDurations', o.filterPrePostExtrapolationDurations, ...
                    'downsampleTo',                 o.downsampleTo, ...
                    'bandpassFreqSpan',             o.bandpassFreqSpan, ...
                    'badChannelDetectionMethod',    cellstr(o.badChannelDetectionMethod), ...
                    'badChannelThreshold',          o.badChannelThreshold, ...
                    'initialEyeComponentThreshold', o.initialEyeComponentThreshold, ...
                    'SOUNDlambda',                  o.SOUNDlambda, ...
                    'leadFieldPath',                o.leadFieldPath, ...
                    'doDecayRemovalPerTrial',       strcmpi(o.doDecayRemovalPerTrial, 'on'), ...
                    'lineNoiseFreq',                o.lineNoiseFreq, ...
                    'lineNoiseNumHarmonics',        o.lineNoiseNumHarmonics, ...
                    'ICAType',                      o.ICAType, ...
                    'brainComponentThreshold',      o.brainComponentThreshold, ...
                    'stimMuscleComponentThreshold', o.stimMuscleComponentThreshold, ...
                    'onOverRejection',              o.onOverRejection, ...
                    'doPostICAArtifactInterpolation', strcmpi(o.doPostICAArtifactInterpolation, 'on'), ...
                    'doDebug',                      strcmpi(o.doDebug, 'on'), ...
                    'doPlotFinalTimtopo',           strcmpi(o.doPlotFinalTimtopo, 'on'), ...
                    'plotXLim',                     o.plotXLim, ...
                    'plotTPOIs',                    o.plotTPOIs, ...
                    'plotChans',                    cellstr(o.plotChans));
                EEG = eeg_checkset( EEG );

                % The orchestrator hands back only the cleaned EEG - its
                % provenance, QC images and intermediate datasets are files in
                % o.outputDir. Fold them into the report, and drop its final
                % .mat when a later Save New Set writes the same dataset as a
                % .set, so the result exists once, in nestapp's own format.
                hOpts                   = struct();
                hOpts.keepIntermediates = opts.aaratepKeepIntermediates;
                hOpts.channelLabels     = channelLabelsOf(EEG);
                [fileReport, hInfo] = aaratepHarvest( ...
                    fileReport, o.outputDir, o.outputFilePrefix, hOpts);
                % Its final .mat and Save New Set's .set hold the same
                % dataset, so only one should survive - but the .set does not
                % exist yet. Deleting here would act on a prediction: a failing
                % Quality Gate or a throwing step between here and the save
                % would leave the result on disk nowhere at all. Queue it and
                % delete once the file has actually come through.
                if ~isempty(hInfo.droppableFinalMat) && savesNewSetLater(spec, si)
                    pendingDeletes{end+1} = hInfo.droppableFinalMat; %#ok<AGROW>
                end
                mdTxt = 'NOT FOUND';
                if hInfo.mdFound; mdTxt = 'read'; end
                sendWorkerLog(opts.logQueue, wLabel, ...
                    'AARATEP: metadata %s, %d QC image(s), %d file(s) removed (%.0f MB)', ...
                    mdTxt, numel(hInfo.figures), numel(hInfo.removed), ...
                    hInfo.bytesFreed / 1e6);

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
                % An unset index reaches pop_newset as '' or NaN and dies inside
                % finputcheck with "argument 'retrieve' must be numeric" - a
                % message that names neither this step nor the fix. Validate up
                % front. (Default is [] = unset; a GUI-cleared cell yields NaN.)
                if ~(isnumeric(setIndex) && isscalar(setIndex) ...
                        && isfinite(setIndex) && setIndex >= 1 ...
                        && setIndex == floor(setIndex))
                    error('nestapp:chooseDataSetNoIndex', ...
                        ['Choose Data Set needs a dataset index (a positive ' ...
                         'integer). Set "Dataset index" on the step to the ' ...
                         'ALLEEG slot you want to switch to.']);
                end
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

            case 'Remove Bad Channels (manual)'
                % pop_topochansel draws only the channels that HAVE
                % coordinates and returns indices into that reduced set. If
                % any channel lacks coordinates the returned indices no longer
                % line up with EEG.chanlocs, and pop_select would remove the
                % wrong channels - so require EVERY channel to have a location,
                % not just some.
                if ~isfield(EEG, 'chanlocs') || isempty(EEG.chanlocs) || ...
                        ~isfield(EEG.chanlocs, 'X') || any(cellfun(@isempty, {EEG.chanlocs.X}))
                    error('nestapp:noChanlocsForPicker', ...
                        ['Remove Bad Channels (manual) draws the electrodes on ' ...
                         'the scalp, so it needs a location for EVERY channel. ' ...
                         'Run "Load Channel Location" first, or remove the ' ...
                         'channels that have no coordinates.']);
                end
                % pop_topochansel RETURNS the selection - no base-workspace
                % side effects to plumb - so the user's choice reaches the
                % pipeline by construction. Pass chanlocs explicitly: called
                % bare it scrapes the caller's workspace for them.
                o = varinToStruct(varin);
                badIdx = pop_topochansel(EEG.chanlocs);
                badIdx = reshape(badIdx, 1, []);

                protect = {};
                if isfield(o, 'protect') && ~isempty(o.protect); protect = o.protect; end
                protect = protect(~cellfun(@isempty, protect));
                if ~isempty(protect) && ~isempty(badIdx)
                    labels  = {EEG.chanlocs.labels};
                    keepIdx = badIdx(ismember(lower(labels(badIdx)), lower(protect)));
                    if ~isempty(keepIdx)
                        nestLog('CHAN', 'Manual: keeping %d protected channel(s): %s', ...
                            numel(keepIdx), strjoin(labels(keepIdx), ', '));
                        badIdx = setdiff(badIdx, keepIdx);
                    end
                end

                if isempty(badIdx)
                    nestLog('CHAN', 'Manual: no channels selected; none removed.');
                else
                    badBefore = {EEG.chanlocs.labels};
                    nestLog('CHAN', 'Manual: removing %d channel(s): %s', ...
                        numel(badIdx), strjoin(badBefore(badIdx), ', '));
                    EEG = pop_select(EEG, 'nochannel', badIdx);
                    EEG = eeg_checkset( EEG );
                    EEG = recordBadChannels(EEG, badBefore);
                end

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

                badBefore = {EEG.chanlocs.labels};
                if sum(importantElects)
                    vars{1,ind2+1} = find(~importantElects);
                    EEG = pop_rejchan(EEG, vars{:});
                else
                    EEG = pop_rejchan(EEG, vars{:});
                end
                % Record what was rejected so "Interpolate Channels" restores
                % only bad channels, not intentionally-removed ones.
                EEG = recordBadChannels(EEG, badBefore);

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
                badBefore = {EEG.chanlocs.labels};
                EEG = clean_artifacts(EEG,vars{:});
                EEG = eeg_checkset( EEG );
                EEG = recordBadChannels(EEG, badBefore);

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
                badBefore = {EEG.chanlocs.labels};
                EEG = pop_clean_rawdata(EEG, vars{:});
                EEG = eeg_checkset( EEG );
                EEG = recordBadChannels(EEG, badBefore);

            case 'Remove Baseline'
                vars = convertContainedStringsToChars(varin);
                % "Time range" (ms) and "Point range" (samples) are alternative
                % ways to specify the baseline window; timerange wins if both
                % are set. "Channels" restricts which channels are corrected.
                % pop_rmbase(EEG, timerange, pointrange, chanlist).
                timerange  = parseRmbaseVec(vars, 'timerange');
                pointrange = parseRmbaseVec(vars, 'pointrange');
                chanlist   = parseRmbaseVec(vars, 'chanlist');
                if isnumeric(timerange) && numel(timerange) == 2
                    timerange(1) = max(timerange(1), EEG.times(1));
                    timerange(2) = min(timerange(2), EEG.times(end));
                end
                EEG = pop_rmbase(EEG, timerange, pointrange, chanlist);
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
                % The return value MUST be assigned. pop_tesa_detrend returns a
                % new struct; calling it bare detrended a copy and discarded it,
                % so this step was a complete no-op that still logged "Linear
                % detrend complete." See test_stepCharacterization.
                %
                % The exponential and double fits need two toolboxes, and
                % tesa_detrend crashes cryptically when either is absent - a
                % message naming neither the real dependency nor this step.
                % Guard up front and fail with something the user can act on.
                % Linear uses only polyfit and needs no toolbox.
                if any(strcmpi(Tdetrend, {'exponential', 'double'}))
                    % fit() does the actual exponential fit.
                    [fitAvailable, fitReason] = ensureCurveFittingFit();
                    if ~fitAvailable
                        error('nestapp:detrendNeedsCurveFitting', ...
                            ['TESA De-Trend "%s" needs the Curve Fitting Toolbox. %s ' ...
                             'Use "linear" instead, or install/license the toolbox.'], ...
                            Tdetrend, fitReason);
                    end
                    % tesa_detrend detects the Curve Fitting Toolbox with
                    % extractfield - itself a Mapping Toolbox function - so
                    % without Mapping the step dies in its own check
                    % ("Undefined function 'extractfield'") even when the fit
                    % it wants to run is available. That is a limitation of
                    % TESA's check, not of the fit; name it plainly.
                    if isempty(which('extractfield'))
                        error('nestapp:detrendNeedsMappingForCheck', ...
                            ['TESA De-Trend "%s" cannot run here: tesa_detrend checks ' ...
                             'for the Curve Fitting Toolbox using extractfield (Mapping ' ...
                             'Toolbox), which is not installed. This is a quirk of TESA''s ' ...
                             'own toolbox check, not of the fit. Use "linear", or install ' ...
                             'the Mapping Toolbox.'], Tdetrend);
                    end
                end
                EEG = pop_tesa_detrend(EEG, Tdetrend, TtimeWin);
                EEG = eeg_checkset( EEG );

            case 'Re-Sample'
                vars = convertContainedStringsToChars(varin);
                EEG = pop_resample(EEG,vars{2:2:end});
                EEG = eeg_checkset( EEG );

            case 'Re-Reference'
                vars = convertContainedStringsToChars(varin);
                ind = find(strcmp(vars,'ref'));
                ref = vars{ind+1};
                % Average reference (ref = '[]') needs no channel labels.
                % A named reference does, so fail with a legible message if
                % chanlocs are missing rather than dereferencing a double.
                hasLocs = isstruct(EEG.chanlocs) && ~isempty(EEG.chanlocs) ...
                    && isfield(EEG.chanlocs,'labels');
                if ~strcmp(ref,'[]') && ~hasLocs
                    error('nestapp:noChanlocs', ...
                        ['Re-Reference: channel locations are missing, so ' ...
                         'reference ''%s'' cannot be resolved. An earlier ' ...
                         'step likely dropped EEG.chanlocs.'], ref);
                end
                if hasLocs && ~ismember(ref,{EEG.chanlocs.labels}) && ~strcmp(ref,'[]')
                    % Only prompt when there is a UI to prompt on. On a
                    % parallel worker uiFigure is [] and there is no display,
                    % so an inputdlg here would block a batch that is already
                    % running. Fail with the information the user would have
                    % needed to answer it instead.
                    if isempty(opts.uiFigure)
                        error('nestapp:refChannelMissing', ...
                            ['Re-Reference: reference channel ''%s'' is not in ' ...
                             'the data. Available channels: %s. Set a valid ' ...
                             'reference on the Re-Reference step and re-run.'], ...
                            ref, strjoin({EEG.chanlocs.labels}, ', '));
                    end
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
                % interpchan: pop_reref interpolates the removed reference
                % channels when this is [] (its own GUI passes exactly that);
                % passing the string 'on' instead reaches isreal('on')==true and
                % indexes EEG.urchanlocs('on'), which crashes. Translate the
                % on/off toggle to what pop_reref expects, and pull it out
                % before stripEmptyVarin so the [] enable value is not deleted.
                indInterp = find(strcmp(vars, 'interpchan'));
                interpOn = false;
                if ~isempty(indInterp)
                    interpOn = strcmpi(vars{indInterp+1}, 'on');
                    vars([indInterp, indInterp+1]) = [];
                end
                vars = stripEmptyVarin(vars);
                if interpOn
                    vars = [vars, {'interpchan', []}];
                end
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
                % Expand the [first last] range to explicit indices, clamped to
                % the channels present. The old inline arithmetic dropped the
                % last channel on small montages and errored (1:[1 64]) on >=64.
                vars{ind+1} = cleanlineChanList(vars{ind+1}, EEG.nbchan);
                EEG = pop_cleanline(EEG, vars{:});
                EEG = eeg_checkset( EEG );

            % ---- TESA 1.2 -------------------------------------------
            % Reachable only when the installed TESA declares 1.2 or later;
            % stepAvailability hides these from the picker and blocks a run
            % that references one otherwise. Argument shapes are upstream's:
            % robustdetrend / robustdemean take positional arguments, the
            % rest take name/value pairs. nestapp exposes times in ms and
            % converts where upstream wants seconds.
            case 'Robust Detrend (TESA)'
                o = varinToStruct(varin);
                EEG = pop_tesa_robustdetrend(EEG, o.timeWindow, ...
                    o.thresholdSD, o.polyOrder, o.excludeWindow);
                EEG = eeg_checkset( EEG );

            case 'Robust Demean (TESA)'
                o = varinToStruct(varin);
                EEG = pop_tesa_robustdemean(EEG, o.timeWindow, ...
                    o.thresholdSD, o.excludeWindow);
                EEG = eeg_checkset( EEG );

            case 'Modified Bandpass Filter (TESA)'
                o = varinToStruct(varin);
                % The Butterworth design requires an even order (it fails deep
                % in an inputParser with @(x)mod(x,2)==0). Catch it here with a
                % message that names the step, as the sibling Frequency Filter
                % (TESA) does, rather than letting the cryptic one surface.
                if ~isempty(o.filtOrder) && mod(o.filtOrder, 2) ~= 0
                    error('nestapp:modBandpassOddOrder', ...
                        ['Modified Bandpass Filter (TESA): filter order must be ' ...
                         'an even number (got %g). Use 4, 6, ...'], o.filtOrder);
                end
                % artifactTimespan and pieceWiseTimeToExtend are SECONDS
                % upstream; the step exposes ms like every other nestapp step.
                EEG = pop_tesa_modifiedbandpassfilter(EEG, ...
                    'lowCutoff',             o.lowCutoff, ...
                    'highCutoff',            o.highCutoff, ...
                    'filterMethod',          o.filterMethod, ...
                    'artifactTimespan',      [o.artifactStartMs, o.artifactEndMs] * 1e-3, ...
                    'pieceWiseTimeToExtend', o.extendMs * 1e-3, ...
                    'filtOrder',             o.filtOrder, ...
                    'filtType',              o.filtType);
                EEG = eeg_checkset( EEG );

            case 'Detect Bad Channels (TESA)'
                o = varinToStruct(varin);
                % artifactTimespan is in ms for this one - upstream differs
                % between these two functions, so do not factor the conversion.
                dbcArgs = { ...
                    'detectionMethod',  o.detectionMethod, ...
                    'replaceMethod',    o.replaceMethod, ...
                    'artifactTimespan', [o.artifactStartMs, o.artifactEndMs]};
                % Pass threshold ONLY when the user set one. Empty lets upstream
                % apply its method-specific default (9 PREP / 20 DDWiener), and
                % fromASR asserts the threshold is empty - so forcing a number
                % broke DDWiener silently and made fromASR throw.
                if isfield(o, 'threshold') && ~isempty(o.threshold)
                    dbcArgs = [dbcArgs, {'threshold', o.threshold}];
                end
                EEG = pop_tesa_detectbadchannels(EEG, dbcArgs{:});
                EEG = eeg_checkset( EEG );

            case 'Detect Bad Channels (RANSAC)'
                o = varinToStruct(varin);
                % clean_channels correlates channels over time windows, so it
                % is a continuous-data operation - not meaningful on epochs.
                if size(EEG.data, 3) > 1
                    error('nestapp:ransacNeedsContinuous', ...
                        ['Detect Bad Channels (RANSAC) works on CONTINUOUS ' ...
                         'data (it correlates channels over time windows). ' ...
                         'Place it before Epoching.']);
                end
                % clean_channels itself requires channel coordinates and errors
                % clearly if most are missing, so no extra guard here.
                %
                % RANSAC correlation is wrecked by slow drift, so we DETECT on a
                % throwaway high-passed copy but REMOVE the flagged channels from
                % the untouched data - the kept signal is never filtered by this
                % step (deferring the real high-pass past decay removal). This
                % mirrors how EEGLAB's own clean_artifacts conditions the data
                % (clean_drifts -> clean_channels) before detecting.
                labelsBefore = {EEG.chanlocs.labels};
                EEGdetect = EEG;
                if o.detectHighpassHz > 0
                    % pop_eegfiltnew (firfilt, core EEGLAB) FIR high-pass; the
                    % copy may ring over the residual TMS decay, but that ringing
                    % is shared across the montage and does not make good
                    % channels look decorrelated - it only informs which
                    % channels are bad.
                    EEGdetect = pop_eegfiltnew(EEGdetect, 'locutoff', o.detectHighpassHz);
                end
                EEGdetect = clean_channels(EEGdetect, o.corrThreshold, ...
                    o.noiseThreshold, o.windowLenS, o.maxBrokenTime, ...
                    o.numSamples, o.subsetFraction);
                % Channels that did not survive detection, as indices into the
                % still-untouched EEG (its channel order matches labelsBefore).
                badIdx = find(~ismember(labelsBefore, {EEGdetect.chanlocs.labels}));
                clear EEGdetect   % free the large filtered copy before pop_select reallocates
                if ~isempty(badIdx)
                    EEG = pop_select(EEG, 'nochannel', badIdx);
                end
                % Record the removals as bad so "Interpolate Channels" can
                % restore them (vs. intentionally-removed channels).
                EEG = recordBadChannels(EEG, labelsBefore);
                EEG = eeg_checkset( EEG );

            case 'Fit Artifact Model (TESA)'
                o = varinToStruct(varin);
                EEG = pop_tesa_fitartifactmodel(EEG, ...
                    'tPulseMs',          o.tPulseMs, ...
                    'tPulseDurationMs',  o.tPulseDurationMs, ...
                    'tSkipMs',           o.tSkipMs, ...
                    'tSelectMs',         o.tSelectMs, ...
                    't_fitmodel_p_ms',   o.tFitPosMs, ...
                    't_fitmodel_n_ms',   o.tFitNegMs, ...
                    't_lincomb_ms',      o.tLinCombMs, ...
                    'epoch_range',       o.epochRangeMs);
                EEG = eeg_checkset( EEG );

            case 'Interactive Channel Reject (TESA)'
                % Takes no options and returns the data with the user's
                % selection applied - it blocks on uiwait until they choose.
                badBefore = {EEG.chanlocs.labels};
                EEG = pop_tesa_interactivechanreject(EEG);
                EEG = eeg_checkset( EEG );
                EEG = recordBadChannels(EEG, badBefore);

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
                fileReport = recordRejectedTrials(fileReport, rejepochs);
                EEG = eeg_checkset( EEG );

            case 'Run ICA (FastICA)'
                vars = convertContainedStringsToChars(varin);
                EEG = runIcaEngine(EEG, 'fastica', vars);

            case 'Run ICA (Infomax)'
                vars = convertContainedStringsToChars(varin);
                % runica takes a numeric 'extended' flag (1 = extended infomax).
                eidx = find(strcmpi(vars, 'extended'), 1);
                if ~isempty(eidx)
                    vars{eidx+1} = double(strcmpi(vars{eidx+1}, 'on'));
                end
                EEG = runIcaEngine(EEG, 'runica', vars);

            case 'Run ICA (Picard)'
                % Picard's params (mode, maxiter) are forwarded by pop_runica to
                % picard(); mode 'standard' = unconstrained (Infomax-like).
                vars = convertContainedStringsToChars(varin);
                EEG = runIcaEngine(EEG, 'picard', vars);

            case 'Label ICA Components'
                vars = convertContainedStringsToChars(varin);
                % Registry key is 'version' (the underlying pop_iclabel
                % argument). Older dispatch looked up 'iclabelVersion'
                % and crashed with "Unable to perform assignment with 0
                % elements" when the lookup returned empty.
                ind = find(strcmpi(vars, 'version'), 1);
                if isempty(ind)
                    iclabelVersion = 'default';
                else
                    iclabelVersion = vars{ind+1};
                end
                EEG = pop_iclabel(EEG, iclabelVersion);
                EEG = eeg_checkset( EEG );

            case 'Flag ICA Components for Rejection'
                vars = convertContainedStringsToChars(varin);
                classNames = {'Brain','Muscle','Eye','Heart','LineNoise', ...
                              'ChannelNoise','Other'};
                threshold = nan(7,2);
                for nflag = 2:2:14
                    t = vars{nflag};
                    cls = nflag/2;
                    if isempty(t) || (numel(t) == 2 && all(isnan(t)))
                        threshold(cls,:) = [NaN NaN];   % disabled for this class
                    elseif numel(t) == 2
                        threshold(cls,:) = t(:)';
                    else
                        % A single number scalar-expands to [x x], and
                        % pop_icflag then tests p > x & p < x, which is never
                        % true - the class is silently NOT flagged. Refuse it.
                        error('nestapp:icflagBadThreshold', ...
                            ['Flag ICA Components: the %s threshold must be a ' ...
                             '[min max] pair (e.g. [0.9 1]), or empty to ' ...
                             'disable. A single number flags nothing.'], ...
                            classNames{cls});
                    end
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
                    nCompDec = size(EEG.icaweights, 1);
                    act2D = computeICAActivation(EEG);
                    act3D = reshape(act2D, nCompDec, size(EEG.data,2), size(EEG.data,3));
                    % Cache activations so pop_subcomp's eeg_getica doesn't recompute them.
                    if isempty(EEG.icaact)
                        EEG.icaact = act3D;
                    end
                    % Per-component variance share using TESA's definition: the
                    % variance of the TRIAL-AVERAGED activation, normalised to
                    % sum to 100% across components (tesa_compselect.m:420-422).
                    % Matching the Remove ICA Components (TESA) path makes
                    % varRemoved comparable across every removal step. The prior
                    % var(act)/var(data) form ignored the spatial weighting of
                    % each component and badly under-counted (e.g. AARATEP round
                    % 2 reported 38 removed comps as only 0.55% variance).
                    compVarsAbs = var(mean(act3D, 3, 'omitnan'), 0, 2)';
                    sv = sum(compVarsAbs);
                    if sv > 0
                        pendingICAStats.compVarPct = double(compVarsAbs / sv * 100);
                    end
                end
                if isfield(EEG,'etc') && isfield(EEG.etc,'ic_classification') && ...
                        isfield(EEG.etc.ic_classification,'ICLabel') && ...
                        isfield(EEG.etc.ic_classification.ICLabel,'classifications')
                    pendingICAStats.iclabelProbs = ...
                        EEG.etc.ic_classification.ICLabel.classifications;
                end
                % Per-component category labels written by custom classifiers
                % (AARATEP muscle, ARTIST decay) so removed ICs are attributed
                % to a category even without ICLabel.
                if isfield(EEG,'etc') && isfield(EEG.etc,'nestappICClass')
                    pendingICAStats.classLabels = EEG.etc.nestappICClass;
                end

                % Render-ready snapshot of the FULL pre-rejection
                % decomposition so the QC figure can show every component
                % and which were rejected - pop_subcomp below deletes the
                % flagged ones. Topos + size + labels only (no .data), so
                % it stays small and serialises across parfor workers.
                latestICASnapshot = pendingICAStats;
                latestICASnapshot.icawinv     = EEG.icawinv;
                latestICASnapshot.chanlocs    = EEG.chanlocs;
                latestICASnapshot.icachansind = EEG.icachansind;
                latestICASnapshot.capturedStep = si;
                latestICASnapshot.capturedName = stepName;
                if isfield(EEG,'etc') && isfield(EEG.etc,'ic_classification') && ...
                        isfield(EEG.etc.ic_classification,'ICLabel') && ...
                        isfield(EEG.etc.ic_classification.ICLabel,'classes')
                    latestICASnapshot.iclabelClasses = ...
                        EEG.etc.ic_classification.ICLabel.classes;
                end

                if ~(isnumeric(EEG.reject.gcompreject) || islogical(EEG.reject.gcompreject))
                    EEG.reject.gcompreject = zeros(1, size(EEG.icaweights, 1));
                end

                % Structural guard: when the flag step marks EVERY component for
                % removal, pop_subcomp is asked to project the data onto an empty
                % surviving set and dies with an opaque matrix-dimension error.
                % Fail here instead, with an interpretable, trackable reason. The
                % trigger is purely structural (all components flagged -> no signal
                % left); the reference/ground note is a diagnostic hint, not a
                % threshold - removing 100% of the decomposition is unrecoverable
                % whatever the cause.
                nCompDecomp = size(EEG.icaweights, 1);
                nFlagged    = sum(logical(EEG.reject.gcompreject));
                if isempty(var_comp) && ~keepcomp && nCompDecomp > 0 && nFlagged >= nCompDecomp
                    error('nestapp:allComponentsFlagged', ...
                        ['All %d independent components were flagged for removal - no ' ...
                         'component would survive, so nothing can be reconstructed. This ' ...
                         'typically means no brain signal was recoverable (e.g. a bad ' ...
                         'reference or ground corrupting the whole montage, which is ' ...
                         'common-mode and so leaves no single channel as an outlier).'], ...
                        nCompDecomp);
                end

                EEG = pop_subcomp( EEG, var_comp, plotag, keepcomp);
                EEG.ICA_Rejected_Comp = ICA_Rejected_Comp;
                EEG = eeg_checkset( EEG );

            case 'Interpolate Channels'
                method = step.params.method;
                % Interpolate ONLY channels that were removed *as bad* by a
                % bad-channel step (recorded in EEG.etc.nestapp.badChannels by
                % recordBadChannels). Channels dropped on purpose - e.g. by
                % "Remove un-needed Channels" or a keep-list - also land in
                % EEG.chaninfo.removedchans, because pop_select records
                % everything it removes; restoring all of removedchans would
                % resurrect those intentionally-dropped channels. We take the
                % channel-location structs (with coordinates) from
                % removedchans but keep only entries whose label is on the
                % bad-channel list.
                %
                % We additionally skip stale / placeholder removedchans
                % entries - those with empty coordinates or a label that is
                % already live. pop_interp would otherwise append them to
                % chanlocs with no matching data row, after which eeg_checkset
                % discards the entire chanlocs struct on the size mismatch and
                % the next step crashes with a cryptic dot-indexing error.
                % See test_interpolateChannelsBadRemovedchans.
                rc = EEG.chaninfo.removedchans;
                if ~isempty(rc)
                    liveLabels = {EEG.chanlocs.labels};
                    if isfield(EEG, 'etc') && isfield(EEG.etc, 'nestapp') && ...
                            isfield(EEG.etc.nestapp, 'badChannels')
                        badLabels = EEG.etc.nestapp.badChannels;
                    else
                        badLabels = {};
                    end
                    isBad   = arrayfun(@(c) ismember(c.labels, badLabels), rc);
                    isValid = arrayfun(@(c) ~isempty(c.X) && ...
                        ~ismember(c.labels, liveLabels), rc);
                    rc = rc(isBad & isValid);
                end
                if isempty(rc)
                    fprintf(['Interpolate Channels: no bad removed ' ...
                        'channels to restore; skipping interpolation.\n']);
                else
                    % No time-range arg: pop_interp discards its 4th argument
                    % on entry (pop_interp.m:61 sets t_range=''), so it never
                    % worked. Call the plain 3-arg form.
                    EEG = pop_interp(EEG, rc, method);
                    interpElecs = [interpElecs; num2cell(rc)]; %#ok<AGROW>
                    EEG.interpElecs = interpElecs;
                end
                EEG.setname = [EEG.setname '_interp'];
                EEG.filename = [EEG.setname '.set'];
                EEG.datfile  = [EEG.setname '.fdt'];
                EEG = eeg_checkset( EEG );
                % Safety net: surface a legible error at the offending step
                % if chanlocs were dropped, instead of letting a downstream
                % step fail on an empty (double) EEG.chanlocs.
                if ~isstruct(EEG.chanlocs) || isempty(EEG.chanlocs)
                    error('nestapp:chanlocsDropped', ...
                        ['Interpolate Channels discarded channel locations ' ...
                         '(chanlocs/data size mismatch). Check ' ...
                         'EEG.chaninfo.removedchans for stale or ' ...
                         'coordinate-less entries.']);
                end

            case 'Find TMS Pulses (TESA)'
                vars = convertContainedStringsToChars(varin);
                ind = find(strcmp(vars,'elec'));
                elec = vars{ind+1};
                vars([ind,ind+1]) = [];
                if iscell(elec)
                    elec = elec{:};
                end
                % tesa_findpulse stores options by case-sensitive dynamic
                % field, so nestapp's lowercase keys must be mapped to its
                % exact spelling or the value is silently dropped and the
                % upstream default label is used instead.
                vars = renameVarinKeys(vars, {'tmslabel','pairlabel'}, ...
                                             {'tmsLabel','pairLabel'});
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
                ind1 = find(strcmp(vars,'types'), 1);
                type = {};
                if ~isempty(ind1); type = vars{ind1+1}; end
                % Require explicit event marker label(s). Epoching around an
                % empty type list silently epochs every event - rarely what is
                % intended, and wrong for an ERP task (oddball / go-no-go) where
                % you epoch around specific stimulus/response codes.
                if isempty(type)
                    error('nestapp:epochingNoTypes', ...
                        ['Epoching needs the event marker label(s) to epoch around ' ...
                         '(e.g. {''S 1'',''S 2''} for an oddball task). Set the ' ...
                         'Epoching step''s "Event label(s)" field to your task triggers.']);
                end
                ind2 = find(strcmp(vars,'timelim'), 1);
                timelim = vars{ind2+1};
                vars([ind1, ind1+1, ind2, ind2+1]) = [];
                vars = stripEmptyVarin(vars);
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
                % pop_tesa_compselect classifies AND removes (its internal
                % pop_subcomp reduces EEG.icawinv to survivors). Grab the
                % full decomposition now so the QC snapshot can show every
                % component and which TESA flagged for removal.
                tesaPreWinv = EEG.icawinv;
                tesaPreChan = EEG.chanlocs;
                tesaPreCind = EEG.icachansind;
                EEG = pop_tesa_compselect( EEG,vars{:});
                EEG = eeg_checkset( EEG );
                tesaSnap = tesaICASnapshot(EEG, tesaPreWinv, ...
                    tesaPreChan, tesaPreCind, si, stepName);
                if isstruct(tesaSnap) && ~isempty(fieldnames(tesaSnap))
                    latestICASnapshot = tesaSnap;   % keep prior snap if this failed
                end

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
                % The PC / "Variance kept" param is stored as a string, so
                % after a UITable edit or a saved pipeline it arrives as the
                % char "{'data', 90}" instead of the cell tesa_sspsir needs
                % (it runs 1:PC internally -> "colon operator with char
                % operands" error). Normalise it back to a cell / number.
                pcIdx = find(strcmpi(vars, 'PC'), 1);
                if ~isempty(pcIdx)
                    vars{pcIdx+1} = parseSspsirPC(vars{pcIdx+1});
                end
                EEG = pop_tesa_sspsir(EEG, vars{:});
                EEG = eeg_checkset( EEG );

            case 'Source-Informed Sensor Cleaning (SOUND)'
                o = varinToStruct(varin);
                EEG = pop_tesa_sound(EEG, 'lambdaValue', o.lambdaValue, 'iter', o.iter);
                EEG = eeg_checkset( EEG );

            case 'Interpolate Missing Data (AR-Blend)'
                ensureAaratepOnPath();
                opts2 = varinToStruct(varin);
                artifactTimespan = [opts2.artifactStartMs, opts2.artifactEndMs] * 1e-3;
                fitDur = opts2.prePostFitMs * 1e-3;
                EEG = c_EEG_ReplaceEpochTimeSegment(EEG, ...
                    'timespanToReplace',   artifactTimespan, ...
                    'method',              'ARExtrapolation', ...
                    'prePostFitDurations', [fitDur, fitDur]);
                EEG = eeg_checkset( EEG );

            case 'Remove Decay Artifact'
                ensureAaratepOnPath();
                % c_TMSEEG_fitAndRemoveDecayArtifact calls fit() on doubles.
                % Make sure it dispatches to the Curve Fitting Toolbox (repairs a
                % mis-set path if it can); otherwise fail with a clear message
                % instead of the cryptic "Undefined function 'fit' for input
                % arguments of type 'double'".
                [fitAvailable, fitReason] = ensureCurveFittingFit();
                if ~fitAvailable
                    error('nestapp:CurveFittingUnavailable', ...
                        ['Remove Decay Artifact requires the Curve Fitting Toolbox ' ...
                         'fit() function. %s'], fitReason);
                end
                opts2 = varinToStruct(varin);
                artifactTimespan = [opts2.artifactStartMs, opts2.artifactEndMs] * 1e-3;
                % Upstream c_TMSEEG_Preprocess_AARATEPPipeline.m line 336:
                %   doDecayRemovalPerTrial = true  -> 'none'
                %   doDecayRemovalPerTrial = false -> 'mean'
                if strcmpi(opts2.perTrial, 'on')
                    trialAggMethod = 'none';
                else
                    trialAggMethod = 'mean';
                end
                EEG = c_TMSEEG_fitAndRemoveDecayArtifact(EEG, ...
                    'artifactTimespan',                 artifactTimespan, ...
                    'trialAggMethod_timeCourseRemoval', trialAggMethod, ...
                    'aggTrimPercent',                   10);
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

                % A manual step must hand back what the user actually chose.
                %
                % This used to call pop_rejmenu(EEG, 1), which is declared with
                % NO output and reports the user's marks by mutating EEG in the
                % BASE workspace through its callback strings. A pipeline runs
                % with `global EEG`, which is a different variable, so nothing
                % the user marked ever came back: only pop_jointprob's automatic
                % rejections were applied and every manual selection was
                % silently discarded.
                %
                % eegplot's 'command' callback is the documented way to harvest
                % a selection: it is evaluated with TMPREJ holding the marked
                % regions. Routing that through a global gets it into this
                % scope, and eegplot2trial converts regions to trial indices.
                autoRej = logical(EEG.reject.rejjp(:)');
                NESTAPP_TMPREJ = [];
                winrej = trial2eegplot(autoRej, zeros(EEG.nbchan, EEG.trials), ...
                                       EEG.pnts, [0.9 0.5 0.5]);
                eegplot(EEG.data, 'srate', EEG.srate, 'winlength', 5, ...
                    'limits', [EEG.xmin EEG.xmax] * 1000, 'winrej', winrej, ...
                    'butlabel', 'Confirm rejections', ...
                    'command', 'global NESTAPP_TMPREJ; NESTAPP_TMPREJ = TMPREJ;');
                uiwait(gcf);

                % The user's marks are authoritative: confirming replaces the
                % automatic set rather than adding to it, so a trial they
                % un-marked is genuinely kept.
                EEG.BadTr = harvestMarkedTrials(NESTAPP_TMPREJ, ...
                                                EEG.pnts, EEG.trials);
                NESTAPP_TMPREJ = [];
                EEG = pop_rejepoch( EEG, EEG.BadTr ,0);
                fileReport = recordRejectedTrials(fileReport, EEG.BadTr);
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
                % vars{:} - pop_tesa_tepextract takes varargin and forwards it
                % as name/value pairs. Passing the cell itself makes
                % varargin{1} a cell, so the pair count is odd and TESA throws
                % "EXAMPLE needs key/value pairs" (tesa_tepextract.m:103).
                EEG = pop_tesa_tepextract( EEG, type, vars{:} );

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
                % vars{:} for the same reason as Extract TEP above: this takes
                % varargin, so vars(:) would arrive as a single cell argument
                % rather than as name/value pairs.
                EEG = pop_tesa_peakanalysis( EEG, input, direction, peak, peakWin, vars{:} );

            case 'TEP Peak Output'
                vars = convertContainedStringsToChars(varin);
                vars = stripEmptyVarin(vars);
                % Capture the peak table upstream returns - the step's whole
                % purpose - and store it on the dataset. It was previously
                % discarded, so the step produced nothing observable. With
                % averageWin or calcType='area' set, output.amp is a windowed
                % measure that is NOT on EEG.ROI, so this is the only place it
                % is preserved (and it round-trips through pop_saveset).
                tepPeakOutput = pop_tesa_peakoutput( EEG, vars{:} );
                if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc)
                    EEG.etc = struct();
                end
                if ~isfield(EEG.etc, 'nestapp') || ~isstruct(EEG.etc.nestapp)
                    EEG.etc.nestapp = struct();
                end
                EEG.etc.nestapp.tepPeakOutput = tepPeakOutput;

            case 'Quality Gate'
                % Pass the running rejection tally as context so the
                % gate's maxRejected*Pct metrics can compare against
                % the original counts logged at Load Data / Epoching.
                qgContext = struct( ...
                    'channels', fileReport.channels, ...
                    'trials',   fileReport.trials);
                gate = qualityGate(EEG, step.params, qgContext);
                gate.stepIndex = si;
                gate.stepName  = stepName;
                % Auto-disambiguate unset / default gate labels with the
                % step index so the dashboard, batch finalizer, and PDF
                % don't collapse multiple "gate" entries into one row.
                if isempty(gate.label) || strcmp(gate.label, 'gate')
                    gate.label = sprintf('gate-%02d', si);
                end
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

            otherwise
                % A name with no case fell straight through this switch and
                % was logged as a completed step having done nothing at all -
                % the quietest possible failure. It is reachable: specFromSaved
                % keeps steps it cannot resolve (with a warning), and a
                % mis-edited case label silently orphans its own step.
                error('nestapp:unknownStep', ...
                    ['Step %d ("%s") has no implementation. It is in the ' ...
                     'pipeline but nothing would run for it. If this pipeline ' ...
                     'was saved by a newer nestapp, the step may not exist ' ...
                     'in this version.'], si, stepName);

        end % switch

        %% Post-step metrics and report update
        if isstruct(EEG) && ~isempty(EEG)
            nChanAfter  = EEG.nbchan;
            nEpochAfter = size(EEG.data, 3);
            if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
                labelsAfter = {EEG.chanlocs.labels};
            else
                labelsAfter = {};
            end
            if strcmp(stepName, 'Load Data')
                fileReport.channels.original = EEG.nbchan;
            end
            fileReport.channels.final = EEG.nbchan;
            if nChanAfter < nChanBefore
                if any(strcmp(stepName, channelRejectionSteps()))
                    fileReport.channels.nRejected = fileReport.channels.nRejected + ...
                        (nChanBefore - nChanAfter);
                    % Record which labels disappeared so the report can name
                    % the rejected channels, not just count them. 'stable'
                    % keeps montage order; force a row for clean concatenation.
                    removedNames = setdiff(labelsBefore, labelsAfter, 'stable');
                    fileReport.channels.rejectedNames = ...
                        [fileReport.channels.rejectedNames, removedNames(:)'];
                    % Split by reason so the report and summary can keep
                    % intentional removals ("Remove un-needed Channels") distinct
                    % from quality-based bad-channel detection.
                    if strcmp(stepName, 'Remove un-needed Channels')
                        fileReport.channels.unneededNames = ...
                            [fileReport.channels.unneededNames, removedNames(:)'];
                    else
                        fileReport.channels.badChannelNames = ...
                            [fileReport.channels.badChannelNames, removedNames(:)'];
                    end
                end
            end
            if any(strcmp(stepName, {'Interpolate Channels','Interpolate Missing Data (TESA)'})) ...
                    && nChanAfter > nChanBefore
                fileReport.channels.nInterpolated = fileReport.channels.nInterpolated + ...
                    (nChanAfter - nChanBefore);
                addedNames = setdiff(labelsAfter, labelsBefore, 'stable');
                fileReport.channels.interpolatedNames = ...
                    [fileReport.channels.interpolatedNames, addedNames(:)'];
            end
            if strcmp(stepName, 'Epoching') && fileReport.trials.original == 0
                fileReport.trials.original    = size(EEG.data, 3);
                fileReport.trials.survivingIdx = 1:fileReport.trials.original;
                fileReport.trials.rejectedIndices = [];
            end
            if size(EEG.data, 3) > 1
                fileReport.trials.final = nEpochAfter;
                if nEpochBefore > 1 && nEpochAfter < nEpochBefore
                    fileReport.trials.rejected = fileReport.trials.rejected + ...
                        (nEpochBefore - nEpochAfter);
                end
            end
            if any(strcmp(stepName, icaDecompositionSteps())) && ~isempty(EEG.icaweights)
                % Open a round per decomposition so the round count is correct
                % even for a round that removes nothing (e.g. AARATEP's 2nd ICA).
                % Every ICA engine must be listed in icaDecompositionSteps -
                % when the per-engine names (Picard/FastICA/Infomax) were
                % missing here, a second decomposition never opened its own
                % round and its removal was folded into the prior round, hiding
                % the whole pass in the report.
                fileReport = openICARound(fileReport, size(EEG.icaweights, 1));
            end
            if strcmp(stepName, 'Remove Flagged ICA Components') && ...
                    isfield(pendingICAStats, 'rejMask')
                rMask = logical(pendingICAStats.rejMask(:)');
                nRej  = sum(rMask);
                if nRej > 0
                    removal = struct( ...
                        'nComponents', numel(rMask), ...
                        'nRejected',   nRej, ...
                        'varRemoved',  NaN, 'varMin', NaN, 'varMax', NaN);
                    if isfield(pendingICAStats, 'compVarPct') && ...
                            numel(pendingICAStats.compVarPct) == numel(rMask)
                        rejPct          = pendingICAStats.compVarPct(rMask);
                        removal.varRemoved = sum(rejPct);
                        removal.varMin     = min(rejPct);
                        removal.varMax     = max(rejPct);
                    end
                    removal.categories = icaCategoriesFromFlags(rMask, pendingICAStats);
                    fileReport         = addICARemoval(fileReport, removal);
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
                removal = struct( ...
                    'nComponents', numel(cl.compClass), ...
                    'nRejected',   nRejTESA, ...
                    'varRemoved',  NaN, 'varMin', NaN, 'varMax', NaN);
                removal.categories.names    = TESA_CATS;
                removal.categories.nRemoved = zeros(1, numel(TESA_CATS));
                removal.categories.varShare = zeros(1, numel(TESA_CATS));
                hasVars = isfield(cl, 'compVars') && numel(cl.compVars) >= numel(cl.compClass);
                if hasVars && nRejTESA > 0
                    rejPct             = double(cl.compVars(rejIdx));
                    removal.varRemoved = sum(rejPct);
                    removal.varMin     = min(rejPct);
                    removal.varMax     = max(rejPct);
                end
                for ci = 1:numel(TESA_CODES)
                    inCat = (cl.compClass == TESA_CODES(ci));
                    removal.categories.nRemoved(ci) = sum(inCat);
                    if hasVars
                        removal.categories.varShare(ci) = sum(cl.compVars(inCat));
                    end
                end
                fileReport = addICARemoval(fileReport, removal);
            end
        else
            nChanAfter  = nChanBefore;
            nEpochAfter = nEpochBefore;
        end % if isstruct(EEG)

        stepRec.name         = stepName;
        stepRec.params       = step.params;   % drives the methods narrative (methodsClause)
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
                % Serialize the PNG export across parallel workers via a lock
                % under the (shared) batch root - concurrent exportgraphics on
                % headless workers can deadlock and wedge the whole batch.
                qcOpts.exportLockDir = opts.batchCtx.batchRoot;
                % Most recent pre-rejection ICA state - lets the ICA panel
                % show all components (rejected ones marked) instead of the
                % post-pop_subcomp survivors sitting in the live EEG.
                qcOpts.icaSnapshot = latestICASnapshot;
                % Hand rejected-trial info to the renderer so the
                % heatmap can show gaps + red bars at the original
                % positions of removed epochs.
                qcOpts.rejectedTrialIdx = fileReport.trials.rejectedIndices;
                qcOpts.originalTrials   = fileReport.trials.original;
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
            % Leave a partial, FAIL-labelled report on disk so a skipped
            % file is still documented (not silently dropped).
            persistFailedReport(fileReport, opts, 'qualityFail', si, stepName, err.message);
            if ~isempty(opts.progressQueue)
                send(opts.progressQueue, struct( ...
                    'fi', opts.fileIndex, 'si', 0, ...
                    'nSteps', nSteps, 'stepName', 'Skipped', 'failed', true));
            end
            rethrow(err);
        end

        warning('An error occurred at file %s at step %d: %s', fileName, si, stepName);

        % Save a metadata-only debug bundle (never raw data) for bug reports.
        if opts.saveErrorBundle
            try
                if ~isempty(opts.batchCtx) && isfield(opts.batchCtx, 'batchRoot')
                    bundleParent = fullfile(opts.batchCtx.batchRoot, 'debug');
                else
                    bundleParent = '';
                end
                saveErrorBundle(bundleParent, struct( ...
                    'err', err, 'EEG', EEG, 'spec', spec, ...
                    'stepName', stepName, 'stepIndex', si, ...
                    'fileName', fileName, 'pipelineName', opts.pipelineName));
            catch bundleErr
                sendWorkerLog(opts.logQueue, wLabel, ...
                    'Error bundle save FAILED (continuing): %s', bundleErr.message);
            end
        end

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
                % Leave a partial, FAIL-labelled report on disk so the
                % failed file is still documented.
                persistFailedReport(fileReport, opts, 'error', si, stepName, err.message);
                % Release the dialog slot first - the success sentinel at the
                % bottom of this function is unreachable from here, and without
                % a release later files' progress messages get dropped.
                if ~isempty(opts.progressQueue)
                    send(opts.progressQueue, struct( ...
                        'fi', opts.fileIndex, 'si', 0, ...
                        'nSteps', nSteps, 'stepName', 'Failed', 'failed', true));
                end
                % The cause is ATTACHED, not just interpolated. Rethrowing
                % with only err.message discarded the inner identifier, so a
                % caller could no longer tell "this step has no
                % implementation" (nestapp:unknownStep) from "this step ran
                % and failed on this data" - a distinction the dispatch
                % coverage check depends on, and which it was reduced to
                % recovering by matching English in the message text.
                outer = MException('nestapp:stepFailed', ...
                    'Step %d (%s) failed: %s', si, stepName, err.message);
                throw(addCause(outer, err));
            end
        end
    end
end % step loop

%% Redundant copies, now that the file has come through
% Queued during the run and only acted on here: every step succeeded, so
% whatever was going to replace these files has been written.
for pd = 1:numel(pendingDeletes)
    try
        delete(pendingDeletes{pd});
        sendWorkerLog(opts.logQueue, wLabel, 'Removed redundant %s', pendingDeletes{pd});
    catch
        % A locked file costs disk, not results.
    end
end

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

function persistFailedReport(report, opts, kind, si, stepName, message)
% PERSISTFAILEDREPORT  Write a partial, FAIL-labelled report to disk so a
%   file abandoned mid-pipeline still leaves a report instead of vanishing.
%   kind: 'qualityFail' (failed a Quality Gate with skip-on-fail) or 'error'
%   (a step threw). Best-effort - reporting problems must never mask the
%   original failure, so everything is wrapped in try/catch.
if isempty(opts.batchCtx)
    return   % no batch context (e.g. ad-hoc run) - nowhere structured to write
end
if ~isfield(report, 'quality') || ~isstruct(report.quality)
    report.quality = struct('figures', {{}}, 'gates', {{}}, 'worstVerdict', 'Fail');
end
report.quality.worstVerdict = 'Fail';        % overall outcome = failed
report.failure = struct('failed', true, 'kind', kind, ...
    'stepIndex', si, 'stepName', stepName, 'message', message);
try
    exportReport(report, opts.batchCtx);
    if opts.autoExportPDF
        exportFileReportPDF(report, opts.batchCtx);
    end
catch writeErr
    warning('nestapp:failedReport', ...
        'Could not write partial failed report for "%s": %s', ...
        report.inputFile, writeErr.message);
end
end

function fileReport = recordRejectedTrials(fileReport, localIdx)
% Map locally-indexed rejected trials back to original-trial numbers
% and append to the cumulative list. Maintains a surviving-trial map
% so chained rejection steps stay correct even when multiple bad-
% epoch passes run in the same pipeline. No-op when rejection ran
% before any Epoching step (no surviving map yet) or when localIdx
% is empty.
if isempty(localIdx), return, end
if ~isfield(fileReport.trials, 'survivingIdx') ...
        || isempty(fileReport.trials.survivingIdx)
    return
end
localIdx = sort(localIdx(:)');
keep = localIdx >= 1 & localIdx <= numel(fileReport.trials.survivingIdx);
localIdx = localIdx(keep);
if isempty(localIdx), return, end
rejectedOriginal = fileReport.trials.survivingIdx(localIdx);
fileReport.trials.rejectedIndices = sort([fileReport.trials.rejectedIndices, ...
                                          rejectedOriginal]);
fileReport.trials.survivingIdx(localIdx) = [];
end

function cats = tesaICACategories()
% TESA component category names, in the order TESA's icaCompClass.compClass codes map to.
cats = {'TMS Muscle','Blink','Eye Move','Muscle','Elec Noise','Sensory','Reject'};
end

function codes = tesaICAClassCodes()
% TESA compClass integer codes, matched positionally to tesaICACategories().
codes = [3, 4, 5, 6, 7, 8, 2];
end

function snap = tesaICASnapshot(EEG, preWinv, preChan, preCind, stepIdx, stepName)
% TESAICASNAPSHOT  Pre-removal ICA snapshot from TESA compselect output.
%   pop_tesa_compselect stores per-component classes in EEG.icaCompClass
%   (full set) and removes the flagged ones via an internal pop_subcomp.
%   preWinv is the decomposition captured BEFORE that call, so it still
%   holds every component and aligns with the full compClass vector. The
%   returned struct matches what renderQualityFigure's buildICAView wants.
snap = struct([]);
if ~isfield(EEG, 'icaCompClass') || ~isstruct(EEG.icaCompClass) || ...
        isempty(fieldnames(EEG.icaCompClass))
    return
end
keys = fieldnames(EEG.icaCompClass);
cl   = EEG.icaCompClass.(keys{end});
if ~isfield(cl, 'compClass') || isempty(cl.compClass)
    return
end
codes = cl.compClass(:)';
n     = numel(codes);
if isempty(preWinv) || size(preWinv, 2) ~= n
    return   % can't align topos with classes - skip rather than mislabel
end

snap = struct();
snap.icawinv      = preWinv;
snap.chanlocs     = preChan;
snap.icachansind  = preCind;
snap.rejMask      = codes > 1;                 % code 1 = Keep
snap.classLabels  = tesaCompClassLabels(codes);
snap.source       = 'TESA';
if isfield(cl, 'compVars') && numel(cl.compVars) >= n
    snap.compVarPct = double(reshape(cl.compVars(1:n), 1, []));
end
snap.capturedStep = stepIdx;
snap.capturedName = stepName;
end

function labels = tesaCompClassLabels(codes)
% TESACOMPCLASSLABELS  Map TESA compClass codes (1..8) to display labels.
%   Codes mirror pop_tesa_compselect; the strings are shown on each QC
%   topo tile next to the component's variance.
map    = {'Keep','Reject','TMS Muscle','Blink','Eye Move','Muscle', ...
          'Elec Noise','Sensory'};
labels = repmat({'Reject'}, 1, numel(codes));
for k = 1:numel(codes)
    c = codes(k);
    if c >= 1 && c <= numel(map)
        labels{k} = map{c};
    end
end
end


function v = worseVerdict(a, b)
% Return the worst severity of two verdicts: Fail > Marginal > Pass > NotChecked.
order = {'NotChecked', 'Pass', 'Marginal', 'Fail'};
ia = find(strcmp(a, order));
ib = find(strcmp(b, order));
if isempty(ia), ia = 0; end
if isempty(ib), ib = 0; end
v = order{max(ia, ib)};
end

function tf = savesNewSetLater(spec, si)
% True when a Save New Set step follows position si. AARATEP's final .mat and
% that step's .set hold the same dataset, so exactly one of them should
% survive - and only this tells us the .set is actually coming.
tf = any(strcmp({spec(si+1:end).name}, 'Save New Set'));
end

function p = savedSetPath(EEG, vars)
% Where Save New Set actually wrote, '' when it did not write at all.
%
% Prefer what pop_saveset stamped on the dataset: it reflects the .set
% extension it appended and any normalisation it applied, so it is the file
% that exists rather than the name we asked for. Fall back to the composed
% savenew value, which is still the destination when a save happened without
% updating the struct.
p = '';
if isstruct(EEG) && isfield(EEG, 'filename') && ~isempty(EEG.filename) ...
        && isfield(EEG, 'filepath') && ~isempty(EEG.filepath) ...
        && endsWith(EEG.filename, '.set', 'IgnoreCase', true)
    p = fullfile(EEG.filepath, EEG.filename);
    return
end

k = find(strcmp(vars, 'savenew'), 1);
if ~isempty(k) && numel(vars) > k && ischar(vars{k+1}) && ~isempty(vars{k+1})
    p = vars{k+1};
    if ~endsWith(p, '.set', 'IgnoreCase', true); p = [p '.set']; end
end
end

function labels = channelLabelsOf(EEG)
% Channel labels of the current dataset, for resolving upstream's channel
% indices to names. Empty when the montage has no labels.
labels = {};
if isstruct(EEG) && isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) ...
        && isfield(EEG.chanlocs, 'labels')
    labels = {EEG.chanlocs.labels};
end
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

function v = parseRmbaseVec(vars, key)
% Extract a numeric-vector param from a key-value varin cell, coercing the
% '[]' placeholder and string forms to a numeric value (or [] when absent).
v = [];
ind = find(strcmpi(vars, key), 1);
if isempty(ind) || ind + 1 > numel(vars)
    return
end
v = vars{ind + 1};
if ischar(v) || isstring(v)
    if strcmp(char(v), '[]')
        v = [];
    else
        v = str2num(char(v)); %#ok<ST2NM>
    end
end
end
