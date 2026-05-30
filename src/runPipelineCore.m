
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [allReports, allSummaries] = runPipelineCore(spec, filePaths, opts)
% RUNPIPELINECORE Execute a typed pipeline spec against a list of files.
%   [allReports, allSummaries] = RUNPIPELINECORE(spec, filePaths, opts)
%
%   spec      - struct array of PipelineStep (name + params struct)
%   filePaths - cell array of full absolute file paths
%   opts      - struct with optional fields:
%     .uiFigure     UIFigure handle (required for interactive steps/dialogs)
%     .pipelineName pipeline name string for EEG.history provenance
%     .statusBar    text component for status line messages
%     .parallel     logical - request parallel execution (default false)
%     .chanLocFile  pre-selected channel location file path (default '')
%
%   See also: processOneFile, initPipelineReport, exportReport

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'pipelineName'), opts.pipelineName = ''; end
if ~isfield(opts, 'statusBar'),    opts.statusBar    = []; end
if ~isfield(opts, 'parallel'),     opts.parallel     = false; end
if ~isfield(opts, 'chanLocFile'),  opts.chanLocFile  = ''; end
if ~isfield(opts, 'outputRoot'),   opts.outputRoot   = ''; end

persistent cachedNestappSrc cachedEeglabGenpath

nFiles = numel(filePaths);
nSteps = numel(spec);

allSummaries = {};
allReports   = {};

if nFiles == 0 || nSteps == 0; return; end

% Dependency check before touching any file.
stepNames = {spec.name};
[depsOk, depsMsg] = checkStepDependencies(stepNames, filePaths);
if ~depsOk
    uialert(opts.uiFigure, depsMsg, 'Missing Dependencies', 'Icon', 'error');
    return
end

% One-time deprecation log: warn the user about Quality Gate params
% renamed in this version. The gate itself silently aliases them.
warnDeprecatedGateParams(spec);

% Unified batch output: every run gets one timestamped folder.
% Destination priority: opts.outputRoot (programmatic override, used
% by tests) > nestapp.outputRoot pref > common parent of inputs.
% Layout pref picks between 'typeBased' and 'perInput'. All writers
% resolve their destinations through outputPaths(batchCtx, kind, stem).
layout   = getpref('nestapp', 'outputLayout', 'typeBased');
batchCtx = buildBatchContext(filePaths, opts.pipelineName, layout, opts.outputRoot);

% Debug log: when the 'debugLog' pref is on, tee the full run trace to a
% file in the batch folder (closed on cleanup, even if the run errors).
if getpref('nestapp', 'debugLog', false)
    lp = nestDebugLog('start', batchCtx.batchRoot);
    debugLogCleanup = onCleanup(@() nestDebugLog('stop'));
    nestLog('CFG', 'Debug log: %s', lp);
end

% Save-on-error bundle (metadata only) - read here, threaded to workers.
opts.saveErrorBundle = getpref('nestapp', 'saveErrorBundle', true);

nestLog('CFG', 'Output root: %s', batchCtx.outputRoot);
nestLog('CFG', 'Batch folder (%s layout): %s', batchCtx.layout, batchCtx.batchRoot);

% Pre-flight overwrite check.
if getpref('nestapp', 'suppressEEGLABDialogs', true)
    warnIfOverwriteFiles(spec, filePaths, opts, batchCtx);
end

% Auto Quality Report: render a QC PNG after every Quality Gate step
% when the pref is on. Prefs are read on the client once and threaded
% through to workers via opts.
autoQualityReport = getpref('nestapp', 'autoQualityReport', false);
qcAttribute       = 'minmax_no_tms';
qcTmsWindow       = [0 25];
qcTmsAutoDetect   = true;
if autoQualityReport
    qcAttribute = getpref('nestapp', 'qualityAttribute', 'minmax_no_tms');
    if ~any(strcmp(qcAttribute, qualityAttributeModes()))
        nestLog('QC', ['Invalid qualityAttribute pref "%s"; ' ...
            'falling back to minmax_no_tms'], qcAttribute);
        qcAttribute = 'minmax_no_tms';
    end

    qcTmsWindow = getpref('nestapp', 'qualityTmsWindow', [0 25]);
    if ~(isnumeric(qcTmsWindow) && numel(qcTmsWindow) == 2 ...
            && qcTmsWindow(2) > qcTmsWindow(1))
        nestLog('QC', 'Invalid qualityTmsWindow pref; falling back to [0 25] ms');
        qcTmsWindow = [0 25];
    end

    qcTmsAutoDetect = logical(getpref('nestapp', 'qualityTmsAutoDetect', true));

    nestLog('QC', ['Auto Quality Report enabled (attribute=%s, ' ...
        'tmsWindow=[%g %g] ms). A QC PNG is captured after every ' ...
        'Quality Gate step.'], ...
        qcAttribute, qcTmsWindow(1), qcTmsWindow(2));
end

% Quality Gate skip-on-fail preference. Off by default; when on, any
% Quality Gate step that emits 'Fail' short-circuits the rest of that
% file's pipeline (other files keep going). Incompatible with batch
% mode (batch verdicts are not known until after the run completes);
% we warn at the post-run step if both are seen together.
skipOnQualityFail = getpref('nestapp', 'skipOnQualityFail', false);
if skipOnQualityFail
    nestLog('QC', 'skipOnQualityFail = true (Quality Gate Fail short-circuits the file)');
end

% Per-file PDF report auto-export pref (Phase 4).
autoExportPDF = getpref('nestapp', 'autoExportPDF', false);
if autoExportPDF
    nestLog('QC', 'autoExportPDF = true (one PDF per file alongside the .mat report)');
end

% Citations for the methods this pipeline uses, derived from its steps.
% Logged once per batch so the references end up in the run log alongside the
% data, where the user will look when writing the methods section.
cites = stepCitations({spec.name});
if ~isempty(cites)
    nestLog('CITE', 'Methods used in this pipeline - please cite:');
    for ci = 1:numel(cites)
        nestLog('CITE', '  %s', cites(ci).reference);
        if ~isempty(cites(ci).doi)
            nestLog('CITE', '  DOI: %s', cites(ci).doi);
        end
    end
    nestLog('CITE', 'See THIRD_PARTY_NOTICES.md for vendored dependencies.');
end

% Parallel guard: requires PCT, no interactive steps, and >1 file.
useParallel = false;
if opts.parallel
    if nFiles <= 1
        parallelSkipMsg(opts.statusBar, ...
            sprintf('Parallel mode skipped: only %d file selected (need >1).', nFiles));
    elseif ~license('test', 'Distrib_Computing_Toolbox')
        parallelSkipMsg(opts.statusBar, ...
            'Parallel mode skipped: Parallel Computing Toolbox not licensed.');
    else
        interactiveSteps = findInteractiveSteps(spec, opts);
        if ~isempty(interactiveSteps)
            parallelSkipMsg(opts.statusBar, ...
                sprintf('Parallel mode skipped: interactive step(s): %s', ...
                strjoin(interactiveSteps, ', ')));
        else
            useParallel = true;
        end
    end
end

% Pool setup before the progress dialog so startup time doesn't inflate
% the first file's apparent duration.
nBars = 1;   % serial uses one slot; parallel uses one slot per worker
if useParallel
    maxWorkers = getpref('nestapp', 'maxParallelWorkers', 4);
    nWorkers   = min(nFiles, maxWorkers);
    nestLog('PAR', 'runPipelineCore: %d files, maxWorkers=%d', nFiles, maxWorkers);

    pool = gcp('nocreate');
    if isempty(pool)
        nestLog('PAR', 'No pool - starting parpool(%d)...', nWorkers);
        t0 = tic; parpool(nWorkers); pool = gcp('nocreate');
        nestLog('PAR', 'parpool ready (%d workers, %.2fs)', pool.NumWorkers, toc(t0));
    elseif pool.NumWorkers ~= nWorkers
        nestLog('PAR', 'Pool size mismatch (%d vs %d) - restarting...', pool.NumWorkers, nWorkers);
        t0 = tic; delete(pool); parpool(nWorkers); pool = gcp('nocreate');
        nestLog('PAR', 'New parpool ready (%d workers, %.2fs)', pool.NumWorkers, toc(t0));
    else
        nestLog('PAR', 'Reusing pool (%d workers)', pool.NumWorkers);
    end

    % Propagate paths to workers only (spmd skips the client - avoids
    % shadowing MATLAB built-ins with EEGLAB subdirectories on the client).
    % genpath(eeglab) is expensive (~30 s); cache it across calls.
    nestappSrc = fileparts(which('runPipelineCore'));
    if isempty(cachedEeglabGenpath) || ~strcmp(cachedNestappSrc, nestappSrc)
        nestLog('PAR', 'Building EEGLAB genpath cache...');
        t0 = tic;
        cachedNestappSrc    = nestappSrc;
        cachedEeglabGenpath = genpath(fileparts(which('eeglab')));
        nestLog('PAR', 'genpath done (%.2fs)', toc(t0));
    end
    eeglabGenpath = cachedEeglabGenpath;
    nestLog('PAR', 'Propagating paths to workers...');
    t0 = tic;
    spmd
        % Suppress MATLAB:dispatcher:nameConflict on each worker for
        % the duration of the addpath. EEGLAB plugins ship functions
        % named gather / labindex / numlabs (same as MATLAB's
        % parallel built-ins); the shadowing is expected and noisy
        % when it fires once per worker per run.
        % NB: spmd disallows anonymous functions, so we can't use
        % onCleanup here - manually re-enable on the way out.
        warning('off', 'MATLAB:dispatcher:nameConflict');
        if ~isempty(nestappSrc),    addpath(nestappSrc);    end
        if ~isempty(eeglabGenpath), addpath(eeglabGenpath); end
        warning('on', 'MATLAB:dispatcher:nameConflict');
    end
    nestLog('PAR', 'spmd done (%.2fs)', toc(t0));
    nCores           = feature('numcores');
    threadsPerWorker = max(1, floor(nCores / pool.NumWorkers));
    nestLog('PAR', 'Workers: %d | CPU cores: %d | BLAS threads/worker: %d', ...
        pool.NumWorkers, nCores, threadsPerWorker);

    nBars = min(nFiles, pool.NumWorkers);
end

% Unified N-bar progress dialog.  Serial uses nBars=1 (one slot cycling
% through each file); parallel uses one slot per worker.  Both modes use the
% same createProgressDlg / updateProgressDlg pair and the same message format.
dlg = createProgressDlg(opts.uiFigure, nBars, nFiles);
dlgCleanup = onCleanup(@() closeDlg(dlg));

reports   = cell(nFiles, 1);
cancelled = false;
failed    = struct('fi', {}, 'name', {}, 'step', {}, 'stepName', {}, ...
                   'message', {}, 'kind', {});

if useParallel
    % DataQueue carries per-step progress, log messages, and file-done
    % sentinels from workers - all routed through updateProgressDlg.
    q = parallel.pool.DataQueue;
    afterEach(q, @(msg) updateProgressDlg(dlg, msg, nFiles, false, []));

    % Strip all UI handles - workers cannot access graphics objects.
    wOpts = opts;
    wOpts.uiFigure       = [];
    wOpts.statusBar      = [];
    wOpts.progressFcn    = [];
    wOpts.onStepError    = [];
    wOpts.onPickChanFile = [];
    wOpts.progressQueue  = q;              % per-step progress + file-done sentinel
    wOpts.logQueue       = q;              % log msgs share the same queue
    wOpts.nWorkers       = pool.NumWorkers; % actual count for BLAS thread cap
    wOpts = applyQCOpts(wOpts, batchCtx, autoQualityReport, qcAttribute, qcTmsWindow, skipOnQualityFail, qcTmsAutoDetect, autoExportPDF);

    nestLog('PAR', 'Submitting %d futures...', nFiles);
    for fi = 1:nFiles
        fOpts           = wOpts;
        fOpts.fileIndex = fi;
        futures(fi) = parfeval(@processOneFile, 2, spec, filePaths{fi}, fOpts); %#ok<AGROW>
    end

    % Track which failed futures we've already captured + drained.
    drained = false(1, nFiles);

    % Poll until all futures finish or user cancels.
    while true
        pause(0.25); drawnow;
        if ~isvalid(dlg.fig) || dlg.fig.UserData.cancelRequested
            nestLog('PAR', 'Cancel requested - cancelling futures');
            cancel(futures);
            cancelled = true;
            t0 = tic;
            while toc(t0) < 30
                termStates = {futures.State};
                if all(strcmp(termStates,'finished') | ...
                       strcmp(termStates,'failed')   | ...
                       strcmp(termStates,'cancelled'))
                    break
                end
                pause(0.25); drawnow;
            end
            break
        end
        states = {futures.State};
        % Drain any newly-errored future immediately so we can record
        % its message via parseFailure and so the unfetched error
        % doesn't keep tripping fetchOutputs later. NB: a parfeval
        % future whose worker threw lands in State == 'finished' with
        % non-empty .Error - NOT State == 'failed'.
        for fi = 1:nFiles
            if ~drained(fi) && futureErrored(futures(fi))
                rec = drainFailedFuture(fi, futures(fi), filePaths{fi}, ~cancelled);
                if ~isempty(rec)
                    failed(end+1) = rec; %#ok<AGROW>
                end
                drained(fi) = true;
            end
        end
        if all(strcmp(states, 'finished') | strcmp(states, 'failed')); break; end
    end
    nestLog('PAR', 'Poll loop exited (cancelled=%d)', cancelled);

    % Errored futures and successful futures both end up in
    % State == 'finished'; the only safe way to fetch outputs is to
    % check for an error first and skip fetchOutputs on those (it
    % rethrows the worker's error on every call).
    for fi = 1:nFiles
        if futureErrored(futures(fi))
            if ~drained(fi)
                rec = drainFailedFuture(fi, futures(fi), filePaths{fi}, ~cancelled);
                if ~isempty(rec), failed(end+1) = rec; end %#ok<AGROW>
                drained(fi) = true;
            end
        elseif strcmp(futures(fi).State, 'finished')
            [reports{fi}, ~] = fetchOutputs(futures(fi));
        end
    end
    delete(futures);

else
    for fi = 1:nFiles
        if dlg.fig.UserData.cancelRequested; cancelled = true; break; end

        fOpts = opts;
        % Wrap progressFcn to produce the same message struct that workers
        % send via DataQueue - so both paths share updateProgressDlg.
        fOpts.progressFcn = @(si, sn) updateProgressDlg(dlg, ...
            struct('fi', fi, 'si', si, 'nSteps', nSteps, 'stepName', sn), ...
            nFiles, true, opts.statusBar);
        fOpts.onStepError = @(si, sn, err) uiconfirm(opts.uiFigure, ...
            sprintf('Error at step %d (%s):\n%s\n\nContinue to next step?', si, sn, err.message), ...
            'Step Failed', 'Options', {'Continue','Abort'}, ...
            'DefaultOption', 'Continue', 'CancelOption', 'Abort');
        fOpts.onPickChanFile = @() pickChanFile(opts.uiFigure);
        fOpts.progressQueue  = [];   % serial uses progressFcn, not DataQueue
        fOpts.fileIndex      = fi;
        fOpts = applyQCOpts(fOpts, batchCtx, autoQualityReport, qcAttribute, qcTmsWindow, skipOnQualityFail, qcTmsAutoDetect, autoExportPDF);

        try
            [reports{fi}, ~] = processOneFile(spec, filePaths{fi}, fOpts);
        catch err
            if strcmp(err.identifier, 'nestapp:cancelled')
                cancelled = true; break;
            end
            [~, fname, fext] = fileparts(filePaths{fi});
            rec = parseFailure(fi, [fname fext], err.message);
            failed(end+1) = rec; %#ok<AGROW>
            logFileFailure('SERIAL', rec);
            continue
        end

        % File-done sentinel: turn the slot green and advance the overall bar.
        % (Parallel mode sends this from within processOneFile via progressQueue.)
        updateProgressDlg(dlg, ...
            struct('fi', fi, 'si', 0, 'nSteps', nSteps, 'stepName', 'Done'), ...
            nFiles, false, []);

        if ~getpref('nestapp', 'hideEEGLABWindow', true)
            eeglab redraw
        end
    end
end

closeDlg(dlg);

nestLog('CFG', 'Batch artifacts saved to: %s', batchCtx.batchRoot);

% Batch-mode Quality Gate verdicts: resolved across all completed
% reports using median + N * MAD cutoffs. Pending verdicts inside
% successful reports become Pass / Marginal / Fail; reports with no
% batch-mode gates are unaffected.
if ~cancelled
    successReports = reports(~cellfun(@isempty, reports));
    if hasPendingBatchGates(successReports)
        if skipOnQualityFail
            nestLog('QC', ['Batch-mode Quality Gates ignore skipOnQualityFail ' ...
                '(verdicts are not known until after the run completes)']);
        end
        successReports = finalizeBatchVerdicts(successReports);
        % Write resolved reports back into the per-file slots so the
        % downstream summary / CSV / save sees the final verdicts.
        slot = 0;
        for fi = 1:nFiles
            if ~isempty(reports{fi})
                slot = slot + 1;
                reports{fi} = successReports{slot};
            end
        end
    end
end

% Post-run failure recovery: if some (but not all) files failed, give the
% user a chance to abandon the whole run before reports are generated.
if ~cancelled && ~isempty(failed)
    nSuccess = sum(~cellfun(@isempty, reports));
    if nSuccess > 0
        decision = promptFailureRecovery(opts.uiFigure, failed, nFiles);
        if strcmp(decision, 'Abandon Run')
            cancelled = true;
        end
    end
end

% Collect summaries for all successfully processed files.
summaries = cell(nFiles, 1);
for fi = 1:nFiles
    if ~isempty(reports{fi})
        [summaries{fi}, ~] = exportReport(reports{fi}, batchCtx);
    end
end
allReports   = reports(~cellfun(@isempty, reports));
allSummaries = summaries(~cellfun(@isempty, summaries));

% Batch-level artifacts: spec snapshot, dashboard PNG, summary CSV.
% Each wrapped in its own try so one failure can't take down the run.
if ~cancelled && ~isempty(allReports)
    writeBatchArtifacts(batchCtx, spec, opts.pipelineName, allReports, failed);
end

if cancelled
    % Discard any partially-completed reports - a cancelled run is not a result.
    error('nestapp:cancelled', 'Pipeline cancelled by user.');
end
if isempty(allReports)
    mode = 'serial';
    if useParallel; mode = 'parallel'; end
    error('nestapp:allFilesFailed', ...
        'All %d files failed in %s mode. Check the console log for details.', nFiles, mode);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unified progress dialog
%
% nBars=1 for serial (one slot cycling through files one at a time);
% nBars=N for parallel (one slot per worker, all active simultaneously).
% updateProgressDlg handles both cases identically.

function dlg = createProgressDlg(parentFig, nBars, nFiles)
PAD    = 12;
figW   = 440;
barW   = figW - 2*PAD;
btnH   = 28;
barH   = 10;
lblH   = 18;
rowH   = lblH + 4 + barH + 8;
titleH = 26;

yBtn        = PAD;
yOverallLbl = yBtn + btnH + PAD;
yOverallBar = yOverallLbl + lblH + 4;
ySlot1      = yOverallBar + barH + PAD;
figH        = ySlot1 + nBars * rowH + PAD + titleH;

if ~isempty(parentFig) && isvalid(parentFig)
    pPos = parentFig.Position;
    % Full-size overlay panel that blocks nestapp interaction while the pipeline runs.
    dlg.overlay = uipanel(parentFig, ...
        'Position',        [0, 0, pPos(3), pPos(4)], ...
        'BackgroundColor', [0.82 0.82 0.84], ...
        'BorderType',      'none');
    cX = (pPos(3) - figW) / 2;
    cY = (pPos(4) - figH) / 2;
    dlg.fig = uipanel(dlg.overlay, ...
        'Position',        [cX, cY, figW, figH], ...
        'BackgroundColor', [0.97 0.97 0.98], ...
        'BorderType',      'none');
    uilabel(dlg.fig, 'Text', 'Running Pipeline', ...
        'FontWeight', 'bold', 'FontSize', 13, ...
        'HorizontalAlignment', 'center', ...
        'Position', [0, figH - titleH, figW, titleH]);
else
    sc   = get(0, 'ScreenSize');
    figX = (sc(3) - figW) / 2;
    figY = (sc(4) - figH) / 2;
    dlg.overlay = [];
    dlg.fig = uifigure('Name', 'Running Pipeline', ...
        'Position', [figX figY figW figH], ...
        'Color',    [0.97 0.97 0.98], ...
        'Resize',   'off');
end

% Headless (no embedding app figure): also stream progress to the command
% window. When the pipeline runs under tests or a batch CLI there is no human
% watching the dialog, so a stall leaves no trace - the console echo (with
% nestLog timestamps) shows the last step reached before a hang.
dlg.streamConsole = isempty(dlg.overlay);

% slotMap(fi)=slot tracks which bar slot is assigned to each file.
% slotAvailable marks which slots are free to accept a new file.
dlg.fig.UserData = struct( ...
    'cancelRequested', false, ...
    'nDone',           0, ...
    'slotMap',         zeros(1, nFiles), ...
    'slotAvailable',   true(1, nBars));

uibutton(dlg.fig, 'push', 'Text', 'Cancel', ...
    'Position', [(figW-100)/2, yBtn, 100, btnH], ...
    'ButtonPushedFcn', @(~,~) setCancelFlag(dlg.fig));

dlg.overallLabel = uilabel(dlg.fig, ...
    'Text',       sprintf('0 / %d files complete', nFiles), ...
    'FontWeight', 'bold', ...
    'FontSize',   13, ...
    'Position',   [PAD, yOverallLbl, barW, lblH]);
uilabel(dlg.fig, 'Text', '', ...
    'BackgroundColor', [0.85 0.85 0.87], ...
    'Position',        [PAD, yOverallBar, barW, barH]);
dlg.overallFill = uilabel(dlg.fig, 'Text', '', ...
    'BackgroundColor', [0.16 0.67 0.47], ...
    'Position',        [PAD, yOverallBar, 0, barH]);

dlg.labels = gobjects(1, nBars);
dlg.fills  = gobjects(1, nBars);
for i = 1:nBars
    yLbl = ySlot1 + (i-1) * rowH;
    yBar = yLbl + lblH + 4;
    dlg.labels(i) = uilabel(dlg.fig, ...
        'Text',     'Idle', ...
        'FontSize', 11, ...
        'Position', [PAD, yLbl, barW, lblH]);
    uilabel(dlg.fig, 'Text', '', ...
        'BackgroundColor', [0.85 0.85 0.87], ...
        'Position',        [PAD, yBar, barW, barH]);
    dlg.fills(i) = uilabel(dlg.fig, 'Text', '', ...
        'BackgroundColor', [0.23 0.51 0.96], ...
        'Position',        [PAD, yBar, 0, barH]);
end

drawnow;
end

function updateProgressDlg(dlg, msg, nFiles, throwOnCancel, statusBar)
% Unified handler for both serial (throwOnCancel=true) and parallel (false).
% msg format - per-step: struct(fi, si, nSteps, stepName)
%              sentinel:  struct(fi, si=0, nSteps, stepName='Done')
%              log line:  struct(log=true, ts, label, text)
%
% Slots are assigned dynamically: a slot is claimed when the file's first
% step message arrives and released when the sentinel arrives.  This avoids
% the mod-based static assignment that breaks when workers finish at
% different speeds.
if nargin < 5; statusBar = []; end

if isfield(msg, 'log')
    fprintf('[%s][%s] %s\n', msg.ts, msg.label, msg.text);
    return
end

if ~isvalid(dlg.fig); return; end

% Streaming Quality Gate verdict (Phase 4). Repaints the slot label
% and recolors the bar fill until the next step start message
% restores the normal step-progress display.
if isfield(msg, 'gateVerdict') && ~isempty(msg.gateVerdict)
    if dlg.streamConsole
        nestLog('PROG', 'File %d \x2014 [%s] %s  (%d/%d)', ...
            msg.fi, upper(msg.gateVerdict), msg.stepName, msg.si, msg.nSteps);
    end
    udGate = dlg.fig.UserData;
    slot = udGate.slotMap(msg.fi);
    if slot > 0
        dlg.labels(slot).Text = sprintf( ...
            'File %d \x2014 [%s] %s  (%d/%d)', ...
            msg.fi, upper(msg.gateVerdict), msg.stepName, ...
            msg.si, msg.nSteps);
        dlg.fills(slot).BackgroundColor = verdictFill(msg.gateVerdict);
    end
    return
end

% In serial mode, flush queued UI events so a Cancel click is registered
% before we read the flag.  Not safe to call drawnow from afterEach handlers.
if throwOnCancel
    drawnow;
    if dlg.fig.UserData.cancelRequested
        error('nestapp:cancelled', 'Pipeline cancelled by user.');
    end
end

ud   = dlg.fig.UserData;
barW = dlg.overallLabel.Position(3);

if msg.si == 0
    % Sentinel: file is fully done on the worker (success or failure).
    slot = ud.slotMap(msg.fi);
    if slot == 0; return; end   % guard against duplicate sentinels
    ud.nDone              = ud.nDone + 1;
    ud.slotAvailable(slot) = true;   % release slot for the next file
    ud.slotMap(msg.fi)    = 0;
    dlg.fig.UserData = ud;
    nDone = ud.nDone;
    isFailed = isfield(msg, 'failed') && msg.failed;
    if dlg.streamConsole
        doneStatus = 'Done';
        if isFailed; doneStatus = 'FAILED'; end
        nestLog('PROG', 'File %d \x2014 %s  (%d / %d files complete)', ...
            msg.fi, doneStatus, nDone, nFiles);
    end
    if isFailed
        dlg.fills(slot).BackgroundColor = [0.85 0.27 0.27];   % red
        dlg.labels(slot).Text = sprintf('File %d \x2014 FAILED', msg.fi);
    else
        dlg.fills(slot).BackgroundColor = [0.16 0.67 0.47];   % green
        dlg.labels(slot).Text = sprintf('File %d \x2014 Done', msg.fi);
    end
    dlg.fills(slot).Position(3)     = barW;
    dlg.overallFill.Position(3) = round(barW * nDone / nFiles);
    dlg.overallLabel.Text       = sprintf('%d / %d files complete', nDone, nFiles);
else
    % Step starting: claim a slot on this file's first message.
    slot = ud.slotMap(msg.fi);
    if slot == 0
        avail = find(ud.slotAvailable, 1);
        if isempty(avail)
            nestLog('PROG', 'No free slot for file %d (step %d) - dropped', msg.fi, msg.si);
            return
        end
        slot = avail;
        ud.slotMap(msg.fi)     = slot;
        ud.slotAvailable(slot) = false;
        dlg.fig.UserData = ud;
    end
    if dlg.streamConsole
        nestLog('PROG', 'File %d \x2014 %s  (%d/%d)', ...
            msg.fi, msg.stepName, msg.si, msg.nSteps);
    end
    dlg.fills(slot).BackgroundColor = [0.23 0.51 0.96];
    dlg.fills(slot).Position(3)     = round(barW * msg.si / msg.nSteps);
    dlg.labels(slot).Text = sprintf('File %d \x2014 %s  (%d/%d)', ...
        msg.fi, msg.stepName, msg.si, msg.nSteps);
    if ~isempty(statusBar) && isvalid(statusBar)
        statusBar.Text = sprintf('  File %d / %d \x2014 %s', msg.fi, nFiles, msg.stepName);
    end
end
% Render the label updates before the step executes.  The drawnow at the top
% of this function flushes cancel clicks but runs before the label is set, so
% without this second call the user sees the previous step's label for the
% entire duration of the current step.  Not safe to call from afterEach handlers.
if throwOnCancel
    drawnow;
end
end

function setCancelFlag(fig)
ud = fig.UserData;
ud.cancelRequested = true;
fig.UserData = ud;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shared helpers

function c = verdictFill(verdict)
% Color the per-slot bar fill briefly while a Quality Gate verdict is
% on display in the progress dialog. Reset by the next step-start
% message back to the normal blue. Absolute-mode gates only ever
% return Pass/Marginal/Fail; batch-mode Pending is resolved
% post-batch and never streamed.
switch verdict
    case 'Pass',     c = [0.20 0.70 0.30];   % green
    case 'Marginal', c = [0.95 0.80 0.20];   % yellow
    case 'Fail',     c = [0.85 0.20 0.20];   % red
    otherwise,       c = [0.70 0.70 0.70];   % gray
end
end

function tf = hasPendingBatchGates(reports)
% Scan completed reports for any Quality Gate left in 'Pending' state.
tf = false;
for k = 1:numel(reports)
    r = reports{k};
    if ~isfield(r, 'quality') || ~isfield(r.quality, 'gates'), continue, end
    for gi = 1:numel(r.quality.gates)
        if strcmpi(r.quality.gates{gi}.verdict, 'Pending')
            tf = true; return
        end
    end
end
end

function opts = applyQCOpts(opts, batchCtx, autoQualityReport, attribute, tmsWindow, ...
        skipOnQualityFail, tmsAutoDetect, autoExportPDF)
% Assign the QC opts onto a worker/serial options struct.
opts.batchCtx          = batchCtx;
opts.autoQualityReport = autoQualityReport;
opts.qcAttribute       = attribute;
opts.qcTmsWindow       = tmsWindow;
opts.skipOnQualityFail = skipOnQualityFail;
opts.qcTmsAutoDetect   = tmsAutoDetect;
opts.autoExportPDF     = autoExportPDF;
end

function parallelSkipMsg(statusBar, msg)
fprintf('[nestapp] %s\n', msg);
if ~isempty(statusBar) && isvalid(statusBar)
    statusBar.Text = ['  ', msg];
end
end

function chFile = pickChanFile(uiFig)
[chName, chPath] = uigetfile('*.*', 'Select a channel location file', ...
    'Parent', uiFig);
if isequal(chName, 0)
    error('nestapp:cancelled', 'Channel location file selection cancelled.');
end
chFile = fullfile(chPath, chName);
end

function tf = futureErrored(f)
% A parfeval future whose worker threw has State == 'finished' with
% non-empty .Error in current MATLAB releases. Older releases used
% State == 'failed'. Check both shapes so we don't rely on the
% release-specific behaviour.
if strcmp(f.State, 'failed')
    tf = true;
    return
end
if strcmp(f.State, 'finished')
    try
        tf = ~isempty(f.Error);
    catch
        tf = false;
    end
    return
end
tf = false;
end

function rec = drainFailedFuture(fi, future, filePath, recordIt)
% Consume a failed parfeval future's error so MATLAB does not surface
% "One or more futures resulted in an error" later. Returns a
% parseFailure record when recordIt is true (genuine pre-cancel
% failure); returns [] when recordIt is false (cancel-induced).
rec = [];
errMsg = '';
try
    errMsg = future.Error.message;
catch
    % rare: future has no .Error - drain anyway.
end
try
    fetchOutputs(future);
catch
    % expected - error is now considered consumed
end
if recordIt
    [~, fname, fext] = fileparts(filePath);
    rec = parseFailure(fi, [fname fext], errMsg);
    logFileFailure('PAR', rec);
end
end

function rec = parseFailure(fi, name, msg)
% Parse a processOneFile error message into a structured failure record.
% Two known shapes (both emitted via error() from processOneFile):
%
%   nestapp:stepFailed
%     "Step 21 (Remove ICA Components (TESA)) failed: <root cause>"
%     -> rec.kind = 'errored'
%
%   nestapp:qualityFail
%     "Step 3 (Quality Gate "post-load") failed: <reason>; <reason>"
%     -> rec.kind = 'skipped' (the file was healthy enough to load
%        but failed user-defined quality thresholds; not really an error)
%
% Greedy capture on the step name because it can itself contain parens.
rec = struct('fi', fi, 'name', name, 'step', '', 'stepName', '', ...
    'message', msg, 'kind', 'errored');

% Try the Quality Gate shape first so its quoted-label form does not
% get mis-parsed by the generic regex below.
qg = regexp(msg, '^Step (\d+) \(Quality Gate "([^"]+)"\) failed:\s*(.*)$', ...
    'tokens', 'once', 'dotall');
if ~isempty(qg)
    rec.step     = qg{1};
    rec.stepName = sprintf('Quality Gate "%s"', qg{2});
    rec.message  = qg{3};
    rec.kind     = 'skipped';
    return
end

tok = regexp(msg, '^Step (\d+) \((.+)\) failed:\s*(.*)$', 'tokens', 'once', 'dotall');
if ~isempty(tok)
    rec.step     = tok{1};
    rec.stepName = tok{2};
    rec.message  = tok{3};
end
end

function logFileFailure(label, rec)
% Loud, basename-tagged failure line on stderr so it shows red in the
% MATLAB command window and is easy to spot among the per-step logs.
if isempty(rec.step)
    fprintf(2, '[%s] *** FILE %d %s FAILED: %s\n', ...
        label, rec.fi, rec.name, oneline(rec.message));
else
    fprintf(2, '[%s] *** FILE %d %s FAILED at step %s (%s)\n      %s\n', ...
        label, rec.fi, rec.name, rec.step, rec.stepName, oneline(rec.message));
end
end

function s = oneline(msg)
s = regexprep(msg, '\s*[\r\n]+\s*', ' | ');
end

function decision = promptFailureRecovery(uiFig, failed, nFiles)
% Post-run dialog: list failed files (grouped by kind) and let the user
% decide whether to accept the partial results or abandon the whole run.
nFailed = numel(failed);
nOk     = nFiles - nFailed;

kinds = arrayfun(@(f) failureKind(f), failed, 'UniformOutput', false);
skipped = failed(strcmp(kinds, 'skipped'));
errored = failed(~strcmp(kinds, 'skipped'));

sections = {};
if ~isempty(skipped)
    sections{end+1} = renderGroup('Skipped at Quality Gate:', skipped);
end
if ~isempty(errored)
    sections{end+1} = renderGroup('Errored:', errored);
end
listText = strjoin(sections, [newline newline]);

msg = sprintf(['%d of %d files failed:\n\n%s\n\n' ...
    'Continue with %d successful files, or abandon the whole run?'], ...
    nFailed, nFiles, listText, nOk);

if isempty(uiFig) || ~isvalid(uiFig)
    fprintf(2, '\n%s\n[No UI figure - continuing with successful files]\n\n', msg);
    decision = 'Continue';
    return
end
decision = uiconfirm(uiFig, msg, 'Some Files Failed', ...
    'Options', {'Continue', 'Abandon Run'}, ...
    'DefaultOption', 1, 'CancelOption', 2, ...
    'Icon', 'warning');
end

function k = failureKind(f)
if isfield(f, 'kind') && ~isempty(f.kind)
    k = f.kind;
else
    k = 'errored';
end
end

function txt = renderGroup(header, group)
MAX_SHOW = 10;
nShow = min(numel(group), MAX_SHOW);
lines = cell(1, nShow);
for k = 1:nShow
    f = group(k);
    if isempty(f.step)
        lines{k} = sprintf('  %s -- %s', f.name, oneline(f.message));
    elseif strcmp(failureKind(f), 'skipped')
        lines{k} = sprintf('  %s -- %s: %s', f.name, f.stepName, oneline(f.message));
    else
        lines{k} = sprintf('  %s -- step %s (%s)', f.name, f.step, f.stepName);
    end
end
body = strjoin(lines, newline);
if numel(group) > MAX_SHOW
    body = [body, sprintf('\n  ... and %d more', numel(group) - MAX_SHOW)];
end
txt = [header, newline, body];
end

function steps = findInteractiveSteps(spec, opts)
ALWAYS_INTERACTIVE = {'Visualize EEG Data', 'Remove Bad Trials', 'Choose Data Set'};
steps = {};
for si = 1:numel(spec)
    name = spec(si).name;
    if any(strcmp(name, ALWAYS_INTERACTIVE))
        steps{end+1} = name; %#ok<AGROW>
    elseif strcmp(name, 'Remove ICA Components (TESA)')
        p = spec(si).params;
        if ~isfield(p, 'compCheck') || strcmp(p.compCheck, 'on')
            steps{end+1} = name; %#ok<AGROW>
        end
    elseif strcmp(name, 'Load Channel Location')
        p            = spec(si).params;
        eachFileDiff = isfield(p, 'eachFilediffPath') && strcmp(p.eachFilediffPath, 'yes');
        needChan     = isfield(p, 'needchanloc') && strcmp(p.needchanloc, 'yes');
        if eachFileDiff || (needChan && isempty(opts.chanLocFile))
            steps{end+1} = name; %#ok<AGROW>
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-flight helpers

function warnIfOverwriteFiles(spec, filePaths, opts, batchCtx)
% Throws 'nestapp:cancelled' if the user declines to overwrite.
% Output files now land under <batchRoot>/data/ (typeBased) or
% <batchRoot>/<stem>/ (perInput); the batch folder is freshly named
% per-run, so legitimate collisions are rare. We still warn when the
% destination .set exists (e.g. if outputRoot was just reused).
saveIdx = find(strcmp({spec.name}, 'Save New Set'), 1);
if isempty(saveIdx); return; end

p = spec(saveIdx).params;
savenew = '';
ifn     = 'yes';
if isfield(p, 'savenew') && ischar(p.savenew) && ~strcmp(p.savenew, '[]') && ~isempty(p.savenew)
    savenew = p.savenew;
end
if isfield(p, 'includeFileName')
    ifn = p.includeFileName;
end
if isempty(savenew); return; end

nFiles   = numel(filePaths);
existing = {};
for fi = 1:nFiles
    [~, fbase] = fileparts(filePaths{fi});
    fbase = replace(replace(fbase, ' ', '_'), '-', '_');
    dataDir = outputPaths(batchCtx, 'data', fbase);
    if strcmpi(ifn, 'yes')
        outName = fullfile(dataDir, [fbase, '_', savenew, '.set']);
    else
        outName = fullfile(dataDir, [savenew, '.set']);
    end
    if exist(outName, 'file')
        [~, dispName] = fileparts(outName);
        existing{end+1} = [dispName, '.set']; %#ok<AGROW>
    end
end
if isempty(existing); return; end

msg = sprintf(['%d output file(s) already exist and will be overwritten:\n\n' ...
    '%s\n\nContinue?'], numel(existing), strjoin(existing, '\n'));
answer = uiconfirm(opts.uiFigure, msg, 'Output Files Exist', ...
    'Options', {'Continue', 'Cancel'}, ...
    'DefaultOption', 2, 'CancelOption', 2, ...
    'Icon', 'warning');
if strcmp(answer, 'Cancel')
    error('nestapp:cancelled', 'Run cancelled by user.');
end
end

function warnDeprecatedGateParams(spec)
% Scan a pipeline spec for Quality Gate steps using renamed params and
% emit one CFG log line per unique mapping so users update their saved
% pipelines. Silent when nothing deprecated is found.
aliases = deprecatedGateAliases();
seen = false(size(aliases, 1), 1);
for si = 1:numel(spec)
    if ~strcmp(spec(si).name, 'Quality Gate'), continue, end
    p = spec(si).params;
    for k = 1:size(aliases, 1)
        if seen(k), continue, end
        old = aliases{k, 1};
        if isfield(p, old) && ~isempty(p.(old)) && p.(old) ~= 0
            nestLog('CFG', ['Quality Gate param "%s" is deprecated. ' ...
                'Rename to "%s" in your saved pipeline (the run will ' ...
                'still honor the old name).'], old, aliases{k, 2});
            seen(k) = true;
        end
    end
end
end

function closeDlg(dlg)
if isfield(dlg, 'overlay') && ~isempty(dlg.overlay) && isvalid(dlg.overlay)
    delete(dlg.overlay);
elseif isfield(dlg, 'fig') && ~isempty(dlg.fig) && isvalid(dlg.fig)
    close(dlg.fig);
end
end

function writeBatchArtifacts(batchCtx, spec, pipelineName, reports, failed)
% Drop the run-level artifacts in <batchRoot>/batch (or _batch for
% perInput layout). Every section is independently try/catched so a
% rendering bug never costs the user their per-file outputs.
batchDir = outputPaths(batchCtx, 'batch');

% 1) Pipeline spec snapshot, so this run is reproducible from the
%    artifacts alone.
try
    specPath = fullfile(batchDir, 'spec.mat');
    save(specPath, 'spec', 'pipelineName');
catch err
    nestLog('CFG', 'Could not save spec.mat: %s', err.message);
end

% 2) Per-file summary CSV.
try
    csvPath = fullfile(batchDir, 'session_summary.csv');
    writeSessionSummaryCsv(csvPath, reports, failed);
catch err
    nestLog('CFG', 'Could not write session_summary.csv: %s', err.message);
end

% 3) Dashboard PNG - only when at least one report carries a Quality
%    Gate (otherwise the dashboard is empty and we save nothing).
try
    if anyReportHasGates(reports)
        pngPath = fullfile(batchDir, 'dashboard.png');
        renderDashboardToFile(reports, pngPath);
    end
catch err
    nestLog('CFG', 'Could not render dashboard PNG: %s', err.message);
end
end

function writeSessionSummaryCsv(csvPath, reports, failed)
% One row per processed file. Failed files are pulled from the
% structured failure log.
rows = cell(0, 6);
for k = 1:numel(reports)
    r = reports{k};
    [~, stem] = fileparts(r.inputFile);
    nSteps = numel(r.steps);
    nErr = 0;
    durS = 0;
    for si = 1:nSteps
        s = r.steps{si};
        if isfield(s, 'duration'), durS = durS + s.duration; end
    end
    verdict = '';
    if isfield(r, 'quality') && isfield(r.quality, 'worstVerdict')
        verdict = r.quality.worstVerdict;
    end
    rows(end+1, :) = {stem, 'ok', nSteps, nErr, durS, verdict}; %#ok<AGROW>
end
for k = 1:numel(failed)
    f = failed(k);
    [~, stem] = fileparts(f.name);
    kind = 'errored';
    if isfield(f, 'kind') && ~isempty(f.kind), kind = f.kind; end
    rows(end+1, :) = {stem, kind, NaN, 1, NaN, ''}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', ...
    {'stem','status','n_steps','n_errors','duration_s','quality_verdict'});
writetable(T, csvPath);
end

function renderDashboardToFile(reports, pngPath)
% Render the Session Quality Dashboard to an offscreen uifigure and
% export it as a PNG using exportapp (uifigure-safe).
fig = uifigure('Visible', 'off', 'Position', [100 100 1400 900]);
cleanup = onCleanup(@() delete(fig));
renderDashboardPanel(fig, reports);
drawnow;
exportapp(fig, pngPath);
end

