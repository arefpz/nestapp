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

% Pre-flight overwrite check.
if getpref('nestapp', 'suppressEEGLABDialogs', true)
    warnIfOverwriteFiles(spec, filePaths, opts);
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
        if ~isempty(nestappSrc),    addpath(nestappSrc);    end
        if ~isempty(eeglabGenpath), addpath(eeglabGenpath); end
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
failed    = struct('fi', {}, 'name', {}, 'step', {}, 'stepName', {}, 'message', {});

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

    nestLog('PAR', 'Submitting %d futures...', nFiles);
    for fi = 1:nFiles
        fOpts           = wOpts;
        fOpts.fileIndex = fi;
        futures(fi) = parfeval(@processOneFile, 2, spec, filePaths{fi}, fOpts); %#ok<AGROW>
    end

    % Poll until all futures finish or user cancels.
    while true
        pause(0.25); drawnow;
        if ~isvalid(dlg.fig) || dlg.fig.UserData.cancelRequested
            nestLog('PAR', 'Cancel requested - cancelling futures');
            cancel(futures);
            cancelled = true;
            % Wait for workers to reach a terminal state before closing the
            % dialog.  cancel() is asynchronous - workers may still be
            % mid-EEGLAB-call and need time to wind down.
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
        if all(strcmp(states, 'finished') | strcmp(states, 'failed')); break; end
    end
    nestLog('PAR', 'Poll loop exited (cancelled=%d)', cancelled);

    finalStates = {futures.State};
    for fi = 1:nFiles
        if strcmp(finalStates{fi}, 'finished')
            [reports{fi}, ~] = fetchOutputs(futures(fi));
        elseif strcmp(finalStates{fi}, 'failed') && ~cancelled
            % Only record genuine pre-cancel failures; cancel-induced ones are expected.
            [~, fname, fext] = fileparts(filePaths{fi});
            rec = parseFailure(fi, [fname fext], futures(fi).Error.message);
            failed(end+1) = rec; %#ok<AGROW>
            logFileFailure('PAR', rec);
        end
    end

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
        [pd, ~, ~]    = fileparts(filePaths{fi});
        [summaries{fi}, ~] = exportReport(reports{fi}, [pd, filesep]);
    end
end
allReports   = reports(~cellfun(@isempty, reports));
allSummaries = summaries(~cellfun(@isempty, summaries));

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

function rec = parseFailure(fi, name, msg)
% Parse a processOneFile error message into a structured failure record.
% Messages from nestapp:stepFailed look like:
%   "Step 21 (Remove ICA Components (TESA)) failed: <root cause>"
% Greedy capture on the step name because it can itself contain parens.
rec = struct('fi', fi, 'name', name, 'step', '', 'stepName', '', 'message', msg);
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
% Post-run dialog: list failed files and let the user decide whether to
% accept the partial results or abandon the whole run.
nFailed = numel(failed);
nOk     = nFiles - nFailed;
MAX_SHOW = 10;
nShow   = min(nFailed, MAX_SHOW);
lines   = cell(1, nShow);
for k = 1:nShow
    f = failed(k);
    if isempty(f.step)
        lines{k} = sprintf('  %s -- %s', f.name, oneline(f.message));
    else
        lines{k} = sprintf('  %s -- step %s (%s)', f.name, f.step, f.stepName);
    end
end
listText = strjoin(lines, newline);
if nFailed > MAX_SHOW
    listText = [listText, sprintf('\n  ... and %d more', nFailed - MAX_SHOW)];
end
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

function warnIfOverwriteFiles(spec, filePaths, opts)
% Throws 'nestapp:cancelled' if the user declines to overwrite.
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
    [fdir, fbase] = fileparts(filePaths{fi});
    if strcmpi(ifn, 'yes')
        stem    = replace(replace(fullfile(fdir, fbase), ' ', '_'), '-', '_');
        outName = [stem, '_', savenew, '.set'];
    else
        outName = fullfile(fdir, [savenew, '.set']);
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

function closeDlg(dlg)
if isfield(dlg, 'overlay') && ~isempty(dlg.overlay) && isvalid(dlg.overlay)
    delete(dlg.overlay);
elseif isfield(dlg, 'fig') && ~isempty(dlg.fig) && isvalid(dlg.fig)
    close(dlg.fig);
end
end

