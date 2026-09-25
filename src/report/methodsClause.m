
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function clause = methodsClause(stepName, params)
% METHODSCLAUSE  One journal-style clause describing a pipeline step.
%   clause = METHODSCLAUSE(stepName, params) returns a single natural-language
%   clause (lower-case start, no trailing period) describing what the step did,
%   with the parameter values a reader needs to judge or reproduce it: the
%   thresholds, windows and criteria that decide what was removed. Everything
%   else is left to the pipeline specification saved with the run, which the
%   paragraph points to.
%
%   Every registry step is covered one of three ways:
%     - a case below, which may return '' when the parameters make the step a
%       no-op (Load Channel Location without coordCheck=remove);
%     - methodsSilentSteps, for steps that never change the data (I/O, event
%       detection, labelling, analysis output, QC gates);
%     - otherwise, a generic clause naming the step and its main settings.
%       Never silence: a step that changes the data and is missing from the
%       paragraph is worse than one described plainly. MethodsProseTest fails
%       if any registry step lands here, so the generic clause is only ever
%       seen for steps this file does not know yet.
%
%   Where a parameter is left empty so that upstream code picks the value, the
%   clause states the value that code actually uses (pop_autorej: 1000 uV and
%   5 SD) rather than omitting it or inventing one.
%
%   The clauses are ordered, repeated and punctuated by pipelineProse. Only
%   established algorithm names appear (FastICA, infomax, ICLabel, SOUND,
%   RANSAC, clean_rawdata/ASR, CleanLine, SSP-SIR, EDM); nestapp's parameter
%   keys never do.
%
%   See also: pipelineProse, methodsSilentSteps, methodsNarrative, stepRegistry

    clause = '';
    p = params;
    if ~isstruct(p); p = struct(); end

    name = canonicalStepName(stepName);   % resolve legacy names from old reports
    switch name
        % ── channels ─────────────────────────────────────────────────────
        case 'Remove un-needed Channels'
            keep = cellstrOf(getf(p, 'channel', {}));
            drop = cellstrOf(getf(p, 'nochannel', {}));
            if ~isempty(keep)
                clause = sprintf('the recording was restricted to %d selected channel%s', ...
                    numel(keep), plural(numel(keep)));
            elseif ~isempty(drop)
                clause = sprintf('channel%s not used for analysis (%s) %s removed', ...
                    plural(numel(drop)), labelList(drop), wasWere(numel(drop)));
            end

        case 'Load Channel Location'
            if strcmpi(txt(p, 'coordCheck', 'off'), 'remove')
                clause = 'channels without valid electrode coordinates were removed';
            end

        case 'Remove Bad Channels (manual)'
            clause = 'bad channels were identified by visual inspection and removed';

        case 'Interactive Channel Reject (TESA)'
            clause = 'bad channels were identified by visual inspection and removed';

        case 'Remove Bad Channels'
            meas = measureName(txt(p, 'measure', 'kurt'));
            clause = sprintf('bad channels were identified by %s (> %g SD) and removed', ...
                meas, num(p, 'threshold', 5));

        case 'Detect Bad Channels (RANSAC)'
            clause = sprintf(['bad channels were identified by RANSAC (correlation ' ...
                'with the predicted signal < %g in more than %g%% of %g-s windows, ' ...
                'or high-frequency noise > %g SD'], num(p, 'corrThreshold', 0.8), ...
                100 * num(p, 'maxBrokenTime', 0.4), num(p, 'windowLenS', 5), ...
                num(p, 'noiseThreshold', 4));
            hp = num(p, 'detectHighpassHz', 1);
            if hp > 0
                clause = sprintf(['%s; detected on a %g-Hz high-passed copy of ' ...
                    'the data) and removed'], clause, hp);
            else
                clause = [clause ') and removed'];
            end

        case 'Detect Bad Channels (TESA)'
            clause = tesaBadChannelClause(p);

        case 'Interpolate Channels'
            clause = sprintf('removed channels were interpolated (%s)', txt(p, 'method', 'spherical'));

        % ── epochs and baseline ──────────────────────────────────────────
        case 'Epoching'
            tl = vec(p, 'timelim');
            if numel(tl) >= 2
                ev = 'each event';
                types = getf(p, 'types', {});
                if ~isempty(types) && any(contains(lower(string(types)), 'tms'))
                    ev = 'the TMS pulse';
                end
                clause = sprintf('the data were segmented into epochs from %s ms relative to %s', ...
                    rng2(tl, 1000), ev);
            end

        case 'Remove Baseline'
            tr = vec(p, 'timerange');
            if numel(tr) >= 2
                if tr(2) <= 0
                    clause = sprintf('the data were baseline-corrected from %s ms', rng2(tr));
                else
                    clause = 'the data were demeaned over the whole epoch';
                end
            elseif ~isempty(vec(p, 'pointrange'))
                clause = 'the data were baseline-corrected over the pre-stimulus interval';
            else
                clause = 'the data were demeaned over the whole epoch';
            end

        case 'Remove Bad Epoch'
            % pop_autorej with empty options uses its own defaults: 1000 uV,
            % 5 SD, and at most 5% of epochs per iteration.
            clause = sprintf(['epochs with extreme amplitudes (beyond +/-1000 uV) or ' ...
                'improbable activity (joint probability beyond %g SD, iterating while ' ...
                'fewer than %g%% of epochs were marked) were rejected'], ...
                num(p, 'startprob', 5), num(p, 'maxrej', 5));

        case 'Remove Bad Trials'
            clause = sprintf(['improbable epochs (joint probability beyond %g SD within ' ...
                'a channel or %g SD across channels) were marked and rejected after ' ...
                'visual confirmation'], num(p, 'localThresh', 5), num(p, 'globalThresh', 5));

        case 'Visualize EEG Data'
            st = vec(p, 'state');
            if numel(st) >= 3 && st(3) == 1
                clause = 'data marked during visual inspection were rejected';
            end

        case 'Automatic Continuous Rejection'
            fl = vec(p, 'freqlimit');
            if numel(fl) < 2; fl = [35 128]; end
            verb = 'removed';
            if strcmpi(txt(p, 'correct', 'remove'), 'blank'); verb = 'blanked'; end
            clause = sprintf(['continuous data segments with abnormal spectral power ' ...
                '(> %g dB in %g-%g Hz, %g-s windows) were %s'], ...
                num(p, 'threshold', 10), fl(1), fl(2), num(p, 'epochlength', 0.5), verb);

        % ── TMS pulse artifact ───────────────────────────────────────────
        case 'Fix TMS Pulse (TESA)'
            clause = sprintf('TMS pulse markers were re-aligned to the pulse artifact on %s', ...
                txt(p, 'elec', 'Cz'));

        case 'Remove TMS Artifacts (TESA)'
            cut = vec(p, 'cutTimesTMS');
            if numel(cut) >= 2
                clause = sprintf('the TMS pulse and early muscle artifact (%s ms) were removed', ...
                    rng2(cut));
            end

        case 'Interpolate Missing Data (TESA)'
            m = txt(p, 'interpolation', 'linear');
            clause = sprintf('the excised window was reconstructed by %s interpolation', m);
            w = vec(p, 'interpWin');
            if strcmpi(m, 'cubic') && numel(w) >= 2
                clause = sprintf('%s (fitted on %g ms before and %g ms after)', clause, w(1), w(2));
            end

        case 'Interpolate Missing Data (AR-Blend)'
            clause = sprintf(['the peri-pulse window (%g to %g ms) was reconstructed by ' ...
                'autoregressive-blended interpolation'], ...
                num(p, 'artifactStartMs', -2), num(p, 'artifactEndMs', 12));

        case 'Remove Decay Artifact'
            per = ''; if strcmpi(txt(p, 'perTrial', 'off'), 'on'); per = ' on a per-trial basis'; end
            clause = sprintf('residual TMS-evoked decay was removed%s (%g to %g ms)', ...
                per, num(p, 'artifactStartMs', -2), num(p, 'artifactEndMs', 12));

        case 'Median Filter 1D'
            w = vec(p, 'timeWin'); if numel(w) < 2; w = [-2 2]; end
            what = 'the TMS pulse';
            if strcmpi(txt(p, 'event_type', 'spike'), 'spike'); what = 'each spike artifact'; end
            clause = sprintf('a %g-point median filter was applied around %s (%s ms)', ...
                num(p, 'mdorder', 5), what, rng2(w));

        case 'Fit Artifact Model (TESA)'
            fitW = vec(p, 'tFitPosMs'); if numel(fitW) < 2; fitW = [0 60]; end
            clause = sprintf(['a parametric model of the TMS-evoked artifact (fitted from ' ...
                '%s ms after the pulse) was subtracted'], rng2(fitW));

        case 'SSP SIR'
            tr = vec(p, 'timeRange'); if numel(tr) < 2; tr = [5 50]; end
            clause = sprintf(['TMS-evoked muscle artifact was suppressed with SSP-SIR ' ...
                '(artifact subspace estimated from %s ms%s)'], rng2(tr), ssPcText(getf(p, 'PC', [])));

        case 'Find Artifacts EDM (TESA)'
            w = vec(p, 'tmsMuscleWin'); if numel(w) < 2; w = [11 50]; end
            clause = sprintf(['artifactual components were identified with the enhanced ' ...
                'deflation method (EDM; TMS-evoked muscle threshold %g over %s ms) and removed'], ...
                num(p, 'tmsMuscleThresh', 10), rng2(w));

        % ── resampling, detrending, filtering ────────────────────────────
        case 'Re-Sample'
            clause = sprintf('the data were downsampled to %g Hz', num(p, 'freq', 1000));

        case 'De-Trend Epoch'
            clause = sprintf('a %s trend was removed from each epoch', polyName(num(p, 'npoly', 1)));

        case 'TESA De-Trend'
            w = vec(p, 'timeWin'); if numel(w) < 2; w = [11 500]; end
            kind = lower(txt(p, 'detrend', 'linear'));
            if strcmp(kind, 'double'); kind = 'double-exponential'; end
            clause = sprintf('a %s trend fitted over %s ms was removed', kind, rng2(w));

        case 'Robust Detrend (TESA)'
            clause = sprintf('a robust %s trend (%s) was removed from each epoch', ...
                polyName(num(p, 'polyOrder', 1)), robustFitText(p));

        case 'Robust Demean (TESA)'
            clause = sprintf('each epoch was robustly demeaned (%s)', robustFitText(p));

        case 'Frequency Filter'
            clause = firClause(p);

        case 'Frequency Filter (TESA)'
            clause = tesaFilterClause(p);

        case 'Modified Bandpass Filter (TESA)'
            clause = modifiedFilterClause(p);

        case 'Frequency Filter (CleanLine)'
            clause = sprintf('line noise was removed with CleanLine (%s Hz)', ...
                joinHz(numVec(p, 'linefreqs', [60 120])));

        % ── continuous cleaning ──────────────────────────────────────────
        case {'Automatic Cleaning Data', 'Clean Artifacts'}
            clause = asrClause(p);

        case 'Source-Informed Sensor Cleaning (SOUND)'
            clause = sprintf(['sensor noise was suppressed using the SOUND algorithm ' ...
                '(lambda = %g, %g iterations)'], num(p, 'lambdaValue', 0.2), num(p, 'iter', 10));

        % ── ICA ──────────────────────────────────────────────────────────
        case {'Run TESA ICA', 'Run ICA (FastICA)'}
            clause = 'the data were decomposed into independent components by FastICA';

        case 'Run ICA (Infomax)'
            if strcmpi(txt(p, 'extended', 'on'), 'on')
                clause = 'the data were decomposed into independent components by extended infomax';
            else
                clause = 'the data were decomposed into independent components by infomax';
            end

        case 'Run ICA (Picard)'
            if strcmpi(txt(p, 'mode', 'standard'), 'ortho')
                clause = 'the data were decomposed into independent components by Picard (Picard-O)';
            else
                clause = 'the data were decomposed into independent components by Picard';
            end

        case 'Flag ICA Components for Rejection'
            conds = iclabelConditions(p);
            if ~isempty(conds)
                clause = sprintf(['components were flagged for removal when ICLabel ' ...
                    'classified them as %s'], listJoin(conds, 'or'));
            end

        case 'Remove Flagged ICA Components'
            comps = numVec(p, 'components', []);
            if isempty(comps)
                clause = 'the flagged components were removed';
            elseif num(p, 'keepcomp', 0) == 1
                clause = sprintf('all components except %s were removed', compList(comps));
            else
                clause = sprintf('independent component%s %s %s removed', ...
                    plural(numel(comps)), compList(comps), wasWere(numel(comps)));
            end

        case 'Remove ICA Components (TESA)'
            dets = tesaCompDetectors(p);
            if ~isempty(dets)
                clause = sprintf('%s components were identified with TESA''s component classifier and removed', ...
                    listJoin(dets, 'and'));
                if strcmpi(txt(p, 'compCheck', 'off'), 'on')
                    clause = strrep(clause, 'and removed', 'and removed after visual review');
                end
            end

        case 'AARATEP Pipeline (whole)'
            % One clause for the whole pipeline: it is a published, named
            % method, so the citation carries the detail and enumerating its
            % internal stages here would only invite them to drift from what
            % upstream actually does.
            clause = ['the data were preprocessed with the AARATEP pipeline ' ...
                '(automated artifact rejection for TMS-EEG)'];

        case 'Re-Reference'
            clause = rerefClause(p);

        % ── user code ────────────────────────────────────────────────────
        case 'Manual Command'
            % Always the same sentence. The description field is a free-text
            % label for the pipeline editor - migration notes, reminders - and
            % is not written for a methods section, so it is never quoted; the
            % command itself is in the pipeline specification.
            clause = 'a custom processing step was applied (see the pipeline specification)';

        otherwise
            if ~any(strcmp(name, methodsSilentSteps()))
                clause = genericClause(name, p);
            end
    end
end

% ── per-step helpers ──────────────────────────────────────────────────────────

function s = tesaBadChannelClause(p)
    method = txt(p, 'detectionMethod', 'PREP_deviation');
    thr    = num(p, 'threshold', []);
    switch method
        case 'PREP_deviation'
            if isempty(thr); thr = 9; end
            how = sprintf('robust amplitude deviation (> %g robust SD)', thr);
        case {'TESA_DDWiener', 'TESA_DDWiener_PerTrial'}
            if isempty(thr); thr = 20; end
            how = sprintf('data-driven Wiener estimation (> %g MAD)', thr);
            if strcmp(method, 'TESA_DDWiener_PerTrial'); how = [how(1:end-1) ', per trial)']; end
        case 'fromASR'
            how = 'clean_rawdata/ASR';
        otherwise
            how = method;
    end
    switch lower(txt(p, 'replaceMethod', 'interpolate'))
        case 'interpolate'; fate = 'interpolated';
        case 'remove';      fate = 'removed';
        case 'nan';         fate = 'set to missing';
        otherwise;          fate = 'marked';
    end
    s = sprintf('bad channels were identified by %s and %s', how, fate);
end

function s = robustFitText(p)
% Shared by Robust Detrend and Robust Demean: the fit window, the excluded
% window, and the outlier criterion.
    w = vec(p, 'timeWindow'); if numel(w) < 2; w = [-1000 1000]; end
    s = sprintf('fitted over %s ms', rng2(w));
    ex = vec(p, 'excludeWindow');
    if numel(ex) >= 2; s = sprintf('%s excluding %s ms', s, rng2(ex)); end
    s = sprintf('%s, with samples beyond %g SD excluded as outliers', s, num(p, 'thresholdSD', 3));
end

function s = firClause(p)
% pop_eegfiltnew: zero-phase FIR unless minphase; revfilt inverts the band.
    lo = num(p, 'locutoff', 0); hi = num(p, 'hicutoff', 0);
    if isempty(lo); lo = 0; end
    if isempty(hi); hi = 0; end
    design = 'zero-phase FIR';
    if num(p, 'minphase', 0) == 1; design = 'causal minimum-phase FIR'; end
    ord = num(p, 'filtorder', 0);
    if ord > 0; design = sprintf('%s, order %g', design, ord); end
    if num(p, 'revfilt', 0) == 1 && lo > 0 && hi > 0
        s = sprintf('the data were band-stop filtered at %g-%g Hz', lo, hi);
    elseif lo > 0 && hi > 0
        s = sprintf('the data were band-pass filtered from %g to %g Hz', lo, hi);
    elseif lo > 0
        s = sprintf('the data were high-pass filtered at %g Hz', lo);
    elseif hi > 0
        s = sprintf('the data were low-pass filtered at %g Hz', hi);
    else
        s = '';
        return
    end
    s = sprintf('%s (%s)', s, design);
end

function s = tesaFilterClause(p)
% Butterworth band-pass / band-stop (TESA pop_tesa_filtbutter): 'high' is the
% high-pass edge, 'low' the low-pass edge.
    type = txt(p, 'type', 'bandpass');
    hi = num(p, 'high', 1); lo = num(p, 'low', 80); ord = ord2(num(p, 'ord', 4));
    switch lower(type)
        case 'bandpass'
            s = sprintf('the data were band-pass filtered from %g to %g Hz (%s Butterworth)', hi, lo, ord);
        case 'bandstop'
            s = sprintf('the data were band-stop filtered at %g-%g Hz (%s Butterworth)', hi, lo, ord);
        case 'highpass'
            s = sprintf('the data were high-pass filtered at %g Hz (%s Butterworth)', hi, ord);
        case 'lowpass'
            s = sprintf('the data were low-pass filtered at %g Hz (%s Butterworth)', lo, ord);
        otherwise
            s = '';
    end
end

function s = modifiedFilterClause(p)
    lo = num(p, 'lowCutoff', 0); hi = num(p, 'highCutoff', 0);
    if isempty(lo); lo = 0; end
    if isempty(hi); hi = 0; end
    stop = strcmpi(txt(p, 'filtType', 'auto'), 'bandstop');
    if lo > 0 && hi > 0 && stop
        edges = sprintf('%g-%g Hz band-stop', lo, hi);
    elseif lo > 0 && hi > 0
        edges = sprintf('%g-%g Hz band-pass', lo, hi);
    elseif lo > 0
        edges = sprintf('%g Hz high-pass', lo);
    elseif hi > 0
        edges = sprintf('%g Hz low-pass', hi);
    else
        s = '';
        return
    end
    if strcmpi(txt(p, 'filterMethod', 'butterworth'), 'butterworth')
        design = sprintf('%s Butterworth', ord2(num(p, 'filtOrder', 4)));
    else
        design = 'FIR';
    end
    s = sprintf(['the data were filtered with a modified %s filter (%s), with the ' ...
        '%g to %g ms pulse window bridged by autoregressive extrapolation to limit ' ...
        'artifact spread'], design, edges, num(p, 'artifactStartMs', -2), ...
        num(p, 'artifactEndMs', 12));
end

function s = asrClause(p)
% clean_rawdata / clean_artifacts. Every criterion can be 'off'.
    parts = {};
    fl = num(p, 'FlatlineCriterion', []);
    if ~isempty(fl); parts{end+1} = sprintf('flat for > %g s', fl); end
    ch = num(p, 'ChannelCriterion', []);
    if ~isempty(ch); parts{end+1} = sprintf('channel correlation < %g', ch); end
    ln = num(p, 'LineNoiseCriterion', []);
    if ~isempty(ln); parts{end+1} = sprintf('line noise > %g SD', ln); end
    bu = num(p, 'BurstCriterion', []);
    if ~isempty(bu); parts{end+1} = sprintf('bursts > %g SD', bu); end
    wc = num(p, 'WindowCriterion', []);
    if ~isempty(wc); parts{end+1} = sprintf('windows with > %g%% bad channels removed', 100 * wc); end
    burstFate = 'repaired';
    if strcmpi(txt(p, 'BurstRejection', 'off'), 'on'); burstFate = 'removed'; end
    s = sprintf(['bad channels were removed and artifact bursts %s using ' ...
        'clean_rawdata/ASR'], burstFate);
    if ~isempty(parts); s = sprintf('%s (%s)', s, strjoin(parts, ', ')); end
end

function conds = iclabelConditions(p)
% Each targeted ICLabel class with its probability range, e.g. 'eye (p >= 0.9)'.
% A Brain range flags components by LOW brain probability.
    conds = {};
    map = {'Eye','eye'; 'Muscle','muscle'; 'Heart','cardiac'; ...
           'LineNoise','line-noise'; 'ChannelNoise','channel-noise'; 'Other','other'};
    for i = 1:size(map, 1)
        r = numVec(p, map{i,1}, [NaN NaN]);
        if numel(r) >= 2 && ~all(isnan(r))
            conds{end+1} = sprintf('%s (%s)', map{i,2}, probRange(r)); %#ok<AGROW>
        end
    end
    b = numVec(p, 'Brain', [NaN NaN]);
    if numel(b) >= 2 && ~all(isnan(b))
        conds{end+1} = sprintf('brain with %s', probRange(b));
    end
end

function s = probRange(r)
    if r(2) >= 1
        s = sprintf('p >= %g', r(1));
    elseif r(1) <= 0
        s = sprintf('p <= %g', r(2));
    else
        s = sprintf('%g <= p <= %g', r(1), r(2));
    end
end

function dets = tesaCompDetectors(p)
% Enabled TESA compselect detectors, each with its criterion.
    dets = {};
    if isOn(p, 'tmsMuscle')
        w = vec(p, 'tmsMuscleWin'); if numel(w) < 2; w = [11 30]; end
        dets{end+1} = sprintf('TMS-evoked muscle (%s ms amplitude ratio > %g)', ...
            rng2(w), num(p, 'tmsMuscleThresh', 8));
    end
    if isOn(p, 'blink')
        dets{end+1} = sprintf('eye-blink (z > %g at %s)', num(p, 'blinkThresh', 2.5), ...
            labelList(cellstrOf(getf(p, 'blinkElecs', {'Fp1','Fp2'}))));
    end
    if isOn(p, 'move')
        dets{end+1} = sprintf('eye-movement (z > %g at %s)', num(p, 'moveThresh', 2), ...
            labelList(cellstrOf(getf(p, 'moveElecs', {'F7','F8'}))));
    end
    if isOn(p, 'muscle')
        f = vec(p, 'muscleFreqIn'); if numel(f) < 2; f = [30 100]; end
        dets{end+1} = sprintf('persistent muscle (%g-%g Hz spectral slope > %g)', ...
            f(1), f(2), num(p, 'muscleThresh', -0.31));
    end
    if isOn(p, 'elecNoise')
        dets{end+1} = sprintf('electrode-noise (z > %g)', num(p, 'elecNoiseThresh', 4));
    end
end

function s = ssPcText(pc)
% SSP-SIR artifact dimension: {'data', 90} or a fixed count, possibly stored
% as the string "{'data', 90}".
    s = '';
    if ischar(pc) || isstring(pc)
        v = regexp(char(pc), '[\d.]+', 'match', 'once');
        if contains(char(pc), 'data') && ~isempty(v)
            s = sprintf(', dimension explaining %s%% of the variance', v);
        elseif ~isempty(v)
            s = sprintf(', %s artifact component%s', v, plural(str2double(v)));
        end
    elseif iscell(pc) && numel(pc) >= 2 && isnumeric(pc{2})
        s = sprintf(', dimension explaining %g%% of the variance', pc{2});
    elseif isnumeric(pc) && isscalar(pc) && isfinite(pc)
        s = sprintf(', %g artifact component%s', pc, plural(pc));
    end
end

function s = rerefClause(p)
% A named channel reference reads as itself; anything empty/'[]'/non-char is the
% average reference.
    ref = getf(p, 'ref', 'Cz');
    if (ischar(ref) || isstring(ref)) && ~isempty(strtrim(char(ref))) ...
            && ~strcmp(strtrim(char(ref)), '[]')
        s = sprintf('the data were re-referenced to %s', char(ref));
    else
        s = 'the data were re-referenced to the common average';
    end
end

function s = genericClause(name, p)
% For a step this file has no case for: name it, with up to three of its
% scalar settings in the registry's own words. Plain, but never silent.
    s = sprintf('the "%s" step was applied', name);
    step = registryStep(name);
    if isempty(step) || ~isfield(step, 'params') || isempty(step.params); return; end
    bits = {};
    for prm = step.params(:)'
        if numel(bits) == 3; break; end
        if ~isfield(p, prm.key); continue; end
        v = p.(prm.key);
        if isnumeric(v) && isscalar(v) && isfinite(v)
            val = sprintf('%g', v);
        elseif (ischar(v) || isstring(v)) && ~isempty(char(v)) && numel(char(v)) <= 20
            val = char(v);
        else
            continue
        end
        unit = ''; if ~isempty(prm.unit); unit = [' ' prm.unit]; end
        bits{end+1} = sprintf('%s %s%s', lower(prm.friendlyName), val, unit); %#ok<AGROW>
    end
    if ~isempty(bits); s = sprintf('%s (%s)', s, strjoin(bits, ', ')); end
end

function step = registryStep(name)
    persistent reg
    if isempty(reg); reg = stepRegistry(); end
    step = reg(strcmp({reg.name}, name));
    if ~isempty(step); step = step(1); end
end

% ── parameter access ──────────────────────────────────────────────────────────
% Saved pipelines and UITable edits deliver values as numbers, chars, strings,
% 'off', or '[]', so every read goes through one of these.

function v = getf(p, key, dflt)
    if isfield(p, key); v = p.(key); else; v = dflt; end
end

function v = num(p, key, dflt)
% A finite numeric scalar, or dflt when absent / empty / 'off' / NaN / not a
% number.
    v = getf(p, key, dflt);
    if ischar(v) || isstring(v)
        v = str2double(v);
    elseif iscell(v) || isempty(v) || ~isnumeric(v) && ~islogical(v)
        v = dflt;
    end
    if ~isempty(v); v = double(v(1)); end
    if ~isempty(v) && isnan(v); v = dflt; end
end

function v = numVec(p, key, dflt)
% A numeric row, or dflt when absent / not numeric / containing NaN. dflt
% itself may be [NaN NaN], which is how ICLabel marks a class as unused.
    v = getf(p, key, dflt);
    if ischar(v) || isstring(v); v = str2num(char(v)); end %#ok<ST2NM>
    if ~isnumeric(v) && ~islogical(v) || isempty(v) || any(isnan(v(:))); v = dflt; end
    v = double(v(:))';
end

function v = vec(p, key)
    v = numVec(p, key, []);
end

function s = txt(p, key, dflt)
% A char value, or dflt when absent / not text / empty / the literal '[]'.
    s = getf(p, key, dflt);
    if isstring(s) && isscalar(s); s = char(s); end
    if ~ischar(s) || isempty(strtrim(s)) || strcmp(strtrim(s), '[]'); s = dflt; end
end

function tf = isOn(p, key)
    tf = strcmpi(txt(p, key, 'off'), 'on');
end

function c = cellstrOf(v)
    if ischar(v) || isstring(v)
        v = strtrim(char(v));
        if isempty(v) || strcmp(v, '[]'); c = {}; else; c = regexp(v, '[^\s,;{}'']+', 'match'); end
    elseif iscell(v)
        c = cellfun(@char, v(:)', 'UniformOutput', false);
        c = c(~cellfun(@isempty, c));
    else
        c = {};
    end
end

% ── formatting helpers ────────────────────────────────────────────────────────

function s = rng2(v, scale)
% Format a [a b] window as "a to b", optionally scaled (e.g. s -> ms).
    if nargin < 2; scale = 1; end
    v = double(v(:))' * scale;
    if numel(v) >= 2
        s = sprintf('%g to %g', v(1), v(2));
    elseif ~isempty(v)
        s = sprintf('%g', v(1));
    else
        s = '';
    end
end

function s = ord2(n)
% Small ordinal as a word ("fourth-order").
    words = {'first-order','second-order','third-order','fourth-order', ...
             'fifth-order','sixth-order'};
    if isempty(n); n = 4; end
    n = round(n);
    if n >= 1 && n <= numel(words); s = words{n}; else; s = sprintf('order-%g', n); end
end

function s = polyName(n)
    switch n
        case 1; s = 'linear';
        case 2; s = 'quadratic';
        otherwise; s = [ord2(n) ' polynomial'];
    end
end

function s = measureName(m)
    switch lower(m)
        case 'kurt'; s = 'kurtosis';
        case 'spec'; s = 'spectral power';
        case 'prob'; s = 'probability';
        otherwise;   s = m;
    end
end

function s = joinHz(v)
% "60 and 120" from [60 120].
    s = listJoin(arrayfun(@(x) sprintf('%g', x), v(:)', 'UniformOutput', false), 'and');
end

function s = labelList(c)
% "TP9 and TP10"; long lists are counted rather than spelled out.
    if numel(c) > 6
        s = sprintf('%d channels', numel(c));
    else
        s = listJoin(c, 'and');
    end
end

function s = compList(v)
    s = listJoin(arrayfun(@(x) sprintf('%g', x), v(:)', 'UniformOutput', false), 'and');
end

function s = wasWere(n)
    if n == 1; s = 'was'; else; s = 'were'; end
end

function s = listJoin(c, conj)
% Oxford-comma join of a cellstr: {a} -> "a"; {a,b} -> "a and b";
% {a,b,c} -> "a, b, and c".
    c = c(:)';
    switch numel(c)
        case 0; s = '';
        case 1; s = c{1};
        case 2; s = sprintf('%s %s %s', c{1}, conj, c{2});
        otherwise
            s = sprintf('%s, %s %s', strjoin(c(1:end-1), ', '), conj, c{end});
    end
end
