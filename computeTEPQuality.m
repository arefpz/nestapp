function quality = computeTEPQuality(EEG, EEGraw, report)
% COMPUTETEPQUALITY  Five-axis TEP quality vector (metric spec v1.0).
%
%   quality = COMPUTETEPQUALITY(EEG, EEGraw, report)
%
%   Returns a struct with one sub-struct per quality axis plus a version
%   field. Each axis struct contains: value, range, interpretation, status
%   ('ok' | 'warning' | 'skipped'), status_detail, and metadata.
%
%   Axes
%     .retention         Axis 1: effective data retention              [0,1]
%     .artifactReduction Axis 2: early-window artifact removed         [0,1]
%     .bgRestoration     Axis 3: background spectrum similarity        [0,1]
%     .reproducibility   Axis 4: clean-window split-half r             [0,1]
%     .aepLikeness       Axis 5: AEP topographic similarity (EXPMT)  [-1,1]
%
%   EEGraw must be captured immediately after the Epoching step, before ICA
%   and artifact rejection. Pass [] if unavailable — Axis 2 will be skipped.
%
%   See also: initPipelineReport, runPipeline, exportReport

METRIC_VERSION    = 'v1.0';
LATERAL_CHANNELS  = {'F7','F8','FT7','FT8','T7','T8','TP7','TP8', ...
                     'P7','P8','FC5','FC6','CP5','CP6'};
VERTEX_CHANNELS   = {'Cz','FCz','Fz','FC1','FC2'};

quality = initAxes(METRIC_VERSION);

%% Guard: must be epoched with at least 2 channels
if isempty(EEG) || ~isstruct(EEG) || ~isfield(EEG,'data') || isempty(EEG.data)
    return
end
if ~isfield(EEG,'trials') || EEG.trials < 2
    return
end
if size(EEG.data, 1) < 2
    return
end

times = EEG.times;  % 1 × nTime, milliseconds

%% Warn if pre-stim baseline is shorter than the 200 ms minimum for Axis 3
preStimMs = abs(times(1));
if preStimMs < 200
    quality.baselineWarning = sprintf( ...
        ['Pre-stimulus baseline is %.0f ms — shorter than the 200 ms ' ...
         'minimum required for Axis 3 (background restoration). ' ...
         'Axis 3 will be skipped.'], preStimMs);
end

%% Axis 1 — Effective Data Retention
quality.retention = computeRetention(EEG, report);

%% Axis 2 — Early-Window Artifact Reduction
quality.artifactReduction = computeArtifactReduction(EEG, EEGraw, times, LATERAL_CHANNELS);

%% Axis 3 — Background Restoration
quality.bgRestoration = computeBgRestoration(EEG, times, preStimMs);

%% Axis 4 — Clean-Window Reproducibility
quality.reproducibility = computeReproducibility(EEG, times, LATERAL_CHANNELS);

%% Axis 5 — AEP-Likeness Flag (experimental)
quality.aepLikeness = computeAEPLikeness(EEG, times, VERTEX_CHANNELS);

end

%% ---- Axis 1 ----------------------------------------------------------------

function ax = computeRetention(EEG, report)
ax = makeAxis('retention', [0 1]);

if ~isstruct(report) || ~isfield(report,'trials') || report.trials.original == 0
    ax.status        = 'skipped';
    ax.status_detail = 'Trial count before pipeline unavailable.';
    return
end

trialRet = EEG.trials / report.trials.original;

% Interpolated channels are reconstructed, not measured — count as rejected.
% EEG.nbchan is always equal to original after interpolation, so we use
% the rejection count from the report rather than the final channel count.
nRejCh     = report.channels.nRejected + report.channels.nInterpolated;
channelRet = max(0, (report.channels.original - nRejCh) / report.channels.original);
retention  = trialRet * channelRet;

ax.value          = retention;
ax.interpretation = sprintf( ...
    ['Retained %.0f%% of trials and %.0f%% of channels ' ...
     '(%d rejected/interpolated; combined: %.0f%%).'], ...
    trialRet*100, channelRet*100, nRejCh, retention*100);
ax.metadata.trialRetention   = trialRet;
ax.metadata.channelRetention = channelRet;
ax.metadata.nChannelsRejected = nRejCh;
end

%% ---- Axis 2 ----------------------------------------------------------------

function ax = computeArtifactReduction(EEG, EEGraw, times, lateralChannels)
ax = makeAxis('artifactReduction', [0 1]);

if isempty(EEGraw) || ~isstruct(EEGraw) || ~isfield(EEGraw,'data') || isempty(EEGraw.data)
    ax.status        = 'skipped';
    ax.status_detail = 'Pre-cleaning epochs (EEGraw) not available.';
    return
end

timeMask = times >= 11 & times <= 30;
if ~any(timeMask)
    ax.status        = 'skipped';
    ax.status_detail = 'Time window 11-30 ms not present in epoch.';
    return
end

% Identify lateral channels present in both datasets
cleanNames = upper({EEG.chanlocs.labels});
rawNames   = upper({EEGraw.chanlocs.labels});
latUpper   = upper(lateralChannels);
commonLat  = intersect(intersect(cleanNames, rawNames), latUpper);

if numel(commonLat) < 2
    % Fall back to top 25% of channels by absolute y-coordinate (left-right axis)
    if isfield(EEG.chanlocs, 'Y') && ~isempty([EEG.chanlocs.Y])
        absY   = abs([EEG.chanlocs.Y]);
        thresh = prctile(absY(absY > 0), 75);
        fbIdx  = find(absY >= thresh);
        if numel(fbIdx) >= 2
            commonLat        = cleanNames(fbIdx);
            ax.status        = 'warning';
            ax.status_detail = 'Standard lateral channel names not found; used top 25% by y-coordinate.';
        end
    end
    if numel(commonLat) < 2
        ax.status        = 'skipped';
        ax.status_detail = 'Fewer than 2 lateral channels identified.';
        return
    end
end

reduction = NaN(1, numel(commonLat));
for ci = 1:numel(commonLat)
    cClean = find(strcmp(cleanNames, commonLat{ci}), 1);
    cRaw   = find(strcmp(rawNames,   commonLat{ci}), 1);
    if isempty(cClean) || isempty(cRaw), continue; end

    avgClean = mean(EEG.data(cClean,  timeMask, :), 3);
    avgRaw   = mean(EEGraw.data(cRaw, timeMask, :), 3);

    varRaw = var(avgRaw(:));
    if varRaw == 0, continue; end
    reduction(ci) = max(0, min(1, 1 - var(avgClean(:)) / varRaw));
end

reduction = reduction(~isnan(reduction));
if isempty(reduction)
    ax.status        = 'skipped';
    ax.status_detail = 'No valid lateral channels after variance check.';
    return
end

ax.value          = median(reduction);
ax.interpretation = sprintf( ...
    ['Reduced early-window (11-30 ms) variance by %.0f%% at lateral channels. ' ...
     'Very low suggests residual artifact; very high may indicate over-cleaning.'], ...
    ax.value * 100);
ax.metadata.nLateralChannels = numel(reduction);
ax.metadata.windowMs         = [11 30];
end

%% ---- Axis 3 ----------------------------------------------------------------

function ax = computeBgRestoration(EEG, times, preStimMs)
ax = makeAxis('bgRestoration', [0 1]);
MIN_PRESTIM_MS = 200;

if preStimMs < MIN_PRESTIM_MS
    ax.status        = 'skipped';
    ax.status_detail = sprintf('Pre-stim baseline %.0f ms < %.0f ms minimum.', ...
        preStimMs, MIN_PRESTIM_MS);
    return
end

preStart  = max(times(1), -500);
preEnd    = -10;
postStart = 500;
postEnd   = min(times(end), 1000);

preMask  = times >= preStart & times <= preEnd;
postMask = times >= postStart & times <= postEnd;

if ~any(preMask)
    ax.status = 'skipped'; ax.status_detail = 'Pre-stim window not present.'; return
end
if ~any(postMask)
    ax.status = 'skipped'; ax.status_detail = 'Post-stim window (500-1000 ms) not present.'; return
end
if postEnd < 1000
    ax.status        = 'warning';
    ax.status_detail = sprintf('Post-stim window shortened to %.0f ms (epoch ends at %.0f ms).', ...
        postEnd, times(end));
end

% Grand average across channels then concatenate trials → single signal
sfreq    = EEG.srate;
preData  = reshape(EEG.data(:, preMask,  :), size(EEG.data,1), []);
postData = reshape(EEG.data(:, postMask, :), size(EEG.data,1), []);
preSig   = mean(preData,  1)';
postSig  = mean(postData, 1)';

winPre  = max(2, min(round(0.5 * sfreq), length(preSig)));
winPost = max(2, min(round(0.5 * sfreq), length(postSig)));

try
    [pxxPre,  f] = pwelch(preSig,  hann(winPre),  round(winPre/2),  [], sfreq);
    [pxxPost, ~] = pwelch(postSig, hann(winPost), round(winPost/2), [], sfreq);
catch ME
    ax.status        = 'skipped';
    ax.status_detail = sprintf('pwelch failed: %s', ME.message);
    return
end

fMask = f >= 5 & f <= 45;
if sum(fMask) < 2
    ax.status = 'skipped'; ax.status_detail = '5-45 Hz band not resolvable at this sample rate.'; return
end

logPre  = log(pxxPre(fMask)  + eps);
logPost = log(pxxPost(fMask) + eps);

if std(logPre) == 0 || std(logPost) == 0
    ax.status = 'skipped'; ax.status_detail = 'Constant PSD in one window.'; return
end

C    = corrcoef(logPre, logPost);
bgR  = max(0, C(1,2));

ax.value          = bgR;
ax.interpretation = sprintf( ...
    ['Background spectrum (5-45 Hz) similarity between pre-stim and ' ...
     '500-1000 ms post-pulse: %.2f. Lower values suggest residual contamination.'], bgR);
ax.metadata.preWinMs  = [preStart preEnd];
ax.metadata.postWinMs = [postStart postEnd];
end

%% ---- Axis 4 ----------------------------------------------------------------

function ax = computeReproducibility(EEG, times, lateralChannels)
ax = makeAxis('reproducibility', [-1 1]);
MIN_TRIALS = 20;
N_SPLITS   = 100;

if EEG.trials < MIN_TRIALS
    ax.status        = 'warning';
    ax.status_detail = sprintf('Only %d trials survived (minimum %d for stable split-half).', ...
        EEG.trials, MIN_TRIALS);
    return
end

timeMask = times >= 30 & times <= 80;
if ~any(timeMask)
    ax.status = 'skipped'; ax.status_detail = 'Time window 30-80 ms not present.'; return
end

chanNames = upper({EEG.chanlocs.labels});
latUpper  = upper(lateralChannels);
isLat     = ismember(chanNames, latUpper);
nonLatIdx = find(~isLat);

if isempty(nonLatIdx)
    nonLatIdx        = 1:EEG.nbchan;
    ax.status        = 'warning';
    ax.status_detail = 'No non-lateral channels found; using all channels.';
end

data = EEG.data(nonLatIdx, :, :);
nCh  = size(data, 1);

rng(42);
rVals = zeros(1, N_SPLITS);
for iter = 1:N_SPLITS
    perm  = randperm(EEG.trials);
    half1 = perm(1:floor(EEG.trials/2));
    half2 = perm(floor(EEG.trials/2)+1:end);

    avg1 = mean(data(:, timeMask, half1), 3);  % nCh × nT
    avg2 = mean(data(:, timeMask, half2), 3);

    rPerCh = NaN(1, nCh);
    for ci = 1:nCh
        if std(avg1(ci,:)) == 0 || std(avg2(ci,:)) == 0, continue; end
        C = corrcoef(avg1(ci,:)', avg2(ci,:)');
        rPerCh(ci) = C(1,2);
    end
    rVals(iter) = mean(rPerCh, 'omitnan');
end

ax.value          = median(rVals, 'omitnan');
ax.interpretation = sprintf( ...
    'Split-half reproducibility of 30-80 ms TEP (median of %d random splits): %.2f. Read alongside trial retention.', ...
    N_SPLITS, ax.value);
ax.metadata.nSplits   = N_SPLITS;
ax.metadata.nNonLatCh = numel(nonLatIdx);
ax.metadata.windowMs  = [30 80];
ax.metadata.nTrials   = EEG.trials;
end

%% ---- Axis 5 ----------------------------------------------------------------

function ax = computeAEPLikeness(EEG, times, vertexChannels)
ax = makeAxis('aepLikeness', [-1 1]);
CAVEAT = 'Cannot fully assess AEP contamination without sham data — interpret with caution.';

timeMask = times >= 90 & times <= 180;
if ~any(timeMask)
    ax.status = 'skipped'; ax.status_detail = ['Window 90-180 ms not present. ' CAVEAT]; return
end

chanNames = upper({EEG.chanlocs.labels});
vertIdx   = find(ismember(chanNames, upper(vertexChannels)));

if isempty(vertIdx)
    ax.status        = 'skipped';
    ax.status_detail = ['Vertex channels (Cz/FCz/Fz/FC1/FC2) not found in montage. ' CAVEAT];
    return
end

grandMean = mean(EEG.data, 3);
topo      = mean(grandMean(:, timeMask), 2);  % ch × 1

aepTemplate          = zeros(EEG.nbchan, 1);
aepTemplate(vertIdx) = 1;

nTopo = norm(topo);
nAEP  = norm(aepTemplate);
if nTopo == 0 || nAEP == 0
    ax.status = 'skipped'; ax.status_detail = ['Zero-norm topography. ' CAVEAT]; return
end

ax.value          = dot(topo / nTopo, aepTemplate / nAEP);
ax.status_detail  = CAVEAT;
ax.interpretation = sprintf( ...
    ['[EXPERIMENTAL] Late TEP (90-180 ms) similarity to canonical AEP topography: %.2f. ' ...
     'High values suggest auditory-dominated response. %s'], ax.value, CAVEAT);
ax.metadata.vertexChannels = vertexChannels;
ax.metadata.windowMs       = [90 180];
end

%% ---- Struct factories ------------------------------------------------------

function quality = initAxes(version)
quality.version           = version;
quality.baselineWarning   = '';
quality.retention         = makeAxis('retention',         [0  1]);
quality.artifactReduction = makeAxis('artifactReduction', [0  1]);
quality.bgRestoration     = makeAxis('bgRestoration',     [0  1]);
quality.reproducibility   = makeAxis('reproducibility',   [-1 1]);
quality.aepLikeness       = makeAxis('aepLikeness',       [-1 1]);
end

function ax = makeAxis(name, range)
ax.name           = name;
ax.value          = NaN;
ax.range          = range;
ax.interpretation = '';
ax.status         = 'ok';
ax.status_detail  = '';
ax.metadata       = struct();
end
