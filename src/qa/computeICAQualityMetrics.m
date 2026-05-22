function metrics = computeICAQualityMetrics(EEG)
% COMPUTEICAQUALITYMETRICS  Per-component ICA quality features.
%   metrics = COMPUTEICAQUALITYMETRICS(EEG) returns a struct array, one
%   entry per ICA component, with quality metrics ported from the TMSEEG
%   tmseeg_multiples_topos triage logic.
%
%   Returns an empty struct array if EEG lacks ICA fields (icaweights,
%   icasphere, icawinv) - safe to call before any ICA step has run.
%
%   metrics(k):
%     .compIdx          k (component index)
%     .emgRatio         pow(30-50 Hz) / pow(0-10 Hz) of the time course
%     .kurtosisTopo     kurtosis(EEG.icawinv(:,k))
%     .eogScore         frontal CBI weighting (0 if chanlocs lacks coords)
%     .tmsWindowVar     max peak-trough of projected component in
%                       a TMS-relevant window (0..150 ms if epoched)
%     .classification   'EMG' | 'Electrode' | 'EOG' | 'OK'
%
%   Classification precedence (first match wins): EMG > Electrode > EOG.
%
%   Reference defaults (from TMSEEG):
%     EMG ratio > 1.1
%     kurtosis  > 15
%     EOG: frontal CBI dominant on a frontal electrode set

EMG_RATIO_THRESH = 1.1;
KURTOSIS_THRESH  = 15;

metrics = struct('compIdx', {}, 'emgRatio', {}, 'kurtosisTopo', {}, ...
    'eogScore', {}, 'tmsWindowVar', {}, 'classification', {});

if ~isfield(EEG, 'icaweights') || isempty(EEG.icaweights) || ...
   ~isfield(EEG, 'icasphere')  || isempty(EEG.icasphere)  || ...
   ~isfield(EEG, 'icawinv')    || isempty(EEG.icawinv)
    return
end

nComp = size(EEG.icawinv, 2);
if nComp == 0
    return
end

% Determine the frontal electrode mask (TMSEEG's frontalelec definition).
[frontMask, frontMask2] = frontalMasks(EEG);
hasFrontal = any(frontMask) && any(frontMask2);

% Component activations - compute once for the whole batch.
icaact = computeICAActivations(EEG);   % nComp x nSamp(*nTrials)
nbchan = EEG.nbchan;

% TMS window for tmsWindowVar (only meaningful when epoched around 0 ms).
if isfield(EEG, 'times') && ~isempty(EEG.times) && size(EEG.data, 3) > 0
    tmsMask = EEG.times >= 0 & EEG.times <= 150;
else
    tmsMask = [];
end

metrics = repmat(struct('compIdx', 0, 'emgRatio', NaN, 'kurtosisTopo', NaN, ...
    'eogScore', 0, 'tmsWindowVar', NaN, 'classification', 'OK'), 1, nComp);

for k = 1:nComp
    metrics(k).compIdx = k;

    % Kurtosis of the topography (single-electrode artifact signature).
    metrics(k).kurtosisTopo = kurtosis(EEG.icawinv(:,k));

    % EMG ratio: pow(30-50) / pow(0-10) from pwelch of the activation.
    metrics(k).emgRatio = emgRatioFromActivation(icaact(k,:), EEG.srate);

    % EOG score: frontal CBI weighting (0..1). High = likely blink/EOG.
    if hasFrontal
        col = abs(EEG.icawinv(:,k));
        norm = sqrt(sum(col.^2));
        if norm > 0
            cbi = col / norm;
            frontalShare = sum(cbi(frontMask));
            frontalShare2 = sum(cbi(frontMask2));
            % Mirrors TMSEEG: requires inner-frontal to dominate outer-frontal.
            if frontalShare > frontalShare2
                metrics(k).eogScore = frontalShare;
            end
        end
    end

    % TMS-window variance: max peak-trough of icawinv * icaact in window.
    if ~isempty(tmsMask) && any(tmsMask)
        if size(EEG.data, 3) > 1
            actMat  = reshape(icaact(k,:), [], size(EEG.data, 2), size(EEG.data, 3));
            meanAct = squeeze(mean(actMat, 3));
        else
            meanAct = icaact(k, 1:size(EEG.data, 2));
        end
        proj = EEG.icawinv(:,k) * meanAct;
        windowed = proj(:, tmsMask);
        metrics(k).tmsWindowVar = max(max(windowed, [], 2) - min(windowed, [], 2));
    end

    % Classification - first match wins.
    if metrics(k).emgRatio > EMG_RATIO_THRESH
        metrics(k).classification = 'EMG';
    elseif metrics(k).kurtosisTopo > KURTOSIS_THRESH
        metrics(k).classification = 'Electrode';
    elseif metrics(k).eogScore > 0.5
        metrics(k).classification = 'EOG';
    else
        metrics(k).classification = 'OK';
    end
end
end

% -- local helpers ---------------------------------------------------------

function [front, frontOuter] = frontalMasks(EEG)
% Build inner/outer frontal masks from chanlocs cartesian coords.
% Matches TMSEEG: theta within +/-30 deg, phi within +/-30 deg = inner frontal.
front      = false(EEG.nbchan, 1);
frontOuter = false(EEG.nbchan, 1);

if ~isfield(EEG, 'chanlocs') || isempty(EEG.chanlocs)
    return
end
hasCoords = isfield(EEG.chanlocs, 'X') && isfield(EEG.chanlocs, 'Y') ...
         && isfield(EEG.chanlocs, 'Z');
if ~hasCoords
    return
end

x = [EEG.chanlocs.X];
y = [EEG.chanlocs.Y];
z = [EEG.chanlocs.Z];
if numel(x) ~= EEG.nbchan || any(cellfun(@isempty, {EEG.chanlocs.X}))
    return
end

[px, py]   = cart2sph(x, y, z);
theta      = px / pi * 180;
phi        = py / pi * 180;
front      = ((abs(theta) <= 30) & (abs(phi) <= 30))';
frontOuter = ((abs(theta) <= 60) & (abs(theta) > 30) & (abs(phi) <= 30))';
end

function act = computeICAActivations(EEG)
% Return nComp x (pnts * trials) component activations. Reuses cached
% EEG.icaact when present, otherwise computes the textbook formula.
if isfield(EEG, 'icaact') && ~isempty(EEG.icaact)
    act = reshape(EEG.icaact, size(EEG.icaact, 1), []);
    return
end
chanInd = 1:EEG.nbchan;
if isfield(EEG, 'icachansind') && ~isempty(EEG.icachansind)
    chanInd = EEG.icachansind;
end
W = EEG.icaweights * EEG.icasphere;
dataMat = reshape(EEG.data(chanInd, :, :), numel(chanInd), []);
act = W * dataMat;
end

function ratio = emgRatioFromActivation(comp, srate)
% Power ratio (30-50 Hz) / (0-10 Hz) via pwelch. Falls back to NaN if
% the signal is too short for a reliable estimate.
comp = double(comp(:));
if numel(comp) < 256
    ratio = NaN;
    return
end
[pxx, f] = pwelch(comp, [], [], [], srate);
lowMask  = f > 0  & f < 10;
highMask = f > 30 & f < 50;
if ~any(lowMask) || ~any(highMask)
    ratio = NaN;
    return
end
lowPow  = mean(pxx(lowMask));
highPow = mean(pxx(highMask));
if lowPow <= 0
    ratio = NaN;
else
    ratio = highPow / lowPow;
end
end
