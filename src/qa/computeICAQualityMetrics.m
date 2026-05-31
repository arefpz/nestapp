function metrics = computeICAQualityMetrics(EEG)
% COMPUTEICAQUALITYMETRICS  Per-component ICA classification + features.
%   metrics = COMPUTEICAQUALITYMETRICS(EEG) returns a struct array, one
%   entry per ICA component.
%
%   The classification comes from whichever classifier the pipeline
%   actually used, in this preference order:
%
%     1. TESA              EEG.icaCompClass.TESA<N>.compClass present
%                            (set by pop_tesa_compselect)
%     2. ICLabel           EEG.etc.ic_classification.ICLabel.classifications
%                            present (set by pop_iclabel)
%     3. Heuristic         neither classifier output is attached - falls
%                            back to fresh EMG-ratio / topo-kurtosis /
%                            frontal-CBI scoring (ported from TMSEEG).
%                            Slow path; only used when ICA was run but
%                            no classifier has tagged the components yet.
%
%   This avoids the divergence between what the visualization shows and
%   what the text report says: both now reflect the classifier that
%   actually drove the rejection decision.
%
%   metrics(k) common fields:
%     .compIdx          k
%     .source           'TESA' | 'ICLabel' | 'Heuristic'
%     .classification   string label (depends on source)
%     .kept             logical, true if the classifier marked this
%                       component as a keeper (Brain / Keep / OK).
%
%   metrics(k) extra fields (filled only on the Heuristic path; NaN /
%   default otherwise so callers can blindly index them):
%     .emgRatio         pow(30-50 Hz) / pow(0-10 Hz)
%     .kurtosisTopo     kurtosis(EEG.icawinv(:,k))
%     .eogScore         frontal CBI weighting
%     .tmsWindowVar     max peak-trough of projected component in
%                       0..150 ms (NaN if not epoched around 0)
%
%   Returns an empty struct array if EEG has no ICA fields at all.

metrics = struct('compIdx', {}, 'source', {}, 'classification', {}, ...
    'kept', {}, 'emgRatio', {}, 'kurtosisTopo', {}, 'eogScore', {}, ...
    'tmsWindowVar', {});

if ~hasICA(EEG)
    return
end

nComp = size(EEG.icawinv, 2);
if nComp == 0
    return
end

if hasTESAClassifications(EEG)
    metrics = fromTESA(EEG, nComp);
elseif hasICLabelClassifications(EEG)
    metrics = fromICLabel(EEG, nComp);
else
    metrics = computeHeuristic(EEG, nComp);
end
end

% -- attached-classifier paths --------------------------------------------

function metrics = fromTESA(EEG, nComp)
% Read the latest TESA round's compClass codes (1..8).
roundKeys = fieldnames(EEG.icaCompClass);
cl        = EEG.icaCompClass.(roundKeys{end});
codes     = cl.compClass(:);
if numel(codes) < nComp
    % Partial classification - pad with "Keep" so downstream indexing works.
    codes(end+1:nComp) = 1;
end

labels = tesaCodeLabels();   % map 1..8 -> string

metrics = repmat(blankMetric('TESA'), 1, nComp);
for k = 1:nComp
    code = codes(k);
    if code >= 1 && code <= numel(labels) && ~isempty(labels{code})
        label = labels{code};
    else
        label = sprintf('Code %d', code);
    end
    metrics(k).compIdx        = k;
    metrics(k).classification = label;
    metrics(k).kept           = (code == 1);
end
end

function metrics = fromICLabel(EEG, nComp)
probs = EEG.etc.ic_classification.ICLabel.classifications;
if size(probs, 1) < nComp
    probs(end+1:nComp, :) = 0;   % treat missing rows as Other (=0 across)
end

labels = {'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'};

metrics = repmat(blankMetric('ICLabel'), 1, nComp);
for k = 1:nComp
    [~, idx] = max(probs(k, :));
    if idx < 1 || idx > numel(labels)
        idx = numel(labels);   % default to 'Other'
    end
    metrics(k).compIdx        = k;
    metrics(k).classification = labels{idx};
    metrics(k).kept           = strcmp(labels{idx}, 'Brain');
end
end

% -- heuristic fallback (only when no classifier output attached) ---------

function metrics = computeHeuristic(EEG, nComp)
EMG_RATIO_THRESH = 1.1;
KURTOSIS_THRESH  = 15;

[frontMask, frontMask2] = frontalMasks(EEG);
hasFrontal = any(frontMask) && any(frontMask2);

icaact = computeICAActivation(EEG);

if isfield(EEG, 'times') && ~isempty(EEG.times) && size(EEG.data, 3) > 0
    tmsMask = EEG.times >= 0 & EEG.times <= 150;
else
    tmsMask = [];
end

metrics = repmat(blankMetric('Heuristic'), 1, nComp);
for k = 1:nComp
    metrics(k).compIdx      = k;
    metrics(k).kurtosisTopo = kurtosis(EEG.icawinv(:, k));
    metrics(k).emgRatio     = emgRatioFromActivation(icaact(k, :), EEG.srate);

    if hasFrontal
        col = abs(EEG.icawinv(:, k));
        nrm = sqrt(sum(col.^2));
        if nrm > 0
            cbi = col / nrm;
            if sum(cbi(frontMask)) > sum(cbi(frontMask2))
                metrics(k).eogScore = sum(cbi(frontMask));
            end
        end
    end

    if ~isempty(tmsMask) && any(tmsMask)
        if size(EEG.data, 3) > 1
            actMat  = reshape(icaact(k,:), [], size(EEG.data, 2), size(EEG.data, 3));
            meanAct = squeeze(mean(actMat, 3));
        else
            meanAct = icaact(k, 1:size(EEG.data, 2));
        end
        proj     = EEG.icawinv(:, k) * meanAct;
        windowed = proj(:, tmsMask);
        metrics(k).tmsWindowVar = max(max(windowed, [], 2) - min(windowed, [], 2));
    end

    % Precedence: EMG > Electrode > EOG > OK.
    if metrics(k).emgRatio > EMG_RATIO_THRESH
        metrics(k).classification = 'EMG';
        metrics(k).kept           = false;
    elseif metrics(k).kurtosisTopo > KURTOSIS_THRESH
        metrics(k).classification = 'Electrode';
        metrics(k).kept           = false;
    elseif metrics(k).eogScore > 0.5
        metrics(k).classification = 'EOG';
        metrics(k).kept           = false;
    else
        metrics(k).classification = 'OK';
        metrics(k).kept           = true;
    end
end
end

% -- small helpers --------------------------------------------------------

function m = blankMetric(source)
m = struct('compIdx', 0, 'source', source, 'classification', '', ...
    'kept', true, 'emgRatio', NaN, 'kurtosisTopo', NaN, ...
    'eogScore', 0, 'tmsWindowVar', NaN);
end

function tf = hasICA(EEG)
tf = isfield(EEG, 'icaweights') && ~isempty(EEG.icaweights) ...
  && isfield(EEG, 'icasphere')  && ~isempty(EEG.icasphere) ...
  && isfield(EEG, 'icawinv')    && ~isempty(EEG.icawinv);
end

function tf = hasTESAClassifications(EEG)
tf = isfield(EEG, 'icaCompClass') && isstruct(EEG.icaCompClass) ...
  && ~isempty(fieldnames(EEG.icaCompClass));
if ~tf, return, end
% Verify the latest round has the field we need.
roundKeys = fieldnames(EEG.icaCompClass);
cl        = EEG.icaCompClass.(roundKeys{end});
tf        = isfield(cl, 'compClass') && ~isempty(cl.compClass);
end

function tf = hasICLabelClassifications(EEG)
tf = isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification') ...
  && isfield(EEG.etc.ic_classification, 'ICLabel') ...
  && isfield(EEG.etc.ic_classification.ICLabel, 'classifications') ...
  && ~isempty(EEG.etc.ic_classification.ICLabel.classifications);
end

function labels = tesaCodeLabels()
% TESA compClass code -> display label.
% Codes 1..8 mirror the codes pop_tesa_compselect emits. Index by code:
labels = cell(1, 8);
labels{1} = 'Keep';
labels{2} = 'Reject';
labels{3} = 'TMS Muscle';
labels{4} = 'Blink';
labels{5} = 'Eye Move';
labels{6} = 'Muscle';
labels{7} = 'Elec Noise';
labels{8} = 'Sensory';
end

function [front, frontOuter] = frontalMasks(EEG)
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

function ratio = emgRatioFromActivation(comp, srate)
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
