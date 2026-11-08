
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function figPath = renderQualityFigure(EEG, outPath, opts)
% RENDERQUALITYFIGURE  4-panel QC montage for a single EEG snapshot.
%   figPath = RENDERQUALITYFIGURE(EEG, outPath, opts) writes a PNG to
%   outPath and returns the resolved path. The parent folder is created
%   if it does not exist.
%
%   opts.panels.attribMatrix   default true
%   opts.panels.icaGrid        default true (auto-collapses if no ICA)
%   opts.panels.butterfly      default true
%   opts.panels.psd            default true
%   opts.size                  [w h] in pixels, default [1600 1200]
%   opts.title                 string for the suptitle (e.g. file basename)
%   opts.stepLabel             e.g. "Step 14 / Remove ICA Components (TESA)"
%   opts.attribute             forwarded to computeAttributeMatrix
%   opts.tmsWindow             [tStart tEnd] in ms, forwarded to
%                              computeAttributeMatrix (default unset -
%                              uses computeAttributeMatrix's default)
%
%   Failure isolation: the caller (processOneFile) wraps this in try/catch
%   so a rendering bug never aborts the pipeline. Internally the function
%   uses onCleanup to close the invisible figure on any error.

if nargin < 3 || ~isstruct(opts), opts = struct(); end
if ~isfield(opts, 'panels'),        opts.panels = struct(); end
if ~isfield(opts.panels, 'attribMatrix'), opts.panels.attribMatrix = true; end
if ~isfield(opts.panels, 'icaGrid'),      opts.panels.icaGrid      = true; end
if ~isfield(opts.panels, 'butterfly'),    opts.panels.butterfly    = true; end
if ~isfield(opts.panels, 'psd'),          opts.panels.psd          = true; end
if ~isfield(opts, 'size'),       opts.size       = [1600 1200];          end
if ~isfield(opts, 'title'),      opts.title      = '';                   end
if ~isfield(opts, 'stepLabel'),  opts.stepLabel  = '';                   end
if ~isfield(opts, 'attribute'),  opts.attribute  = 'minmax_no_tms';      end
if ~isfield(opts, 'maxICs'),     opts.maxICs     = 25;                   end
if ~isfield(opts, 'icaSnapshot'),opts.icaSnapshot = struct([]);          end
if ~isfield(opts, 'exportLockDir'), opts.exportLockDir = '';             end

parentDir = fileparts(outPath);
if ~isempty(parentDir) && ~exist(parentDir, 'dir')
    mkdir(parentDir);
end

fig = outputFigure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 opts.size(1) opts.size(2)], ...
    'PaperPositionMode', 'auto');
cleanup = onCleanup(@() closeFigSafely(fig));

% Compute ICA metrics once and reuse for the topo panel and the
% suptitle. On the heuristic path this avoids running pwelch per
% component twice; on TESA / ICLabel paths it avoids the field-copy
% work too.
% Prefer the pre-rejection snapshot (all components, real reject flags)
% when the pipeline captured one; otherwise fall back to the live EEG.
icaView = buildICAView(EEG, opts);

t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel 1: channel x trial attribute heatmap ---
if opts.panels.attribMatrix
    nexttile(t, 1);
    drawAttributeMatrix(EEG, opts);
end

% --- Panel 2: ICA component topo grid (or placeholder) ---
if opts.panels.icaGrid
    nexttile(t, 2);
    drawICAGrid(icaView, opts);
end

% --- Panel 3: butterfly ---
if opts.panels.butterfly
    nexttile(t, 3);
    drawButterfly(EEG);
end

% --- Panel 4: PSD per channel ---
if opts.panels.psd
    nexttile(t, 4);
    drawPSD(EEG);
end

title(t, buildSuperTitle(EEG, opts, icaView.metrics), 'Interpreter', 'none', ...
    'FontWeight', 'bold');

% Serialize the rasterization across parallel workers: concurrent
% exportgraphics on headless workers can deadlock on shared OS graphics
% resources, which is what wedges one random file in a big batch (see
% acquireQCExportLock). The figure build above stays fully parallel; only
% this export is mutually exclusive. No lock dir (tests / interactive /
% single process) -> no locking. The lock releases as this function returns.
if ~isempty(opts.exportLockDir)
    exportLock = acquireQCExportLock(opts.exportLockDir); %#ok<NASGU>
end
exportgraphics(fig, outPath, 'Resolution', 150);
figPath = outPath;
end

% -- panel helpers ---------------------------------------------------------

function drawAttributeMatrix(EEG, opts)
attrOpts = struct('attribute', opts.attribute);
if isfield(opts, 'tmsWindow') && ~isempty(opts.tmsWindow)
    attrOpts.tmsWindow = opts.tmsWindow;
end
try
    [SM, summary] = computeAttributeMatrix(EEG, attrOpts);
catch err
    text(0.5, 0.5, sprintf('Channel x Trial Quality Map unavailable\n%s', err.message), ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis off
    return
end
% Reorder rows into an anterior-posterior montage: left-temporal,
% left para-sagittal, medial, right para-sagittal, right-temporal
% (top to bottom). Within each band, sorted from most anterior to
% most posterior. Falls back to identity ordering when chanlocs
% aren't usable.
nbchan = summary.nbchan;
order  = montageOrder(EEG, nbchan);
SM     = SM(order, :);
summary.flatChanMask = summary.flatChanMask(order);
summary.satChanMask  = summary.satChanMask(order);

% If processOneFile told us which original trials were rejected,
% reconstruct a full-width matrix with NaN columns at the rejected
% positions so the X axis still shows the original trial timeline.
[plotSM, rejectedX] = embedRejectedTrials(SM, opts);
imagesc(plotSM, 'AlphaData', ~isnan(plotSM));
colormap(gca, 'parula');
cb = colorbar;
cb.Label.String = 'Noise score (log; brighter = noisier)';
xlabel('Trial');
ylabel('Channel');
title('Channel x Trial Quality Map');
subtitle(attributeDisplayName(opts.attribute));

% Y-axis: cap at ~16 ticks and label each with the channel name from
% chanlocs when available (Fp1, Cz, ...) instead of a bare index.
% tickIdx are positions in the *reordered* SM, so look up the
% original chanlocs index via order().
step    = max(1, floor(nbchan/16));
tickIdx = 1:step:nbchan;
yticks(tickIdx);
yticklabels(channelLabels(EEG, order(tickIdx)));
% Disable the tex interpreter so labels like 'EXT_1' don't render
% with a subscript.
set(gca, 'TickLabelInterpreter', 'none');
hold on
plotWidth = size(plotSM, 2);
for k = 1:nbchan
    if summary.flatChanMask(k)
        plot([0.5 plotWidth+0.5], [k k], 'Color', [0.9 0.2 0.2], 'LineWidth', 1);
    elseif summary.satChanMask(k)
        plot([0.5 plotWidth+0.5], [k k], 'Color', [0.85 0.2 0.85], 'LineWidth', 1);
    end
end
% Vertical red bars at rejected-trial positions so you can see WHEN
% in the run the bad epochs were dropped. The heatmap columns at
% those positions are NaN gaps; the bar makes them visible.
for x = rejectedX
    plot([x x], [0.5 nbchan+0.5], 'Color', [0.85 0.2 0.2], 'LineWidth', 1.5);
end
hold off
end

function [plotSM, rejectedX] = embedRejectedTrials(SM, opts)
% Reconstruct the SM matrix on the full original-trial axis when we
% know which trials were rejected upstream. Surviving columns sit at
% their original positions; rejected columns become NaN. Returns the
% extended matrix plus the list of rejected x-positions to mark with
% a vertical bar. When no rejection metadata is available (or the
% surviving count doesn't match SM's width), returns SM unchanged.
plotSM    = SM;
rejectedX = [];
if ~isfield(opts, 'rejectedTrialIdx') || ~isfield(opts, 'originalTrials')
    return
end
nOrig = opts.originalTrials;
rej   = opts.rejectedTrialIdx;
if ~isnumeric(nOrig) || isempty(nOrig) || nOrig <= 0, return, end
if isempty(rej), return, end
rej   = unique(rej(rej >= 1 & rej <= nOrig));
keptX = setdiff(1:nOrig, rej);
if numel(keptX) ~= size(SM, 2)
    % Trial counts don't match - bail out rather than draw a
    % misaligned heatmap.
    return
end
plotSM = nan(size(SM, 1), nOrig);
plotSM(:, keptX) = SM;
rejectedX = rej(:)';
end

function drawICAGrid(view, opts)
metrics = view.metrics;
if isempty(metrics)
    text(0.5, 0.5, 'ICA not yet computed', ...
        'HorizontalAlignment', 'center', 'FontSize', 14, ...
        'Units', 'normalized');
    axis off
    title('ICA components');
    return
end

src = view.source;
% "Captured at step N" makes clear these are the components as they were
% when ICA rejection ran, not the post-cleaning survivors.
if ~isempty(view.capturedStep)
    capNote = sprintf('  -  ICs captured at step %d', view.capturedStep);
else
    capNote = '';   % no snapshot: panel reflects the live (post-clean) EEG
end

% topoplot needs chanlocs; fall back to a textual summary if absent.
hasChanlocs = ~isempty(view.chanlocs) && isfield(view.chanlocs, 'X') ...
           && ~isempty(view.chanlocs(1).X) && ~isempty(which('topoplot'));

nTotal = numel(metrics);
nKept  = sum([metrics.kept]);
nRej   = nTotal - nKept;

if ~hasChanlocs
    summary = classificationCounts(metrics);
    txt = sprintf(['ICA topos unavailable (no chanlocs / no topoplot)\n' ...
        '%d components (%d kept / %d rejected) - classified by %s\n%s'], ...
        nTotal, nKept, nRej, src, summary);
    text(0.5, 0.5, txt, 'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis off
    title(sprintf('ICA components (%s)%s', src, capNote));
    return
end

% Pick which components to draw: largest by variance, but always keep the
% rejected ones so a small artifact IC is never hidden by the cap.
order = pickDisplayOrder(view, opts.maxICs);
nShow = numel(order);
nCol  = ceil(sqrt(nShow));
nRow  = ceil(nShow / nCol);

% Take over the current tile with a nested grid layout.
parentAx = gca;
parentPos = parentAx.Position;
delete(parentAx);
ax0 = axes('Position', parentPos);
% Short single-line title - green/red borders carry the kept/rejected
% meaning, so the title stays clear of the figure's main title.
title(ax0, sprintf('ICA components (%s)%s', src, capNote), 'FontSize', 9);
axis(ax0, 'off');

innerWidth  = parentPos(3) / nCol;
innerHeight = parentPos(4) / nRow * 0.95;
yTop = parentPos(2) + parentPos(4) * 0.95;
fig  = ancestor(ax0, 'figure');

for i = 1:nShow
    k = order(i);
    [rIdx, cIdx] = ind2sub([nRow nCol], i);
    px = parentPos(1) + (cIdx-1) * innerWidth;
    py = yTop          - rIdx     * innerHeight;
    tilePos = [px py innerWidth*0.95 innerHeight*0.82];
    ax = axes('Position', tilePos);
    try
        % Through topoplotInto, not bare topoplot: topoplot sets the FIGURE
        % colormap, which turned the quality heatmap above (parula) turbo.
        topoplotInto(ax, view.icawinv(:,k), view.chanlocs, 'electrodes', 'off');
    catch
        text(0.5, 0.5, '?', 'HorizontalAlignment','center','Units','normalized','Parent',ax);
        axis(ax,'off');
    end
    % topoplot calls "axis off", which hides the axes box - so the border
    % can't live on the axes. Draw it as a figure annotation rectangle over
    % the tile: green = accepted (kept), red = rejected.
    if metrics(k).kept
        borderColor = [0.15 0.65 0.20];   % green - accepted
    else
        borderColor = [0.85 0.20 0.20];   % red - rejected
    end
    annotation(fig, 'rectangle', tilePos, 'Color', borderColor, 'LineWidth', 2);
    % Tile title: index, class, and the component's % variance (when known).
    pv = view.compSize(k);
    if isnan(pv), vtxt = ''; else, vtxt = sprintf('  %.1f%%', pv); end
    title(ax, sprintf('%d %s%s', k, metrics(k).classification, vtxt), 'FontSize', 8);
end
end

% -- ICA view construction (snapshot-preferred) ----------------------------

function view = buildICAView(EEG, opts)
% Assemble what drawICAGrid / the suptitle need. Preference order:
%   1. opts.icaSnapshot  - the full pre-rejection decomposition captured by
%      processOneFile just before pop_subcomp (all comps, real reject flags,
%      per-comp variance, capture step). This is what makes the panel show
%      every component and which were removed.
%   2. live EEG          - post-rejection survivors (legacy behaviour) when
%      no snapshot was captured (e.g. a gate before any ICA removal).
view = struct('metrics', [], 'icawinv', [], 'chanlocs', [], ...
    'source', '', 'compSize', [], 'capturedStep', []);

snap = opts.icaSnapshot;
if isValidSnapshot(snap)
    view.icawinv  = snap.icawinv;
    view.chanlocs = snap.chanlocs;
    view.source   = snapshotSource(snap);
    if isfield(snap, 'capturedStep'), view.capturedStep = snap.capturedStep; end
    [view.metrics, view.compSize] = metricsFromSnapshot(snap);
    return
end

m = computeICAQualityMetrics(EEG);
view.metrics = m;
if ~isempty(m) && isfield(EEG, 'icawinv') && ~isempty(EEG.icawinv)
    view.icawinv  = EEG.icawinv;
    if isfield(EEG, 'chanlocs'), view.chanlocs = EEG.chanlocs; end
    view.source   = m(1).source;
    view.compSize = nan(1, numel(m));   % no size info - keep natural order
end
end

function tf = isValidSnapshot(snap)
tf = isstruct(snap) && ~isempty(fieldnames(snap)) ...
  && isfield(snap, 'icawinv') && ~isempty(snap.icawinv) ...
  && isfield(snap, 'rejMask') && ~isempty(snap.rejMask);
end

function src = snapshotSource(snap)
if isfield(snap, 'source') && ~isempty(snap.source)
    src = snap.source;                       % caller-declared (e.g. 'TESA')
elseif isfield(snap, 'iclabelProbs') && ~isempty(snap.iclabelProbs)
    src = 'ICLabel';
elseif isfield(snap, 'classLabels') && ~isempty(snap.classLabels)
    src = 'Flags';
else
    src = 'ICA';
end
end

function [m, compSize] = metricsFromSnapshot(snap)
n   = size(snap.icawinv, 2);
rej = logical(reshape(snap.rejMask, 1, []));
if numel(rej) < n, rej(end+1:n) = false; end
rej = rej(1:n);

labels = snapshotLabels(snap, n);
src    = snapshotSource(snap);

if isfield(snap, 'compVarPct') && numel(snap.compVarPct) >= n
    compSize = double(reshape(snap.compVarPct(1:n), 1, []));
else
    compSize = nan(1, n);
end

m = repmat(struct('compIdx', 0, 'source', src, 'classification', '', ...
    'kept', true, 'compSize', NaN), 1, n);
for k = 1:n
    m(k).compIdx        = k;
    m(k).classification = labels{k};
    m(k).kept           = ~rej(k);
    m(k).compSize       = compSize(k);
end
end

function labels = snapshotLabels(snap, n)
% Per-component class label: ICLabel argmax if probabilities were captured,
% else the custom flag labels (AARATEP muscle / ARTIST decay), else blank.
labels = repmat({''}, 1, n);
if isfield(snap, 'iclabelProbs') && ~isempty(snap.iclabelProbs)
    probs = snap.iclabelProbs;
    if size(probs, 1) < n, probs(end+1:n, :) = 0; end
    classes = {'Brain','Muscle','Eye','Heart','Line Noise','Channel Noise','Other'};
    if isfield(snap, 'iclabelClasses') && ~isempty(snap.iclabelClasses)
        classes = snap.iclabelClasses;
    end
    for k = 1:n
        [~, idx] = max(probs(k, :));
        if idx < 1 || idx > numel(classes), idx = numel(classes); end
        labels{k} = char(classes{idx});
    end
elseif isfield(snap, 'classLabels') && ~isempty(snap.classLabels)
    cl = snap.classLabels;
    for k = 1:min(n, numel(cl))
        if ~isempty(cl{k}), labels{k} = char(cl{k}); end
    end
end
end

function order = pickDisplayOrder(view, maxICs)
% Largest-variance-first order, capped at maxICs, but every rejected
% component is always included (even if it pushes past the cap).
n  = numel(view.metrics);
sz = view.compSize;
if isempty(sz) || all(isnan(sz))
    bySize = 1:n;
else
    [~, bySize] = sort(sz, 'descend', 'MissingPlacement', 'last');
end
if n <= maxICs
    order = bySize;
    return
end

isRej = ~[view.metrics.kept];
rejBySize  = bySize(isRej(bySize));
keptBySize = bySize(~isRej(bySize));
nKeptShow  = max(0, maxICs - numel(rejBySize));
chosen     = [rejBySize, keptBySize(1:min(nKeptShow, numel(keptBySize)))];

% Lay the chosen set out largest-first for a tidy grid.
if isempty(sz) || all(isnan(sz))
    order = chosen;
else
    [~, ord2] = sort(view.compSize(chosen), 'descend', 'MissingPlacement', 'last');
    order = chosen(ord2);
end
end

function drawButterfly(EEG)
if isempty(EEG.data)
    text(0.5, 0.5, 'No data', 'HorizontalAlignment','center','Units','normalized');
    axis off; title('Butterfly'); return
end
nTrials = max(size(EEG.data, 3), 1);
if nTrials > 1
    grand = mean(EEG.data, 3);
else
    grand = EEG.data(:, :, 1);
end
if isfield(EEG, 'times') && numel(EEG.times) == size(grand, 2)
    xt = EEG.times;
    xlab = 'Time (ms)';
else
    xt = (1:size(grand, 2)) / EEG.srate * 1000;
    xlab = 'Time (ms, from start)';
end
plot(xt, grand', 'LineWidth', 0.5);
xlabel(xlab);
ylabel('Amplitude (uV)');
title(sprintf('Butterfly (mean over %d trial%s)', nTrials, plural(nTrials)));
grid on
% Tight x-limits matching the data range.
if numel(xt) >= 2
    xlim([xt(1) xt(end)]);
end
end

function drawPSD(EEG)
if isempty(EEG.data) || size(EEG.data, 2) < 256
    text(0.5, 0.5, 'PSD unavailable (segment too short)', ...
        'HorizontalAlignment','center','Units','normalized');
    axis off; title('PSD per channel'); return
end
nbchan = size(EEG.data, 1);
data2D = reshape(EEG.data, nbchan, []);    % view; never transposed (see below)
% Per-channel so we never materialise a full (samples x channels) copy, and so
% qaWelchPsd bounds the spectrum (decimated bandwidth + fixed resolution). This
% replaces an unbounded pwelch that, on raw 5 kHz multi-trial data, built ~10^5
% bins and wedged exportgraphics on a parallel worker.
[p, f, fsEff] = qaWelchPsd(data2D(1, :), EEG.srate);   % first channel sizes the output
pxx = zeros(numel(p), nbchan);
pxx(:, 1) = p;
for c = 2:nbchan
    pxx(:, c) = qaWelchPsd(data2D(c, :), EEG.srate);
end
% Bounded point count (~F_max/res points/channel) keeps the transparent overlay
% cheap to rasterize, so we keep the original per-channel alpha.
loglog(f, pxx, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.4]);
hold on
loglog(f, mean(pxx, 2), 'k', 'LineWidth', 1.5);
hold off
xlabel('Frequency (Hz)');
ylabel('Power');
title('PSD per channel (mean bold)');
grid on
if numel(f) >= 2
    xlim([max(f(2), 0.5) f(end)]);
end
% Note the working rate when the QA PSD was computed on downsampled data, so the
% panel is not mistaken for the raw spectrum (lower effective rate / bandwidth).
if fsEff < EEG.srate
    text(0.99, 0.02, sprintf('PSD computed at %g Hz (downsampled from %g Hz)', fsEff, EEG.srate), ...
        'Units','normalized', 'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
        'FontSize', 7, 'Color', [0.45 0.45 0.45], 'Interpreter','none');
end
end

% -- misc ------------------------------------------------------------------

function name = attributeDisplayName(mode)
% Plain-English label for the attribute mode, shown as the heatmap
% subtitle so users do not have to recognize the internal token.
switch mode
    case 'minmax'
        name = 'Peak-to-peak amplitude (full epoch)';
    case 'minmax_no_tms'
        name = 'Peak-to-peak amplitude (TMS pulse window excluded)';
    case 'highfreq'
        name = 'High-frequency activity (muscle / movement)';
    otherwise
        name = mode;   % unknown mode - fall back to the raw token
end
end

function s = classificationCounts(metrics)
% Group classifications by label and produce "Label: N" lines.
labels = {metrics.classification};
[u, ~, ic] = unique(labels);
counts = accumarray(ic, 1);
parts  = cell(1, numel(u));
for k = 1:numel(u)
    parts{k} = sprintf('%s: %d', u{k}, counts(k));
end
s = strjoin(parts, ', ');
end


function s = buildSuperTitle(EEG, opts, metrics)
parts = {};
if ~isempty(opts.stepLabel), parts{end+1} = opts.stepLabel; end
if ~isempty(opts.title),     parts{end+1} = opts.title;     end
nbchan  = fieldOr(EEG, 'nbchan', size(EEG.data,1));
nTrials = fieldOr(EEG, 'trials', max(size(EEG.data,3),1));
srate   = fieldOr(EEG, 'srate', NaN);
parts{end+1} = sprintf('nbchan=%d trials=%d srate=%g Hz', nbchan, nTrials, srate);

% Append ICA kept/rejected counts when available - source aware.
if ~isempty(metrics)
    nKept = sum([metrics.kept]);
    nRej  = numel(metrics) - nKept;
    parts{end+1} = sprintf('ICA (%s): %d kept / %d rejected', ...
        metrics(1).source, nKept, nRej);
end
s = strjoin(parts, '  |  ');
end

function order = montageOrder(EEG, nbchan)
% Return a 1:nbchan permutation that arranges channels by montage:
%   1. Left  temporal   (top)
%   2. Left  para-sagittal
%   3. Medial / midline
%   4. Right para-sagittal
%   5. Right temporal   (bottom)
% Within each band, sorted from most anterior to most posterior
% (high X to low X) per the EEGLAB convention X+ = nose, Y+ = left.
% Bins are defined by |Y| / max(|Y|), so the same thresholds work for
% any coordinate scale (unit sphere, mm, etc.). Falls back to the
% identity ordering when chanlocs lack X/Y, when X/Y are missing for
% any channel, or when all channels are on the midline.
order = 1:nbchan;
if ~isfield(EEG, 'chanlocs') || isempty(EEG.chanlocs), return, end
if numel(EEG.chanlocs) < nbchan, return, end
if ~all(isfield(EEG.chanlocs, {'X', 'Y'})), return, end

cl = EEG.chanlocs(1:nbchan);
if any(cellfun(@isempty, {cl.X})) || any(cellfun(@isempty, {cl.Y}))
    return
end
X = [cl.X];
Y = [cl.Y];
maxAbsY = max(abs(Y));
if maxAbsY == 0, return, end

absNormY = abs(Y) / maxAbsY;
MEDIAL_THRESH   = 0.20;   % |Y|/max(|Y|) below this counts as midline
TEMPORAL_THRESH = 0.60;   % above this counts as far-lateral / temporal

isMedial    = absNormY <  MEDIAL_THRESH;
isTemporal  = absNormY >= TEMPORAL_THRESH;
isPara      = ~isMedial & ~isTemporal;

leftTemp  = find(isTemporal & Y > 0);
leftPara  = find(isPara     & Y > 0);
medial    = find(isMedial);
rightPara = find(isPara     & Y < 0);
rightTemp = find(isTemporal & Y < 0);

order = [ ...
    sortByX(leftTemp,  X), ...
    sortByX(leftPara,  X), ...
    sortByX(medial,    X), ...
    sortByX(rightPara, X), ...
    sortByX(rightTemp, X)];

% Append any channel a band missed (e.g. Y == 0 on the right boundary)
% in their original order so SM rows still match 1:nbchan.
missing = setdiff(1:nbchan, order, 'stable');
if ~isempty(missing)
    order = [order, missing];
end
end

function s = sortByX(idx, X)
if isempty(idx), s = []; return, end
[~, ord] = sort(X(idx), 'descend');
s = reshape(idx(ord), 1, []);
end

function labels = channelLabels(EEG, idx)
% Best-effort channel labels for the given indices, using chanlocs
% when available and falling back to "ch N" otherwise. Handles every
% partial-chanlocs case (no field, empty struct, missing labels
% field, empty per-channel label) so the Y axis always renders
% something readable regardless of the input montage.
n = numel(idx);
labels = cell(1, n);
haveChanlocs = isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) ...
            && isfield(EEG.chanlocs, 'labels');
for k = 1:n
    i = idx(k);
    if haveChanlocs && i <= numel(EEG.chanlocs) && ~isempty(EEG.chanlocs(i).labels)
        labels{k} = char(EEG.chanlocs(i).labels);
    else
        labels{k} = sprintf('ch %d', i);
    end
end
end


function closeFigSafely(fig)
if ishandle(fig)
    try
        close(fig);
    catch
        % ignore - rendering already failed
    end
end
end
