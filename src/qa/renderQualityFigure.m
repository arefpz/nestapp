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

parentDir = fileparts(outPath);
if ~isempty(parentDir) && ~exist(parentDir, 'dir')
    mkdir(parentDir);
end

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 opts.size(1) opts.size(2)], ...
    'PaperPositionMode', 'auto');
cleanup = onCleanup(@() closeFigSafely(fig));

% Compute ICA metrics once and reuse for the topo panel and the
% suptitle. On the heuristic path this avoids running pwelch per
% component twice; on TESA / ICLabel paths it avoids the field-copy
% work too.
icaMetrics = computeICAQualityMetrics(EEG);

t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel 1: channel x trial attribute heatmap ---
if opts.panels.attribMatrix
    nexttile(t, 1);
    drawAttributeMatrix(EEG, opts);
end

% --- Panel 2: ICA component topo grid (or placeholder) ---
if opts.panels.icaGrid
    nexttile(t, 2);
    drawICAGrid(EEG, icaMetrics);
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

title(t, buildSuperTitle(EEG, opts, icaMetrics), 'Interpreter', 'none', ...
    'FontWeight', 'bold');

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

function drawICAGrid(EEG, metrics)
if isempty(metrics)
    text(0.5, 0.5, 'ICA not yet computed', ...
        'HorizontalAlignment', 'center', 'FontSize', 14, ...
        'Units', 'normalized');
    axis off
    title('ICA components');
    return
end

% spectopo / topoplot need chanlocs; fall back to a textual summary if absent.
hasChanlocs = isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) ...
           && isfield(EEG.chanlocs, 'X') && ~isempty(EEG.chanlocs(1).X) ...
           && ~isempty(which('topoplot'));

nComp = numel(metrics);
nCol  = ceil(sqrt(nComp));
nRow  = ceil(nComp / nCol);

nKept = sum([metrics.kept]);
nRej  = nComp - nKept;
src   = metrics(1).source;

if ~hasChanlocs
    summary = classificationCounts(metrics);
    txt = sprintf(['ICA topos unavailable (no chanlocs / no topoplot)\n' ...
        '%d components (%d kept / %d rejected) - classified by %s\n%s'], ...
        nComp, nKept, nRej, src, summary);
    text(0.5, 0.5, txt, 'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis off
    title(sprintf('ICA components (%s)', src));
    return
end

% Take over the current tile with a nested grid layout.
parentAx = gca;
parentPos = parentAx.Position;
delete(parentAx);
ax0 = axes('Position', parentPos);
title(ax0, sprintf('ICA components (%s)  -  border colored by class; gray = kept', src));
axis(ax0, 'off');

innerWidth  = parentPos(3) / nCol;
innerHeight = parentPos(4) / nRow * 0.95;   % leave space for tile title
yTop = parentPos(2) + parentPos(4) * 0.95;

for k = 1:nComp
    [rIdx, cIdx] = ind2sub([nRow nCol], k);
    px = parentPos(1) + (cIdx-1) * innerWidth;
    py = yTop          - rIdx     * innerHeight;
    ax = axes('Position', [px py innerWidth*0.95 innerHeight*0.85]);
    try
        topoplot(EEG.icawinv(:,k), EEG.chanlocs, 'electrodes', 'off');
    catch
        text(0.5, 0.5, '?', 'HorizontalAlignment','center','Units','normalized','Parent',ax);
        axis(ax,'off');
    end
    borderColor = classColor(metrics(k).classification, metrics(k).kept);
    set(ax, 'XColor', borderColor, 'YColor', borderColor, ...
            'LineWidth', 2, 'Box', 'on', 'XTick', [], 'YTick', []);
    title(ax, sprintf('%d %s', k, metrics(k).classification), 'FontSize', 8);
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
data2D = reshape(EEG.data, nbchan, []);
[pxx, f] = pwelch(data2D', [], [], [], EEG.srate);
loglog(f, pxx, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.4]);
hold on
loglog(f, mean(pxx, 2), 'k', 'LineWidth', 1.5);
hold off
xlabel('Frequency (Hz)');
ylabel('Power');
title('PSD per channel (mean bold)');
grid on
if numel(f) >= 2
    xlim([max(f(2), 0.5) min(f(end), EEG.srate/2)]);
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

function c = classColor(cls, kept)
% Border color per component classification. Source-agnostic - merges
% TESA, ICLabel, and Heuristic category names into a small palette so
% the visualization reads the same regardless of which classifier ran.
%   kept = true  -> always gray, regardless of label (Brain / Keep / OK).
%   kept = false -> color by artifact family:
%                     muscle/EMG          -> red
%                     electrode noise     -> orange
%                     eye / blink / EOG   -> yellow
%                     other / unknown     -> magenta
GRAY    = [0.70 0.70 0.70];
RED     = [0.85 0.20 0.20];
ORANGE  = [0.95 0.50 0.10];
YELLOW  = [0.90 0.80 0.10];
MAGENTA = [0.80 0.30 0.80];

if nargin >= 2 && kept
    c = GRAY;
    return
end

switch cls
    case {'EMG', 'Muscle', 'TMS Muscle'},                       c = RED;
    case {'Electrode', 'Elec Noise', 'Channel Noise',  ...
          'Line Noise'},                                        c = ORANGE;
    case {'EOG', 'Blink', 'Eye', 'Eye Move'},                   c = YELLOW;
    case {'Heart', 'Sensory', 'Reject', 'Other'},               c = MAGENTA;
    case {'OK', 'Keep', 'Brain'},                               c = GRAY;
    otherwise,                                                  c = MAGENTA;
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

function s = plural(n)
if n == 1, s = ''; else, s = 's'; end
end

function s = buildSuperTitle(EEG, opts, metrics)
parts = {};
if ~isempty(opts.stepLabel), parts{end+1} = opts.stepLabel; end
if ~isempty(opts.title),     parts{end+1} = opts.title;     end
nbchan  = getOr(EEG, 'nbchan', size(EEG.data,1));
nTrials = getOr(EEG, 'trials', max(size(EEG.data,3),1));
srate   = getOr(EEG, 'srate', NaN);
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

function v = getOr(s, field, default)
if isfield(s, field) && ~isempty(s.(field))
    v = s.(field);
else
    v = default;
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
