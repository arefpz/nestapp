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

t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel 1: channel x trial attribute heatmap ---
if opts.panels.attribMatrix
    nexttile(t, 1);
    drawAttributeMatrix(EEG, opts);
end

% --- Panel 2: ICA component topo grid (or placeholder) ---
if opts.panels.icaGrid
    nexttile(t, 2);
    drawICAGrid(EEG);
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

title(t, buildSuperTitle(EEG, opts), 'Interpreter', 'none', ...
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
    text(0.5, 0.5, sprintf('Attribute matrix unavailable\n%s', err.message), ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis off
    return
end
imagesc(SM);
colormap(gca, 'parula');
cb = colorbar;
cb.Label.String = 'log attribute';
xlabel('Trial');
ylabel('Channel');
title(sprintf('Attribute matrix (%s)', opts.attribute), 'Interpreter', 'none');

% Mark flagged channels with text on the y-axis.
nbchan = summary.nbchan;
yticks(1:max(1, floor(nbchan/16)):nbchan);
hold on
for k = 1:nbchan
    if summary.flatChanMask(k)
        plot([0.5 size(SM,2)+0.5], [k k], 'Color', [0.9 0.2 0.2], 'LineWidth', 1);
    elseif summary.satChanMask(k)
        plot([0.5 size(SM,2)+0.5], [k k], 'Color', [0.85 0.2 0.85], 'LineWidth', 1);
    end
end
hold off
end

function drawICAGrid(EEG)
metrics = computeICAQualityMetrics(EEG);
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

if ~hasChanlocs
    txt = sprintf('ICA topos unavailable (no chanlocs / no topoplot)\n%d components\nEMG: %d  Electrode: %d  EOG: %d', ...
        nComp, ...
        sum(strcmp({metrics.classification}, 'EMG')), ...
        sum(strcmp({metrics.classification}, 'Electrode')), ...
        sum(strcmp({metrics.classification}, 'EOG')));
    text(0.5, 0.5, txt, 'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis off
    title('ICA components');
    return
end

% Take over the current tile with a nested grid layout.
parentAx = gca;
parentPos = parentAx.Position;
delete(parentAx);
ax0 = axes('Position', parentPos);
title(ax0, 'ICA components (border: red=EMG, orange=Elec, yellow=EOG)');
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
    borderColor = classColor(metrics(k).classification);
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

function c = classColor(cls)
switch cls
    case 'EMG',       c = [0.85 0.20 0.20];
    case 'Electrode', c = [0.95 0.50 0.10];
    case 'EOG',       c = [0.90 0.80 0.10];
    otherwise,        c = [0.70 0.70 0.70];
end
end

function s = plural(n)
if n == 1, s = ''; else, s = 's'; end
end

function s = buildSuperTitle(EEG, opts)
parts = {};
if ~isempty(opts.stepLabel), parts{end+1} = opts.stepLabel; end
if ~isempty(opts.title),     parts{end+1} = opts.title;     end
nbchan  = getOr(EEG, 'nbchan', size(EEG.data,1));
nTrials = getOr(EEG, 'trials', max(size(EEG.data,3),1));
srate   = getOr(EEG, 'srate', NaN);
parts{end+1} = sprintf('nbchan=%d trials=%d srate=%g Hz', nbchan, nTrials, srate);

% Append ICA classification counts when available.
metrics = computeICAQualityMetrics(EEG);
if ~isempty(metrics)
    nEMG  = sum(strcmp({metrics.classification}, 'EMG'));
    nElec = sum(strcmp({metrics.classification}, 'Electrode'));
    nEOG  = sum(strcmp({metrics.classification}, 'EOG'));
    parts{end+1} = sprintf('ICA: %d EMG / %d Elec / %d EOG', nEMG, nElec, nEOG);
end
s = strjoin(parts, '  |  ');
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
