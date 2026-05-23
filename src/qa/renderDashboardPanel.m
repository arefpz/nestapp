function renderDashboardPanel(parent, reports, opts)
% RENDERDASHBOARDPANEL  Paint the Session Quality Dashboard.
%   RENDERDASHBOARDPANEL(parent, reports)
%   RENDERDASHBOARDPANEL(parent, reports, opts)
%
%   parent  : uipanel or uifigure that the dashboard renders into.
%             Any existing children are deleted first so the same
%             function can refresh the live view in place and also
%             paint an offscreen figure for PNG export.
%   reports : cell array of pipeline report structs.
%   opts    : optional struct
%     .title           override header text (default 'Session Quality Overview')
%     .onRefresh       function handle, called with no args when the
%                        in-app Refresh button is clicked. If empty the
%                        Refresh button is hidden (PNG export path).
%     .onExport        function handle for Export Dashboard PNG.
%                        If empty the button is hidden.
%     .onFailedRowClick function handle taking (filename) when a row in
%                        the failed-files table is selected; used to
%                        jump to that file's text report. Empty = no
%                        callback wired.

if nargin < 3 || ~isstruct(opts), opts = struct(); end
if ~isfield(opts, 'title'),           opts.title = 'Session Quality Overview'; end
if ~isfield(opts, 'onRefresh'),       opts.onRefresh = []; end
if ~isfield(opts, 'onExport'),        opts.onExport = []; end
if ~isfield(opts, 'onFailedRowClick'),opts.onFailedRowClick = []; end

% Clear any existing children so this works as both first-paint and refresh.
delete(allchild(parent));

verdicts = aggregateGateVerdicts(reports);
metrics  = aggregateMetricDistributions(reports);

drawHeader(parent, opts.title, verdicts);
drawHeatmap(parent, verdicts);
drawFailedTable(parent, reports, opts.onFailedRowClick);
drawHistograms(parent, metrics);
drawButtons(parent, opts.onRefresh, opts.onExport);
end

% -- header ----------------------------------------------------------------

function drawHeader(parent, titleText, verdicts)
nFiles = numel(verdicts.files);
if nFiles == 0
    sub = 'No reports yet - run a pipeline with a Quality Gate.';
else
    sub = sprintf('%d files: %d Pass / %d Marginal / %d Fail / %d Pending', ...
        nFiles, verdicts.counts.Pass, verdicts.counts.Marginal, ...
        verdicts.counts.Fail, verdicts.counts.Pending);
end
uilabel(parent, 'Text', titleText, 'FontWeight', 'bold', 'FontSize', 14, ...
    'Position', [10, parentTop(parent) - 25, parentWidth(parent) - 20, 22]);
uilabel(parent, 'Text', sub, ...
    'Position', [10, parentTop(parent) - 48, parentWidth(parent) - 20, 22], ...
    'FontColor', [0.4 0.4 0.4]);
end

% -- verdict heatmap -------------------------------------------------------

function drawHeatmap(parent, v)
W = parentWidth(parent);
H = parentTop(parent);
% Layout: left half, vertical band starting below header.
pos = [10, round(H * 0.42), round(W * 0.48), round(H * 0.42)];
ax = uiaxes(parent, 'Position', pos);
% Filenames and gate labels contain underscores; default 'tex' interpreter
% would render them as subscripts (e.g. rtmsct_214_1_SPL).
ax.TickLabelInterpreter = 'none';
title(ax, 'Verdict heatmap (files x gates)', 'Interpreter', 'none');

if isempty(v.verdicts)
    text(ax, 0.5, 0.5, 'No verdicts to plot', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis(ax, 'off');
    return
end

imagesc(ax, v.verdicts);
colormap(ax, verdictColormap());
caxis(ax, [0 4]);

ax.YTick = 1:numel(v.files);
ax.YTickLabel = v.files;
ax.XTick = 1:numel(v.gates);
ax.XTickLabel = v.gates;
ax.XTickLabelRotation = 30;
ax.YDir = 'reverse';
xlabel(ax, 'Gate');
ylabel(ax, 'File');
end

function cmap = verdictColormap()
% Index by code: 0 NotChecked, 1 Pass, 2 Marginal, 3 Fail, 4 Pending.
cmap = [ ...
    0.85 0.85 0.85;   % 0 gray
    0.20 0.70 0.30;   % 1 green
    0.95 0.80 0.20;   % 2 yellow
    0.85 0.20 0.20;   % 3 red
    0.30 0.45 0.85];  % 4 blue (Pending)
end

% -- failed-files table ----------------------------------------------------

function drawFailedTable(parent, reports, onRowClick)
W = parentWidth(parent);
H = parentTop(parent);
pos = [round(W * 0.50), round(H * 0.42), round(W * 0.48), round(H * 0.42)];

rows = collectFailures(reports);
uilabel(parent, 'Text', 'Failed / Marginal files', ...
    'FontWeight', 'bold', ...
    'Position', [pos(1), pos(2) + pos(4) + 2, pos(3), 18]);

t = uitable(parent, 'Position', pos, ...
    'ColumnName', {'File', 'Gate', 'Verdict', 'Reasons'}, ...
    'ColumnWidth', {120, 90, 70, 'auto'}, ...
    'Data', rows);

if ~isempty(onRowClick)
    t.CellSelectionCallback = @(src, ev) handleSelect(src, ev, onRowClick);
end
end

function rows = collectFailures(reports)
% One row per file. Roll up every Marginal/Fail gate's label and
% reasons so a file is never split across multiple rows: the user gets
% the full picture of why a file was flagged in one place.
rows = {};
for ri = 1:numel(reports)
    r = reports{ri};
    if ~isstruct(r) || ~isfield(r, 'quality') ...
            || ~isfield(r.quality, 'gates'), continue, end
    [~, name] = fileparts(r.inputFile);

    flaggedGates = {};
    allReasons   = {};
    worstSeen    = 'Pass';
    for gi = 1:numel(r.quality.gates)
        g = r.quality.gates{gi};
        if ~any(strcmp(g.verdict, {'Marginal', 'Fail'})), continue, end
        flaggedGates{end+1} = g.label; %#ok<AGROW>
        worstSeen = worstVerdict(worstSeen, g.verdict);
        if isfield(g, 'reasons') && ~isempty(g.reasons)
            for k = 1:numel(g.reasons)
                allReasons{end+1} = sprintf('[%s] %s', g.label, g.reasons{k}); %#ok<AGROW>
            end
        end
    end
    if isempty(flaggedGates), continue, end

    gateStr   = strjoin(flaggedGates, ', ');
    reasonStr = strjoin(allReasons,   '; ');
    rows(end+1, :) = {name, gateStr, worstSeen, reasonStr}; %#ok<AGROW>
end
end

function v = worstVerdict(a, b)
order = {'Pass', 'Marginal', 'Fail'};
ia = find(strcmp(a, order)); if isempty(ia), ia = 0; end
ib = find(strcmp(b, order)); if isempty(ib), ib = 0; end
v = order{max(ia, ib)};
end

function handleSelect(src, ev, cb)
if isempty(ev.Indices), return, end
row = ev.Indices(1);
if row < 1 || row > size(src.Data, 1), return, end
fileName = src.Data{row, 1};
try
    cb(fileName);
catch
    % Swallow callback errors - dashboard must not crash the app.
end
end

% -- per-metric histograms -------------------------------------------------

function drawHistograms(parent, metrics)
W = parentWidth(parent);
H = parentTop(parent);
yTop = round(H * 0.36);
yBot = 50;
hPanel = yTop - yBot;
if hPanel < 60 || isempty(metrics)
    if isempty(metrics)
        uilabel(parent, 'Text', 'No enabled metrics to plot', ...
            'Position', [10, yTop - 22, W - 20, 20], ...
            'FontColor', [0.5 0.5 0.5]);
    end
    return
end

uilabel(parent, 'Text', 'Per-metric distributions', 'FontWeight', 'bold', ...
    'Position', [10, yTop, W - 20, 18]);

nCols = min(4, numel(metrics));
nRows = ceil(numel(metrics) / nCols);
cellW = floor((W - 20) / nCols) - 5;
cellH = floor((hPanel - 22) / nRows) - 5;

for k = 1:numel(metrics)
    col = mod(k - 1, nCols);
    row = floor((k - 1) / nCols);
    x = 10 + col * (cellW + 5);
    y = yTop - 22 - (row + 1) * cellH - row * 5;
    ax = uiaxes(parent, 'Position', [x, y, cellW, cellH]);
    ax.TickLabelInterpreter = 'none';
    drawOneHistogram(ax, metrics(k));
end
end

function drawOneHistogram(ax, m)
vals = m.values;
if isempty(vals)
    text(ax, 0.5, 0.5, m.displayName, 'HorizontalAlignment', 'center', ...
        'Units', 'normalized');
    axis(ax, 'off');
    return
end

nBins = max(5, min(20, ceil(sqrt(numel(vals)))));
histogram(ax, vals, nBins, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none');
title(ax, m.displayName, 'FontSize', 9, 'Interpreter', 'none');
ax.FontSize = 8;

hold(ax, 'on');
yl = ylim(ax);
for t = m.absThresholds
    plot(ax, [t t], yl, 'r-', 'LineWidth', 1.2);
end
for t = m.batchCutoffs
    plot(ax, [t t], yl, 'm--', 'LineWidth', 1.2);
end
hold(ax, 'off');
end

% -- buttons ---------------------------------------------------------------

function drawButtons(parent, onRefresh, onExport)
W = parentWidth(parent);

x = W - 200;
if ~isempty(onRefresh)
    uibutton(parent, 'Text', 'Refresh', ...
        'Position', [x, 10, 90, 28], ...
        'ButtonPushedFcn', @(~,~) onRefresh());
    x = x + 100;
end
if ~isempty(onExport)
    uibutton(parent, 'Text', 'Export PNG', ...
        'Position', [x, 10, 90, 28], ...
        'ButtonPushedFcn', @(~,~) onExport());
end
end

% -- layout helpers --------------------------------------------------------

function w = parentWidth(parent)
pos = getPositionPx(parent);
w = pos(3);
end

function h = parentTop(parent)
pos = getPositionPx(parent);
h = pos(4);
end

function p = getPositionPx(parent)
% Returns the inner [width height] in pixels. uipanel and uifigure
% both expose .Position as [x y w h].
if isprop(parent, 'InnerPosition')
    p = parent.InnerPosition;
else
    p = parent.Position;
end
end
