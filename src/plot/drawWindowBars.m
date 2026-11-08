% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function info = drawWindowBars(parent, res, opts)
% DRAWWINDOWBARS  One panel per window of interest, bars by group with intervals.
%   info = DRAWWINDOWBARS(parent, res, opts) draws the quantification view: for
%   each window, a bar per group carrying the group mean of the per-subject
%   window measure, with a confidence interval on it.
%
%   info returns .est (1-by-nWindows cell of the per-window estimates) and
%   .axes (every axes made).
%
%   This is the figure the waveform cannot give you. A curve with a band shows
%   that two groups differ somewhere; a bar per window says by how much, in one
%   number per component, with an interval that states how well that number is
%   pinned. It is also the shape a reader of a paper expects a TEP component
%   comparison to arrive in.
%
%   THE MEASURE IS PER SUBJECT, THEN AVERAGED. Each subject's own curve is
%   measured over the window, and the group mean is taken across those subject
%   values - never the window measure of the group-mean curve. For a mean the
%   two happen to coincide; for a PEAK they do not, and the difference is not
%   small: a peak found on an average is attenuated by every millisecond of
%   latency jitter between subjects, so it systematically understates a
%   component that is present in everyone at slightly different times.
%
%   The arithmetic is computeWindowMeasures', the same function tepWindowTable
%   and the Measures CSV use, so a bar here IS the number that leaves in the
%   export. And the interval comes from curveInterval, the same function the
%   waveform band comes from, so a bar and the band above it cannot disagree
%   about the design - paired intervals stay Cousineau-Morey corrected here too,
%   which a hand-rolled std/sqrt(n) in this file would have quietly lost.
%
%   Each panel autoscales SEPARATELY. The comparison this figure exists to make
%   is between groups within a window; components differ from each other by an
%   order of magnitude, so one shared axis would flatten P30 into the baseline
%   to leave headroom for N100. The zero line is drawn because the measures are
%   signed and a bar chart with a hidden zero is a misleading one.
%
%   opts:
%     .windows  window structs; omit for defaultTEPComponentDefs. Supplying an
%               EMPTY set means no windows, and draws nothing
%     .measure  'mean' (default) | 'peak' | 'area'
%     .mode     'TEP' | 'GMFP' | 'LMFP', for the axis label
%     .level    confidence level. Unset follows the level the run was
%               computed at (res.info.level), NOT a literal here - see below
%     .colors   group colours, default groupColors(nGroups)
%     .barPatterns  one texture per group (see hatchBars), default none
%     .axesFcn  @(parent, position) -> axes, as drawTEPTopo
%
%   THE LEVEL IS NOT DEFAULTED HERE. It used to be a literal 0.95, which
%   matched what groupCurves had actually used only by coincidence: the two
%   numbers were independent, so a run at another level would have drawn these
%   intervals at 0.95 while the status line and the exported figure's footer
%   both stated the run's level. Unset now follows res.info.level, and falls
%   through to curveInterval's own default for a result predating that field -
%   so the number lives in exactly one place at each layer.
%
%   See also: computeWindowMeasures, curveInterval, exploreMeasures, drawTEPTopo

if nargin < 3; opts = struct(); end
nG = numel(res.groups);

% .windows is settled before fillDefaults for the same reason drawTEPTopo does
% it: omitted means "use the standard components", supplied-empty means "there
% are none", and a defaults struct cannot tell those apart.
if ~isfield(opts, 'windows'); opts.windows = defaultTEPComponentDefs(); end
opts = fillDefaults(opts, struct( ...
    'measure', 'mean', 'mode', 'TEP', 'level', [], ...
    'colors', groupColors(max(nG, 1)), 'barPatterns', {{}}, ...
    'axesFcn', @(p, pos) uiaxes(p, 'Position', pos)));

% [] all the way down to curveInterval, which owns the default. An explicit
% draw option wins; failing that the run's own level; failing that whatever
% curveInterval decides.
level = opts.level;
if isempty(level) && isfield(res, 'info') && isfield(res.info, 'level')
    level = res.info.level;
end

info = struct('est', {{}}, 'axes', gobjects(0));
if nG == 0 || isempty(res.time); return; end

w  = opts.windows;
nW = numel(w);
if nW == 0; return; end

% ── per-subject measures, then the group estimate, for every window ──────
est = cell(1, nW);
for k = 1:nW
    perGroup = cell(1, nG);
    for g = 1:nG
        perGroup{g} = subjectMeasures(res.groups(g).curves, res.time, ...
                                      w(k), opts.measure);
    end
    % Each subject contributes ONE number, so this is curveInterval over a
    % one-sample "curve" - which is exactly what it is, and keeps the paired
    % normalisation and Morey correction identical to the waveform's.
    est{k} = curveInterval(perGroup, res.design, level);
end
info.est = est;

% ── layout ──────────────────────────────────────────────────────────────
P = parentRect(parent);
s = chromeScale(P);
PAD    = max(4,  round(12 * s));
LEFT   = max(34, round(52 * s));   % y tick labels and the axis label
BOTTOM = max(30, round(46 * s));   % group names under the bars
TOP    = max(20, round(30 * s));   % the panel title

panelW = (P(3) - 2*PAD - LEFT) / nW;
panelH = max(P(4) - PAD - BOTTOM - TOP, 60);

axesMade = gobjects(1, nW);
for k = 1:nW
    ax = opts.axesFcn(parent, ...
        [PAD + LEFT + (k-1)*panelW, PAD + BOTTOM, panelW - round(14 * s), panelH]);
    drawOnePanel(ax, est{k}, w(k), res, opts, k == 1);
    axesMade(k) = ax;
end
info.axes = axesMade;
end

% ── helpers ─────────────────────────────────────────────────────────────────

function v = subjectMeasures(curves, time, win, measure)
% One scalar per subject (one row of curves), measured the way the CSV measures.
n = size(curves, 1);
v = nan(n, 1);
pol = windowPolarity(win);
for i = 1:n
    m = computeWindowMeasures(curves(i, :), time, win.winStart, win.winEnd, pol);
    switch lower(char(measure))
        case 'peak', v(i) = m.peakAmp;
        case 'area', v(i) = m.area;
        otherwise,   v(i) = m.mean;
    end
end
end

function drawOnePanel(ax, e, win, res, opts, isFirst)
nG  = numel(e);
mu  = arrayfun(@(x) x.mean, e);
lo  = arrayfun(@(x) x.lo,   e);
hi  = arrayfun(@(x) x.hi,   e);

cla(ax, 'reset');
hold(ax, 'on');

BAR_WIDTH = 0.68;
b = bar(ax, 1:nG, mu, BAR_WIDTH, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = opts.colors(1:nG, :);

% Asymmetric on purpose: a paired interval is not symmetric about the mean once
% the Morey correction is applied, and drawing +/- one number would misreport it.
errorbar(ax, 1:nG, mu, mu - lo, hi - mu, 'LineStyle', 'none', ...
         'Color', [0.25 0.27 0.30], 'LineWidth', 0.9, 'CapSize', 5, ...
         'HandleVisibility', 'off');

% The measures are signed, so zero has to be visible or the bars mislead.
plot(ax, [0.4, nG + 0.6], [0 0], 'Color', [0 0 0 0.45], 'LineWidth', 0.75, ...
     'HandleVisibility', 'off');

hold(ax, 'off');
xlim(ax, [0.4, nG + 0.6]);
% Textures last, once the limits are final: they are laid out in pixels.
if ~isempty(opts.barPatterns)
    hatchBars(ax, 1:nG, mu, BAR_WIDTH, opts.barPatterns(1:nG), opts.colors(1:nG, :));
end
ax.XTick      = 1:nG;
ax.XTickLabel = {res.groups.name};
ax.XAxis.FontSize = 8;
title(ax, sprintf('%s  %g-%g', win.name, win.winStart, win.winEnd), 'FontSize', 9);
box(ax, 'off');
if isFirst
    ylabel(ax, measureLabel(opts.measure, opts.mode));
end
end

function s = measureLabel(measure, mode)
switch lower(char(measure))
    case 'peak', s = 'Peak amplitude (\muV)';
    case 'area', s = 'Area (\muV\cdotms)';
    otherwise
        if any(strcmpi(mode, {'GMFP', 'LMFP'}))
            s = sprintf('Mean %s (\\muV)', upper(char(mode)));
        else
            s = 'Mean amplitude (\muV)';
        end
end
end

function r = parentRect(parent)
p = parent.Position;
r = [0 0 p(3) p(4)];
end

function s = chromeScale(P)
% Same reference as drawTEPTopo: the pixel constants are right at roughly the
% in-app panel size, and shrink with a smaller canvas rather than eating it.
s = min([P(3) / 900, P(4) / 500, 1]);
s = max(s, 0.45);
end
