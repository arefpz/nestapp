% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function files = themeGallery(outDir)
% THEMEGALLERY  Render the same figures in every candidate palette, to choose by eye.
%   files = THEMEGALLERY(outDir) writes one PNG per group palette (a TEP
%   overlay with 3 groups, one with 6, and the window bars) plus one page of
%   topographies per diverging map, and returns their paths. Needs EEGLAB for
%   the maps.
%
%   The figures are drawn by nestapp's own drawers on a synthetic cohort with
%   TEP-shaped curves, so what is judged is what a user would see. Tiles carry
%   only the palette's name: accessibility numbers live in auditPalettes and are
%   for developers, not for this page.
%
%   See also: candidatePalettes, auditPalettes

if nargin < 1 || isempty(outDir); outDir = fullfile(tempdir, 'nestapp_theme_gallery'); end
if ~isfolder(outDir); mkdir(outDir); end
[groupSets, divergingSets] = candidatePalettes();

res3 = syntheticResult(3);
res6 = syntheticResult(6);
files = {};

for s = 1:numel(groupSets)
    C = groupSets(s).colors;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1100 430]);
    t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, groupSets(s).name, 'FontWeight', 'bold', 'FontSize', 14);
    drawTEPOverlay(nexttile(t), res3, struct('colors', cycle(C, 3)));
    drawTEPOverlay(nexttile(t), res6, struct('colors', cycle(C, 6)));
    curves = print(fig, '-RGBImage', '-r0');
    delete(fig);

    % The bars plot lays out its own axes across a whole parent.
    fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
                 'Position', [50 50 1000 430]);
    drawWindowBars(fig, res6, struct('colors', cycle(C, 6), ...
        'windows', defaultTEPComponentDefs(), ...
        'axesFcn', @(p, pos) axes('Parent', p, 'Units', 'pixels', 'Position', pos)));
    bars = print(fig, '-RGBImage', '-r0');
    delete(fig);

    files{end+1} = fullfile(outDir, sprintf('group_%02d_%s.png', s, slug(groupSets(s).name))); %#ok<AGROW>
    imwrite([curves, imresize(bars, [size(curves, 1) NaN])], files{end});
end

[vals, locs] = syntheticMaps();
nD = numel(divergingSets); nK = size(vals, 2);
LABEL_W = 170; CELL = 210;
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
             'Position', [50 50 LABEL_W + CELL * nK, CELL * nD]);
for d = 1:nD
    y = CELL * (nD - d);
    annotation(fig, 'textbox', 'Units', 'pixels', 'Position', [8 y LABEL_W - 16 CELL], ...
        'String', divergingSets(d).name, 'EdgeColor', 'none', 'FontWeight', 'bold', ...
        'FontSize', 12, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    for k = 1:nK
        ax = axes('Parent', fig, 'Units', 'pixels', ...
                  'Position', [LABEL_W + CELL * (k - 1) + 10, y + 10, CELL - 20, CELL - 20]);
        drawScalpTopo(ax, vals(:, k), locs, struct('colorbar', false, 'clim', [-1 1]));
        colormap(ax, divergingSets(d).map);
    end
end
files{end+1} = fullfile(outDir, 'topo_maps.png');
exportgraphics(fig, files{end}, 'Resolution', 110);
delete(fig);
end

% -- synthetic data --------------------------------------------------------------

function C = cycle(C, n)
C = C(mod(0:n - 1, size(C, 1)) + 1, :);
end

function s = slug(name)
s = lower(regexprep(name, '[^A-Za-z0-9]+', '_'));
end

function res = syntheticResult(nGroups)
% A cohort whose groups have recognisably TEP-shaped curves that differ a
% little in amplitude and latency, so overlapping lines and bands look real.
[cache, entries] = fakeReducedCache('subjects', 8, 'sessions', nGroups, 'repeats', 0);
t = cache(1).time;
tep = @(t, a, lag) a * ( 0.6 * g(t, 30 + lag, 6) - 1.0 * g(t, 45 + lag, 8) ...
        + 1.2 * g(t, 60 + lag, 10) - 2.6 * g(t, 110 + lag, 22) + 2.0 * g(t, 190 + lag, 35));
rng(3, 'twister');   % explicit generator: EEGLAB's runica can leave the legacy one active
for i = 1:numel(cache)
    grp  = sscanf(entries(i).group, 'g%d');
    subj = sscanf(entries(i).subject, 's%d');
    amp  = (1 + 0.12 * (grp - 1)) * (0.8 + 0.05 * subj);
    lag  = 3 * (grp - 1) + randn;
    base = tep(t, amp, lag) + 0.15 * randn(size(t));
    cache(i).trialAvg = repmat(base, numel(cache(i).labels), 1) .* ...
                        linspace(1, 0.6, numel(cache(i).labels))';
end
res = groupCurves(cache, entries, struct('roi', {cache(1).labels(1:2)}));
end

function y = g(t, mu, sd)
y = exp(-0.5 * ((t - mu) / sd) .^ 2);
end

function [vals, locs] = syntheticMaps()
% Four plausible scalp patterns on a standard 32-channel layout.
labels = {'Fp1','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','T7','C3', ...
          'Cz','C4','T8','CP5','CP1','CP2','CP6','P7','P3','Pz','P4','P8','PO3', ...
          'PO4','O1','Oz','O2','AF3','AF4'};
locs = struct('labels', labels);
locs = pop_chanedit(locs, 'lookup', fullfile(fileparts(which('eeglab')), 'plugins', ...
                    'dipfit', 'standard_BEM', 'elec', 'standard_1005.elc'));
xyz = [[locs.X]' [locs.Y]' [locs.Z]'];
xyz = xyz ./ vecnorm(xyz, 2, 2);
src = @(p) exp(-sum((xyz - p ./ norm(p)) .^ 2, 2) / 0.35);
vals = [src([0.5 0.5 0.7]) - src([-0.6 -0.4 0.6]), ...         % left frontal dipole
        src([0 0 1]) - 0.7 * src([0.9 0 0.2]) - 0.7 * src([-0.9 0 0.2]), ...
        -src([0.2 0 1]) + 0.8 * src([-0.8 0 0.4]), ...         % N100-like
        src([-0.3 0.4 0.8]) - src([0.3 -0.6 0.7])];
vals = vals ./ max(abs(vals));
end
