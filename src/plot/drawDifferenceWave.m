% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function drawDifferenceWave(ax, res, opts)
% DRAWDIFFERENCEWAVE  One group minus the other, with an interval on it.
%   DRAWDIFFERENCEWAVE(ax, res, opts) draws the contrast between exactly two
%   groups from a groupCurves result.
%
%   opts:
%     .colors     n-by-3, first row is the difference colour (for scripts)
%     .groupStyle a groupStyles struct: its .contrast colour and .lineWidth.
%                 NOT a group's colour - the difference is neither group, and
%                 drawing it in group 1's colour made it read as that group
%     .xlim       time limits, default [-50 300]
%     .direction  'second - first' (default) or 'first - second'
%     .level      confidence level, default the one groupCurves used
%
%   WHICH WAY ROUND IS A CHOICE, not a consequence of the order the groups were
%   added. A difference wave inverted is the same result reported upside down,
%   and which sign a reader finds natural depends on the question - "what did
%   the intervention add" is not the same sentence as "what did it take away".
%   Left alone it stays second-minus-first, so nothing that already reads this
%   plot changes; the setting exists so the insertion order does not decide it.
%
%   Flipping is exact negation of the stored estimate, not a recomputation:
%   the bounds are mean -/+ t*sem, so negating the mean and swapping lo with hi
%   gives precisely the interval the other contrast would have produced. The
%   paired normalisation inside .sem is symmetric under negation, so there is
%   nothing left to redo - which is why this is a DRAW option and needs no trip
%   back through groupCurves.
%
%   The estimate arrives as res.contrast, computed by groupCurves at the run's
%   confidence level, exactly as drawTEPOverlay renders res.est. Deriving it
%   here instead put a statistic inside a drawing function and ignored the
%   level the caller asked for.
%
%   A shaded band that excludes zero is not a significance test and is not
%   labelled as one - the app reports estimates, not p-values. It is the
%   interval, and the reader draws their own conclusion.
%
%   See also: differenceInterval, groupCurves, drawTEPOverlay

if nargin < 3; opts = struct(); end
if numel(res.groups) ~= 2
    error('nestapp:differenceNeedsTwoGroups', ...
          'A difference wave needs exactly two groups; got %d.', numel(res.groups));
end

opts = fillDefaults(opts, struct('xlim', [-50 300], ...
                                 'showBands', false, 'windows', [], ...
                                 'direction', 'second - first', 'level', []));
sty = resolveGroupStyle(opts, {res.groups.name});
ink = sty.contrast;
if isfield(opts, 'colors') && ~isempty(opts.colors); ink = opts.colors(1, :); end

A = res.groups(1);
B = res.groups(2);
% groupCurves computed this at the run's confidence level; deriving it here
% would silently ignore that level, as this function used to.
% Re-derived before the flip, so .note carries the drawn level into the title
% and from there into the exported figure's footer.
est = intervalAtLevel(res.contrast, opts.level);

if matchesChoice(opts.direction, 'first - second')
    % Negate, and swap the bounds with it: lo = mean - t*sem becomes the upper
    % bound of the negated estimate. Swapping the group handles too keeps the
    % legend, title and picture describing the same subtraction.
    est.mean = -est.mean;
    [est.lo, est.hi] = deal(-est.hi, -est.lo);
    [A, B] = deal(B, A);
end

d  = est.mean;
lo = est.lo;
hi = est.hi;

time = res.time;
cla(ax, 'reset');
hold(ax, 'on');

fill(ax, [time, fliplr(time)], [lo, fliplr(hi)], ink, ...
     'FaceAlpha', 0.20, 'LineStyle', 'none', 'HandleVisibility', 'off');
plot(ax, time, d, 'Color', ink, 'LineWidth', sty.lineWidth, ...
     'DisplayName', sprintf('%s - %s (n=%d)', B.name, A.name, est.n));

yl = ylim(ax);
plot(ax, [0 0], yl, 'Color', [0 0 0 0.35], 'LineWidth', 0.75, 'HandleVisibility', 'off');
plot(ax, opts.xlim, [0 0], 'Color', [0 0 0 0.55], 'LineWidth', 0.9, 'HandleVisibility', 'off');
ylim(ax, yl);
hold(ax, 'off');

% Shading goes on last, so the y limits the bands span are the final ones.
% shadeTimeWindows stacks its patches to the bottom, so no curve is covered.
if opts.showBands
    shadeTimeWindows(ax, opts.windows);
end

xlim(ax, opts.xlim);
xlabel(ax, 'Time (ms)');
ylabel(ax, 'Difference (\muV)');
title(ax, sprintf('%s minus %s  (%s)', B.name, A.name, est.note));
box(ax, 'off');
legend(ax, 'Location', 'northeast', 'Box', 'off');
end

