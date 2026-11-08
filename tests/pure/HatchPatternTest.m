% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef HatchPatternTest < NestappTestCase
% HATCHPATTERNTEST  The geometry of bar textures, without drawing anything.
%
%   MATLAB has no hatch property, so hatchSegments computes the lines; a
%   texture that leaks out of its bar or runs at the wrong angle would be
%   visible on every printed figure.

    methods (Test)


        function solidHasNoTexture(tc)
            [l, d] = hatchSegments([0 40 0 80], 'solid', 6);
            tc.verifyEmpty(l); tc.verifyEmpty(d);
        end

        function everyHatchSegmentStaysInsideItsBar(tc, pattern)
            rect = [10 50 -30 70];
            [l, d] = hatchSegments(rect, pattern, 6);
            pts = [l(:, ~isnan(l(1, :))), d];
            tc.verifyNotEmpty(pts);
            tol = 1e-9;
            tc.verifyTrue(all(pts(1, :) >= 10 - tol & pts(1, :) <= 50 + tol & ...
                              pts(2, :) >= -30 - tol & pts(2, :) <= 70 + tol));
        end

        function diagonalLinesRunAt45Degrees(tc)
            l = hatchSegments([0 60 0 60], 'diagonal', 6);
            seg = l(:, 1:2);
            tc.verifyEqual(atan2d(diff(seg(2, :)), diff(seg(1, :))), 45, 'AbsTol', 1e-9);
        end

        function hatchLinesAreEvenlySpaced(tc)
            l = hatchSegments([0 60 0 30], 'horizontal', 6);
            ys = unique(round(l(2, ~isnan(l(2, :))), 9));
            tc.verifyEqual(diff(ys(:)), repmat(6, numel(ys) - 1, 1), 'AbsTol', 1e-9);
        end

        function crosshatchIsBothDiagonals(tc)
            a = hatchSegments([0 60 0 60], 'diagonal', 6);
            b = hatchSegments([0 60 0 60], 'backdiagonal', 6);
            c = hatchSegments([0 60 0 60], 'crosshatch', 6);
            tc.verifyEqual(size(c, 2), size(a, 2) + size(b, 2));
        end

        function aZeroHeightBarGetsNothing(tc)
            [l, d] = hatchSegments([0 40 5 5], 'crosshatch', 6);
            tc.verifyEmpty(l); tc.verifyEmpty(d);
        end

        function anUnknownPatternIsRefused(tc)
            tc.verifyError(@() hatchSegments([0 1 0 1], 'tartan', 6), 'nestapp:badPattern');
        end
    end

    properties (TestParameter)
        pattern = {'diagonal', 'backdiagonal', 'crosshatch', 'horizontal', 'vertical', 'dots'}
    end
end
