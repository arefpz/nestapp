% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function themes = plotThemes()
% PLOTTHEMES  The colour themes Explore plots can be drawn in. The first is the default.
%   themes = PLOTTHEMES() returns a struct array, one element per theme:
%     .name         shown in the Colours dialog and stored in preferences
%     .colors       n x 3 RGB, in the order groups receive them
%     .lineStyles   line styles; colour themes use them only once the colours
%                   run out, the print theme from the first group
%     .barPatterns  bar textures (see hatchBars), cycled the same way
%     .styleFirst   true when styles and textures vary from the FIRST group
%                   rather than only after the colours run out
%     .lineWidth    width of the group mean lines
%
%   The palettes were chosen by eye from tools/themeGallery and checked with
%   tools/auditPalettes (contrast on white; separation under simulated colour
%   vision deficiencies). Each palette is cut where that separation falls off,
%   so a group past the cut is told apart by line style and texture rather
%   than by a colour too close to an earlier one:
%     Okabe-Ito  7 - the amber that replaces yellow sits too near the orange
%                    for a colour-blind reader as an 8th colour
%     Vibrant    6 - Tol's grey is left out; he reserves it for missing data
%     Muted      9
%     Soft, Jewel, Earth, Slate  8 each - nestapp palettes, each a mood (a box
%                    of OKLCH chroma and hue) searched for colours that stay
%                    apart under colour vision deficiencies; varied lightness
%                    is what keeps them safe
%     Greyscale  greys alone cannot separate groups, so every group also gets
%                its own line style and bar texture, and lines are heavier so
%                dotted ones survive printing
%
%   Sources carry no licence terms: Okabe & Ito (2008) and Paul Tol's SRON
%   schemes are published for free use; Soft, Jewel, Earth and Slate are
%   nestapp's own.
%
%   See also: groupStyles, groupColors, hatchBars

COLOUR_STYLES = {'-', '--', ':', '-.'};
COLOUR_PATTERNS = {'solid', 'diagonal', 'backdiagonal', 'crosshatch', 'horizontal', 'dots'};

themes = [ ...
    theme('Okabe-Ito', hex('#262626', '#E69F00', '#0072B2', '#009E73', '#D55E00', ...
                           '#56B4E9', '#CC79A7'), ...
          COLOUR_STYLES, COLOUR_PATTERNS, false, 1.75)
    theme('Vibrant', hex('#0077BB', '#EE7733', '#009988', '#EE3377', '#33BBEE', '#CC3311'), ...
          COLOUR_STYLES, COLOUR_PATTERNS, false, 1.75)
    theme('Muted', hex('#332288', '#CC6677', '#117733', '#DDCC77', '#88CCEE', '#882255', ...
                       '#44AA99', '#999933', '#AA4499'), ...
          COLOUR_STYLES, COLOUR_PATTERNS, false, 1.75)
    theme('Soft', hex('#3B638C', '#AA8544', '#6C352D', '#499B9F', '#81632F', '#65836B', ...
                      '#573E55', '#8279B8'), ...
          COLOUR_STYLES, COLOUR_PATTERNS, false, 1.75)
    theme('Jewel', hex('#314398', '#BB6802', '#760A03', '#6079D3', '#99335A', '#681053', ...
                       '#319751', '#A05A9A'), ...
          COLOUR_STYLES, COLOUR_PATTERNS, false, 1.75)
    theme('Earth', hex('#884C28', '#4F91BD', '#014555', '#AC7980', '#8D8E46', '#376C8F', ...
                       '#543232', '#5A7854'), ...
          COLOUR_STYLES, COLOUR_PATTERNS, false, 1.75)
    theme('Slate', hex('#1C3249', '#BC7D25', '#8186C5', '#58290A', '#087158', '#605692', ...
                       '#3C9F8E', '#825618'), ...
          COLOUR_STYLES, COLOUR_PATTERNS, false, 1.75)
    theme('Greyscale', repmat([0.10; 0.45; 0.65; 0.30; 0.55], 1, 3), ...
          COLOUR_STYLES, COLOUR_PATTERNS, true, 2.25)];
end

function t = theme(name, colors, lineStyles, barPatterns, styleFirst, lineWidth)
t = struct('name', name, 'colors', colors, 'lineStyles', {lineStyles}, ...
           'barPatterns', {barPatterns}, 'styleFirst', styleFirst, 'lineWidth', lineWidth);
end

function rgb = hex(varargin)
rgb = zeros(numel(varargin), 3);
for k = 1:numel(varargin)
    h = varargin{k}(2:end);
    rgb(k, :) = [hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))] / 255;
end
end
