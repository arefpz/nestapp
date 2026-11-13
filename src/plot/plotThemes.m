% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function themes = plotThemes()
% PLOTTHEMES  The colour themes Explore plots can be drawn in. The first is the default.
%   themes = PLOTTHEMES() returns a struct array, one element per theme:
%     .name         shown in plot Options and stored in preferences
%     .colors       n x 3 RGB, in the order groups receive them
%     .lineStyles   line styles; colour themes use them only once the colours
%                   run out, the print theme from the first group
%     .barPatterns  bar textures (see hatchBars), cycled the same way
%     .styleFirst   true when styles and textures vary from the FIRST group
%     .lineWidth    width of the group mean lines
%     .contrast     colour for lines that are no group's - a difference wave
%
%   Each palette is cut where colour-blind separation falls off (see
%   tools/auditPalettes; chosen by eye in tools/themeGallery). Okabe-Ito drops
%   the amber that would collide with its orange; Vibrant leaves out Tol's
%   missing-data grey. Greyscale varies line style and texture from the first
%   group and draws heavier lines, so it survives printing. Sources carry no
%   licence terms: Okabe & Ito (2008), Paul Tol's SRON schemes, and nestapp's
%   own Soft, Jewel, Earth and Slate.
%
%   See also: groupStyles, groupColors, hatchBars

themes = [ ...
    theme('Okabe-Ito', hexToRgb('#262626', '#E69F00', '#0072B2', '#009E73', ...
                                '#D55E00', '#56B4E9', '#CC79A7'))
    theme('Vibrant',   hexToRgb('#0077BB', '#EE7733', '#009988', '#EE3377', ...
                                '#33BBEE', '#CC3311'))
    theme('Muted',     hexToRgb('#332288', '#CC6677', '#117733', '#DDCC77', ...
                                '#88CCEE', '#882255', '#44AA99', '#999933', '#AA4499'))
    theme('Soft',      hexToRgb('#3B638C', '#AA8544', '#6C352D', '#499B9F', ...
                                '#81632F', '#65836B', '#573E55', '#8279B8'))
    theme('Jewel',     hexToRgb('#314398', '#BB6802', '#760A03', '#6079D3', ...
                                '#99335A', '#681053', '#319751', '#A05A9A'))
    theme('Earth',     hexToRgb('#884C28', '#4F91BD', '#014555', '#AC7980', ...
                                '#8D8E46', '#376C8F', '#543232', '#5A7854'))
    theme('Slate',     hexToRgb('#1C3249', '#BC7D25', '#8186C5', '#58290A', ...
                                '#087158', '#605692', '#3C9F8E', '#825618'))
    theme('Greyscale', repmat([0.10; 0.45; 0.65; 0.30; 0.55], 1, 3), true, 2.25)];
end

function t = theme(name, colors, styleFirst, lineWidth)
if nargin < 3; styleFirst = false; end
if nargin < 4; lineWidth  = 1.75;  end
t = struct('name', name, 'colors', colors, ...
           'lineStyles', {{'-', '--', ':', '-.'}}, ...
           'barPatterns', {{'solid', 'diagonal', 'backdiagonal', 'crosshatch', ...
                            'horizontal', 'dots'}}, ...
           'styleFirst', styleFirst, 'lineWidth', lineWidth, ...
           'contrast', [0.15 0.15 0.15]);
end
