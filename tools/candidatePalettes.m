% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [groupSets, divergingSets] = candidatePalettes()
% CANDIDATEPALETTES  Every palette under consideration for the plot themes.
%   [groupSets, divergingSets] = CANDIDATEPALETTES()
%
%   groupSets      struct array .name .family .colors (n x 3 RGB in [0,1])
%   divergingSets  struct array .name .family .map (256 x 3)
%
%   A developer tool for choosing themes (tools/themeGallery, tools/auditPalettes),
%   not shipped behaviour. Limited to sources with no licence terms to honour:
%     - Okabe & Ito (2008) and Paul Tol's SRON schemes, both published for free use
%     - MATLAB's own named palettes, read at run time with orderedcolors (R2023b+),
%       never copied
%     - nestapp originals, designed here in OKLCH
%
%   See also: themeGallery, auditPalettes

groupSets = struct('name', {}, 'family', {}, 'colors', {});

% -- accessible ---------------------------------------------------------------
groupSets(end+1) = entry('Okabe-Ito (current)', 'accessible', groupColors(8));
groupSets(end+1) = entry('Okabe-Ito, amber for yellow', 'accessible', hex( ...
    '#262626', '#E69F00', '#0072B2', '#009E73', '#D55E00', '#56B4E9', '#CC79A7', '#B8860B'));
groupSets(end+1) = entry('Tol bright', 'accessible', hex( ...
    '#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'));
groupSets(end+1) = entry('Tol vibrant', 'accessible', hex( ...
    '#0077BB', '#EE7733', '#009988', '#EE3377', '#33BBEE', '#CC3311', '#BBBBBB'));
groupSets(end+1) = entry('Tol muted', 'accessible', hex( ...
    '#332288', '#CC6677', '#117733', '#DDCC77', '#88CCEE', '#882255', '#44AA99', ...
    '#999933', '#AA4499'));
groupSets(end+1) = entry('Tol high contrast', 'accessible', hex( ...
    '#004488', '#DDAA33', '#BB5566'));

% -- MATLAB's own, read at run time -------------------------------------------
if exist('orderedcolors', 'file')
    for nm = {'gem', 'gem12', 'glow', 'sail', 'reef', 'meadow', 'dye', 'earth'}
        try
            groupSets(end+1) = entry(['MATLAB ' nm{1}], 'matlab', orderedcolors(nm{1})); %#ok<AGROW>
        catch
            % A name missing from this release is simply not offered.
        end
    end
end

% -- nestapp originals --------------------------------------------------------
% Each has a mood - a box of OKLCH lightness, chroma and hue - and within it
% the colours were picked greedily (scratch search, 2026-11-08): each next one
% maximises its worst CIEDE2000 distance to those already chosen under normal
% and simulated protan/deutan/tritan vision, at >= 3:1 contrast on white.
% Constant-lightness drafts looked even but collapsed for colour-blind
% readers; varying lightness is what keeps them apart.
groupSets(end+1) = entry('Soft', 'original', hex( ...         % C 0.05-0.10
    '#3B638C', '#AA8544', '#6C352D', '#499B9F', '#81632F', '#65836B', '#573E55', '#8279B8'));
groupSets(end+1) = entry('Jewel', 'original', hex( ...        % C 0.11-0.17
    '#314398', '#BB6802', '#760A03', '#6079D3', '#99335A', '#681053', '#319751', '#A05A9A'));
groupSets(end+1) = entry('Earth', 'original', hex( ...        % warm hues + slate blue
    '#884C28', '#4F91BD', '#014555', '#AC7980', '#8D8E46', '#376C8F', '#543232', '#5A7854'));
groupSets(end+1) = entry('Slate', 'original', hex( ...        % cool, one warm accent
    '#1C3249', '#BC7D25', '#8186C5', '#58290A', '#087158', '#605692', '#3C9F8E', '#825618'));

% -- print -------------------------------------------------------------------
groupSets(end+1) = entry('Greyscale', 'print', repmat([0.10; 0.45; 0.65; 0.30; 0.55], 1, 3));

% -- diverging maps for topoplots ----------------------------------------------
divergingSets = struct('name', {}, 'family', {}, 'map', {});
divergingSets(end+1) = struct('name', 'Red-blue (current)', 'family', 'current', ...
                              'map', divergingColormap());
divergingSets(end+1) = struct('name', 'Soft red-blue', 'family', 'original', ...
    'map', divergingOklch([0.48 0.11 255], [0.50 0.13 25]));
divergingSets(end+1) = struct('name', 'Teal-orange', 'family', 'original', ...
    'map', divergingOklch([0.45 0.09 200], [0.60 0.14 55]));
divergingSets(end+1) = struct('name', 'Purple-green', 'family', 'original', ...
    'map', divergingOklch([0.42 0.13 300], [0.50 0.12 145]));
end

% -- helpers -------------------------------------------------------------------

function e = entry(name, family, colors)
e = struct('name', name, 'family', family, 'colors', colors);
end

function rgb = hex(varargin)
rgb = zeros(numel(varargin), 3);
for k = 1:numel(varargin)
    h = varargin{k}(2:end);
    rgb(k, :) = [hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))] / 255;
end
end

function map = divergingOklch(cold, warm)
% Cold end -> white -> warm end, interpolated in OKLab so lightness rises and
% falls evenly on both limbs.
n = 256;
t = linspace(-1, 1, n)';
toLab = @(p) [p(1), p(2) * cosd(p(3)), p(2) * sind(p(3))];
c = toLab(cold); w = toLab(warm); white = [1 0 0];
lab = zeros(n, 3);
for i = 1:n
    s = abs(t(i));
    if t(i) < 0, e = c; else, e = w; end
    lab(i, :) = white * (1 - s) + e * s;
end
map = oklabToSrgb(lab);
end

function rgb = oklabToSrgb(lab)
% OKLab (Ottosson 2020) to gamma-encoded sRGB, clipped to gamut.
M2 = [1  0.3963377774  0.2158037573
      1 -0.1055613458 -0.0638541728
      1 -0.0894841775 -1.2914855480];
M1 = [ 4.0767416621 -3.3077115913  0.2309699292
      -1.2684380046  2.6097574011 -0.3413193965
      -0.0041960863 -0.7034186147  1.7076147010];
lms = (lab * M2') .^ 3;
lin = lms * M1';
lin = min(max(lin, 0), 1);
rgb = 1.055 * lin .^ (1 / 2.4) - 0.055;
small = lin <= 0.0031308;
rgb(small) = 12.92 * lin(small);
rgb = min(max(rgb, 0), 1);
end
