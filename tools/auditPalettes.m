% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function T = auditPalettes(prefixes)
% AUDITPALETTES  Developer report: how separable and legible each palette is.
%   T = AUDITPALETTES() audits every group palette in candidatePalettes and
%   returns (and prints) a table. T = AUDITPALETTES(prefixes) sets the group
%   counts to evaluate, default [3 4 6 8].
%
%   For each palette, and each prefix of its first n colours (groups are drawn
%   in palette order, so the first few colours matter most):
%     minContrast  lowest WCAG contrast of any colour on white. Lines need
%                  >= 3:1 to read; below that a colour fades into the page
%     dE_<n>       smallest CIEDE2000 distance between any two of the first n
%                  colours, under normal vision (N) and simulated protan (P),
%                  deutan (D) and tritan (T) vision (Machado et al. 2009, full
%                  severity). Roughly: < 10 easily confused, 10-20 distinct
%                  side by side, > 20 distinct at a glance
%     greyGap      smallest lightness step between any two of the first 4
%                  colours once printed in greyscale (CIELAB L*)
%
%   Internal only: these numbers inform which themes ship and in what order.
%   They are never shown in the app.
%
%   See also: candidatePalettes, themeGallery

if nargin < 1; prefixes = [3 4 6 8]; end
sets = candidatePalettes();

rows = cell(numel(sets), 1);
for s = 1:numel(sets)
    C = sets(s).colors;
    r = struct('palette', string(sets(s).name), 'family', string(sets(s).family), ...
               'n', size(C, 1), 'minContrast', min(contrastOnWhite(C)));
    for n = prefixes
        if n > size(C, 1)
            r.(sprintf('dE%d_N', n)) = NaN; r.(sprintf('dE%d_worstCVD', n)) = NaN;
            continue
        end
        P = C(1:n, :);
        r.(sprintf('dE%d_N', n)) = minPairDE(P);
        r.(sprintf('dE%d_worstCVD', n)) = min([minPairDE(simulateCvd(P, 'protan')), ...
                                               minPairDE(simulateCvd(P, 'deutan')), ...
                                               minPairDE(simulateCvd(P, 'tritan'))]);
    end
    lab = srgbToLab(C(1:min(4, end), :));
    r.greyGap = min(diff(sort(lab(:, 1))));
    rows{s} = r;
end
T = struct2table([rows{:}]);
T = varfun(@(x) roundIfNumeric(x), T, 'OutputFormat', 'table');
T.Properties.VariableNames = regexprep(T.Properties.VariableNames, '^Fun_', '');
if nargout == 0; disp(T); end
end

% -- metrics -------------------------------------------------------------------

function x = roundIfNumeric(x)
if isnumeric(x); x = round(x, 1); end
end

function c = contrastOnWhite(rgb)
% WCAG 2 contrast ratio of each colour against white.
c = 1.05 ./ (relLuminance(rgb) + 0.05);
end

function Y = relLuminance(rgb)
Y = srgbToLinear(rgb) * [0.2126; 0.7152; 0.0722];
end

function lin = srgbToLinear(rgb)
lin = ((rgb + 0.055) / 1.055) .^ 2.4;
small = rgb <= 0.04045;
lin(small) = rgb(small) / 12.92;
end

function rgb = linearToSrgb(lin)
lin = min(max(lin, 0), 1);
rgb = 1.055 * lin .^ (1 / 2.4) - 0.055;
small = lin <= 0.0031308;
rgb(small) = 12.92 * lin(small);
end

function out = simulateCvd(rgb, kind)
% Machado, Oliveira & Fernandes (2009), severity 1.0, applied in linear RGB.
switch kind
    case 'protan'
        M = [0.152286 1.052583 -0.204868; 0.114503 0.786281 0.099216; -0.003882 -0.048116 1.051998];
    case 'deutan'
        M = [0.367322 0.860646 -0.227968; 0.280085 0.672501 0.047413; -0.011820 0.042940 0.968881];
    case 'tritan'
        M = [1.255528 -0.076749 -0.178779; -0.078411 0.930809 0.147602; 0.004733 0.691367 0.303900];
end
out = linearToSrgb(srgbToLinear(rgb) * M');
end

function lab = srgbToLab(rgb)
% sRGB -> CIELAB (D65).
xyz = srgbToLinear(rgb) * [0.4124564 0.3575761 0.1804375
                           0.2126729 0.7151522 0.0721750
                           0.0193339 0.1191920 0.9503041]';
xyz = xyz ./ [0.95047 1.00000 1.08883];
f = xyz .^ (1/3);
small = xyz <= (6/29)^3;
f(small) = xyz(small) / (3 * (6/29)^2) + 4/29;
lab = [116 * f(:, 2) - 16, 500 * (f(:, 1) - f(:, 2)), 200 * (f(:, 2) - f(:, 3))];
end

function d = minPairDE(rgb)
lab = srgbToLab(rgb);
d = Inf;
for i = 1:size(lab, 1) - 1
    for j = i + 1:size(lab, 1)
        d = min(d, ciede2000(lab(i, :), lab(j, :)));
    end
end
end

function dE = ciede2000(lab1, lab2)
% CIEDE2000 colour difference (Sharma, Wu & Dalal 2005), kL = kC = kH = 1.
L1 = lab1(1); a1 = lab1(2); b1 = lab1(3);
L2 = lab2(1); a2 = lab2(2); b2 = lab2(3);
Cbar = (hypot(a1, b1) + hypot(a2, b2)) / 2;
G = 0.5 * (1 - sqrt(Cbar^7 / (Cbar^7 + 25^7)));
a1p = (1 + G) * a1; a2p = (1 + G) * a2;
C1p = hypot(a1p, b1); C2p = hypot(a2p, b2);
h1p = mod(atan2d(b1, a1p), 360); h2p = mod(atan2d(b2, a2p), 360);
dLp = L2 - L1; dCp = C2p - C1p;
if C1p * C2p == 0
    dhp = 0;
else
    dhp = h2p - h1p;
    if dhp > 180, dhp = dhp - 360; elseif dhp < -180, dhp = dhp + 360; end
end
dHp = 2 * sqrt(C1p * C2p) * sind(dhp / 2);
Lbp = (L1 + L2) / 2; Cbp = (C1p + C2p) / 2;
if C1p * C2p == 0
    hbp = h1p + h2p;
elseif abs(h1p - h2p) <= 180
    hbp = (h1p + h2p) / 2;
elseif h1p + h2p < 360
    hbp = (h1p + h2p + 360) / 2;
else
    hbp = (h1p + h2p - 360) / 2;
end
T  = 1 - 0.17 * cosd(hbp - 30) + 0.24 * cosd(2 * hbp) + 0.32 * cosd(3 * hbp + 6) ...
     - 0.20 * cosd(4 * hbp - 63);
dTheta = 30 * exp(-((hbp - 275) / 25)^2);
RC = 2 * sqrt(Cbp^7 / (Cbp^7 + 25^7));
SL = 1 + 0.015 * (Lbp - 50)^2 / sqrt(20 + (Lbp - 50)^2);
SC = 1 + 0.045 * Cbp;
SH = 1 + 0.015 * Cbp * T;
RT = -sind(2 * dTheta) * RC;
dE = sqrt((dLp / SL)^2 + (dCp / SC)^2 + (dHp / SH)^2 + RT * (dCp / SC) * (dHp / SH));
end
