% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function Y = relativeLuminance(rgb)
% RELATIVELUMINANCE  WCAG relative luminance of each sRGB colour (rows of rgb).
%   Y = RELATIVELUMINANCE(rgb) linearises gamma-encoded sRGB and weights it,
%   giving 0 for black and 1 for white. Contrast against white is then
%   1.05 / (Y + 0.05).
%
%   One definition for the hatch ink choice (hatchBars) and the palette audit
%   (tools/auditPalettes), so what the audit measures is what the app uses.
%
%   See also: hatchBars, auditPalettes

lin = ((rgb + 0.055) / 1.055) .^ 2.4;
small = rgb <= 0.04045;
lin(small) = rgb(small) / 12.92;
Y = lin * [0.2126; 0.7152; 0.0722];
end
