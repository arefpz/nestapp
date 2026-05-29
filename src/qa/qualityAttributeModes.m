
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function modes = qualityAttributeModes()
% QUALITYATTRIBUTEMODES  Valid values for the qualityAttribute preference.
%   Single source of truth - referenced by computeAttributeMatrix
%   (validation), runPipelineCore (validation), and the
%   attributeDisplayName lookup in renderQualityFigure.
modes = {'minmax', 'minmax_no_tms', 'highfreq'};
end
