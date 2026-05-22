function modes = qualityAttributeModes()
% QUALITYATTRIBUTEMODES  Valid values for the qualityAttribute preference.
%   Single source of truth - referenced by computeAttributeMatrix
%   (validation), runPipelineCore (validation), and the
%   attributeDisplayName lookup in renderQualityFigure.
modes = {'minmax', 'minmax_no_tms', 'highfreq'};
end
