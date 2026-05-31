
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function vars = stripEmptyVarin(vars)
% STRIPEMPTYVARIN  Remove key-value pairs where the value is an "empty" sentinel.
%   vars = STRIPEMPTYVARIN(vars)
%
%   Handles the old '[]' string sentinel, numeric [], NaN scalar, and empty cell -
%   all representations of "not set" used across old and new typed-model param formats.
isEmptyVal = cellfun(@(v) (ischar(v) && strcmp(v,'[]')) || ...
                          (isnumeric(v) && isempty(v)) || ...
                          (isnumeric(v) && isscalar(v) && isnan(v)) || ...
                          (iscell(v) && isempty(v)) || ...
                          (iscell(v) && isscalar(v) && ischar(v{1}) && strcmp(v{1},'[]')), ...
                          vars(2:2:end));
toRemove  = find(isEmptyVal);
removeIdx = sort([2*toRemove-1, 2*toRemove]);
vars(removeIdx) = [];
end
