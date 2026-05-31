
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function cats = mergeCategories(cats, add)
% MERGECATEGORIES  Union-merge one ICA category tally into another by name.
%   cats = MERGECATEGORIES(cats, add)
%
%   Both are structs with fields names{1xC}, nRemoved(1xC), varShare(1xC).
%   Categories from `add` are added to `cats` by matching name; new names are
%   appended. Empty categories (zero count and zero variance) are skipped so
%   unused schemes (e.g. the default ICLabel 7 when a TESA round is recorded)
%   don't add zero rows.
%
%   See also: addICARemoval, recomputeICATotals, icaCategoriesFromFlags

for i = 1:numel(add.names)
    if add.nRemoved(i) == 0 && add.varShare(i) == 0
        continue
    end
    j = find(strcmp(cats.names, add.names{i}), 1);
    if isempty(j)
        cats.names{end+1}    = add.names{i};   %#ok<AGROW>
        cats.nRemoved(end+1) = add.nRemoved(i); %#ok<AGROW>
        cats.varShare(end+1) = add.varShare(i); %#ok<AGROW>
    else
        cats.nRemoved(j) = cats.nRemoved(j) + add.nRemoved(i);
        cats.varShare(j) = cats.varShare(j) + add.varShare(i);
    end
end
end
