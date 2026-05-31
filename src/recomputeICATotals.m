
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function ica = recomputeICATotals(ica)
% RECOMPUTEICATOTALS  Derive top-level ICA report fields from the per-round records.
%   ica = RECOMPUTEICATOTALS(ica)
%
%   With ICA recorded one round per decomposition (see openICARound /
%   addICARemoval), the top-level summary fields are derived from the rounds:
%     nComponents - the first (original) decomposition's component count
%     nRejected   - total components removed across all rounds
%     nKept       - components surviving the final round (nComp - nRej of the
%                   last round); equals nComponents - nRejected because each
%                   round re-decomposes the previous round's residual
%     categories  - union of every round's per-category tally
%     varRemoved / varMin / varMax - from the first round only (variance is not
%                   additive across different ICA bases)
%
%   See also: openICARound, addICARemoval, mergeCategories, buildReportText

if isempty(ica.rounds)
    return
end

ica.nComponents = ica.rounds{1}.nComponents;

total = 0;
cats  = struct('names', {{}}, 'nRemoved', [], 'varShare', []);
for i = 1:numel(ica.rounds)
    total = total + ica.rounds{i}.nRejected;
    cats  = mergeCategories(cats, ica.rounds{i}.categories);
end
ica.nRejected = total;

last      = ica.rounds{end};
ica.nKept = last.nComponents - last.nRejected;

ica.categories = cats;
ica.varRemoved = ica.rounds{1}.varRemoved;
ica.varMin     = ica.rounds{1}.varMin;
ica.varMax     = ica.rounds{1}.varMax;
end
