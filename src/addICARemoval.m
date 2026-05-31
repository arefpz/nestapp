
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function fileReport = addICARemoval(fileReport, removal)
% ADDICAREMOVAL  Fold a component removal into the current ICA round.
%   fileReport = ADDICAREMOVAL(fileReport, removal)
%
%   Every ICA removal path (ICLabel / AARATEP / ARTIST / manual via "Remove
%   Flagged ICA Components", and TESA compselect) calls this so all of them
%   report with the same per-round + per-category shape. The removal is
%   attributed to the most recent round opened by openICARound.
%
%   removal fields:
%     nRejected    - components removed
%     varRemoved   - % data variance removed (NaN if unavailable)
%     varMin/varMax- per-component variance range (NaN if unavailable)
%     categories   - struct(names, nRemoved, varShare)
%     nComponents  - (optional) decomposition size, used only to size a
%                    defensively-opened round when no Run ICA preceded it
%
%   See also: openICARound, recomputeICATotals, mergeCategories

if isempty(fileReport.ica.rounds)
    % Defensive: a removal with no preceding Run ICA round. Open one so totals
    % stay sane, sized to the removal's decomposition if it reported one.
    nComp = removal.nRejected;
    if isfield(removal, 'nComponents'); nComp = removal.nComponents; end
    fileReport = openICARound(fileReport, nComp);
end

k   = numel(fileReport.ica.rounds);
rnd = fileReport.ica.rounds{k};

rnd.nRejected = rnd.nRejected + removal.nRejected;
if isnan(rnd.varRemoved)   % first removal this round carries the variance figures
    rnd.varRemoved = removal.varRemoved;
    rnd.varMin     = removal.varMin;
    rnd.varMax     = removal.varMax;
end
rnd.categories = mergeCategories(rnd.categories, removal.categories);

fileReport.ica.rounds{k} = rnd;
fileReport.ica = recomputeICATotals(fileReport.ica);
end
