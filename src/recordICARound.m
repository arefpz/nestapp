
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function fileReport = recordICARound(fileReport, rnd)
% RECORDICAROUND  Append one ICA component-removal round to a report.
%   fileReport = RECORDICAROUND(fileReport, rnd)
%
%   Single entry point that every ICA removal path (TESA compselect, ICLabel,
%   AARATEP muscle, ARTIST decay, manual) feeds, so all of them report with the
%   same shape: a detected count, N rounds, and a per-category breakdown.
%
%   rnd fields:
%     roundNum     - 1-based index of this removal round
%     nComponents  - components in the decomposition this round acted on
%     nRejected    - components removed this round
%     varRemoved   - % data variance removed this round (NaN if unavailable)
%     varMin/varMax- per-component variance range (NaN if unavailable)
%     categories   - struct(names{1xC}, nRemoved(1xC), varShare(1xC))
%
%   Effects on fileReport.ica:
%     - appends rnd to .rounds
%     - accumulates .nRejected, recomputes .nKept (against .nComponents, which
%       is set once at Run ICA)
%     - merges the round's categories into .categories by name (union); empty
%       categories are skipped so unused schemes don't clutter the report
%     - sets top-level .varRemoved/.varMin/.varMax from the first round only
%       (variance across different ICA bases is not additive)
%
%   See also: processOneFile, initPipelineReport, buildReportText

fileReport.ica.rounds{end+1} = rnd;
fileReport.ica.nRejected = fileReport.ica.nRejected + rnd.nRejected;
fileReport.ica.nKept     = fileReport.ica.nComponents - fileReport.ica.nRejected;

if isscalar(fileReport.ica.rounds)
    fileReport.ica.varRemoved = rnd.varRemoved;
    fileReport.ica.varMin     = rnd.varMin;
    fileReport.ica.varMax     = rnd.varMax;
end

fileReport.ica.categories = mergeCategories(fileReport.ica.categories, rnd.categories);
end

function cats = mergeCategories(cats, add)
% Union-merge add into cats by category name. Skip empty categories so a
% scheme that flagged nothing (e.g. the default ICLabel 7 when a TESA run is
% recorded) doesn't append zero rows.
for i = 1:numel(add.names)
    if add.nRemoved(i) == 0 && add.varShare(i) == 0
        continue
    end
    j = find(strcmp(cats.names, add.names{i}), 1);
    if isempty(j)
        cats.names{end+1}    = add.names{i};
        cats.nRemoved(end+1) = add.nRemoved(i);
        cats.varShare(end+1) = add.varShare(i);
    else
        cats.nRemoved(j) = cats.nRemoved(j) + add.nRemoved(i);
        cats.varShare(j) = cats.varShare(j) + add.varShare(i);
    end
end
end
