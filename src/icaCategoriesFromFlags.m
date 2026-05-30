
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function cats = icaCategoriesFromFlags(rMask, p)
% ICACATEGORIESFROMFLAGS  Per-category breakdown of removed ICA components.
%   cats = ICACATEGORIESFROMFLAGS(rMask, p)
%
%   Builds the category tally for the components removed by "Remove Flagged
%   ICA Components". Each rejected IC (rMask true) is attributed to exactly one
%   category by priority:
%     1) a custom-classifier label (AARATEP 'Muscle', ARTIST 'Decay'), from
%        p.classLabels (1xnComp cellstr written by markICClass), if present
%     2) else the ICLabel argmax category, from p.iclabelProbs, if ICLabel ran
%     3) else 'Manual'
%
%   Inputs
%     rMask - 1xnComp logical, true for removed components
%     p     - pendingICAStats struct, optional fields: classLabels (cellstr),
%             iclabelProbs (nComp x 7), compVarPct (1xnComp % data variance)
%
%   Output
%     cats - struct(names{1xC}, nRemoved(1xC), varShare(1xC)) listing only the
%            categories that actually appear, for recordICARound to merge.
%
%   See also: markICClass, recordICARound, processOneFile

ICLABEL_NAMES = {'Brain','Muscle','Eye','Heart','Line Noise','Ch Noise','Other'};
rMask = logical(rMask(:)');
n     = numel(rMask);
rej   = find(rMask);
catOf = repmat({''}, 1, n);

if isfield(p, 'classLabels') && numel(p.classLabels) == n
    for c = rej
        if ~isempty(p.classLabels{c}); catOf{c} = p.classLabels{c}; end
    end
end
if isfield(p, 'iclabelProbs') && size(p.iclabelProbs, 1) >= n
    [~, best] = max(p.iclabelProbs, [], 2);
    for c = rej
        if isempty(catOf{c}); catOf{c} = ICLABEL_NAMES{best(c)}; end
    end
end
for c = rej
    if isempty(catOf{c}); catOf{c} = 'Manual'; end
end

hasVar = isfield(p, 'compVarPct') && numel(p.compVarPct) == n;
names = {}; nRemoved = []; varShare = [];
for c = rej
    nm  = catOf{c};
    vsc = 0;
    if hasVar; vsc = p.compVarPct(c); end
    j = find(strcmp(names, nm), 1);
    if isempty(j)
        names{end+1}    = nm;   %#ok<AGROW>
        nRemoved(end+1) = 1;    %#ok<AGROW>
        varShare(end+1) = vsc;  %#ok<AGROW>
    else
        nRemoved(j) = nRemoved(j) + 1;
        varShare(j) = varShare(j) + vsc;
    end
end
cats = struct('names', {names}, 'nRemoved', nRemoved, 'varShare', varShare);
end
