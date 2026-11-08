% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function s = capSummary(labelSets)
% CAPSUMMARY  Which files share a cap, and what each cap rule would keep.
%   s = CAPSUMMARY(labelSets) takes a cell array of cellstr label lists, one
%   per file, and returns:
%     .which         1xN, the cap file i is on, numbered in order of first
%                    appearance
%     .nVariants     how many distinct caps there are
%     .counts        files per cap
%     .modal         1xN logical, files on the cap most files carry (ties go
%                    to the cap seen first)
%     .modalLabels   that cap's labels, in its first file's order
%     .sharedLabels  the modal labels every file carries
%     .dropped       the modal labels some file lacks
%
%   Order and case do not matter: the same electrodes in a different order, or
%   spelled 'Fp1' and 'FP1', are the same cap.
%
%   One definition for groupCurves, which applies the cap rule, and the app,
%   which has to describe both rules BEFORE computing anything so the user can
%   choose. Two copies would let the prompt's numbers drift from the result's.
%
%   See also: groupCurves, groupedCacheEntries, electrodeAvailability

s = struct('which', zeros(1, 0), 'nVariants', 0, 'counts', zeros(0, 1), ...
           'modal', false(1, 0), 'modalLabels', {{}}, 'sharedLabels', {{}}, ...
           'dropped', {{}});
if isempty(labelSets); return; end

sigs = cellfun(@(L) strjoin(sort(lower(cellstr(L))), ','), labelSets, ...
               'UniformOutput', false);
[~, first, which] = unique(sigs, 'stable');
s.which     = which(:)';
s.nVariants = numel(first);
s.counts    = accumarray(s.which(:), 1);

[~, winner]   = max(s.counts);      % max takes the first on a tie
s.modal       = (s.which == winner);
s.modalLabels = labelSets{first(winner)};

% One representative per cap: files on the same cap carry the same labels.
inAll          = electrodeAvailability(s.modalLabels, labelSets(first));
s.sharedLabels = s.modalLabels(inAll);
s.dropped      = s.modalLabels(~inAll);
end
