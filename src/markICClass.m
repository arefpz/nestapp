
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function EEG = markICClass(EEG, mask, label)
% MARKICCLASS  Tag flagged ICs with a category label for the run report.
%   EEG = MARKICCLASS(EEG, mask, label)
%
%   Records `label` (e.g. 'Muscle', 'Decay') for the ICs selected by the
%   logical `mask` on EEG.etc.nestappICClass, a 1xnComp cellstr. When the
%   "Remove Flagged ICA Components" step later removes flagged ICs, it reads
%   these labels so the report can attribute removed components to a category
%   even when the flags came from a custom classifier (AARATEP / ARTIST)
%   rather than ICLabel. A later classifier wins for an IC flagged twice.
%
%   See also: aaratepMuscleClassifier, artistFlagDecayICs, recordICARound

mask  = logical(mask(:)');
nComp = numel(mask);
if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc)
    EEG.etc = struct();
end
if ~isfield(EEG.etc, 'nestappICClass') || numel(EEG.etc.nestappICClass) ~= nComp
    EEG.etc.nestappICClass = repmat({''}, 1, nComp);
end
for i = find(mask)
    EEG.etc.nestappICClass{i} = label;
end
end
