
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function fileReport = openICARound(fileReport, nComponents)
% OPENICAROUND  Begin a new ICA decomposition round on a PipelineReport.
%   fileReport = OPENICAROUND(fileReport, nComponents)
%
%   Called once per ICA decomposition (Run ICA / Run TESA ICA) so the recorded
%   round count reflects the number of decompositions - even a round that ends
%   up removing nothing is still a round. Component removals fold into the most
%   recent round via addICARemoval.
%
%   See also: addICARemoval, recomputeICATotals, processOneFile

rnd = struct( ...
    'roundNum',    numel(fileReport.ica.rounds) + 1, ...
    'nComponents', nComponents, ...
    'nRejected',   0, ...
    'varRemoved',  NaN, 'varMin', NaN, 'varMax', NaN);
rnd.categories = struct('names', {{}}, 'nRemoved', [], 'varShare', []);

fileReport.ica.rounds{end+1} = rnd;
fileReport.ica = recomputeICATotals(fileReport.ica);
end
