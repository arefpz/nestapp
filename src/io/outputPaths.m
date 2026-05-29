
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function dirPath = outputPaths(batchCtx, kind, stem)
% OUTPUTPATHS  Resolve the directory for one artifact category.
%   dirPath = OUTPUTPATHS(batchCtx, kind)
%   dirPath = OUTPUTPATHS(batchCtx, kind, stem)
%
%   The single authority on where each artifact lives inside a batch
%   folder. Filename composition stays in each writer (e.g.
%   reportArtifactName for *.pdf / *.mat); this function only owns the
%   *directory* and creates it on demand.
%
%   batchCtx : struct returned by buildBatchContext.
%   kind     : 'data'    .set output destination
%              'reports' .pdf / .mat / per-file session log
%              'qc'      per-checkpoint PNG folder for one input file
%              'batch'   run-level artifacts (dashboard, spec, csv,
%                          cross-file log)
%   stem     : input file basename (required for 'data', 'reports', 'qc';
%              unused for 'batch').
%
%   Layout map:
%
%     typeBased
%       data    -> <batchRoot>/data/
%       reports -> <batchRoot>/reports/
%       qc      -> <batchRoot>/qc/<stem>/
%       batch   -> <batchRoot>/batch/
%
%     perInput
%       data    -> <batchRoot>/<stem>/
%       reports -> <batchRoot>/<stem>/
%       qc      -> <batchRoot>/<stem>/qc/
%       batch   -> <batchRoot>/_batch/
%
%   See also: buildBatchContext, commonResultsRoot.

if nargin < 3, stem = ''; end

needsStem = any(strcmp(kind, {'data', 'reports', 'qc'}));
if needsStem && isempty(stem)
    error('nestapp:outputPaths:missingStem', ...
        'outputPaths(kind="%s") requires a non-empty stem.', kind);
end

switch batchCtx.layout
    case 'typeBased'
        dirPath = typeBasedPath(batchCtx.batchRoot, kind, stem);
    case 'perInput'
        dirPath = perInputPath(batchCtx.batchRoot, kind, stem);
    otherwise
        error('nestapp:outputPaths:badLayout', ...
            'Unknown layout "%s".', batchCtx.layout);
end

if ~exist(dirPath, 'dir')
    mkdir(dirPath);
end
end

% -- per-layout dispatchers -----------------------------------------------

function p = typeBasedPath(batchRoot, kind, stem)
switch kind
    case 'data',    p = fullfile(batchRoot, 'data');
    case 'reports', p = fullfile(batchRoot, 'reports');
    case 'qc',      p = fullfile(batchRoot, 'qc', stem);
    case 'batch',   p = fullfile(batchRoot, 'batch');
    otherwise,      error('nestapp:outputPaths:badKind', ...
                          'Unknown kind "%s".', kind);
end
end

function p = perInputPath(batchRoot, kind, stem)
switch kind
    case 'data',    p = fullfile(batchRoot, stem);
    case 'reports', p = fullfile(batchRoot, stem);
    case 'qc',      p = fullfile(batchRoot, stem, 'qc');
    case 'batch',   p = fullfile(batchRoot, '_batch');
    otherwise,      error('nestapp:outputPaths:badKind', ...
                          'Unknown kind "%s".', kind);
end
end
