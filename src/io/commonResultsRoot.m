
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function root = commonResultsRoot(filePaths)
% COMMONRESULTSROOT  Common parent folder for a list of files.
%   root = COMMONRESULTSROOT(filePaths) returns the directory that holds
%   every file in filePaths when they all share one parent, or the
%   parent of the first file when they don't. Used as the default
%   anchor for the batch output folder when the user hasn't set
%   nestapp.outputRoot in Preferences.
parents = cellfun(@(p) fileparts(p), filePaths, 'UniformOutput', false);
if numel(unique(parents)) == 1
    root = parents{1};
else
    root = fileparts(filePaths{1});
end
end
