function cleanup = hideFromPath(fnName)
% HIDEFROMPATH  Temporarily remove a function from the MATLAB path.
%
%   Syntax
%     cleanup = hideFromPath(fnName)
%
%   Inputs
%     fnName  char - a function name (e.g. 'tesa_peakanalysis').
%
%   Outputs
%     cleanup  onCleanup - restores the path when it goes out of scope.
%
%   Removes every path folder that provides fnName so which(fnName) returns
%   empty for the duration of a test - letting "behaviour when this plugin is
%   absent" tests run deterministically on a fully-equipped machine instead
%   of skipping. The caller should ASSERT that which(fnName) is now empty so
%   the test fails loudly if the function could not be hidden (rather than
%   silently passing). The onCleanup restores the path even if the test errors.
%
%   Example
%     c = hideFromPath('tesa_peakanalysis');
%     testCase.assertEmpty(which('tesa_peakanalysis'), 'could not hide TESA');
%     ... exercise the no-TESA path ...
%
%   See also: which, rmpath, onCleanup

builtinPrefix = 'built-in';
hits = which(fnName, '-all');
if ~iscell(hits); hits = {hits}; end

removed = {};
for i = 1:numel(hits)
    p = hits{i};
    if isempty(p) || startsWith(p, builtinPrefix)
        continue
    end
    d = fileparts(p);
    if isempty(d); continue; end
    onPath = any(strcmp(strsplit(path, pathsep), d));
    if onPath
        rmpath(d);
        removed{end+1} = d; %#ok<AGROW>
    end
end

cleanup = onCleanup(@() restore(removed));
end

function restore(dirs)
for i = 1:numel(dirs)
    if isfolder(dirs{i})
        addpath(dirs{i});
    end
end
end
