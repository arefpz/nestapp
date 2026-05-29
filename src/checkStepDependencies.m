
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [ok, msg] = checkStepDependencies(stepNames, filePaths)
% CHECKSTEPDEPENDENCIES  Verify required plugins are on the MATLAB path.
%
%   [ok, msg] = checkStepDependencies(stepNames, filePaths)
%
%   Inputs
%     stepNames  - cell array of step name strings (app.SelectedListBox.Items)
%     filePaths  - cell array of full file paths being processed (app.file).
%                  Used to filter format-specific loaders (e.g. bva-io only
%                  checked when .vhdr files are selected).
%
%   Outputs
%     ok   - true if all dependencies are satisfied
%     msg  - formatted message listing missing plugins ('' when ok=true)
%
%   Checks are done via which(), so the function must be called after
%   EEGLAB has been added to the path (i.e. after app startup).

if nargin < 2
    filePaths = {};
end

% Vendored AARATEP helpers ship with nestapp under third_party/ but are
% only added to the path lazily during step dispatch. Add them now - only
% when an AARATEP step is actually selected - so the which() probes below
% see the bundled functions instead of reporting them as missing plugins.
% Gating avoids a ~280-file genpath walk on every non-AARATEP pre-flight.
aaratepSteps = {'Interpolate Missing Data (AR-Blend)', ...
                'Remove Decay Artifact', ...
                'Flag ICA Components (AARATEP Muscle)'};
if any(ismember(stepNames, aaratepSteps))
    try
        ensureAaratepOnPath();
    catch
        % If the vendored tree is genuinely absent, the which() checks
        % below report the AARATEP steps as missing with the bundled note.
    end
end

% Build extension set from file paths for format-specific dep filtering.
[~,~,extList] = cellfun(@fileparts, filePaths, 'UniformOutput', false);
exts = unique(lower(extList));

steps    = stepRegistry();
nameList = {steps.name};

% missing: containers.Map keyed by plugin name
missing = containers.Map('KeyType','char','ValueType','any');

for i = 1:numel(stepNames)
    idx = find(strcmp(nameList, stepNames{i}), 1);
    if isempty(idx); continue; end
    reqs = steps(idx).requires;
    if isempty(reqs); continue; end
    for j = 1:numel(reqs)
        r = reqs(j);
        % Skip format-specific loaders when no files are selected or no file of
        % that format is present in the selection.
        if ~isempty(r.fileExt) && (isempty(filePaths) || ~any(strcmpi(exts, r.fileExt)))
            continue
        end
        if isempty(which(r.fn))
            if ~isKey(missing, r.plugin)
                missing(r.plugin) = struct( ...
                    'installNote', r.installNote, ...
                    'steps',       {{}}, ...
                    'fns',         {{}});
            end
            entry = missing(r.plugin);
            entry.steps{end+1} = stepNames{i};
            entry.fns{end+1}   = r.fn;
            missing(r.plugin)  = entry;
        end
    end
end

if isempty(missing)
    ok  = true;
    msg = '';
    return
end

ok    = false;
lines = {'Missing plugins - install before running:', ''};
pluginNames = keys(missing);
for i = 1:numel(pluginNames)
    plugin = pluginNames{i};
    entry  = missing(plugin);
    usedBy = strjoin(unique(entry.steps), ', ');
    lines{end+1} = sprintf('  %s', plugin);                          %#ok<AGROW>
    lines{end+1} = sprintf('    Steps:   %s', usedBy);               %#ok<AGROW>
    lines{end+1} = sprintf('    Install: %s', entry.installNote);    %#ok<AGROW>
    lines{end+1} = '';                                                %#ok<AGROW>
end
msg = strjoin(lines, newline);
end
