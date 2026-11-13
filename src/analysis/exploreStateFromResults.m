% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [state, report] = exploreStateFromResults(source)
% EXPLORESTATEFROMRESULTS  Read a saved Results file back into Explore state.
%   [state, report] = EXPLORESTATEFROMRESULTS(source)
%
%   source may be:
%     - the struct exploreResults built, or
%     - the struct load() returns for a .mat containing one, under ANY variable
%       name. The Results exit saves it as `out` and assigns it as `tepResults`,
%       and a user is free to rename it, so the variable name is not part of the
%       format.
%
%   state:
%     .entries    1-by-N struct (path, subject, group, subjectConfident) -
%                 the files table, in exploreDataset's shape
%     .roi        cellstr of electrode labels
%     .windows    windows of interest
%     .design     'paired' | 'unpaired'
%     .mode       'TEP' | 'GMFP' | 'LMFP'
%     .plot       plot name to select, '' when the file does not name one
%     .plotParams name/params struct array of per-plot settings
%     .capRule    'exclude' | 'shared' for mixed caps, '' when none was needed
%
%   report:
%     .ok       false when the source is not a results file at all
%     .missing  paths in the file that are not on disk now
%     .notes    cellstr of things the caller should tell the user
%
%   THE RESULTS FILE IS THE SESSION FORMAT. It already carries the files table,
%   the ROI, the windows, the design and the plot name, so a separate session
%   file would be a second artifact holding the same state - the "four confusing
%   exports" problem this rework exists to remove. What was missing was never a
%   way to write the state out; it was a way to read it back.
%
%   Paths are NOT resolved to the stored curves. A saved result holds group
%   averages, which is enough to redraw the figure and nothing else: change the
%   ROI or move one recording between groups and the averages are wrong. Resuming
%   WORK means recomputing from the files, so this returns the paths and reports
%   which of them have gone, and the caller reloads. Restoring a picture that
%   silently stops responding to the rail would be worse than refusing.
%
%   Nothing here is fatal except "this is not a results file". A file listing
%   recordings that have since moved is still worth opening for the groups, ROI
%   and windows it defines - that is the expensive human judgement the file
%   exists to preserve - so the missing paths are reported and the rest loads.
%
%   See also: exploreResults, exploreDataset, loadReducedSets, groupCurves

state  = emptyState();
report = struct('ok', false, 'missing', {{}}, 'notes', {{}});

out = findResultsStruct(source);
if isempty(out)
    report.notes{end+1} = ['That file does not contain a nestapp analysis. ' ...
        'Open a .mat saved by the Results exit.'];
    return
end
report.ok = true;

state.entries = normaliseEntries(fieldOr(out, 'files', []));
state.roi     = cellstr(fieldOr(out, 'roi', {}));
state.windows = fieldOr(out, 'windows', []);
state.design  = charOr(fieldOr(out, 'design', 'unpaired'), 'unpaired');
state.mode    = charOr(fieldOr(out, 'mode', 'TEP'), 'TEP');

prov = fieldOr(out, 'provenance', struct());
state.plot = charOr(fieldOr(prov, 'plot', ''), '');
state.plotParams = fieldOr(out, 'plotParams', struct('name', {}, 'params', {}));

% The cap rule the analysis was computed under, so reopening reproduces it
% without asking again. A file from before the rule was recorded that still
% reports mixed caps was computed under 'exclude' - the only rule there was.
montage = fieldOr(prov, 'montage', struct());
state.capRule = lower(charOr(fieldOr(montage, 'rule', ''), ''));
if isempty(state.capRule) && isequal(fieldOr(montage, 'uniform', true), false)
    state.capRule = 'exclude';
end
if ~any(strcmp(state.capRule, {'', 'exclude', 'shared'}))
    state.capRule = '';
end

if isempty(state.entries)
    report.notes{end+1} = 'The file names no recordings, so no groups were restored.';
end

% Which recordings have moved since the file was written.
present = arrayfun(@(e) isfile(e.path), state.entries);
if any(~present)
    report.missing = {state.entries(~present).path};
    state.entries  = state.entries(present);
    report.notes{end+1} = sprintf( ...
        ['%d of %d recording%s could not be found and %s left out. The groups, ' ...
         'ROI and windows were restored; re-add the files to include them.'], ...
        sum(~present), numel(present), plural(numel(present)), ...
        isAre(sum(~present)));
end

% Provenance is a note, never a refusal - a result written by another version is
% still the analysis someone wants back.
saved = charOr(fieldOr(prov, 'nestapp', ''), '');
if ~isempty(saved) && ~strcmp(saved, nestappVersion())
    report.notes{end+1} = sprintf( ...
        'Saved by nestapp %s; this is %s.', saved, nestappVersion());
end
end

% ── helpers ─────────────────────────────────────────────────────────────────

function s = emptyState()
s = struct('entries', struct('path', {}, 'subject', {}, 'group', {}, ...
                             'subjectConfident', {}), ...
           'roi', {{}}, 'windows', [], 'design', 'unpaired', 'mode', 'TEP', ...
           'plot', '', 'plotParams', struct('name', {}, 'params', {}), ...
           'capRule', '');
end

function out = findResultsStruct(source)
% A results struct, or the first variable in a loaded file that is one. Keyed on
% the fields that make it an analysis rather than on a variable name, so a file
% saved as `out`, assigned as `tepResults`, or renamed by the user all open.
out = [];
if ~isstruct(source) || isempty(source); return; end
if looksLikeResults(source); out = source; return; end

names = fieldnames(source);
for k = 1:numel(names)
    v = source.(names{k});
    if isstruct(v) && isscalar(v) && looksLikeResults(v)
        out = v;
        return
    end
end
end

function tf = looksLikeResults(s)
tf = all(isfield(s, {'files', 'windows', 'design'}));
end

function e = normaliseEntries(files)
% exploreDataset's shape exactly, so the restored table is indistinguishable
% from one built by adding folders.
e = struct('path', {}, 'subject', {}, 'group', {}, 'subjectConfident', {});
if isempty(files) || ~isstruct(files); return; end
for k = 1:numel(files)
    e(k).path    = charOr(fieldOr(files(k), 'path', ''), '');
    e(k).subject = charOr(fieldOr(files(k), 'subject', ''), '');
    e(k).group   = charOr(fieldOr(files(k), 'group', ''), '');
    % Absent means confident: a file written before the flag was saved, or one
    % assembled by hand, states ids deliberately. Only a recorded false is a
    % guess, so the benefit of the doubt goes the way that does not invent
    % warnings about ids nobody guessed.
    e(k).subjectConfident = logical(fieldOr(files(k), 'subjectConfident', true));
end
end


function c = charOr(v, default)
if ischar(v) && isrow(v)
    c = v;
elseif isstring(v) && isscalar(v)
    c = char(v);
else
    c = default;
end
end
