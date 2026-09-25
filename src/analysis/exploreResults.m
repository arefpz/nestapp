% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function out = exploreResults(res, entries, opts)
% EXPLORERESULTS  Everything a figure was built from, in one struct.
%   out = EXPLORERESULTS(res, entries, opts) packages a groupCurves result for
%   the Results exit: saved as a .mat, or dropped into the base workspace, for
%   someone who wants to keep working in MATLAB.
%
%   out fields:
%     .time        1xT common time base (ms)
%     .groups      per group: name, subjects, curves (nSubj x T), chanMeans,
%                  nFiles, nSubjects, and .files - the individual recordings
%                  behind the subject collapse (.name .subject .curve)
%     .est         the interval per group, as drawn
%     .contrast    the two-group difference estimate, when there is one
%     .channels    the montage the curves were computed on, plus chanlocs
%     .roi         the ROI electrodes
%     .windows     the windows of interest
%     .measures    the per-subject table (exploreMeasures), so the small
%                  tabular form travels with the curves it came from
%     .design      'paired' or 'unpaired'
%     .files       path, subject and group per input file
%     .excluded    files dropped for being on a different cap, and the
%                  subjects dropped from a paired estimate
%     .provenance  app version, EEGLAB version, timestamp, and the montage
%                  report - what is needed to say where a figure came from
%
%   Why a struct and not a CSV. This is curves at sampling rate: 63 channels x
%   ~2000 samples per group, which is millions of cells with no units and no
%   structure once flattened. The measures table is the CSV; this is the object.
%
%   Why provenance is embedded rather than a sidecar file: a sidecar gets
%   separated from what it describes the first time someone copies the result
%   somewhere. The pipeline reports already record their own provenance for the
%   same reason.
%
%   opts:
%     .roi      ROI electrode labels used
%     .windows  windows of interest used
%     .mode       'TEP' | 'GMFP' | 'LMFP'
%     .plot       name of the plot the figure showed, if any
%     .plotParams per-plot settings, name/params struct array
%
%   THIS FILE IS THE SESSION FORMAT. It carries the files table, the ROI, the
%   windows, the design, the plot and its settings, which is everything the tab
%   needs to come back - so File > Load Analysis reads it and there is no
%   separate session artifact holding a second copy of the same state.
%
%   See also: groupCurves, exploreMeasures, exploreDataset, nestappVersion

if nargin < 2; entries = struct('path', {}, 'subject', {}, 'group', {}); end
if nargin < 3; opts = struct(); end
opts = fillDefaults(opts, struct('roi', {{}}, 'windows', [], ...
                                 'mode', 'TEP', 'plot', '', ...
                                 'plotParams', struct('name', {}, 'params', {})));
if isempty(opts.windows); opts.windows = defaultTEPComponentDefs(); end

out          = struct();
out.time     = res.time;
out.groups   = res.groups;
out.est      = res.est;
out.contrast = fieldOr(res, 'contrast', struct([]));
out.design   = res.design;
out.roi      = opts.roi;
out.windows  = opts.windows;
out.mode     = opts.mode;
% The per-plot settings, so reopening this file restores the picture and not
% just the data behind it. Everything else the Explore tab holds was already
% here; without these the file is one field short of being a session.
out.plotParams = opts.plotParams;

out.channels = struct('labels', {res.channelLabels}, 'chanlocs', res.chanlocs);

resForMeasures      = res;
resForMeasures.mode = opts.mode;
out.measures        = exploreMeasures(resForMeasures, opts.windows);

if isempty(entries)
    out.files = struct('path', {}, 'subject', {}, 'group', {}, ...
                       'subjectConfident', {});
else
    % subjectConfident travels too. It is false where a subject rule could
    % not be applied because the file lacked the chosen filename part - drop
    % it and reopening an analysis loses which ids were never really read.
    out.files = struct('path', {entries.path}, 'subject', {entries.subject}, ...
                       'group', {entries.group}, ...
                       'subjectConfident', {entries.subjectConfident});
end

montage = struct();
if isfield(res, 'info') && isfield(res.info, 'montage')
    montage = res.info.montage;
end
out.excluded = struct( ...
    'montageFiles', {fieldOr(montage, 'excluded', {})}, ...
    'unpairedSubjects', {fieldOr(res, 'dropped', {})});

out.provenance = struct( ...
    'nestapp',   nestappVersion(), ...
    'eeglab',    eeglabVersionOrUnknown(), ...
    'matlab',    version(), ...
    'created',   char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'plot',      opts.plot, ...
    'timeBase',  fieldOr(res, 'info', struct()), ...
    'montage',   montage);
end

% ── helpers ─────────────────────────────────────────────────────────────────


function v = eeglabVersionOrUnknown()
% A result may legitimately be built with no EEGLAB on the path (the curve layer
% needs none), so this records what it can rather than failing.
v = 'unknown';
try
    if exist('eeg_getversion', 'file') == 2
        v = eeg_getversion();
    end
catch
end
if isempty(v); v = 'unknown'; end
end
