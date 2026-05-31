function buildTemplates()
% BUILDTEMPLATES  Regenerate all built-in pipeline template .mat files.
%
%   Run this from the MATLAB command window (after run_nestapp, or from
%   within src/) whenever template definitions or stepRegistry defaults
%   change.  The generated .mat files in src/templates/ are committed to
%   version control and loaded at runtime - no override logic runs in the
%   app itself.
%
%   See also: stepRegistry, nestapp

addpath(fileparts(mfilename('fullpath')));  % ensure src/ is on path
reg    = stepRegistry();
outDir = fullfile(fileparts(mfilename('fullpath')), 'templates');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% 1 - TMS-EEG / TEP (TESA)
% Two-round FastICA pipeline per Rogasch et al. 2017, restructured to match
% the TESA User Manual step order precisely:
%   - Bad channels removed before epoching (manual step 5)
%   - Full-epoch demean before any ICA (manual step 7)
%   - Bad trial removal before Round 1 ICA (manual step 11)
%   - TMS artifact re-cut before each ICA round (manual steps 12, 17)
%   - Cut extended to 15 ms + re-interpolated after Round 1 (manual steps 14-15)
%   - Bandpass + 60 Hz bandstop placed between ICA rounds (manual step 16)
%   - Final TMS interpolation after Round 2 ICA (manual step 19)
%   - Re-reference after channel interpolation (manual step 21)
steps = { ...
    'Load Data', 'Load Channel Location', 'Remove un-needed Channels', ...
    'Find TMS Pulses (TESA)', 'Remove Bad Channels', ...
    'Epoching', 'Remove Baseline', ...
    'Remove TMS Artifacts (TESA)', 'Interpolate Missing Data (TESA)', 'Re-Sample', ...
    'Remove Bad Epoch', ...
    'Remove TMS Artifacts (TESA)', ...
    'Run TESA ICA', 'Remove ICA Components (TESA)', ...
    'Remove TMS Artifacts (TESA)', 'Interpolate Missing Data (TESA)', ...
    'Frequency Filter (TESA)', 'Frequency Filter (TESA)', ...
    'Remove TMS Artifacts (TESA)', ...
    'Run TESA ICA', 'Remove ICA Components (TESA)', ...
    'Interpolate Missing Data (TESA)', ...
    'Interpolate Channels', 'Re-Reference', ...
    'Remove Baseline', 'Save New Set'};
ovs = emptyOvs(steps);
ovs = setOv(ovs, steps, 'Epoching', 'timelim', [-1, 1]);
% Demean over full epoch before ICA (manual step 7: subtract mean of entire epoch).
ovs = setOv(ovs, steps, 'Remove Baseline', 'timerange', [-1000 1000], 1);
% TMS artifact cut windows: default [-2 10] ms for occurrences 1-2,
% extended to [-2 15] ms after Round 1 (occurrences 3-4).
ovs = setOv(ovs, steps, 'Remove TMS Artifacts (TESA)', 'cutTimesTMS', [-2 15], 3);
ovs = setOv(ovs, steps, 'Remove TMS Artifacts (TESA)', 'cutTimesTMS', [-2 15], 4);
% Interpolation: narrow window pre-downsample (~5 samples at 5 kHz),
% wider window post-downsample (5 samples at 1 kHz).
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpolation', 'cubic', 1);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpWin',     [1 1],   1);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpolation', 'cubic', 2);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpWin',     [5 5],   2);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpolation', 'cubic', 3);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpWin',     [5 5],   3);
% Filters between ICA rounds (manual step 16): occ 1 = bandpass, occ 2 = bandstop.
% Bandpass uses defaults (1-80 Hz, order 4).
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'type', 'bandstop', 2);
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'high', 58,         2);
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'low',  62,         2);
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'ord',  2,          2);
% Re-reference after channel interpolation so the montage is complete (manual step 21).
ovs = setOv(ovs, steps, 'Re-Reference', 'ref', '[]');
% Round 2 ICA: all artifact detectors on; Round 1 default (tmsMuscle only) is correct.
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'blink',     'on', 2);
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'move',      'on', 2);
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'muscle',    'on', 2);
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'elecNoise', 'on', 2);
ovs = setOv(ovs, steps, 'Save New Set', 'savenew', 'tesa');
saveMat(reg, steps, ovs, 'TMS-EEG / TEP (TESA)', fullfile(outDir, '1_tesa_tep.mat'));

%% 2 - Resting-State EEG
% Continuous-data pipeline per PREP (Bigdely-Shamlo 2015) with structural
% improvements from Delorme 2023 (Sci Rep, doi:10.1038/s41598-023-27528-0).
steps = { ...
    'Load Data', 'Load Channel Location', 'Remove un-needed Channels', ...
    'Frequency Filter', 'Frequency Filter (CleanLine)', 'Automatic Cleaning Data', ...
    'Re-Reference', 'Run ICA', 'Label ICA Components', ...
    'Flag ICA Components for Rejection', 'Remove Flagged ICA Components', ...
    'Interpolate Channels', 'Save New Set'};
ovs = emptyOvs(steps);
ovs = setOv(ovs, steps, 'Frequency Filter', 'locutoff', 0.5);
ovs = setOv(ovs, steps, 'Frequency Filter', 'hicutoff', 40);
ovs = setOv(ovs, steps, 'Automatic Cleaning Data', 'FlatlineCriterion', 4);
ovs = setOv(ovs, steps, 'Automatic Cleaning Data', 'ChannelCriterion',  0.85);
ovs = setOv(ovs, steps, 'Re-Reference', 'ref', '[]');
% 0.8 is the practical threshold for resting data; 0.9 default flags nothing.
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Muscle', [0.8, 1]);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Eye',    [0.8, 1]);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Heart',  [0.9, 1]);
ovs = setOv(ovs, steps, 'Save New Set', 'savenew', 'resting');
saveMat(reg, steps, ovs, 'Resting-State EEG', fullfile(outDir, '2_resting_state.mat'));

%% 3 - Minimal (Delorme 2023)
% Minimal pipeline per Delorme 2023 ("EEG is better left alone", Sci Rep).
% HPF 0.5 Hz only; no LPF; no explicit re-reference.
steps = { ...
    'Load Data', 'Load Channel Location', 'Remove un-needed Channels', ...
    'Frequency Filter', 'Automatic Cleaning Data', 'Run ICA', ...
    'Label ICA Components', 'Flag ICA Components for Rejection', ...
    'Remove Flagged ICA Components', 'Interpolate Channels', 'Save New Set'};
ovs = emptyOvs(steps);
ovs = setOv(ovs, steps, 'Frequency Filter', 'locutoff', 0.5);
ovs = setOv(ovs, steps, 'Frequency Filter', 'hicutoff', 0);   % 0 = no LPF
ovs = setOv(ovs, steps, 'Automatic Cleaning Data', 'FlatlineCriterion', 4);
ovs = setOv(ovs, steps, 'Automatic Cleaning Data', 'ChannelCriterion',  0.85);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Muscle', [0.8, 1]);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Eye',    [0.8, 1]);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Heart',  [0.9, 1]);
ovs = setOv(ovs, steps, 'Save New Set', 'savenew', 'minimal');
saveMat(reg, steps, ovs, 'Minimal (Delorme 2023)', fullfile(outDir, '3_minimal.mat'));

fprintf('buildTemplates: done - %s\n', outDir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local helpers

function saveMat(reg, steps, ovs, templateName, outPath)
% Build a v3 pipeline spec and save it in the same format as a user-saved pipeline.
    n    = numel(steps);
    spec = repmat(struct('name', '', 'params', struct()), 1, n);
    for i = 1:n
        s        = makePipelineStep(steps{i}, reg);
        ovFields = fieldnames(ovs{i});
        for fi = 1:numel(ovFields)
            key = ovFields{fi};
            if ~isfield(s.params, key)
                error('buildTemplates:badKey', ...
                    'Key "%s" is not a param of step "%s".', key, steps{i});
            end
            s.params.(key) = ovs{i}.(key);
        end
        spec(i) = s;
    end
    pipelineName = templateName;
    save(outPath, 'pipelineName', 'spec');
    fprintf('  %s\n', outPath);
end

function ovs = emptyOvs(steps)
    ovs = repmat({struct()}, 1, numel(steps));
end

function ovs = setOv(ovs, steps, stepName, key, value, occurrence)
% Set one override field, identified by step name and optional occurrence
% index (for steps that appear more than once, e.g. two ICA rounds).
    if nargin < 6; occurrence = 1; end
    idx = find(strcmp(steps, stepName));
    if isempty(idx)
        error('buildTemplates:badStep', 'Step "%s" not found.', stepName);
    end
    if occurrence > numel(idx)
        error('buildTemplates:badOccurrence', ...
              'Step "%s" has %d occurrence(s); requested %d.', ...
              stepName, numel(idx), occurrence);
    end
    ovs{idx(occurrence)}.(key) = value;
end
