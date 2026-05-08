function [results, warnings] = batchTEPExtract(filePaths, roiElectrodes, varargin)
% BATCHTEPEXTRACT  Extract TEP peaks across multiple .set files.
%
%   [results, warnings] = batchTEPExtract(filePaths, roiElectrodes)
%   [results, warnings] = batchTEPExtract(__, 'compDefs', compDefs)
%   [results, warnings] = batchTEPExtract(__, 'smoothWin', 5)
%   [results, warnings] = batchTEPExtract(__, 'csvPath', 'out.csv')
%   [results, warnings] = batchTEPExtract(__, 'progressFcn', @(i,n) ...)
%   [results, warnings] = batchTEPExtract(__, 'loadFcn', @myLoader)
%
%   Inputs
%     filePaths     - 1×N cell array of full paths to .set files
%     roiElectrodes - cell array of electrode label strings for the ROI
%
%   Optional name-value arguments
%     'compDefs'    struct array (tepPeakFinder format). Default: canonical 6-component.
%     'smoothWin'   moving-average window in samples (default 5, same as plotTEP).
%     'csvPath'     write results to this CSV path. Empty = no file written.
%     'progressFcn' function handle @(iFile, nFiles) called before each file.
%     'loadFcn'     function handle @(fullPath) -> EEG. Default: pop_loadset.
%                   Use for testing (pass synthetic EEG structs without EEGLAB I/O).
%
%   Outputs
%     results  - MATLAB table, one row per file × per component (long format).
%                Columns: file, roi_electrodes, component, polarity, found,
%                latency_ms, amplitude_uv, win_start_ms, win_end_ms,
%                nom_latency_ms, n_trials, n_channels_roi.
%     warnings - cell array of strings, one per skipped or degraded file.
%
%   CSV format: long format (one row per file × component). Researchers pivot
%   to wide in R/SPSS as needed; merge with subject metadata via the 'file' column.
%
%   Example
%     files   = {'sub01_cleaned.set', 'sub02_cleaned.set'};
%     files   = fullfile('/data', files);
%     roi     = {'FC1', 'FC3', 'C1', 'C3'};
%     [T, w]  = batchTEPExtract(files, roi, 'csvPath', 'tep_peaks.csv');
%
%   See also: tepPeakFinder, plotTEP

% ── parse inputs ──────────────────────────────────────────────────────────
p = inputParser;
p.addRequired('filePaths',      @(x) iscell(x) && ~isempty(x));
p.addRequired('roiElectrodes',  @(x) iscell(x) && ~isempty(x));
p.addParameter('compDefs',    [], @(x) isempty(x) || isstruct(x));
p.addParameter('smoothWin',   5,  @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('csvPath',     '', @ischar);
p.addParameter('progressFcn', [], @(x) isempty(x) || isa(x,'function_handle'));
p.addParameter('loadFcn',     [], @(x) isempty(x) || isa(x,'function_handle'));
p.parse(filePaths, roiElectrodes, varargin{:});
opts = p.Results;

% ── guard: TESA required ───────────────────────────────────────────────────
if isempty(opts.loadFcn) && isempty(which('tesa_peakanalysis'))
    error('batchTEPExtract:noTESA', ...
        ['TESA toolbox not found. Add TESA to the MATLAB path before ' ...
         'calling batchTEPExtract.']);
end

% ── resolve canonical component definitions ────────────────────────────────
% Resolve early so buildNaNRows can use them without calling tepPeakFinder.
compDefs = opts.compDefs;
if isempty(compDefs)
    compDefs = defaultCompDefs();
end

% ── per-file loop ──────────────────────────────────────────────────────────
nFiles  = numel(opts.filePaths);
rows    = {};
warnings = {};

for fi = 1:nFiles
    if ~isempty(opts.progressFcn)
        opts.progressFcn(fi, nFiles);
    end

    fp = opts.filePaths{fi};
    [fdir, fname, fext] = fileparts(fp); %#ok<ASGLU> fdir/fext referenced inside evalc string

    %% Load
    try
        if ~isempty(opts.loadFcn)
            EEG = opts.loadFcn(fp);
        else
            evalc("EEG = pop_loadset('filename', [fname fext], 'filepath', fdir)");
        end
    catch ME
        warnings{end+1} = sprintf('%s: load failed — %s', fname, ME.message); %#ok<AGROW>
        rows = appendNaNRows(rows, fname, opts.roiElectrodes, compDefs);
        continue
    end

    %% Guard: must be epoched
    if ~isstruct(EEG) || ~isfield(EEG,'trials') || EEG.trials < 2
        warnings{end+1} = sprintf('%s: skipped — not epoched (trials=%d)', ...
            fname, getfield_safe(EEG, 'trials', 0)); %#ok<AGROW>
        rows = appendNaNRows(rows, fname, opts.roiElectrodes, compDefs);
        continue
    end

    %% Resolve ROI channels
    allLabels  = {EEG.chanlocs.labels};
    roiPresent = intersect(opts.roiElectrodes, allLabels, 'stable');
    nRoiFound  = numel(roiPresent);
    if nRoiFound == 0
        warnings{end+1} = sprintf('%s: skipped — none of the requested ROI electrodes found', fname); %#ok<AGROW>
        rows = appendNaNRows(rows, fname, opts.roiElectrodes, compDefs);
        continue
    end
    if nRoiFound < numel(opts.roiElectrodes)
        warnings{end+1} = sprintf('%s: partial ROI — %d of %d requested electrodes found', ...
            fname, nRoiFound, numel(opts.roiElectrodes)); %#ok<AGROW>
    end
    roiIdx = find(ismember(allLabels, roiPresent));

    %% Grand-mean waveform (ROI average, trial average, smoothed)
    roiData  = mean(EEG.data(roiIdx, :, :), 1);          % 1 × T × nTrials
    waveform = mean(roiData, 3);                           % 1 × T
    waveform = smoothdata(waveform, 'movmean', opts.smoothWin);

    %% Peak detection
    try
        peaks = tepPeakFinder(waveform, EEG.times, compDefs);
    catch ME
        warnings{end+1} = sprintf('%s: peak detection failed — %s', fname, ME.message); %#ok<AGROW>
        % Pass nRoiFound so the NaN rows still record how many ROI channels were found.
        rows = appendNaNRows(rows, fname, opts.roiElectrodes, compDefs, nRoiFound);
        continue
    end

    %% Append one row per component
    roiStr = strjoin(roiPresent, ',');
    for ci = 1:numel(peaks)
        pk = peaks(ci);
        cd = compDefs(ci);
        rows{end+1} = { ...
            fname,          ...
            roiStr,         ...
            pk.name,        ...
            pk.polarity,    ...
            double(pk.found), ...
            pk.latencyMs,   ...
            pk.amplitudeUV, ...
            cd.winStart,    ...
            cd.winEnd,      ...
            cd.nomLatency,  ...
            double(EEG.trials), ...
            double(nRoiFound) ...
        }; %#ok<AGROW>
    end
end

% ── assemble table ─────────────────────────────────────────────────────────
COL_NAMES = {'file','roi_electrodes','component','polarity','found', ...
    'latency_ms','amplitude_uv','win_start_ms','win_end_ms', ...
    'nom_latency_ms','n_trials','n_channels_roi'};

if isempty(rows)
    results = cell2table(cell(0, numel(COL_NAMES)), 'VariableNames', COL_NAMES);
else
    results = cell2table(vertcat(rows{:}), 'VariableNames', COL_NAMES);
end

% ── write CSV ──────────────────────────────────────────────────────────────
if ~isempty(opts.csvPath)
    writeCSV(results, opts.csvPath);
end
end

%% ── local helpers ────────────────────────────────────────────────────────

function rows = appendNaNRows(rows, fname, roiElectrodes, compDefs, nRoiFound)
% Append one NaN-filled row per component for a file that was skipped or failed.
% nRoiFound: how many ROI channels were actually found before the failure (default 0).
if nargin < 5; nRoiFound = 0; end
roiStr = strjoin(roiElectrodes, ',');
for ci = 1:numel(compDefs)
    cd = compDefs(ci);
    rows{end+1} = { ...
        fname, roiStr, cd.name, cd.polarity, ...
        0, NaN, NaN, cd.winStart, cd.winEnd, cd.nomLatency, 0, nRoiFound ...
    }; %#ok<AGROW>
end
end

function defs = defaultCompDefs()
% Canonical 6-component windows following Beck et al. 2024 (Hum Brain Mapp, 45:e70048).
defs = struct( ...
    'name',       {'N15',  'P30',  'N45',  'P60',  'N100', 'P180'}, ...
    'polarity',   {'neg',  'pos',  'neg',  'pos',  'neg',  'pos'}, ...
    'nomLatency', {15,     30,     45,     60,     100,    180}, ...
    'winStart',   {10,     20,     40,     50,     70,     150}, ...
    'winEnd',     {20,     40,     55,     70,     150,    240});
end

function v = getfield_safe(s, field, default)
% Return s.(field) if it exists, else default.
if isstruct(s) && isfield(s, field)
    v = s.(field);
else
    v = default;
end
end

function writeCSV(T, csvPath)
% Write results table to CSV using manual fprintf (matching ExportReportsCSV pattern).
fid = fopen(csvPath, 'w');
if fid == -1
    error('batchTEPExtract:csvWriteFailed', 'Cannot open file for writing: %s', csvPath);
end
cleanup = onCleanup(@() fclose(fid));
% Header
fprintf(fid, '%s\n', strjoin(T.Properties.VariableNames, ','));
% Rows — roi_electrodes may contain commas, so quote it
for ri = 1:height(T)
    fprintf(fid, '%s,"%s",%s,%s,%d,%.4f,%.4f,%.1f,%.1f,%.1f,%d,%d\n', ...
        T.file{ri}, T.roi_electrodes{ri}, T.component{ri}, T.polarity{ri}, ...
        T.found(ri), T.latency_ms(ri), T.amplitude_uv(ri), ...
        T.win_start_ms(ri), T.win_end_ms(ri), T.nom_latency_ms(ri), ...
        T.n_trials(ri), T.n_channels_roi(ri));
end
end
