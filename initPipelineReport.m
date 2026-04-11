function report = initPipelineReport(inputFile)
% INITPIPELINEREPORT  Create a fresh PipelineReport struct for one EEG file.
%
%   report = INITPIPELINEREPORT(inputFile) returns a struct that accumulates
%   processing metrics as each pipeline step runs. Pass the full path to the
%   file being processed. Call EXPORTREPORT after the step loop completes.
%
%   Fields
%     inputFile             - full path to the file being processed
%     processedAt           - datenum when processing started
%     steps                 - cell array of step records (one per completed step)
%     channels.original     - EEG.nbchan at Load Data
%     channels.nRejected    - cumulative channels removed (bad channel steps)
%     channels.nInterpolated- cumulative channels interpolated
%     channels.final        - EEG.nbchan after last step
%     trials.original       - epoch count at first Epoching step
%     trials.rejected       - cumulative rejected epochs
%     trials.final          - epoch count after last step
%     ica.nComponents       - components identified (set after Run ICA / Run TESA ICA)
%     ica.nRejected         - cumulative components removed
%     ica.nKept             - nComponents - nRejected
%     ica.varRemoved        - total % data variance removed (NaN if unavailable)
%     ica.varMin / varMax   - per-component variance range (NaN if unavailable)
%     ica.categories        - ICLabel breakdown (populated when ICLabel was run)
%     teps                  - struct for TEP results (reserved for M5)
%
%   Each element of .steps is a struct with fields:
%     name, chansBefore, chansAfter, trialsBefore, trialsAfter, duration, timestamp
%
%   See also: exportReport, runPipeline

report.inputFile   = inputFile;
report.processedAt = now;   % datenum — compatible with exportReport's datestr() calls

report.steps = {};

report.channels.original      = 0;
report.channels.nRejected     = 0;
report.channels.nInterpolated = 0;
report.channels.final         = 0;

report.trials.original = 0;
report.trials.rejected = 0;
report.trials.final    = 0;

report.ica.nComponents = 0;
report.ica.nRejected   = 0;
report.ica.nKept       = 0;
report.ica.varRemoved  = NaN;
report.ica.varMin      = NaN;
report.ica.varMax      = NaN;
report.ica.rounds      = {}; % per-round structs for multi-round TESA

ICA_CATEGORIES = {'Brain','Muscle','Eye','Heart','Line Noise','Ch Noise','Other'};
report.ica.categories.names    = ICA_CATEGORIES;
report.ica.categories.nRemoved = zeros(1, 7);
report.ica.categories.varShare = zeros(1, 7);

report.teps = struct();
end
