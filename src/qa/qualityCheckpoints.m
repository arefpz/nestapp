function names = qualityCheckpoints()
% QUALITYCHECKPOINTS  Step names after which a quality PNG is captured.
%   names = QUALITYCHECKPOINTS() returns a curated cell array of step
%   names. When the autoQualityReport preference is on, processOneFile
%   renders a quality figure after any step whose name is in this list.
%
%   Edit this file to add or remove checkpoint moments. A unit test
%   (test_qualityCheckpoints) asserts every entry exists in stepRegistry.
names = { ...
    'Load Data', ...                        % raw
    'Load Channel Location', ...            % after coords
    'Remove Bad Channels', ...              % after bad-channel removal
    'Clean Artifacts', ...
    'Automatic Cleaning Data', ...
    'Epoching', ...                          % epoched but un-ICA'd
    'Remove ICA Components (TESA)', ...     % each round (TESA recipes call twice)
    'Remove Flagged ICA Components', ...
    'Interpolate Channels', ...             % near-final
    'Re-Reference', ...                      % final-ish
    'Quality Gate'};                         % capture EEG state at every gate
end
