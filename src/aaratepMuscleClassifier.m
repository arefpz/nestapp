function EEG = aaratepMuscleClassifier(EEG, opts)
% AARATEPMUSCLECLASSIFIER  Flag ICs as TMS-induced muscle (AARATEP).
%   EEG = AARATEPMUSCLECLASSIFIER(EEG, opts) marks ICs whose mean
%   rectified activation in a short post-pulse window is much larger than
%   the rest of the epoch, indicating localized TMS-evoked muscle bursts.
%   Flagged ICs are written to EEG.reject.gcompreject so the existing
%   "Remove Flagged ICA Components" step can remove them.
%
%   Algorithm is lifted verbatim from
%   c_TMSEEG_Preprocess_AARATEPPipeline.m (lines 400-412) and credited
%   under the AARATEPPipeline MIT license (see THIRD_PARTY_NOTICES.md).
%
%   Reference:
%     Cline C.C. et al. (2021). Advanced Artifact Removal for Automated
%     TMS-EEG Data Processing. IEEE NER. doi:10.1109/NER49283.2021.9441147.
%
%   Inputs:
%     EEG  - EEGLAB struct with ICA decomposition.
%     opts - struct with fields:
%              winStartMs        (default 11)
%              winEndMs          (default 30)
%              muscleThreshold   (default 8)  ratio of window-mean to epoch-mean
%
%   Output:
%     EEG with EEG.reject.gcompreject updated.

    arguments
        EEG  struct
        opts.winStartMs      (1,1) double = 11
        opts.winEndMs        (1,1) double = 30
        opts.muscleThreshold (1,1) double = 8
    end

    if ~isfield(EEG, 'icaweights') || isempty(EEG.icaweights)
        error('aaratepMuscleClassifier:NoICA', ...
            'ICA must be run before AARATEP muscle classification.');
    end

    icaact = eeg_getica(EEG);
    numComps = size(EEG.icaweights, 1);

    indicesInWin = EEG.times >= opts.winStartMs & EEG.times < opts.winEndMs;
    if ~any(indicesInWin)
        error('aaratepMuscleClassifier:EmptyWindow', ...
            'Window [%g, %g] ms does not overlap EEG.times.', ...
            opts.winStartMs, opts.winEndMs);
    end

    tmsMuscleRatio = nan(numComps, 1);
    for iC = 1:numComps
        muscleScore       = abs(mean(icaact(iC, :, :), 3));
        winScore          = mean(muscleScore(:, indicesInWin), 2);
        tmsMuscleRatio(iC) = winScore / mean(muscleScore);
    end

    toReject = tmsMuscleRatio > opts.muscleThreshold;

    if ~isfield(EEG, 'reject') || ~isfield(EEG.reject, 'gcompreject') || ...
            numel(EEG.reject.gcompreject) ~= numComps
        EEG.reject.gcompreject = false(1, numComps);
    end
    EEG.reject.gcompreject = EEG.reject.gcompreject(:)' | toReject(:)';

    nestLog('AARATEP', 'Muscle classifier: %d / %d ICs flagged (ratio > %g in [%g, %g] ms).', ...
        sum(toReject), numComps, opts.muscleThreshold, opts.winStartMs, opts.winEndMs);
end
