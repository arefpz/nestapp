
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function EEG = artistRejectBadTrials(EEG, opts)
% ARTISTREJECTBADTRIALS  Reject bad trials per ARTIST stage 2 (Wu 2018).
%   EEG = ARTISTREJECTBADTRIALS(EEG, opts) marks epochs as bad if the
%   z-scored absolute amplitude exceeds opts.zThreshold for more than
%   opts.epochInterpFraction of channels, excluding the immediate
%   post-pulse window [opts.excludeStartMs, opts.excludeEndMs].
%
%   Reference:
%     Wu W. et al. (2018). ARTIST: A fully automated artifact rejection
%     algorithm for single-pulse TMS-EEG data. Hum Brain Mapp 39(4):1607.
%     doi:10.1002/hbm.23938. §2.2.1 (stage 2).
%
%   Inputs:
%     EEG   - EEGLAB struct, must be epoched (size(EEG.data,3) > 1).
%     opts  - struct with fields:
%               zThreshold          (default 3)
%               epochInterpFraction (default 0.20)
%               excludeStartMs      (default 0)
%               excludeEndMs        (default 50)
%
%   Output:
%     EEG with bad epochs removed via pop_rejepoch.

    arguments
        EEG  struct
        opts.zThreshold          (1,1) double = 3
        opts.epochInterpFraction (1,1) double = 0.20
        opts.excludeStartMs      (1,1) double = 0
        opts.excludeEndMs        (1,1) double = 50
    end

    if size(EEG.data, 3) < 2
        warning('artistRejectBadTrials:NotEpoched', ...
            'Data is not epoched; skipping bad-trial rejection.');
        return
    end

    [nChan, nTime, nTrial] = size(EEG.data);

    % Build analysis-window mask excluding immediate post-pulse samples.
    inExclude = EEG.times >= opts.excludeStartMs & EEG.times < opts.excludeEndMs;
    analysisMask = ~inExclude;
    if ~any(analysisMask)
        warning('artistRejectBadTrials:EmptyAnalysisWindow', ...
            'Exclusion window covers the whole epoch; using full epoch.');
        analysisMask = true(1, nTime);
    end

    % Per-channel: peak-to-peak amplitude in the analysis window for each trial.
    amp = nan(nChan, nTrial);
    for iChan = 1:nChan
        chData = squeeze(EEG.data(iChan, analysisMask, :));   % time x trial
        amp(iChan, :) = max(chData, [], 1) - min(chData, [], 1);
    end

    % Z-score each channel across trials, then count channels that exceed
    % the threshold per trial.
    z = (amp - mean(amp, 2, 'omitnan')) ./ std(amp, 0, 2, 'omitnan');
    isOutlier = abs(z) > opts.zThreshold;             % nChan x nTrial
    outlierFraction = mean(isOutlier, 1, 'omitnan');  % 1 x nTrial

    badTrials = find(outlierFraction > opts.epochInterpFraction);

    if isempty(badTrials)
        nestLog('ARTIST', 'Bad-trial rejection: no trials marked.');
        return
    end

    nestLog('ARTIST', 'Bad-trial rejection: removing %d / %d epochs (z>%.1f on >%d%% chan).', ...
        numel(badTrials), nTrial, opts.zThreshold, round(100*opts.epochInterpFraction));

    EEG = pop_rejepoch(EEG, badTrials, 0);
end
