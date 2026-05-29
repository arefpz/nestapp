
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function EEG = artistFlagDecayICs(EEG, opts)
% ARTISTFLAGDECAYICS  Flag ICs dominated by post-pulse decay (ARTIST stage 1).
%   EEG = ARTISTFLAGDECAYICS(EEG, opts) marks ICs whose mean rectified
%   activation in the early post-TMS window exceeds a magnitude threshold.
%   Flagged ICs are written to EEG.reject.gcompreject so the existing
%   "Remove Flagged ICA Components" step can remove them.
%
%   Reference:
%     Wu W. et al. (2018). ARTIST: A fully automated artifact rejection
%     algorithm for single-pulse TMS-EEG data. Hum Brain Mapp 39(4):1607.
%     doi:10.1002/hbm.23938. §2.2.1: ICs with "mean magnitude above a
%     certain threshold (30 µV by default) within the first 50 ms after
%     the TMS pulse" are removed.
%
%   Inputs:
%     EEG  - EEGLAB struct with EEG.icaweights / icasphere / icawinv set.
%     opts - struct with fields:
%              winStartMs         (default 0)
%              winEndMs           (default 50)
%              magnitudeThreshold (default 30) in microvolts
%
%   Output:
%     EEG with EEG.reject.gcompreject updated.

    arguments
        EEG  struct
        opts.winStartMs         (1,1) double = 0
        opts.winEndMs           (1,1) double = 50
        opts.magnitudeThreshold (1,1) double = 30
    end

    if ~isfield(EEG, 'icaweights') || isempty(EEG.icaweights)
        error('artistFlagDecayICs:NoICA', ...
            'ICA must be run before decay-IC flagging.');
    end

    icaact = eeg_getica(EEG);   % nComp x nTime x nTrial
    numComps = size(icaact, 1);

    inWin = EEG.times >= opts.winStartMs & EEG.times < opts.winEndMs;
    if ~any(inWin)
        error('artistFlagDecayICs:EmptyWindow', ...
            'Window [%g, %g] ms does not overlap EEG.times.', ...
            opts.winStartMs, opts.winEndMs);
    end

    % Trial-averaged rectified activation in the post-pulse window.
    winMag = nan(numComps, 1);
    for iC = 1:numComps
        trace = abs(mean(icaact(iC, :, :), 3));
        winMag(iC) = mean(trace(inWin));
    end

    flagged = winMag > opts.magnitudeThreshold;

    if ~isfield(EEG, 'reject') || ~isfield(EEG.reject, 'gcompreject') || ...
            numel(EEG.reject.gcompreject) ~= numComps
        EEG.reject.gcompreject = false(1, numComps);
    end
    EEG.reject.gcompreject = EEG.reject.gcompreject(:)' | flagged(:)';

    nestLog('ARTIST', 'Decay-IC flag: %d / %d ICs flagged (mean |act| > %g µV in [%g, %g] ms).', ...
        sum(flagged), numComps, opts.magnitudeThreshold, opts.winStartMs, opts.winEndMs);
end
