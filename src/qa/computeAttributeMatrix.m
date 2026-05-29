
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [SM, summary] = computeAttributeMatrix(EEG, opts)
% COMPUTEATTRIBUTEMATRIX  Per-(channel, trial) noise / amplitude attribute.
%   [SM, summary] = COMPUTEATTRIBUTEMATRIX(EEG, opts) returns an
%   nbchan x nTrials matrix of a log-scale quality attribute, ported from
%   the Get_SM function in TMSEEG-5.0 (tmseeg_rm_ch_tr_1.m).
%
%   opts.attribute  : 'minmax' | 'minmax_no_tms' | 'highfreq'
%                       (default 'minmax_no_tms')
%   opts.tmsWindow  : [tStart tEnd] in ms, excluded in the no-TMS mode
%                       (default [0 25])
%   opts.freqBand   : [fmin fmax] Hz bandpass before scoring
%                       (default [1 80])
%
%   Modes
%     'minmax'         log(max - min) across the full epoch
%     'minmax_no_tms'  log(max - min) excluding the TMS pulse window
%     'highfreq'       log(sum |first derivative|) - muscle / movement score
%
%   Output
%     SM       nbchan x nTrials log-attribute matrix; flat or saturated
%              channels are returned as NaN rows.
%     summary  struct with:
%                .perChanMedian   (nbchan x 1)  NaN for flagged channels
%                .perTrialMedian  (1 x nTrials)
%                .flatChanMask    (nbchan x 1 logical) var(channel) < eps
%                .satChanMask     (nbchan x 1 logical) max|data| > 250 uV
%                .nbchan, .nTrials, .attribute
%
%   Reference: Frehlich, Mei, Garcia Dominguez, Farzan, Schwartzmann (TMSEEG)

if nargin < 2 || ~isstruct(opts), opts = struct(); end
if ~isfield(opts, 'attribute'),  opts.attribute = 'minmax_no_tms'; end
if ~isfield(opts, 'tmsWindow'),  opts.tmsWindow = [0 25];          end
if ~isfield(opts, 'freqBand'),   opts.freqBand  = [1 80];          end

validAttrs = qualityAttributeModes();
if ~any(strcmp(opts.attribute, validAttrs))
    error('computeAttributeMatrix:badAttribute', ...
        'Unknown attribute "%s". Valid: %s', opts.attribute, strjoin(validAttrs, ', '));
end

% Coerce expected fields.
data    = EEG.data;
times   = EEG.times;
srate   = EEG.srate;
nbchan  = size(data, 1);
nPnts   = size(data, 2);
nTrials = size(data, 3);
if nTrials == 0
    % Continuous data: treat the whole record as one "trial" for the
    % attribute computation. Less informative but avoids crashing on
    % pre-epoched files placed at the Load Data checkpoint.
    data    = reshape(data, nbchan, nPnts, 1);
    nTrials = 1;
end

SAT_THRESHOLD_UV = 250;  % matches the "saturated electrode" heuristic in QC literature

% Flag flat (zero variance) and saturated channels up-front. Excluded from
% the bandpass loop and returned as NaN rows so the heatmap renders them
% distinctly from merely-noisy channels.
chanVar       = var(reshape(data, nbchan, []), 0, 2);
flatChanMask  = chanVar(:) < eps;
chanPeak      = max(abs(reshape(data, nbchan, [])), [], 2);
satChanMask   = chanPeak(:) > SAT_THRESHOLD_UV;
skipMask      = flatChanMask | satChanMask;

% Pre-compute time indices (mirrors TMSEEG t_st/t_end and p_st/p_end).
t_st  = find(times >= opts.tmsWindow(1) - 1e6, 1, 'first');   % effectively start
t_end = find(times <= opts.tmsWindow(2) + 1e6, 1, 'last');    % effectively end
if isempty(t_st),  t_st  = 1;     end
if isempty(t_end), t_end = nPnts; end

p_st  = find(times < opts.tmsWindow(1), 1, 'last');
p_end = find(times > opts.tmsWindow(2), 1, 'first');
if isempty(p_st),  p_st  = 0;        end
if isempty(p_end), p_end = nPnts + 1; end

% Bandpass filter design (matches TMSEEG: 2nd-order Butterworth via SOS).
ord = 2;
nyq = srate / 2;
fLo = max(opts.freqBand(1), 0.01);
fHi = min(opts.freqBand(2), nyq - 0.1);
[z, p, k] = butter(ord, [fLo fHi] / nyq, 'bandpass');
[sos, g]  = zp2sos(z, p, k);

% Time window for the scoring (full epoch or TMS-excluded).
switch opts.attribute
    case 'minmax'
        timeIdx = t_st:t_end;
    case {'minmax_no_tms', 'highfreq'}
        timeIdx = [t_st:p_st, p_end:t_end];
        timeIdx = timeIdx(timeIdx >= 1 & timeIdx <= nPnts);
end

if isempty(timeIdx)
    error('computeAttributeMatrix:emptyWindow', ...
        ['No time samples available after applying TMS window [%g %g] ms. ', ...
        'Either widen tmsWindow or pick a different attribute.'], ...
        opts.tmsWindow(1), opts.tmsWindow(2));
end

% Score every (channel, trial). Channels in skipMask are left as NaN.
SM = nan(nbchan, nTrials);
goodCh = find(~skipMask);
if ~isempty(goodCh)
    for trl = 1:nTrials
        seg = double(data(goodCh, timeIdx, trl))';   % nSamp x nGoodCh
        filt = filtfilt(sos, g, seg);                 % nSamp x nGoodCh
        switch opts.attribute
            case {'minmax', 'minmax_no_tms'}
                trialVal = log(max(filt, [], 1)' - min(filt, [], 1)');
            case 'highfreq'
                trialVal = log(sum(abs(diff(filt, 1, 1))', 2));
        end
        SM(goodCh, trl) = trialVal;
    end
end

summary.perChanMedian  = nan(nbchan, 1);
summary.perTrialMedian = nan(1, nTrials);
if any(~skipMask)
    summary.perChanMedian(~skipMask) = median(SM(~skipMask, :), 2, 'omitnan');
    summary.perTrialMedian           = median(SM(~skipMask, :), 1, 'omitnan');
end
summary.flatChanMask = flatChanMask;
summary.satChanMask  = satChanMask;
summary.nbchan       = nbchan;
summary.nTrials      = nTrials;
summary.attribute    = opts.attribute;
end
