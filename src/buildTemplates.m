
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
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
% Final pre-save baseline correction: TESA / Rogasch 2017 convention
% uses a short pre-stimulus window so each TEP trace is aligned to
% its own pre-pulse baseline before peak analysis.
ovs = setOv(ovs, steps, 'Remove Baseline', 'timerange', [-500 -10],   2);
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

%% 4 - TMS-EEG / TEP (TESA + Quality Gates)
% Same step order as template 1, with four Quality Gate checkpoints
% inserted at the most diagnostic handoffs. All gates run in absolute
% mode with reasonable starting thresholds for a typical 100-150 pulse
% single-pulse TMS protocol. Override the per-metric WarnAt fields to
% tighten the Marginal cutoff without moving the Fail boundary; leave
% any threshold at 0 to disable it.
%
% Skip-on-fail honors setpref('nestapp','skipOnQualityFail',true) so a
% failing file aborts its own pipeline but leaves the rest of the batch
% to run.
%
%   QG1  After Epoching (post Find TMS Pulses + Remove Bad Channels) -
%        trigger sanity (minTriggers) plus a cap on how many channels
%        Remove Bad Channels was allowed to throw away. SM / ICA
%        metrics still need a fit ICA decomposition.
%   QG2  After 'Remove Bad Epoch' (epoched, post-rejection) - enough
%        trials survived, channel rank is intact, no run-away flat
%        channels.
%   QG3  After Round 1 ICA cleanup + interpolation - TMS-evoked muscle
%        was actually removed (EMG fraction low) and the trial pool is
%        still reasonable.
%   QG4  After Re-Reference (final, pre-save) - comprehensive QC:
%        bad-channel / bad-trial fractions and residual electrode
%        artifact comps.
steps = { ...
    'Load Data', 'Load Channel Location', 'Remove un-needed Channels', ...
    'Find TMS Pulses (TESA)', ...
    'Remove Bad Channels', ...
    'Epoching', ...
    'Quality Gate', ...
    'Remove Baseline', ...
    'Remove TMS Artifacts (TESA)', 'Interpolate Missing Data (TESA)', 'Re-Sample', ...
    'Remove Bad Epoch', ...
    'Quality Gate', ...
    'Remove TMS Artifacts (TESA)', ...
    'Run TESA ICA', 'Remove ICA Components (TESA)', ...
    'Remove TMS Artifacts (TESA)', 'Interpolate Missing Data (TESA)', ...
    'Quality Gate', ...
    'Frequency Filter (TESA)', 'Frequency Filter (TESA)', ...
    'Remove TMS Artifacts (TESA)', ...
    'Run TESA ICA', 'Remove ICA Components (TESA)', ...
    'Interpolate Missing Data (TESA)', ...
    'Interpolate Channels', 'Re-Reference', ...
    'Quality Gate', ...
    'Remove Baseline', 'Save New Set'};
ovs = emptyOvs(steps);

% --- Inherit every override from template 1 ---------------------------
ovs = setOv(ovs, steps, 'Epoching', 'timelim', [-1, 1]);
ovs = setOv(ovs, steps, 'Remove Baseline', 'timerange', [-1000 1000], 1);
ovs = setOv(ovs, steps, 'Remove Baseline', 'timerange', [-500 -10],   2);
ovs = setOv(ovs, steps, 'Remove TMS Artifacts (TESA)', 'cutTimesTMS', [-2 15], 3);
ovs = setOv(ovs, steps, 'Remove TMS Artifacts (TESA)', 'cutTimesTMS', [-2 15], 4);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpolation', 'cubic', 1);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpWin',     [1 1],   1);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpolation', 'cubic', 2);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpWin',     [5 5],   2);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpolation', 'cubic', 3);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (TESA)', 'interpWin',     [5 5],   3);
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'type', 'bandstop', 2);
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'high', 58,         2);
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'low',  62,         2);
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'ord',  2,          2);
ovs = setOv(ovs, steps, 'Re-Reference', 'ref', '[]');
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'blink',     'on', 2);
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'move',      'on', 2);
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'muscle',    'on', 2);
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'elecNoise', 'on', 2);
ovs = setOv(ovs, steps, 'Save New Set', 'savenew', 'tesa_qc');

% --- QG1: post-epoching -----------------------------------------------
% Runs after Find TMS Pulses + Remove Bad Channels + Epoching, so
% trigger count and channels-removed are both meaningful here.
% Trigger fail < 50 pulses (file is broken); warn < 70 (below most
% TMS protocols). Channels: fail if removeBadChannels threw away
% > 10% of channels (warn cutoff = 0.8*10 = 8%). Trial-rejection is
% checked on QG2 instead - that gate runs after Remove Bad Epoch
% where rejectedTrialPct actually has signal.
ovs = setOv(ovs, steps, 'Quality Gate', 'gateLabel',           'post-triggers', 1);
ovs = setOv(ovs, steps, 'Quality Gate', 'minTriggers',         50,              1);
ovs = setOv(ovs, steps, 'Quality Gate', 'minTriggersWarnAt',   70,              1);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxRejectedChanPct',  10,              1);

% --- QG2: post-epoch-rejection ----------------------------------------
% Epoched, post-bad-trial removal. Fail if fewer than 50 trials remain
% (insufficient for averaging), warn at 80. Cap cumulative trial
% rejection at 15% (default warn cutoff 12%). Rank should stay close
% to the post-bad-channel-removal count; 0.9 catches gross rank
% deficiency. Flat-channel count > 2 here is a smell - removeBadChannels
% should have caught those already.
ovs = setOv(ovs, steps, 'Quality Gate', 'gateLabel',           'post-bad-epoch', 2);
ovs = setOv(ovs, steps, 'Quality Gate', 'minTrials',           50,               2);
ovs = setOv(ovs, steps, 'Quality Gate', 'minTrialsWarnAt',     80,               2);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxRejectedTrialPct', 15,               2);
ovs = setOv(ovs, steps, 'Quality Gate', 'minRankRatio',        0.9,              2);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxFlatChans',        5,                2);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxFlatChansWarnAt',  2,                2);

% --- QG3: post-Round 1 ICA + interpolation ----------------------------
% Round 1 targets TMS-evoked muscle. Fail if EMG / muscle fraction
% remains > 40% of ICA components; warn at 20%. Cumulative trial
% rejection > 25% means Round 1 lost too many trials (Round 2 still
% has another chance), warn at 10%.
ovs = setOv(ovs, steps, 'Quality Gate', 'gateLabel',                 'post-round1-ica', 3);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxEMGFraction',            0.40,              3);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxEMGFractionWarnAt',      0.20,              3);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxRejectedTrialPct',       25,                3);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxRejectedTrialPctWarnAt', 10,                3);

% --- QG4: final, pre-save ---------------------------------------------
% Comprehensive QC after both ICA rounds, filters, and re-referencing.
% Tighter than QG3 because no further cleanup follows. Cumulative
% rejection counts: > 10% trials or > 15% channels lost = Fail.
ovs = setOv(ovs, steps, 'Quality Gate', 'gateLabel',                  'final', 4);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxRejectedTrialPct',         10,     4);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxRejectedTrialPctWarnAt',   5,      4);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxRejectedChanPct',          15,     4);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxRejectedChanPctWarnAt',    5,      4);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxElectrodeCount',           3,      4);
ovs = setOv(ovs, steps, 'Quality Gate', 'maxElectrodeCountWarnAt',     1,      4);
saveMat(reg, steps, ovs, 'TMS-EEG / TEP (TESA + Quality Gates)', fullfile(outDir, '4_tesa_tep_qc.mat'));

%% 5 - TMS-EEG / ARTIST (paper-faithful, TESA compselect fallback)
% ARTIST (Wu et al. 2018, Hum Brain Mapp 39(4):1607,
% doi:10.1002/hbm.23938). Three-stage TMS-EEG pipeline:
%   Stage 1: cut + interpolate TMS pulse, downsample, baseline, bandpass.
%   Stage 2: z-score trial rejection + RANSAC bad-channel rejection.
%   Stage 3: two ICA rounds. Round 1 catches decay components (mean |act|
%            in [0, 50] ms > 30 uV); Round 2 with PCA dim-reduction +
%            FLDA classifier in the paper. We substitute TESA compselect
%            because the trained FLDA weights are not publicly distributed.
% Numeric values trace to Wu 2018 §2.2.1 - §2.2.3.
% Order per Wu 2018 §2.2.2 ("Preprocessing": epoch + baseline + filter)
% followed by §2.2.1 stages (pulse cut, downsample, stage 2 reject,
% stage 3 ICA). Filtering comes BEFORE bad-trial / bad-channel rejection.
% Two Remove Baseline occurrences match template 1's pattern: an early
% full-epoch demean for ICA stability, then a final pre-save TEP baseline
% [-300, -100] ms per Wu 2018 §2.2.3.
% Notch is TESA Butterworth bandstop 58-62 Hz (paper-faithful substitute
% for FIR notch; safe on epoched data, unlike CleanLine).
% Stage-3 IC classification: TESA pop_tesa_compselect substitutes for
% ARTIST's 23-feature FLDA, which is not publicly distributed.
% Known gap: paper specifies PCA -> 99.9% variance for round-2 ICA; the
% Run ICA step does not yet expose a pca parameter. See THIRD_PARTY_NOTICES.md.
steps = { ...
    'Load Data', 'Load Channel Location', 'Remove un-needed Channels', ...
    'Find TMS Pulses (TESA)', ...
    'Epoching', ...
    'Remove Baseline', ...
    'Remove TMS Artifacts (TESA)', 'Interpolate Missing Data (TESA)', ...
    'Re-Sample', ...
    'Frequency Filter', ...
    'Frequency Filter (TESA)', ...
    'Reject Bad Trials (ARTIST)', ...
    'Remove Bad Channels (ARTIST)', ...
    'Run ICA', ...
    'Flag ICA Components (ARTIST Decay)', ...
    'Remove Flagged ICA Components', ...
    'Run ICA', ...
    'Remove ICA Components (TESA)', ...
    'Interpolate Channels', 'Re-Reference', ...
    'Remove Baseline', 'Save New Set'};
ovs = emptyOvs(steps);
ovs = setOv(ovs, steps, 'Epoching',            'timelim',   [-0.5, 1.5]);
% First baseline: full-epoch demean for ICA stability (matches template 1).
ovs = setOv(ovs, steps, 'Remove Baseline',     'timerange', [-500, 1500], 1);
ovs = setOv(ovs, steps, 'Remove TMS Artifacts (TESA)', 'cutTimesTMS', [0, 10]);
ovs = setOv(ovs, steps, 'Re-Sample',           'freq',      1000);
ovs = setOv(ovs, steps, 'Frequency Filter',    'locutoff',  1);
ovs = setOv(ovs, steps, 'Frequency Filter',    'hicutoff',  100);
% 60 Hz notch via TESA bandstop (paper-faithful substitute for FIR notch).
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'type', 'bandstop');
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'high', 58);
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'low',  62);
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'ord',  2);
% Round 1 ICA: infomax, no PCA. Round 2: infomax (TESA compselect needs runica).
ovs = setOv(ovs, steps, 'Run ICA', 'icatype', 'runica', 1);
ovs = setOv(ovs, steps, 'Run ICA', 'icatype', 'runica', 2);
% Stage-3 compselect: TMS-EEG-aware fallback for the missing FLDA. All
% detectors on, mirroring the round-2 settings of template 1.
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'blink',     'on');
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'move',      'on');
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'muscle',    'on');
ovs = setOv(ovs, steps, 'Remove ICA Components (TESA)', 'elecNoise', 'on');
ovs = setOv(ovs, steps, 'Re-Reference',        'ref',       '[]');
% Final pre-save baseline per Wu 2018 §2.2.3.
ovs = setOv(ovs, steps, 'Remove Baseline',     'timerange', [-300, -100], 2);
ovs = setOv(ovs, steps, 'Save New Set',        'savenew',   'artist');
saveMat(reg, steps, ovs, ...
    'TMS-EEG / ARTIST (paper-faithful, TESA compselect fallback)', ...
    fullfile(outDir, '5_artist.mat'));

%% 6 - TMS-EEG / AARATEP (Cline 2021)
% AARATEP (Cline et al. 2021, IEEE NER, doi:10.1109/NER49283.2021.9441147).
% Helpers vendored from chriscline/AARATEPPipeline under MIT.
% SOUND uses TESA's default lead field (subject-individual anatomy is
% intentionally out of scope - override leadFieldPath in the SOUND step
% to use your own). Numeric values trace to
% third_party/aaratep/c_TMSEEG_Preprocess_AARATEPPipeline.m.
% Order traced verbatim from c_TMSEEG_Preprocess_AARATEPPipeline.m:
%   - pop_reref([]) is applied BEFORE the early eye-IC ICA (line 244)
%     so the early ICA decomposition is on average-referenced data
%   - the early eye-IC pass calls c_TMSEEG_runICLabel with
%     muscleComponentThreshold = NaN (line 254): ONLY Eye gets flagged
%   - round-2 IC rejection ORs muscle flags onto ICLabel flags before a
%     single pop_subcomp - so Label/Flag (ICLabel) must come BEFORE the
%     AARATEP muscle classifier (otherwise pop_icflag overwrites the
%     muscle flags via gcompreject assignment)
%   - upstream uses c_EEG_filter_butterworth bandstop (NOT CleanLine);
%     we use Frequency Filter (TESA) bandstop 58-62 Hz as the equivalent
%     (same algorithm class: Butterworth + filtfilt zero-phase).
% Known gap: upstream Remove Bad Channels uses a {PREP_deviation,
%   TESA_DDWiener_PerTrial} ensemble via c_TMSEEG_detectBadChannels;
%   nestapp's Remove Bad Channels uses pop_rejchan kurtosis with
%   threshold 10. Documented in THIRD_PARTY_NOTICES.md.
steps = { ...
    'Load Data', 'Load Channel Location', 'Remove un-needed Channels', ...
    'Find TMS Pulses (TESA)', ...
    'Epoching', ...
    'Interpolate Missing Data (AR-Blend)', ...
    'Re-Sample', ...
    'Remove Baseline', ...
    'Frequency Filter', ...
    'Remove Bad Channels', ...
    'Re-Reference', ...
    'Run ICA', ...
    'Label ICA Components', ...
    'Flag ICA Components for Rejection', ...
    'Remove Flagged ICA Components', ...
    'Remove Recording Noise (SOUND)', ...
    'Remove Decay Artifact', ...
    'Interpolate Missing Data (AR-Blend)', ...
    'Frequency Filter (TESA)', ...
    'Run ICA', ...
    'Label ICA Components', ...
    'Flag ICA Components for Rejection', ...
    'Flag ICA Components (AARATEP Muscle)', ...
    'Remove Flagged ICA Components', ...
    'Frequency Filter', ...
    'Re-Reference', ...
    'Save New Set'};
ovs = emptyOvs(steps);
% Epoching: example value from README; users can override.
ovs = setOv(ovs, steps, 'Epoching',  'timelim',   [-1, 2]);
% AR-Blend interp: artifactTimespan = [-0.002, 0.012] s -> [-2, 12] ms.
% prePostFitDurations = [20, 20] ms hardcoded upstream.
ovs = setOv(ovs, steps, 'Interpolate Missing Data (AR-Blend)', 'artifactStartMs', -2, 1);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (AR-Blend)', 'artifactEndMs',   12, 1);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (AR-Blend)', 'prePostFitMs',    20, 1);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (AR-Blend)', 'artifactStartMs', -2, 2);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (AR-Blend)', 'artifactEndMs',   12, 2);
ovs = setOv(ovs, steps, 'Interpolate Missing Data (AR-Blend)', 'prePostFitMs',    20, 2);
ovs = setOv(ovs, steps, 'Re-Sample',                     'freq',                    1000);
ovs = setOv(ovs, steps, 'Remove Baseline',               'timerange',               [-500, -10]);
% Two-stage filter: HPF at 1 Hz first, LPF at 200 Hz second.
ovs = setOv(ovs, steps, 'Frequency Filter', 'locutoff', 1,   1);
ovs = setOv(ovs, steps, 'Frequency Filter', 'hicutoff', 0,   1);
ovs = setOv(ovs, steps, 'Frequency Filter', 'locutoff', 0,   2);
ovs = setOv(ovs, steps, 'Frequency Filter', 'hicutoff', 200, 2);
% Bad-channel detection: threshold 10 per badChannelThreshold = 10.
% Note: nestapp's Remove Bad Channels uses pop_rejchan kurtosis; upstream
% uses an ensemble. Same threshold value, different algorithm. Documented.
ovs = setOv(ovs, steps, 'Remove Bad Channels',           'threshold',               10);
% Early average reference before the early eye-IC ICA (upstream line 244).
ovs = setOv(ovs, steps, 'Re-Reference', 'ref', '[]', 1);
% Both ICA rounds: fastica.
ovs = setOv(ovs, steps, 'Run ICA', 'icatype', 'fastica', 1);
ovs = setOv(ovs, steps, 'Run ICA', 'icatype', 'fastica', 2);
% Early eye-IC pass: c_TMSEEG_runICLabel called with eyeThreshold=0.9 and
% muscle/brain/other = NaN (line 254). Only Eye is flagged. Explicitly
% disable Muscle (registry default is [0.9, 1]) and all others.
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Brain',        [NaN, NaN], 1);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Muscle',       [NaN, NaN], 1);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Eye',          [0.9, 1],   1);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Heart',        [NaN, NaN], 1);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'LineNoise',    [NaN, NaN], 1);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'ChannelNoise', [NaN, NaN], 1);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Other',        [NaN, NaN], 1);
% SOUND: lambda = 10^-1.5 ~ 0.0316; iter stays at default 10.
ovs = setOv(ovs, steps, 'Remove Recording Noise (SOUND)', 'lambdaValue', 10^-1.5);
% Decay removal: artifactTimespan = [-2, 12] ms, doDecayRemovalPerTrial = true.
ovs = setOv(ovs, steps, 'Remove Decay Artifact', 'artifactStartMs', -2);
ovs = setOv(ovs, steps, 'Remove Decay Artifact', 'artifactEndMs',   12);
ovs = setOv(ovs, steps, 'Remove Decay Artifact', 'perTrial',        'on');
% Line noise: 60 Hz Butterworth bandstop (paper-faithful substitute
% for upstream's c_EEG_filter_butterworth bandstop). Tight 58-62 Hz band
% per AARATEP's typical practice.
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'type', 'bandstop');
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'high', 58);
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'low',  62);
ovs = setOv(ovs, steps, 'Frequency Filter (TESA)', 'ord',  2);
% Round-2 IC rejection: pop_icflag writes gcompreject from ICLabel
% probabilities FIRST, then the AARATEP muscle classifier ORs its
% flags into gcompreject. Reverse order would let pop_icflag clobber
% the muscle flags. Brain < 0.3 = reject (brainComponentThreshold from
% upstream). Eye/Heart/Line/ChannelNoise use standard 0.8 ICLabel
% thresholds (c_TMSEEG_runICLabel's defaults for the categories upstream
% does not explicitly override).
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Brain',       [0, 0.3], 2);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Eye',         [0.8, 1], 2);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'Heart',       [0.8, 1], 2);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'LineNoise',   [0.8, 1], 2);
ovs = setOv(ovs, steps, 'Flag ICA Components for Rejection', 'ChannelNoise',[0.8, 1], 2);
% Final average reference (upstream line 498).
ovs = setOv(ovs, steps, 'Re-Reference',  'ref',     '[]', 2);
ovs = setOv(ovs, steps, 'Save New Set',  'savenew', 'aaratep');
saveMat(reg, steps, ovs, 'TMS-EEG / AARATEP (Cline 2021)', ...
    fullfile(outDir, '6_aaratep.mat'));

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
