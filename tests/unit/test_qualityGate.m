classdef test_qualityGate < matlab.unittest.TestCase
% TEST_QUALITYGATE  Unit tests for src/qa/qualityGate.m

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            tc.addTeardown(@() rmpath(genpath(fullfile(root, 'src'))));
        end
    end

    methods (Static, Access = private)
        function EEG = makeEEG(nbchan, nTrials, nPnts, srate)
            EEG.data   = randn(nbchan, nPnts, nTrials);
            EEG.times  = (0:nPnts-1) / srate * 1000 - 200;
            EEG.srate  = srate;
            EEG.nbchan = nbchan;
            EEG.trials = nTrials;
            EEG.pnts   = nPnts;
            EEG.event  = struct('type', {}, 'latency', {});
        end

        function EEG = withEvents(EEG, n)
            for k = 1:n
                EEG.event(k).type    = 'TMS';
                EEG.event(k).latency = k * 100;
            end
        end
    end

    methods (Test)
        % -- disabled / default behavior ---------------------------------

        function disabled_gate_passes(tc)
            EEG = test_qualityGate.makeEEG(16, 40, 500, 1000);
            gate = qualityGate(EEG, struct());
            tc.verifyEqual(gate.verdict, 'Pass');
            tc.verifyEmpty(gate.reasons);
            tc.verifyEqual(gate.mode, 'absolute');
        end

        function records_label_and_thresholds(tc)
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);
            params = struct('gateLabel', 'after-load', 'expectedChans', 8);
            gate = qualityGate(EEG, params);
            tc.verifyEqual(gate.label, 'after-load');
            tc.verifyEqual(gate.thresholds.expectedChans, 8);
            tc.verifyEqual(gate.verdict, 'Pass');
        end

        % -- exact-match checks ------------------------------------------

        function expectedChans_mismatch_fails(tc)
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);
            gate = qualityGate(EEG, struct('expectedChans', 64));
            tc.verifyEqual(gate.verdict, 'Fail');
            tc.verifyTrue(any(contains(gate.reasons, 'nbchan')));
        end

        function expectedSrate_mismatch_fails(tc)
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);
            gate = qualityGate(EEG, struct('expectedSrate', 500));
            tc.verifyEqual(gate.verdict, 'Fail');
            tc.verifyTrue(any(contains(gate.reasons, 'srate')));
        end

        % -- min/max thresholds + marginal slack -------------------------

        function minTriggers_fail_below_slack(tc)
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);
            EEG = test_qualityGate.withEvents(EEG, 50);   % below 0.8 * 100 = 80
            gate = qualityGate(EEG, struct('minTriggers', 100, 'marginalSlack', 0.8));
            tc.verifyEqual(gate.verdict, 'Fail');
        end

        function minTriggers_marginal_above_threshold_with_warnAt(tc)
            % WarnAt sits above threshold for min checks. With
            % minTriggers = 100, minTriggersWarnAt = 120, a value of 110
            % is in the marginal band [100, 120).
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);
            EEG = test_qualityGate.withEvents(EEG, 110);
            gate = qualityGate(EEG, struct( ...
                'minTriggers',       100, ...
                'minTriggersWarnAt', 120));
            tc.verifyEqual(gate.verdict, 'Marginal');
        end

        function minTriggers_pass_at_or_above(tc)
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);
            EEG = test_qualityGate.withEvents(EEG, 150);
            gate = qualityGate(EEG, struct('minTriggers', 100));
            tc.verifyEqual(gate.verdict, 'Pass');
        end

        function maxFlatChans_counts_zero_var_channels(tc)
            EEG = test_qualityGate.makeEEG(8, 30, 500, 1000);
            EEG.data(3, :, :) = 0;
            EEG.data(7, :, :) = 0;
            gate = qualityGate(EEG, struct('maxFlatChans', 1));
            tc.verifyEqual(gate.verdict, 'Fail');
            tc.verifyEqual(gate.metrics.nFlatChans, 2);
        end

        function maxSatChans_counts_saturated_channels(tc)
            EEG = test_qualityGate.makeEEG(8, 30, 500, 1000);
            EEG.data(1, :, :) = 500;  % > 250 uV
            % With threshold 5 and value 1, 1 < 0.8*5=4 -> Pass.
            gate = qualityGate(EEG, struct('maxSatChans', 5));
            tc.verifyEqual(gate.metrics.nSatChans, 1);
            tc.verifyEqual(gate.verdict, 'Pass');
            % With threshold 1 and value 1, 1 > 0.8 -> Marginal (right at edge).
            gate = qualityGate(EEG, struct('maxSatChans', 1));
            tc.verifyEqual(gate.verdict, 'Marginal');
        end

        function minRankRatio_catches_rank_deficient_data(tc)
            EEG = test_qualityGate.makeEEG(8, 1, 500, 1000);
            % Force rank deficiency: make 4 channels identical.
            EEG.data(2, :) = EEG.data(1, :);
            EEG.data(3, :) = EEG.data(1, :);
            EEG.data(4, :) = EEG.data(1, :);
            gate = qualityGate(EEG, struct('minRankRatio', 0.9));
            tc.verifyEqual(gate.verdict, 'Fail');
        end

        function minTrials_catches_low_trial_count(tc)
            EEG = test_qualityGate.makeEEG(8, 5, 500, 1000);
            gate = qualityGate(EEG, struct('minTrials', 30, 'marginalSlack', 0.8));
            tc.verifyEqual(gate.verdict, 'Fail');
        end

        % -- worst-of logic ----------------------------------------------

        function fail_dominates_marginal(tc)
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);
            EEG = test_qualityGate.withEvents(EEG, 85);   % marginal for minTriggers=100
            gate = qualityGate(EEG, struct( ...
                'minTriggers', 100, ...
                'expectedChans', 64, ...   % fail (8 != 64)
                'marginalSlack', 0.8));
            tc.verifyEqual(gate.verdict, 'Fail');
            tc.verifyGreaterThanOrEqual(numel(gate.reasons), 2);
        end

        % -- batch mode --------------------------------------------------

        function batch_mode_returns_pending(tc)
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);
            gate = qualityGate(EEG, struct( ...
                'thresholdMode',      'batch', ...
                'maxOutlierTrialPct', 10));
            tc.verifyEqual(gate.verdict, 'Pending');
            tc.verifyEmpty(gate.reasons);
            tc.verifyTrue(~isnan(gate.metrics.pctOutlierTrials));
        end

        function batch_mode_collects_metrics_regardless_of_threshold(tc)
            % In batch mode the threshold is irrelevant - we still
            % collect every metric whose enabling toggle is on.
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);
            EEG = test_qualityGate.withEvents(EEG, 50);
            gate = qualityGate(EEG, struct( ...
                'thresholdMode', 'batch', ...
                'minTriggers', 1));   % any non-zero enables collection
            tc.verifyEqual(gate.metrics.nTriggers, 50);
            tc.verifyEqual(gate.verdict, 'Pending');
        end

        % -- ICA-based checks --------------------------------------------

        function emgFraction_uses_attached_TESA_classifier(tc)
            EEG = test_qualityGate.makeEEG(16, 1, 2000, 1000);
            EEG.icaweights = eye(4, 16);
            EEG.icasphere  = eye(16);
            EEG.icawinv    = randn(16, 4);
            EEG.icachansind = 1:16;
            EEG.icaCompClass.TESA1.compClass = [1 3 3 1];   % 2 of 4 are TMS Muscle
            EEG.icaCompClass.TESA1.compVars  = ones(1, 4);
            gate = qualityGate(EEG, struct('maxEMGFraction', 0.3));
            tc.verifyEqual(gate.metrics.emgFraction, 0.5);
            tc.verifyEqual(gate.verdict, 'Fail');
        end

        function electrodeCount_uses_ICLabel(tc)
            EEG = test_qualityGate.makeEEG(16, 1, 2000, 1000);
            EEG.icaweights = eye(5, 16);
            EEG.icasphere  = eye(16);
            EEG.icawinv    = randn(16, 5);
            EEG.icachansind = 1:16;
            % 3 of 5 components classified as Channel Noise (electrode artifact)
            probs = zeros(5, 7);
            probs(1, 1) = 1;   % Brain
            probs(2, 6) = 1;   % Channel Noise
            probs(3, 6) = 1;
            probs(4, 6) = 1;
            probs(5, 1) = 1;
            EEG.etc.ic_classification.ICLabel.classifications = probs;
            gate = qualityGate(EEG, struct('maxElectrodeCount', 2));
            tc.verifyEqual(gate.metrics.electrodeCount, 3);
            tc.verifyEqual(gate.verdict, 'Fail');
        end

        % -- gate.metrics shape ------------------------------------------

        % -- WarnAt overrides (Phase 4) ----------------------------------

        function maxWarnAt_overrides_slack(tc)
            % maxFlatChans = 10, slack = 0.8 -> warn cutoff at 8.
            % Override with WarnAt = 3 so anything > 3 (and <= 10) is
            % Marginal, anything > 10 is Fail.
            EEG = test_qualityGate.makeEEG(8, 30, 500, 1000);
            for k = 1:5
                EEG.data(k, :, :) = 0;   % 5 flat channels
            end
            gate = qualityGate(EEG, struct( ...
                'maxFlatChans',       10, ...
                'marginalSlack',      0.8, ...
                'maxFlatChansWarnAt', 3));
            tc.verifyEqual(gate.metrics.nFlatChans, 5);
            tc.verifyEqual(gate.verdict, 'Marginal');
        end

        function maxWarnAt_zero_falls_back_to_slack(tc)
            % Same fixture as above; WarnAt = 0 -> slack semantics:
            % nFlatChans = 5 < slack * 10 = 8 -> Pass.
            EEG = test_qualityGate.makeEEG(8, 30, 500, 1000);
            for k = 1:5
                EEG.data(k, :, :) = 0;
            end
            gate = qualityGate(EEG, struct( ...
                'maxFlatChans',       10, ...
                'marginalSlack',      0.8, ...
                'maxFlatChansWarnAt', 0));
            tc.verifyEqual(gate.verdict, 'Pass');
        end

        function minWarnAt_defines_marginal_upper_bound(tc)
            % WarnAt is the upper edge of the marginal band for a min
            % check: anything below threshold fails, anything between
            % threshold and WarnAt is Marginal, anything at or above
            % WarnAt passes.
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);

            % 40 < threshold(50) -> Fail
            EEG = test_qualityGate.withEvents(EEG, 40);
            gate = qualityGate(EEG, struct( ...
                'minTriggers',       50, ...
                'minTriggersWarnAt', 100));
            tc.verifyEqual(gate.verdict, 'Fail');

            % 70 in [50, 100) -> Marginal
            EEG.event = struct('type', {}, 'latency', {});
            EEG = test_qualityGate.withEvents(EEG, 70);
            gate = qualityGate(EEG, struct( ...
                'minTriggers',       50, ...
                'minTriggersWarnAt', 100));
            tc.verifyEqual(gate.verdict, 'Marginal');

            % 120 >= warnAt(100) -> Pass
            EEG.event = struct('type', {}, 'latency', {});
            EEG = test_qualityGate.withEvents(EEG, 120);
            gate = qualityGate(EEG, struct( ...
                'minTriggers',       50, ...
                'minTriggersWarnAt', 100));
            tc.verifyEqual(gate.verdict, 'Pass');
        end

        function minWarnAt_zero_means_no_marginal_band(tc)
            % Without an explicit WarnAt, min checks have no marginal
            % band: anything below threshold fails, anything else passes.
            % (Min metrics like rankRatio cap at 1.0, so a slack-derived
            % upper cutoff would falsely mark perfect data as Marginal.)
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);
            EEG = test_qualityGate.withEvents(EEG, 99);
            gate = qualityGate(EEG, struct( ...
                'minTriggers',   100, ...
                'marginalSlack', 0.8));
            tc.verifyEqual(gate.verdict, 'Fail');
        end

        % -- rejected-pct metrics (Phase 5) ------------------------------

        function maxRejectedChanPct_fails_above_threshold(tc)
            % Regression for the bug where maxBadChanPct never fired
            % even with 8/63 channels removed. The new maxRejectedChanPct
            % reads the running tally from the context, not statistical
            % outliers of the surviving channels.
            EEG = test_qualityGate.makeEEG(55, 10, 500, 1000); % 55 left now
            ctx = struct( ...
                'channels', struct('original', 63, 'nRejected', 8), ...
                'trials',   struct('original', 0,  'rejected',  0));
            gate = qualityGate(EEG, struct('maxRejectedChanPct', 10), ctx);
            % 8/63 = 12.7% > 10 -> Fail.
            tc.verifyEqual(gate.verdict, 'Fail');
            tc.verifyTrue(any(contains(gate.reasons, '% rejected channels')));
        end

        function maxRejectedTrialPct_fails_above_threshold(tc)
            EEG = test_qualityGate.makeEEG(8, 59, 500, 1000);
            ctx = struct( ...
                'channels', struct('original', 8,  'nRejected', 0), ...
                'trials',   struct('original', 80, 'rejected', 21));
            gate = qualityGate(EEG, struct('maxRejectedTrialPct', 15), ctx);
            % 21/80 = 26.25% > 15 -> Fail.
            tc.verifyEqual(gate.verdict, 'Fail');
            tc.verifyTrue(any(contains(gate.reasons, '% rejected trials')));
        end

        function maxRejectedChanPct_warnAt_yields_marginal(tc)
            % 8/63 = 12.7%. Threshold 15 (Fail above), WarnAt 10 (Marginal
            % above) -> Marginal.
            EEG = test_qualityGate.makeEEG(55, 10, 500, 1000);
            ctx = struct( ...
                'channels', struct('original', 63, 'nRejected', 8), ...
                'trials',   struct('original', 0,  'rejected',  0));
            gate = qualityGate(EEG, struct( ...
                'maxRejectedChanPct',       15, ...
                'maxRejectedChanPctWarnAt', 10), ctx);
            tc.verifyEqual(gate.verdict, 'Marginal');
        end

        function maxRejected_passes_when_under_threshold(tc)
            EEG = test_qualityGate.makeEEG(60, 10, 500, 1000);
            ctx = struct( ...
                'channels', struct('original', 63, 'nRejected', 3), ...
                'trials',   struct('original', 80, 'rejected',  4));
            gate = qualityGate(EEG, struct( ...
                'maxRejectedChanPct',  10, ...
                'maxRejectedTrialPct', 15), ctx);
            % 3/63 = 4.8%, 4/80 = 5% - both under their thresholds.
            tc.verifyEqual(gate.verdict, 'Pass');
        end

        function maxRejected_NaN_when_context_missing(tc)
            % No context -> metric is NaN -> check skipped silently.
            EEG = test_qualityGate.makeEEG(8, 10, 500, 1000);
            gate = qualityGate(EEG, struct('maxRejectedChanPct', 10));
            tc.verifyEqual(gate.verdict, 'Pass');
            tc.verifyTrue(isnan(gate.metrics.rejectedChanPct));
        end

        % -- deprecated alias behavior ----------------------------------

        function maxBadChanPct_alias_maps_to_maxOutlier(tc)
            % Saved pipelines using the old key must still work. With one
            % flat channel and seven normal ones, pctOutlierChans is 0%
            % (the flat channel is excluded), so a permissive 0.5%
            % outlier threshold still passes the gate either way - what
            % we are testing here is that the *threshold* gets carried
            % across into the new name.
            EEG = test_qualityGate.makeEEG(8, 30, 500, 1000);
            gate = qualityGate(EEG, struct('maxBadChanPct', 25));
            tc.verifyEqual(gate.thresholds.maxOutlierChanPct, 25);
            % The deprecated field is NOT mirrored in thresholds (it's
            % been folded into the new one).
            tc.verifyFalse(isfield(gate.thresholds, 'maxBadChanPct'));
        end

        function maxBadTrialPct_alias_maps_to_maxOutlier(tc)
            EEG = test_qualityGate.makeEEG(8, 30, 500, 1000);
            gate = qualityGate(EEG, struct('maxBadTrialPct', 12));
            tc.verifyEqual(gate.thresholds.maxOutlierTrialPct, 12);
        end

        function metrics_always_has_cheap_fields(tc)
            EEG = test_qualityGate.makeEEG(8, 5, 500, 1000);
            gate = qualityGate(EEG, struct());
            tc.verifyTrue(isfield(gate.metrics, 'nbchan'));
            tc.verifyTrue(isfield(gate.metrics, 'srate'));
            tc.verifyTrue(isfield(gate.metrics, 'nTriggers'));
            tc.verifyTrue(isfield(gate.metrics, 'nTrials'));
            tc.verifyTrue(isfield(gate.metrics, 'rankRatio'));
        end
    end
end
