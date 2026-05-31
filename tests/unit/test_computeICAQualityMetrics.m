
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef test_computeICAQualityMetrics < matlab.unittest.TestCase
% TEST_COMPUTEICAQUALITYMETRICS  Unit tests for src/qa/computeICAQualityMetrics.m

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            tc.addTeardown(@() rmpath(genpath(fullfile(root, 'src'))));
        end
    end

    methods (Static, Access = private)
        function EEG = makeEEGWithICA(nbchan, nComp, nPnts, srate)
            % Build a minimal EEG with usable ICA fields. The caller is
            % expected to overwrite icawinv / data columns to shape specific
            % topographies and activations for individual tests.
            EEG.srate      = srate;
            EEG.nbchan     = nbchan;
            EEG.pnts       = nPnts;
            EEG.trials     = 1;
            EEG.times      = (0:nPnts-1) / srate * 1000;
            EEG.data       = randn(nbchan, nPnts);
            EEG.icaweights = eye(nComp, nbchan);
            EEG.icasphere  = eye(nbchan);
            EEG.icawinv    = randn(nbchan, nComp);
            EEG.icachansind = 1:nbchan;

            % Cart coords are needed for the frontal mask; place a few at the
            % "front" so the EOG path can fire when topographies favor them.
            EEG.chanlocs = struct('X', {}, 'Y', {}, 'Z', {}, 'labels', {});
            for c = 1:nbchan
                EEG.chanlocs(c).X = 1;
                EEG.chanlocs(c).Y = (c - nbchan/2) * 0.1;  % spread left-right
                EEG.chanlocs(c).Z = 0;
                EEG.chanlocs(c).labels = sprintf('Ch%d', c);
            end
        end
    end

    methods (Test)
        function empty_when_no_ICA(tc)
            EEG = struct('data', randn(8, 1000), 'srate', 1000, 'nbchan', 8, ...
                'pnts', 1000, 'trials', 1, 'times', (0:999), ...
                'chanlocs', struct());
            metrics = computeICAQualityMetrics(EEG);
            tc.verifyEmpty(metrics);
        end

        function spike_topography_classified_as_electrode(tc)
            % Seed rng so the synthetic noise doesn't accidentally produce
            % a high-frequency activation that triggers EMG before
            % Electrode (precedence: EMG > Electrode in the classifier).
            rng(42);
            EEG = test_computeICAQualityMetrics.makeEEGWithICA(32, 4, 2000, 1000);
            % Make component 2 a "fried electrode": one channel huge,
            % others tiny -> high topo kurtosis.
            EEG.icawinv(:, 2) = 0.01 * randn(32, 1);
            EEG.icawinv(10, 2) = 50;
            % Keep the activation tame and low-frequency so EMG ratio
            % stays below 1.1.
            t = (0:1999) / 1000;
            EEG.data = repmat(0.5 * sin(2*pi*2*t), 32, 1) ...
                     + 0.05 * randn(32, 2000);
            metrics = computeICAQualityMetrics(EEG);
            tc.verifyGreaterThan(metrics(2).kurtosisTopo, 15);
            tc.verifyEqual(metrics(2).classification, 'Electrode');
        end

        function high_frequency_component_classified_as_EMG(tc)
            EEG = test_computeICAQualityMetrics.makeEEGWithICA(16, 4, 4000, 1000);
            % Inject a strong 40 Hz component into channel data so the
            % default-identity icaweights produces a high-freq activation
            % at component 1.
            t = (0:3999) / 1000;
            EEG.data(1, :) = 5 * sin(2*pi*40*t) + 0.1*randn(1, 4000);
            % Make every other channel quiet.
            EEG.data(2:end, :) = 0.05 * randn(15, 4000);
            % Keep topographies smooth so kurtosis stays low.
            EEG.icawinv = randn(16, 4) * 0.5;
            metrics = computeICAQualityMetrics(EEG);
            tc.verifyGreaterThan(metrics(1).emgRatio, 1.1);
            tc.verifyEqual(metrics(1).classification, 'EMG');
        end

        function smooth_lowfreq_component_classified_as_OK(tc)
            EEG = test_computeICAQualityMetrics.makeEEGWithICA(16, 4, 4000, 1000);
            t = (0:3999) / 1000;
            % 2 Hz dominant; nothing in 30-50 Hz band.
            EEG.data(1, :) = 5 * sin(2*pi*2*t);
            EEG.data(2:end, :) = 0.05 * randn(15, 4000);
            % Smooth, peripherally-concentrated topography so neither
            % kurtosis (Electrode) nor frontal CBI (EOG) fires.
            EEG.icawinv = zeros(16, 4);
            EEG.icawinv(:, 1) = sin((1:16)' * 0.4);   % distributed, low kurtosis
            EEG.icawinv([1 2 15 16], 1) = 2;          % weight peripheral channels
            metrics = computeICAQualityMetrics(EEG);
            tc.verifyLessThan(metrics(1).kurtosisTopo, 15);
            tc.verifyLessThan(metrics(1).emgRatio,    1.1);
            tc.verifyLessThan(metrics(1).eogScore,    0.5);
            tc.verifyEqual(metrics(1).classification, 'OK');
        end

        function metrics_length_matches_components(tc)
            EEG = test_computeICAQualityMetrics.makeEEGWithICA(16, 7, 2000, 1000);
            metrics = computeICAQualityMetrics(EEG);
            tc.verifyLength(metrics, 7);
            tc.verifyEqual([metrics.compIdx], 1:7);
        end

        function heuristic_source_when_no_classifier_attached(tc)
            EEG = test_computeICAQualityMetrics.makeEEGWithICA(8, 3, 2000, 1000);
            metrics = computeICAQualityMetrics(EEG);
            tc.verifyEqual(metrics(1).source, 'Heuristic');
        end

        function TESA_classifications_take_precedence(tc)
            % When EEG.icaCompClass is attached the heuristic path must
            % NOT run; labels must come from the TESA compClass codes.
            EEG = test_computeICAQualityMetrics.makeEEGWithICA(16, 5, 2000, 1000);
            EEG.icaCompClass.TESA1.compClass = [1 3 1 6 1];   % Keep, TMS Muscle, Keep, Muscle, Keep
            EEG.icaCompClass.TESA1.compVars  = ones(1, 5);

            metrics = computeICAQualityMetrics(EEG);

            tc.verifyEqual({metrics.source}, repmat({'TESA'}, 1, 5));
            tc.verifyEqual({metrics.classification}, ...
                {'Keep','TMS Muscle','Keep','Muscle','Keep'});
            tc.verifyEqual([metrics.kept], [true false true false true]);
            % Heuristic fields should be NaN / defaults (not computed).
            tc.verifyTrue(all(isnan([metrics.emgRatio])));
        end

        function TESA_uses_latest_round_when_multiple_present(tc)
            % If both TESA1 and TESA2 exist, the latest round's
            % classification should win (matches processOneFile.m logic).
            EEG = test_computeICAQualityMetrics.makeEEGWithICA(16, 3, 2000, 1000);
            EEG.icaCompClass.TESA1.compClass = [1 1 1];
            EEG.icaCompClass.TESA1.compVars  = ones(1, 3);
            EEG.icaCompClass.TESA2.compClass = [4 1 6];       % Blink, Keep, Muscle
            EEG.icaCompClass.TESA2.compVars  = ones(1, 3);

            metrics = computeICAQualityMetrics(EEG);
            tc.verifyEqual({metrics.classification}, {'Blink','Keep','Muscle'});
        end

        function ICLabel_classifications_used_when_present(tc)
            % EEG.etc.ic_classification.ICLabel.classifications is a
            % nComp x 7 probability matrix - argmax per row picks the
            % display label. Brain -> kept; everything else -> rejected.
            EEG = test_computeICAQualityMetrics.makeEEGWithICA(16, 4, 2000, 1000);
            probs = zeros(4, 7);
            probs(1, 1) = 0.9;   % Brain
            probs(2, 2) = 0.8;   % Muscle
            probs(3, 3) = 0.7;   % Eye
            probs(4, 5) = 0.6;   % Line Noise
            EEG.etc.ic_classification.ICLabel.classifications = probs;

            metrics = computeICAQualityMetrics(EEG);

            tc.verifyEqual({metrics.source}, repmat({'ICLabel'}, 1, 4));
            tc.verifyEqual({metrics.classification}, {'Brain','Muscle','Eye','Line Noise'});
            tc.verifyEqual([metrics.kept], [true false false false]);
        end

        function TESA_wins_over_ICLabel_when_both_present(tc)
            EEG = test_computeICAQualityMetrics.makeEEGWithICA(16, 2, 2000, 1000);
            EEG.icaCompClass.TESA1.compClass = [4 1];
            EEG.icaCompClass.TESA1.compVars  = ones(1, 2);
            EEG.etc.ic_classification.ICLabel.classifications = ...
                [1 0 0 0 0 0 0; 0 1 0 0 0 0 0];

            metrics = computeICAQualityMetrics(EEG);
            tc.verifyEqual({metrics.source}, {'TESA', 'TESA'});
        end
    end
end
