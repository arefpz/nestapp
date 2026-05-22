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
            EEG = test_computeICAQualityMetrics.makeEEGWithICA(32, 4, 2000, 1000);
            % Make component 2 a "fried electrode": one channel huge,
            % others tiny -> high topo kurtosis.
            EEG.icawinv(:, 2) = 0.01 * randn(32, 1);
            EEG.icawinv(10, 2) = 50;
            % Keep the activation tame so EMG doesn't dominate.
            EEG.data(:, :) = randn(32, 2000) * 0.5;
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
    end
end
