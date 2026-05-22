classdef test_computeAttributeMatrix < matlab.unittest.TestCase
% TEST_COMPUTEATTRIBUTEMATRIX  Unit tests for src/qa/computeAttributeMatrix.m

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
            EEG.times  = (0:nPnts-1) / srate * 1000 - 200;  % -200..+800 ms range
            EEG.srate  = srate;
            EEG.nbchan = nbchan;
            EEG.trials = nTrials;
            EEG.pnts   = nPnts;
        end
    end

    methods (Test)
        function default_shape_matches_data(tc)
            EEG = test_computeAttributeMatrix.makeEEG(16, 40, 500, 1000);
            [SM, summary] = computeAttributeMatrix(EEG);
            tc.verifySize(SM, [16 40]);
            tc.verifyEqual(summary.nbchan,  16);
            tc.verifyEqual(summary.nTrials, 40);
            tc.verifyEqual(summary.attribute, 'minmax_no_tms');
        end

        function inflated_channel_has_top_median(tc)
            EEG = test_computeAttributeMatrix.makeEEG(16, 30, 500, 1000);
            EEG.data(5,:,:) = EEG.data(5,:,:) * 5;
            [~, summary] = computeAttributeMatrix(EEG);
            [~, worstChan] = max(summary.perChanMedian);
            tc.verifyEqual(worstChan, 5, ...
                'Channel 5 was inflated 5x and should top the per-channel median.');
        end

        function inflated_trial_has_top_median(tc)
            EEG = test_computeAttributeMatrix.makeEEG(16, 30, 500, 1000);
            EEG.data(:,:,7) = EEG.data(:,:,7) * 5;
            [~, summary] = computeAttributeMatrix(EEG);
            [~, worstTrial] = max(summary.perTrialMedian);
            tc.verifyEqual(worstTrial, 7);
        end

        function flat_channel_returns_NaN_row(tc)
            EEG = test_computeAttributeMatrix.makeEEG(8, 20, 500, 1000);
            EEG.data(3,:,:) = 0;
            [SM, summary] = computeAttributeMatrix(EEG);
            tc.verifyTrue(summary.flatChanMask(3));
            tc.verifyTrue(all(isnan(SM(3,:))));
            tc.verifyFalse(any(summary.flatChanMask([1 2 4:8])));
        end

        function saturated_channel_is_flagged(tc)
            EEG = test_computeAttributeMatrix.makeEEG(8, 20, 500, 1000);
            EEG.data(2,:,:) = 500;   % > 250 uV threshold
            [SM, summary] = computeAttributeMatrix(EEG);
            tc.verifyTrue(summary.satChanMask(2));
            tc.verifyTrue(all(isnan(SM(2,:))));
        end

        function all_three_modes_return_same_shape(tc)
            EEG = test_computeAttributeMatrix.makeEEG(8, 20, 500, 1000);
            modes = {'minmax', 'minmax_no_tms', 'highfreq'};
            for k = 1:numel(modes)
                SM = computeAttributeMatrix(EEG, struct('attribute', modes{k}));
                tc.verifySize(SM, [8 20], ...
                    sprintf('Mode %s produced wrong shape', modes{k}));
            end
        end

        function invalid_attribute_errors_cleanly(tc)
            EEG = test_computeAttributeMatrix.makeEEG(4, 4, 200, 1000);
            tc.verifyError(@() computeAttributeMatrix(EEG, struct('attribute','bogus')), ...
                'computeAttributeMatrix:badAttribute');
        end

        function continuous_data_treated_as_one_trial(tc)
            EEG = test_computeAttributeMatrix.makeEEG(8, 0, 1000, 1000);
            EEG.data = randn(8, 1000);   % 2D continuous
            [SM, summary] = computeAttributeMatrix(EEG);
            tc.verifySize(SM, [8 1]);
            tc.verifyEqual(summary.nTrials, 1);
        end
    end
end
