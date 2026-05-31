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

        function tmsWindow_is_respected(tc)
            % Synthetic EEG with a large transient at t = 30..45 ms.
            % With tmsWindow=[0 25] the transient is INCLUDED in scoring
            % (score is high). With tmsWindow=[0 50] the transient is
            % EXCLUDED (score is much lower). Proves the parameter is
            % being applied to the scoring window.
            srate = 1000;
            nPnts = 500;
            EEG.data   = 0.1 * randn(4, nPnts, 10);
            EEG.times  = (0:nPnts-1) / srate * 1000 - 100;  % -100..+399 ms
            EEG.srate  = srate;
            EEG.nbchan = 4;
            EEG.trials = 10;
            EEG.pnts   = nPnts;

            % Inject a 30..45 ms transient on every trial, channel 1.
            mask = EEG.times >= 30 & EEG.times <= 45;
            EEG.data(1, mask, :) = 50;   % 50 uV bump

            sumNarrow = computeAttributeMatrix(EEG, struct('tmsWindow', [0 25]));
            sumWide   = computeAttributeMatrix(EEG, struct('tmsWindow', [0 50]));

            % Channel 1's median should be higher when the bump is in-window.
            tc.verifyGreaterThan(median(sumNarrow(1,:), 'omitnan'), ...
                                 median(sumWide(1,:),   'omitnan'), ...
                'Narrow tmsWindow should leave the 30-45 ms transient in the scoring window');
        end

        function default_tmsWindow_is_25(tc)
            % Same setup as above but verifies the default behaves like
            % the narrow case (i.e., default is now [0 25], not [0 50]).
            srate = 1000;
            nPnts = 500;
            EEG.data   = 0.1 * randn(4, nPnts, 10);
            EEG.times  = (0:nPnts-1) / srate * 1000 - 100;
            EEG.srate  = srate;
            EEG.nbchan = 4;
            EEG.trials = 10;
            EEG.pnts   = nPnts;
            mask = EEG.times >= 30 & EEG.times <= 45;
            EEG.data(1, mask, :) = 50;

            sumDefault = computeAttributeMatrix(EEG);                       % default
            sumWide    = computeAttributeMatrix(EEG, struct('tmsWindow', [0 50]));

            tc.verifyGreaterThan(median(sumDefault(1,:), 'omitnan'), ...
                                 median(sumWide(1,:),    'omitnan'), ...
                'Default tmsWindow should be the narrower [0 25] ms');
        end
    end
end
