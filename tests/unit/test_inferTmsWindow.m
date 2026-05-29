
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef test_inferTmsWindow < matlab.unittest.TestCase
% TEST_INFERTMSWINDOW  Unit tests for src/qa/inferTmsWindow.m

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            tc.addTeardown(@() rmpath(genpath(fullfile(root, 'src'))));
        end
    end

    methods (Static, Access = private)
        function EEG = makeEpoched(srate, nPnts, nTrials)
            % EEG.times spans -200..(nPnts/srate*1000 - 200) ms; t=0 lands
            % at sample (200/1000)*srate + 1. With srate=1000, nPnts=500,
            % t=0 is at sample 201.
            EEG.srate  = srate;
            EEG.pnts   = nPnts;
            EEG.trials = nTrials;
            EEG.nbchan = 4;
            EEG.times  = (0:nPnts-1) / srate * 1000 - 200;
            EEG.data   = randn(4, nPnts, nTrials);
            EEG.event  = struct('type', {}, 'latency', {});
        end
    end

    methods (Test)
        function fallback_when_no_events(tc)
            EEG = test_inferTmsWindow.makeEpoched(1000, 500, 10);
            win = inferTmsWindow(EEG, [0 50]);
            tc.verifyEqual(win, [0 50]);
        end

        function fallback_when_no_matching_label(tc)
            EEG = test_inferTmsWindow.makeEpoched(1000, 500, 10);
            EEG.event(1).type    = 'OtherEvent';
            EEG.event(1).latency = 201;
            win = inferTmsWindow(EEG, [5 25]);
            tc.verifyEqual(win, [5 25]);
        end

        function epoched_event_at_t0(tc)
            EEG = test_inferTmsWindow.makeEpoched(1000, 500, 10);
            EEG.event(1).type    = 'TMS';
            EEG.event(1).latency = 201;   % corresponds to t = 0 ms
            win = inferTmsWindow(EEG, [99 99]);
            tc.verifyEqual(win, [0 25]);
        end

        function epoched_event_offset_picks_event_time(tc)
            EEG = test_inferTmsWindow.makeEpoched(1000, 500, 10);
            EEG.event(1).type    = 'TMS';
            EEG.event(1).latency = 206;   % 5 samples after t = 0 (= +5 ms)
            win = inferTmsWindow(EEG, [99 99]);
            tc.verifyEqual(win, [5 5 + 25]);
        end

        function custom_labels_override_defaults(tc)
            EEG = test_inferTmsWindow.makeEpoched(1000, 500, 10);
            EEG.event(1).type    = 'StimOnset';
            EEG.event(1).latency = 211;   % t = 10 ms
            % Default label set doesn't include 'StimOnset' -> fallback.
            tc.verifyEqual(inferTmsWindow(EEG, [99 99]), [99 99]);
            % With custom label list it matches.
            win = inferTmsWindow(EEG, [99 99], struct('labels', {{'StimOnset'}}));
            tc.verifyEqual(win, [10 35]);
        end

        function decayMs_param_respected(tc)
            EEG = test_inferTmsWindow.makeEpoched(1000, 500, 10);
            EEG.event(1).type    = 'TMS';
            EEG.event(1).latency = 201;
            win = inferTmsWindow(EEG, [0 0], struct('decayMs', 50));
            tc.verifyEqual(win, [0 50]);
        end

        function numeric_event_type_converted(tc)
            EEG = test_inferTmsWindow.makeEpoched(1000, 500, 10);
            EEG.event(1).type    = 128;
            EEG.event(1).latency = 201;
            % '128' is not in default labels but the conversion still works
            % when caller supplies a custom label set.
            win = inferTmsWindow(EEG, [99 99], struct('labels', {{'128'}}));
            tc.verifyEqual(win, [0 25]);
        end
    end
end
