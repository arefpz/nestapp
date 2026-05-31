classdef test_renderQualityFigure < matlab.unittest.TestCase
% TEST_RENDERQUALITYFIGURE  Smoke test for src/qa/renderQualityFigure.m
%   Asserts the renderer produces a PNG without error for both an epoched
%   ICA'd EEG and a minimal continuous EEG. The actual pixel contents are
%   not verified - the goal is to confirm the function path is exercised
%   end to end without unhandled failures.

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            tc.addTeardown(@() rmpath(genpath(fullfile(root, 'src'))));
        end
    end

    methods (Static, Access = private)
        function EEG = makeEpochedEEG(nbchan, nTrials, nPnts, srate)
            EEG.data    = randn(nbchan, nPnts, nTrials);
            EEG.times   = (0:nPnts-1) / srate * 1000 - 200;
            EEG.srate   = srate;
            EEG.nbchan  = nbchan;
            EEG.trials  = nTrials;
            EEG.pnts    = nPnts;
            EEG.chanlocs = struct('X',{},'Y',{},'Z',{},'labels',{});
        end
    end

    methods (Test)
        function epoched_no_ICA_renders_PNG(tc)
            EEG = test_renderQualityFigure.makeEpochedEEG(16, 30, 600, 1000);
            outPath = fullfile(tempdir, ['qc_test_', char(matlab.lang.internal.uuid()), '.png']);
            tc.addTeardown(@() safeDelete(outPath));
            renderQualityFigure(EEG, outPath, struct());
            tc.verifyTrue(exist(outPath, 'file') == 2);
            d = dir(outPath);
            tc.verifyGreaterThan(d.bytes, 10000, 'PNG should be > 10 KB');
        end

        function continuous_no_ICA_renders_PNG(tc)
            EEG.data    = randn(8, 2000);
            EEG.times   = (0:1999) / 1000 * 1000;
            EEG.srate   = 1000;
            EEG.nbchan  = 8;
            EEG.trials  = 1;
            EEG.pnts    = 2000;
            EEG.chanlocs = struct('X',{},'Y',{},'Z',{},'labels',{});
            outPath = fullfile(tempdir, ['qc_cont_', char(matlab.lang.internal.uuid()), '.png']);
            tc.addTeardown(@() safeDelete(outPath));
            renderQualityFigure(EEG, outPath, struct());
            tc.verifyTrue(exist(outPath, 'file') == 2);
            d = dir(outPath);
            tc.verifyGreaterThan(d.bytes, 5000);
        end

        function creates_missing_parent_dir(tc)
            EEG = test_renderQualityFigure.makeEpochedEEG(8, 10, 500, 1000);
            nestedDir = fullfile(tempdir, ['qc_nested_', char(matlab.lang.internal.uuid())]);
            outPath   = fullfile(nestedDir, 'foo.png');
            tc.addTeardown(@() safeDeleteDir(nestedDir));
            renderQualityFigure(EEG, outPath, struct());
            tc.verifyTrue(exist(outPath, 'file') == 2);
        end
    end
end

function safeDelete(p)
if exist(p, 'file'), delete(p); end
end

function safeDeleteDir(d)
if exist(d, 'dir'), rmdir(d, 's'); end
end
