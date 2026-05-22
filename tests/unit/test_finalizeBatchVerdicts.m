classdef test_finalizeBatchVerdicts < matlab.unittest.TestCase
% TEST_FINALIZEBATCHVERDICTS  Unit tests for src/qa/finalizeBatchVerdicts.m

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            tc.addTeardown(@() rmpath(genpath(fullfile(root, 'src'))));
        end
    end

    methods (Static, Access = private)
        function r = makeReport(quality)
            r = struct('inputFile', 'fake.set', ...
                'processedAt', datetime('now'), ...
                'quality', quality);
        end

        function g = pendingGate(label, stepIdx, metrics, sigmas, slack)
            if nargin < 4, sigmas = 3; end
            if nargin < 5, slack  = 0.8; end
            % maxBadTrialPct = 1 opts the pctBadTrials field into the
            % batch comparison (value is irrelevant in batch mode).
            g = struct( ...
                'label',      label, ...
                'mode',       'batch', ...
                'verdict',    'Pending', ...
                'reasons',    {{}}, ...
                'metrics',    metrics, ...
                'thresholds', struct( ...
                    'outlierSigmas', sigmas, 'marginalSlack', slack, ...
                    'maxBadTrialPct', 1), ...
                'stepIndex',  stepIdx);
        end
    end

    methods (Test)
        function single_gate_outlier_is_flagged(tc)
            % Five files with pctBadTrials = [1 2 3 4 50]. Outlier at 50.
            vals = [1 2 3 4 50];
            reports = cell(1, 5);
            for k = 1:5
                m = struct('pctBadTrials', vals(k), 'nFlatChans', NaN, ...
                    'nSatChans', NaN, 'pctBadChans', NaN, ...
                    'emgFraction', NaN, 'electrodeCount', NaN);
                q = struct('figures', {{}}, ...
                    'gates', {{ test_finalizeBatchVerdicts.pendingGate('g1', 5, m) }}, ...
                    'worstVerdict', 'Pending');
                reports{k} = test_finalizeBatchVerdicts.makeReport(q);
            end
            reports = finalizeBatchVerdicts(reports);
            verdicts = cellfun(@(r) r.quality.gates{1}.verdict, reports, ...
                'UniformOutput', false);
            tc.verifyEqual(verdicts{5}, 'Fail');
            tc.verifyTrue(all(strcmp(verdicts(1:4), 'Pass')));
        end

        function tight_batch_all_pass(tc)
            % Spread is small relative to the values themselves: median=5,
            % MAD=2, cutoff at sig=3 is 5+3*1.4826*2 ~= 13.9; slack*cutoff
            % ~= 11.1; max value 7 < 11.1 -> all Pass.
            vals = [3 4 5 6 7];
            reports = makeBatch(vals);
            reports = finalizeBatchVerdicts(reports);
            verdicts = cellfun(@(r) r.quality.gates{1}.verdict, reports, ...
                'UniformOutput', false);
            tc.verifyTrue(all(strcmp(verdicts, 'Pass')), ...
                sprintf('Expected all Pass, got: %s', strjoin(verdicts, ', ')));
        end

        function outlierSigmas_controls_strictness(tc)
            % vals=[1 2 3 4 5] -> median=3, MAD=1.
            % sig=20: cutoff = 3 + 20*1.4826 = 32.6; slack*cutoff = 26 -> 5 is Pass
            % sig=1:  cutoff = 3 + 1.4826    = 4.48; 5 > 4.48 -> Fail.
            vals = [1 2 3 4 5];
            reports1 = makeBatch(vals, 'sig', 20);
            reports1 = finalizeBatchVerdicts(reports1);
            reports2 = makeBatch(vals, 'sig', 1);
            reports2 = finalizeBatchVerdicts(reports2);

            tc.verifyEqual(reports1{5}.quality.gates{1}.verdict, 'Pass');
            tc.verifyEqual(reports2{5}.quality.gates{1}.verdict, 'Fail');
        end

        function absolute_mode_gates_untouched(tc)
            m = struct('pctBadTrials', 5, 'nFlatChans', NaN, ...
                'nSatChans', NaN, 'pctBadChans', NaN, ...
                'emgFraction', NaN, 'electrodeCount', NaN);
            absGate = struct( ...
                'label', 'abs', 'mode', 'absolute', ...
                'verdict', 'Pass', 'reasons', {{}}, ...
                'metrics', m, 'thresholds', struct(), 'stepIndex', 3);
            q = struct('figures', {{}}, ...
                'gates', {{absGate}}, 'worstVerdict', 'Pass');
            reports = { test_finalizeBatchVerdicts.makeReport(q) };
            reports = finalizeBatchVerdicts(reports);
            tc.verifyEqual(reports{1}.quality.gates{1}.verdict, 'Pass');
            tc.verifyEqual(reports{1}.quality.gates{1}.mode, 'absolute');
        end

        function worstVerdict_recomputed_after_resolve(tc)
            vals = [1 2 3 4 100];
            reports = makeBatch(vals);
            reports = finalizeBatchVerdicts(reports);
            tc.verifyEqual(reports{5}.quality.worstVerdict, 'Fail');
            tc.verifyEqual(reports{1}.quality.worstVerdict, 'Pass');
        end
    end
end

function reports = makeBatch(vals, varargin)
% Build N reports with one pending gate each carrying pctBadTrials=vals(k).
p = struct('sig', 3, 'slack', 0.8);
for i = 1:2:numel(varargin)
    p.(varargin{i}) = varargin{i+1};
end
reports = cell(1, numel(vals));
for k = 1:numel(vals)
    m = struct('pctBadTrials', vals(k), 'nFlatChans', NaN, ...
        'nSatChans', NaN, 'pctBadChans', NaN, ...
        'emgFraction', NaN, 'electrodeCount', NaN);
    g = struct( ...
        'label',      'g1', ...
        'mode',       'batch', ...
        'verdict',    'Pending', ...
        'reasons',    {{}}, ...
        'metrics',    m, ...
        'thresholds', struct( ...
            'outlierSigmas', p.sig, 'marginalSlack', p.slack, ...
            'maxBadTrialPct', 1), ...
        'stepIndex',  5);
    q = struct('figures', {{}}, 'gates', {{g}}, 'worstVerdict', 'Pending');
    reports{k} = struct('inputFile', 'fake.set', ...
        'processedAt', datetime('now'), 'quality', q);
end
end
