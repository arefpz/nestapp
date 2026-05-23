classdef test_aggregateMetricDistributions < matlab.unittest.TestCase
% TEST_AGGREGATEMETRICDISTRIBUTIONS  Unit tests for src/qa/aggregateMetricDistributions.m

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            tc.addTeardown(@() rmpath(genpath(fullfile(root, 'src'))));
        end
    end

    methods (Static, Access = private)
        function r = reportWithGates(gates)
            r = struct( ...
                'inputFile',   'fake.set', ...
                'processedAt', datetime('now'), ...
                'quality',     struct('figures', {{}}, ...
                                      'gates', {gates}, ...
                                      'worstVerdict', 'Pass'));
        end

        function g = gate(mode, thresholds, metrics, varargin)
            g = struct( ...
                'label',      'g', ...
                'mode',       mode, ...
                'verdict',    'Pass', ...
                'reasons',    {{}}, ...
                'metrics',    metrics, ...
                'thresholds', thresholds, ...
                'stepIndex',  1);
            for k = 1:2:numel(varargin)
                g.(varargin{k}) = varargin{k+1};
            end
        end
    end

    methods (Test)
        function disabled_metrics_absent_from_output(tc)
            mk = @test_aggregateMetricDistributions.reportWithGates;
            mg = @test_aggregateMetricDistributions.gate;
            % All thresholds disabled (= 0) -> nothing to plot.
            r = mk({mg('absolute', struct('maxFlatChans', 0), ...
                struct('nFlatChans', 3))});
            d = aggregateMetricDistributions({r});
            tc.verifyEmpty(d);
        end

        function single_enabled_metric_appears_once(tc)
            mk = @test_aggregateMetricDistributions.reportWithGates;
            mg = @test_aggregateMetricDistributions.gate;
            r = mk({mg('absolute', struct('maxFlatChans', 5), ...
                struct('nFlatChans', 2))});
            d = aggregateMetricDistributions({r});
            tc.verifyLength(d, 1);
            tc.verifyEqual(d(1).name, 'nFlatChans');
            tc.verifyEqual(d(1).values, 2);
            tc.verifyEqual(d(1).absThresholds, 5);
            tc.verifyEmpty(d(1).batchCutoffs);
            tc.verifyEqual(d(1).mode, 'absolute');
        end

        function mixed_modes_record_both_thresholds(tc)
            mk = @test_aggregateMetricDistributions.reportWithGates;
            mg = @test_aggregateMetricDistributions.gate;

            absGate = mg('absolute', struct('maxBadTrialPct', 10), ...
                struct('pctBadTrials', 3));
            batchGate = mg('batch', struct('maxBadTrialPct', 1), ...
                struct('pctBadTrials', 8), ...
                'batchCutoffs', struct('pctBadTrials', 12));

            d = aggregateMetricDistributions({mk({absGate, batchGate})});
            tc.verifyLength(d, 1);
            tc.verifyEqual(d.name, 'pctBadTrials');
            tc.verifyEqual(sort(d.values), [3 8]);
            tc.verifyEqual(d.absThresholds, 10);
            tc.verifyEqual(d.batchCutoffs, 12);
            tc.verifyEqual(d.mode, 'mixed');
        end

        function NaN_metric_values_omitted(tc)
            mk = @test_aggregateMetricDistributions.reportWithGates;
            mg = @test_aggregateMetricDistributions.gate;
            r1 = mk({mg('absolute', struct('maxEMGFraction', 0.3), ...
                struct('emgFraction', NaN))});
            r2 = mk({mg('absolute', struct('maxEMGFraction', 0.3), ...
                struct('emgFraction', 0.15))});
            d = aggregateMetricDistributions({r1, r2});
            tc.verifyLength(d, 1);
            tc.verifyEqual(d.values, 0.15);   % only the non-NaN entry
        end

        function min_style_metric_picked_up_via_minTriggers(tc)
            mk = @test_aggregateMetricDistributions.reportWithGates;
            mg = @test_aggregateMetricDistributions.gate;
            r = mk({mg('absolute', struct('minTriggers', 100), ...
                struct('nTriggers', 85))});
            d = aggregateMetricDistributions({r});
            tc.verifyEqual(d.name, 'nTriggers');
            tc.verifyEqual(d.values, 85);
            tc.verifyEqual(d.absThresholds, 100);
        end
    end
end
