classdef test_aggregateGateVerdicts < matlab.unittest.TestCase
% TEST_AGGREGATEGATEVERDICTS  Unit tests for src/qa/aggregateGateVerdicts.m

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            tc.addTeardown(@() rmpath(genpath(fullfile(root, 'src'))));
        end
    end

    methods (Static, Access = private)
        function r = makeReport(name, gates)
            r = struct( ...
                'inputFile',   [name '.set'], ...
                'processedAt', datetime('now'), ...
                'quality',     struct('figures', {{}}, ...
                                      'gates', {gates}, ...
                                      'worstVerdict', 'Pass'));
        end

        function g = makeGate(label, verdict, stepIdx)
            if nargin < 3, stepIdx = 1; end
            g = struct( ...
                'label',      label, ...
                'mode',       'absolute', ...
                'verdict',    verdict, ...
                'reasons',    {{}}, ...
                'metrics',    struct(), ...
                'thresholds', struct(), ...
                'stepIndex',  stepIdx);
        end
    end

    methods (Test)
        function empty_reports_returns_empty(tc)
            s = aggregateGateVerdicts({});
            tc.verifyEmpty(s.files);
            tc.verifyEmpty(s.gates);
            tc.verifyEmpty(s.verdicts);
            tc.verifyEqual(s.counts.Pass, 0);
        end

        function three_files_two_gates_each(tc)
            mk = @test_aggregateGateVerdicts.makeReport;
            mg = @test_aggregateGateVerdicts.makeGate;
            r1 = mk('a', {mg('raw','Pass',1), mg('post-ICA','Fail',5)});
            r2 = mk('b', {mg('raw','Marginal',1), mg('post-ICA','Pass',5)});
            r3 = mk('c', {mg('raw','Pass',1), mg('post-ICA','Pass',5)});

            s = aggregateGateVerdicts({r1, r2, r3});
            tc.verifyEqual(s.files, {'a','b','c'});
            tc.verifyEqual(s.gates, {'raw','post-ICA'});
            tc.verifyEqual(s.verdicts, [1 3; 2 1; 1 1]);
            tc.verifyEqual(s.counts.Pass,     4);
            tc.verifyEqual(s.counts.Marginal, 1);
            tc.verifyEqual(s.counts.Fail,     1);
        end

        function report_without_gates_yields_NotChecked_row(tc)
            mk = @test_aggregateGateVerdicts.makeReport;
            mg = @test_aggregateGateVerdicts.makeGate;
            r1 = mk('with',    {mg('raw','Pass',1)});
            r2 = mk('without', {});
            r2.quality = rmfield(r2.quality, 'gates');

            s = aggregateGateVerdicts({r1, r2});
            tc.verifyEqual(s.verdicts, [1; 0]);
            tc.verifyEqual(s.gates, {'raw'});
        end

        function same_label_different_stepIndex_collapses(tc)
            mk = @test_aggregateGateVerdicts.makeReport;
            mg = @test_aggregateGateVerdicts.makeGate;
            % Two gates both labeled 'check' at different positions.
            % Worst-of wins for the single column.
            r = mk('a', {mg('check','Pass',2), mg('check','Fail',9)});
            s = aggregateGateVerdicts({r});
            tc.verifyEqual(s.gates, {'check'});
            tc.verifyEqual(s.verdicts, 3);
        end

        function counts_match_matrix(tc)
            mk = @test_aggregateGateVerdicts.makeReport;
            mg = @test_aggregateGateVerdicts.makeGate;
            r1 = mk('a', {mg('g1','Pass',1), mg('g2','Marginal',2)});
            r2 = mk('b', {mg('g1','Fail',1), mg('g2','Fail',2)});
            r3 = mk('c', {mg('g1','Pass',1)});

            s = aggregateGateVerdicts({r1, r2, r3});
            % 3 Pass: r1.g1, r3.g1, plus r3.g2 is NotChecked
            tc.verifyEqual(s.counts.Pass,     2);
            tc.verifyEqual(s.counts.Marginal, 1);
            tc.verifyEqual(s.counts.Fail,     2);
            tc.verifyEqual(sum(s.verdicts(:) == 0), 1);
        end
    end
end
