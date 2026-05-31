classdef test_specFromSaved < matlab.unittest.TestCase
% TEST_SPECFROMSAVED  Unit tests for src/specFromSaved.m.
%   New format only (data.spec). Old-format migration was removed.

    properties
        registry
    end

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fullfile(root, 'src'));
            tc.addTeardown(@rmpath, fullfile(root, 'src'));
        end
    end

    methods (TestMethodSetup)
        function loadRegistry(tc)
            tc.registry = stepRegistry();
        end
    end

    %% New format (data.spec present)
    methods (Test)
        function new_format_single_step_passthrough(tc)
            s = makePipelineStep('Re-Sample', tc.registry);
            data.spec = s;
            [spec, warns] = specFromSaved(data, tc.registry);
            tc.verifyEqual(numel(spec), 1);
            tc.verifyEqual(spec(1).name, 'Re-Sample');
            tc.verifyEmpty(warns);
        end

        function new_format_multi_step_passthrough(tc)
            s1 = makePipelineStep('Load Data',  tc.registry);
            s2 = makePipelineStep('Re-Sample',  tc.registry);
            s3 = makePipelineStep('Run ICA',    tc.registry);
            data.spec = [s1, s2, s3];
            [spec, warns] = specFromSaved(data, tc.registry);
            tc.verifyEqual(numel(spec), 3);
            tc.verifyEqual({spec.name}, {'Load Data','Re-Sample','Run ICA'});
            tc.verifyEmpty(warns);
        end

        function new_format_unknown_step_produces_warning(tc)
            data.spec = struct('name', 'Deleted Step', 'params', struct());
            [spec, warns] = specFromSaved(data, tc.registry);
            tc.verifyEqual(numel(spec), 1);
            tc.verifyEqual(numel(warns), 1);
            tc.verifyTrue(contains(warns{1}, 'Deleted Step'));
        end

        function new_format_params_preserved(tc)
            s = makePipelineStep('Re-Sample', tc.registry);
            s.params.freq = 512;
            data.spec = s;
            [spec, ~] = specFromSaved(data, tc.registry);
            tc.verifyEqual(spec(1).params.freq, 512);
        end
    end

    %% Missing spec field
    methods (Test)
        function no_spec_field_returns_empty_with_warning(tc)
            data = struct('somethingElse', 1);
            [spec, warns] = specFromSaved(data, tc.registry);
            tc.verifyEqual(numel(spec), 0);
            tc.verifyNotEmpty(warns);
        end
    end
end
