classdef test_makePipelineStep < matlab.unittest.TestCase
% TEST_MAKEPIPELINESTEP  Unit tests for src/makePipelineStep.m.
%   Verifies step construction, error handling, and registry integrity.

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

    %% Basic construction
    methods (Test)
        function returns_struct_with_name_and_params(tc)
            s = makePipelineStep('Re-Sample', tc.registry);
            tc.verifyEqual(s.name, 'Re-Sample');
            tc.verifyTrue(isstruct(s.params));
        end

    end

    %% Error handling
    methods (Test)
        function unknown_step_throws_nestapp_error(tc)
            tc.verifyError( ...
                @() makePipelineStep('Not A Real Step', tc.registry), ...
                'nestapp:unknownStep');
        end

        function empty_name_throws(tc)
            tc.verifyError( ...
                @() makePipelineStep('', tc.registry), ...
                'nestapp:unknownStep');
        end
    end

    %% Registry completeness
    methods (Test)
        function all_registry_steps_constructible(tc)
            % Every step in the registry must construct without error.
            for k = 1:numel(tc.registry)
                stepName = tc.registry(k).name;
                s = makePipelineStep(stepName, tc.registry);
                tc.verifyEqual(s.name, stepName, ...
                    sprintf('name mismatch for step %d (%s)', k, stepName));
                tc.verifyTrue(isstruct(s.params), ...
                    sprintf('params not a struct for step %s', stepName));
            end
        end

        function params_fields_match_registry_defaults(tc)
            % step.params must have exactly the fields from registry.defaults.
            for k = 1:numel(tc.registry)
                reg  = tc.registry(k);
                s    = makePipelineStep(reg.name, tc.registry);
                expected = fieldnames(reg.defaults);
                actual   = fieldnames(s.params);
                tc.verifyEqual(sort(actual), sort(expected), ...
                    sprintf('%s: params fields differ from registry defaults', reg.name));
            end
        end

        function registry_has_required_steps(tc)
            names = {tc.registry.name};
            required = {'Load Data','Re-Sample','Run ICA','Epoching', ...
                        'Remove Bad Channels','Label ICA Components'};
            for k = 1:numel(required)
                tc.verifyTrue(ismember(required{k}, names), ...
                    sprintf('Registry missing required step: %s', required{k}));
            end
        end

        function all_params_have_type_field(tc)
            % Phase A invariant: every param entry must carry a type string.
            validTypes = {'scalar','integer','vector','logical','string','stringlist'};
            for k = 1:numel(tc.registry)
                params = tc.registry(k).params;
                for p = 1:numel(params)
                    tc.verifyNotEmpty(params(p).type, ...
                        sprintf('%s.%s: type is empty', tc.registry(k).name, params(p).key));
                    tc.verifyTrue(ismember(params(p).type, validTypes), ...
                        sprintf('%s.%s: invalid type "%s"', ...
                            tc.registry(k).name, params(p).key, params(p).type));
                end
            end
        end
    end
end
