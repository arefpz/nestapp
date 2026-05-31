classdef test_paramsToVarin < matlab.unittest.TestCase
% TEST_PARAMSTOVARIN  Unit tests for src/paramsToVarin.m.
%
%   paramsToVarin is the bridge between the typed spec model and the
%   flat key-value varin cells consumed by runPipelineCore's switch/case.
%   Several case branches use positional access (vars{2}, vars{2:2:end}),
%   so field order is a load-bearing contract, not an implementation detail.

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fullfile(root, 'src'));
            tc.addTeardown(@rmpath, fullfile(root, 'src'));
        end
    end

    %% Output shape
    methods (Test)
        function empty_params_gives_empty_cell(tc)
            varin = paramsToVarin(struct());
            tc.verifyEqual(varin, cell(1, 0));
        end

        function single_field_gives_two_element_cell(tc)
            p.freq = 1000;
            varin = paramsToVarin(p);
            tc.verifyEqual(numel(varin), 2);
        end

        function n_fields_gives_2n_element_cell(tc)
            p.a = 1;
            p.b = 2;
            p.c = 3;
            varin = paramsToVarin(p);
            tc.verifyEqual(numel(varin), 6);
        end
    end

    %% Key-value interleaving
    methods (Test)
        function keys_at_odd_positions(tc)
            p.freq = 1000;
            p.fc   = 0.9;
            p.df   = 0.2;
            varin  = paramsToVarin(p);
            tc.verifyEqual(varin{1}, 'freq');
            tc.verifyEqual(varin{3}, 'fc');
            tc.verifyEqual(varin{5}, 'df');
        end

        function values_at_even_positions(tc)
            p.freq = 1000;
            p.fc   = 0.9;
            p.df   = 0.2;
            varin  = paramsToVarin(p);
            tc.verifyEqual(varin{2}, 1000);
            tc.verifyEqual(varin{4}, 0.9);
            tc.verifyEqual(varin{6}, 0.2);
        end
    end

    %% Field order matches registry order (positional-access invariant)
    methods (Test)
        function resample_field_order_matches_registry(tc)
            % Re-Sample's switch case uses vars{2:2:end} positional access.
            % The value at position 2 must be the first param's value (freq).
            % If field order were random, the wrong value would be passed.
            reg = stepRegistry();
            s   = makePipelineStep('Re-Sample', reg);
            s.params.freq = 512;
            varin = paramsToVarin(s.params);

            % Position 1 = first key, position 2 = first value
            tc.verifyEqual(varin{1}, 'freq', ...
                'freq must be first key — positional access depends on this');
            tc.verifyEqual(varin{2}, 512, ...
                'freq value must be at position 2');

            % vars{2:2:end} must give all values in registry order
            allValues = varin(2:2:end);
            tc.verifyEqual(allValues{1}, 512);
        end

        function vector_value_preserved_exactly(tc)
            p.cutTimesTMS = [-2 10];
            varin = paramsToVarin(p);
            tc.verifyEqual(varin{2}, [-2 10]);
        end

        function string_value_preserved_exactly(tc)
            p.icatype = 'fastica';
            varin = paramsToVarin(p);
            tc.verifyEqual(varin{2}, 'fastica');
        end

        function cell_value_preserved_exactly(tc)
            p.types = {'TMS', 'sham'};
            varin = paramsToVarin(p);
            tc.verifyEqual(varin{2}, {'TMS', 'sham'});
        end
    end

    %% Round-trip: makePipelineStep -> paramsToVarin
    methods (Test)
        function roundtrip_preserves_all_registry_keys(tc)
            reg  = stepRegistry();
            step = makePipelineStep('Epoching', reg);
            varin = paramsToVarin(step.params);

            % Every key from the registry must appear in varin
            regIdx    = find(strcmp({reg.name}, 'Epoching'), 1);
            regKeys   = {reg(regIdx).params.key};
            varinKeys = varin(1:2:end);
            for k = 1:numel(regKeys)
                tc.verifyTrue(ismember(regKeys{k}, varinKeys), ...
                    sprintf('Key %s missing from varin', regKeys{k}));
            end
        end
    end
end
