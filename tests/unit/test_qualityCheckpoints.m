classdef test_qualityCheckpoints < matlab.unittest.TestCase
% TEST_QUALITYCHECKPOINTS  Drift guard for src/qa/qualityCheckpoints.m.
%   Asserts every name returned by qualityCheckpoints() exists in
%   stepRegistry(). Fails loudly if a registry step is renamed without
%   updating the checkpoint list.

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            tc.addTeardown(@() rmpath(genpath(fullfile(root, 'src'))));
        end
    end

    methods (Test)
        function every_checkpoint_exists_in_registry(tc)
            checkpoints = qualityCheckpoints();
            registry    = stepRegistry();
            registryNames = {registry.name};

            for k = 1:numel(checkpoints)
                tc.verifyTrue(any(strcmp(checkpoints{k}, registryNames)), ...
                    sprintf(['qualityCheckpoints lists "%s" but no step ', ...
                    'with that name exists in stepRegistry. ', ...
                    'Either rename the entry or remove it.'], checkpoints{k}));
            end
        end

        function returns_nonempty_cellstr(tc)
            checkpoints = qualityCheckpoints();
            tc.verifyClass(checkpoints, 'cell');
            tc.verifyGreaterThan(numel(checkpoints), 0);
            for k = 1:numel(checkpoints)
                tc.verifyClass(checkpoints{k}, 'char');
            end
        end
    end
end
