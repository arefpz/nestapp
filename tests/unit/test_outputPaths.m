
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef test_outputPaths < matlab.unittest.TestCase
% TEST_OUTPUTPATHS  Unit tests for src/io/outputPaths.m and buildBatchContext.m

    properties (Access = private)
        tmpRoot char
        savedPref char     % outputRoot pref value at test start (restored on teardown)
        hadPref logical
    end

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            tc.addTeardown(@() rmpath(genpath(fullfile(root, 'src'))));
        end
    end

    methods (TestMethodSetup)
        function makeTmpRoot(tc)
            tc.tmpRoot = fullfile(tempdir, ['nestapp_op_', char(matlab.lang.internal.uuid())]);
            mkdir(tc.tmpRoot);
            tc.addTeardown(@() rmdir(tc.tmpRoot, 's'));

            tc.hadPref   = ispref('nestapp', 'outputRoot');
            if tc.hadPref
                tc.savedPref = getpref('nestapp', 'outputRoot');
            else
                tc.savedPref = '';
            end
            tc.addTeardown(@() tc.restoreOutputRootPref());
        end
    end

    methods (Access = private)
        function restoreOutputRootPref(tc)
            if tc.hadPref
                setpref('nestapp', 'outputRoot', tc.savedPref);
            elseif ispref('nestapp', 'outputRoot')
                rmpref('nestapp', 'outputRoot');
            end
        end

        function paths = makeFilePaths(tc, n)
            paths = cell(1, n);
            for k = 1:n
                paths{k} = fullfile(tc.tmpRoot, sprintf('sub_%d.vhdr', k));
            end
        end
    end

    methods (Test)

        % -- buildBatchContext ------------------------------------------

        function batchContext_uses_pref_when_set(tc)
            setpref('nestapp', 'outputRoot', tc.tmpRoot);
            ctx = buildBatchContext(tc.makeFilePaths(2), 'TMS-EEG / TEP (TESA)', 'typeBased');
            tc.verifyEqual(ctx.outputRoot, tc.tmpRoot);
            tc.verifyTrue(startsWith(ctx.batchId, char(ctx.startedAt, 'yyyyMMdd_HHmmss')));
            tc.verifyTrue(endsWith(ctx.batchId, ctx.pipelineSlug));
            tc.verifyTrue(exist(ctx.batchRoot, 'dir') == 7);
        end

        function batchContext_falls_back_to_common_parent(tc)
            if ispref('nestapp', 'outputRoot'); rmpref('nestapp', 'outputRoot'); end
            paths = tc.makeFilePaths(3);
            ctx = buildBatchContext(paths, 'Resting-State EEG', 'perInput');
            tc.verifyEqual(ctx.outputRoot, tc.tmpRoot);
            tc.verifyEqual(ctx.layout, 'perInput');
        end

        function pipelineSlug_sanitizes_punctuation(tc)
            % Pass tmpRoot as the explicit outputRoot override so the
            % batch folder lands in tempdir even when the user has an
            % nestapp.outputRoot pref set in their MATLAB session.
            ctx = buildBatchContext(tc.makeFilePaths(1), ...
                'TMS-EEG / TEP (TESA + Quality Gates)', 'typeBased', tc.tmpRoot);
            % Lowercase + [a-z0-9_]+, no leading/trailing underscores.
            tc.verifyTrue(~isempty(regexp(ctx.pipelineSlug, '^[a-z0-9_]+$', 'once')));
            tc.verifyFalse(startsWith(ctx.pipelineSlug, '_'));
            tc.verifyFalse(endsWith(ctx.pipelineSlug, '_'));
        end

        function pipelineSlug_defaults_when_empty(tc)
            ctx = buildBatchContext(tc.makeFilePaths(1), '', 'typeBased', tc.tmpRoot);
            tc.verifyEqual(ctx.pipelineSlug, 'pipeline');
        end

        function layout_falls_back_to_typeBased_when_unknown(tc)
            ctx = buildBatchContext(tc.makeFilePaths(1), 'p', 'bogus', tc.tmpRoot);
            tc.verifyEqual(ctx.layout, 'typeBased');
        end

        function explicit_outputRoot_override_beats_pref(tc)
            % Regression: tests must not pollute the repo root when the
            % user has nestapp.outputRoot pointed there. The override
            % argument must take precedence over the pref.
            poisonedRoot = fullfile(tempdir, ...
                ['nestapp_poison_', char(matlab.lang.internal.uuid())]);
            mkdir(poisonedRoot);
            tc.addTeardown(@() rmdir(poisonedRoot, 's'));
            setpref('nestapp', 'outputRoot', poisonedRoot);

            ctx = buildBatchContext(tc.makeFilePaths(1), 'p', 'typeBased', tc.tmpRoot);
            tc.verifyEqual(ctx.outputRoot, tc.tmpRoot, ...
                'Override must beat the outputRoot pref.');
            tc.verifyTrue(startsWith(ctx.batchRoot, tc.tmpRoot), ...
                'Batch folder must land under the override root, not the pref root.');
        end

        % -- outputPaths typeBased --------------------------------------

        function typeBased_layout_paths(tc)
            setpref('nestapp', 'outputRoot', tc.tmpRoot);
            ctx = buildBatchContext(tc.makeFilePaths(1), 'p', 'typeBased');

            d1 = outputPaths(ctx, 'data',    'sub_1');
            d2 = outputPaths(ctx, 'reports', 'sub_1');
            d3 = outputPaths(ctx, 'qc',      'sub_1');
            d4 = outputPaths(ctx, 'batch');

            tc.verifyEqual(d1, fullfile(ctx.batchRoot, 'data'));
            tc.verifyEqual(d2, fullfile(ctx.batchRoot, 'reports'));
            tc.verifyEqual(d3, fullfile(ctx.batchRoot, 'qc', 'sub_1'));
            tc.verifyEqual(d4, fullfile(ctx.batchRoot, 'batch'));

            % mkdir-on-demand: every directory we asked for now exists.
            for d = {d1, d2, d3, d4}
                tc.verifyTrue(exist(d{1}, 'dir') == 7);
            end
        end

        % -- outputPaths perInput ---------------------------------------

        function perInput_layout_paths(tc)
            setpref('nestapp', 'outputRoot', tc.tmpRoot);
            ctx = buildBatchContext(tc.makeFilePaths(1), 'p', 'perInput');

            d1 = outputPaths(ctx, 'data',    'sub_1');
            d2 = outputPaths(ctx, 'reports', 'sub_1');
            d3 = outputPaths(ctx, 'qc',      'sub_1');
            d4 = outputPaths(ctx, 'batch');

            tc.verifyEqual(d1, fullfile(ctx.batchRoot, 'sub_1'));
            tc.verifyEqual(d2, fullfile(ctx.batchRoot, 'sub_1'));
            tc.verifyEqual(d3, fullfile(ctx.batchRoot, 'sub_1', 'qc'));
            tc.verifyEqual(d4, fullfile(ctx.batchRoot, '_batch'));

            for d = {d1, d2, d3, d4}
                tc.verifyTrue(exist(d{1}, 'dir') == 7);
            end
        end

        % -- error paths ------------------------------------------------

        function missing_stem_errors_for_perfile_kinds(tc)
            setpref('nestapp', 'outputRoot', tc.tmpRoot);
            ctx = buildBatchContext(tc.makeFilePaths(1), 'p', 'typeBased');
            tc.verifyError(@() outputPaths(ctx, 'data', ''), ...
                'nestapp:outputPaths:missingStem');
            tc.verifyError(@() outputPaths(ctx, 'reports'), ...
                'nestapp:outputPaths:missingStem');
            tc.verifyError(@() outputPaths(ctx, 'qc'), ...
                'nestapp:outputPaths:missingStem');
        end

        function unknown_kind_errors(tc)
            setpref('nestapp', 'outputRoot', tc.tmpRoot);
            ctx = buildBatchContext(tc.makeFilePaths(1), 'p', 'typeBased');
            tc.verifyError(@() outputPaths(ctx, 'bogus'), ...
                'nestapp:outputPaths:badKind');
        end
    end
end
