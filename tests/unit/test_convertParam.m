classdef test_convertParam < matlab.unittest.TestCase
% TEST_CONVERTPARAM  Unit tests for src/convertParam.m.
%   Tests every type branch with representative inputs including edge cases.

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fullfile(root, 'src'));
            tc.addTeardown(@rmpath, fullfile(root, 'src'));
        end
    end

    %% scalar / integer
    methods (Test)
        function scalar_numeric_passthrough(tc)
            tc.verifyEqual(convertParam(3.14, 'scalar'), 3.14);
        end

        function scalar_from_string(tc)
            tc.verifyEqual(convertParam('250', 'scalar'), 250);
        end

        function scalar_from_cell_scalar_string(tc)
            tc.verifyEqual(convertParam({'42'}, 'scalar'), 42);
        end

        function scalar_from_numeric_vector_takes_first(tc)
            tc.verifyEqual(convertParam([5 10 15], 'scalar'), 5);
        end

        function integer_numeric_passthrough(tc)
            tc.verifyEqual(convertParam(1000, 'integer'), 1000);
        end

        function integer_from_string(tc)
            tc.verifyEqual(convertParam('64', 'integer'), 64);
        end
    end

    %% vector
    methods (Test)
        function vector_numeric_row_passthrough(tc)
            tc.verifyEqual(convertParam([1 2 3], 'vector'), [1 2 3]);
        end

        function vector_numeric_col_becomes_row(tc)
            tc.verifyEqual(convertParam([1;2;3], 'vector'), [1 2 3]);
        end

        function vector_from_char_bracketed(tc)
            tc.verifyEqual(convertParam('[-500 -10]', 'vector'), [-500 -10]);
        end

        function vector_from_char_space_separated(tc)
            tc.verifyEqual(convertParam('0 500', 'vector'), [0 500]);
        end

        function vector_from_numeric_cell(tc)
            tc.verifyEqual(convertParam({'-500'; '-10'}, 'vector'), [-500 -10]);
        end

        function vector_from_numeric_cell_row(tc)
            tc.verifyEqual(convertParam({1, 2, 3}, 'vector'), [1 2 3]);
        end

        function vector_non_numeric_cell_passthrough(tc)
            v = {'Fp1', 'Fp2'};
            result = convertParam(v, 'vector');
            % Non-numeric cell should pass through unchanged
            tc.verifyEqual(result, v);
        end

        function vector_non_numeric_char_passthrough(tc)
            result = convertParam('some_channel', 'vector');
            % Non-parseable char passes through as char
            tc.verifyClass(result, 'char');
        end
    end

    %% string
    methods (Test)
        function string_char_passthrough(tc)
            tc.verifyEqual(convertParam('fastica', 'string'), 'fastica');
        end

        function string_from_empty_numeric(tc)
            tc.verifyEqual(convertParam([], 'string'), '');
        end

        function string_from_string_object(tc)
            tc.verifyEqual(convertParam(string('runica'), 'string'), 'runica');
        end

        function string_from_numeric_scalar(tc)
            tc.verifyEqual(convertParam(42, 'string'), '42');
        end
    end

    %% stringlist
    methods (Test)
        function stringlist_cell_passthrough(tc)
            v = {'Fp1', 'Fp2', 'Fp3'};
            tc.verifyEqual(convertParam(v, 'stringlist'), v);
        end

        function stringlist_from_comma_char(tc)
            result = convertParam('Fp1, Fp2, Fp3', 'stringlist');
            tc.verifyEqual(result, {'Fp1', 'Fp2', 'Fp3'});
        end

        function stringlist_from_comma_char_no_spaces(tc)
            result = convertParam('TP9,TP10', 'stringlist');
            tc.verifyEqual(result, {'TP9', 'TP10'});
        end

        function stringlist_from_empty_numeric(tc)
            tc.verifyEqual(convertParam([], 'stringlist'), {});
        end

        function stringlist_from_bracket_sentinel(tc)
            % Regression: old MATLAB template .mat files store "not set"
            % stringlist params as the char '[]'. Must convert to {} not {'[]'}.
            tc.verifyEqual(convertParam('[]', 'stringlist'), {});
        end

        function stringlist_single_string_wrapped(tc)
            result = convertParam(42, 'stringlist');
            tc.verifyClass(result, 'cell');
            tc.verifyEqual(numel(result), 1);
        end

        function stringlist_from_multi_element_string_array(tc)
            % Regression: multi-element string arrays from old template .mat
            % files caused char() to produce a char matrix, crashing strsplit.
            raw = string(["Fp1", "Fp2", "Fp3"]);
            result = convertParam(raw, 'stringlist');
            tc.verifyClass(result, 'cell');
            tc.verifyEqual(result, {'Fp1', 'Fp2', 'Fp3'});
        end
    end

    %% fallback
    methods (Test)
        function unknown_type_passthrough(tc)
            val = struct('x', 1, 'y', 2);
            tc.verifyEqual(convertParam(val, 'some_future_type'), val);
        end
    end
end
