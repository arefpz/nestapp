% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef MethodsProseTest < NestappTestCase
% METHODSPROSETEST  The methods paragraph covers every step and reads correctly.
%
%   The paragraph is built from whatever pipeline a user assembles, so the
%   tests are written against the REGISTRY rather than against any one
%   pipeline: every step, with its defaults and with each parameter pushed to
%   awkward values, must produce a clause that reads as prose. That is what
%   "works beyond our own data" means here - a pipeline nobody has tried is
%   still made of these steps and these parameters.
%
%   What went wrong before, and is pinned here:
%     - 23 data-changing steps fell through to an empty clause and silently
%       vanished (RANSAC, Robust Detrend, Manual Command among them);
%     - repeats were dropped as duplicates (a second ICA, the final
%       re-reference);
%     - the TMS window verb looked at one end only;
%     - a notch filter was described as a band-pass.

    properties (Constant)
        % Silent only for some settings, and silent with their defaults.
        CONDITIONAL = {'Load Channel Location', 'Visualize EEG Data'};
    end

    methods (Test)
        function everyStepIsCoveredOrDeclaredSilent(tc)
            reg    = stepRegistry();
            silent = methodsSilentSteps();
            for s = reg(:)'
                c = methodsClause(s.name, s.defaults);
                if any(strcmp(s.name, silent)) || any(strcmp(s.name, tc.CONDITIONAL))
                    tc.verifyEmpty(c, sprintf('%s should add no sentence with its defaults', s.name));
                else
                    tc.verifyNotEmpty(c, sprintf('%s produces no methods sentence', s.name));
                    tc.verifyFalse(startsWith(c, 'the "'), ...
                        sprintf('%s falls through to the generic clause - give it a case', s.name));
                end
            end
        end

        function silentListNamesRealSteps(tc)
            names = {stepRegistry().name};
            for n = methodsSilentSteps()
                tc.verifyTrue(any(strcmp(n{1}, names)), sprintf('%s is not a registry step', n{1}));
            end
        end

        function clausesReadAsProseForAnyParameterValue(tc)
            % Every parameter of every step, one at a time, set to each value
            % a saved pipeline or a table edit can deliver.
            odd = {[], '', '[]', 'off', 'on', 0, 1, -1, 1e6, [1 2], [NaN NaN], {'a', 'b'}, "str"};
            reg = stepRegistry();
            for s = reg(:)'
                keys = {};
                if isstruct(s.params) && ~isempty(s.params); keys = {s.params.key}; end
                tc.assertProse(s.name, s.defaults, keys);
                for k = keys
                    for v = odd
                        p = s.defaults; p.(k{1}) = v{1};
                        tc.assertProse(s.name, p, keys);
                    end
                end
            end
        end

        function unknownStepGetsAGenericClauseNotSilence(tc)
            c = methodsClause('Some Future Step', struct('x', 3));
            tc.verifyEqual(c, 'the "Some Future Step" step was applied');
        end

        function manualCommandNeverQuotesItsDescription(tc)
            % The description is an editor label (it can hold notes such as
            % "was a step, now a direct call"), so it must not reach the text.
            c = methodsClause('Manual Command', struct('command', 'EEG = foo(EEG);', ...
                'description', 'AARATEP peak-amplitude IC flag (15 uV) - was a step, now a direct call'));
            tc.verifyEqual(c, 'a custom processing step was applied (see the pipeline specification)');
        end

        function thresholdsAreStated(tc)
            tc.verifySubstring(methodsClause('Remove Bad Epoch', struct('startprob', [], 'maxrej', [])), ...
                'joint probability beyond 5 SD');
            c = methodsClause('Flag ICA Components for Rejection', struct( ...
                'Eye', [0.4 1], 'Muscle', [0.2 1], 'Brain', [0 0.3], 'Heart', [NaN NaN]));
            tc.verifySubstring(c, 'eye (p >= 0.4)');
            tc.verifySubstring(c, 'brain with p <= 0.3');
            tc.verifySubstring(c, 'flagged');
            tc.verifySubstring(methodsClause('Remove ICA Components (TESA)', ...
                struct('tmsMuscle', 'on', 'tmsMuscleThresh', 8, 'tmsMuscleWin', [11 30])), ...
                '11 to 30 ms amplitude ratio > 8');
            tc.verifySubstring(methodsClause('Detect Bad Channels (RANSAC)', ...
                stepDefaults('Detect Bad Channels (RANSAC)')), 'correlation with the predicted signal < 0.8');
        end

        function notchIsBandStop(tc)
            c = methodsClause('Frequency Filter', struct('locutoff', 58, 'hicutoff', 62, 'revfilt', 1));
            tc.verifySubstring(c, 'band-stop filtered at 58-62 Hz');
        end

        function modifiedFilterNamesItsDesign(tc)
            c = methodsClause('Modified Bandpass Filter (TESA)', struct('lowCutoff', 1, ...
                'highCutoff', [], 'filterMethod', 'eegfiltnew'));
            tc.verifyFalse(contains(c, 'Butterworth'));
            tc.verifySubstring(c, 'FIR');
        end

        % ── repeats ──────────────────────────────────────────────────────
        function identicalRepeatReadsAgain(tc)
            b = prose({'Run TESA ICA', struct()}, {'Re-Reference', struct('ref', [])}, ...
                      {'Run TESA ICA', struct()}, {'Re-Reference', struct('ref', [])});
            tc.verifySubstring(b, 'by FastICA again');
            tc.verifySubstring(b, 'common average again');
        end

        function immediateRepeatIsStatedOnce(tc)
            b = prose({'Re-Sample', struct('freq', 1000)}, {'Re-Sample', struct('freq', 1000)});
            tc.verifyEqual(count(b, 'downsampled'), 1);
            tc.verifyFalse(contains(b, 'again'));
        end

        function changedRepeatReadsLater(tc)
            b = prose({'Remove Baseline', struct('timerange', [-500 -10])}, ...
                      {'Re-Sample', struct('freq', 1000)}, ...
                      {'Remove Baseline', struct('timerange', [-1000 -2])});
            tc.verifySubstring(b, 'Later, the data were baseline-corrected from -1000 to -2 ms');
        end

        function tmsWindowVerbUsesBothEnds(tc)
            cut = @(w) {'Remove TMS Artifacts (TESA)', struct('cutTimesTMS', w)};
            tc.verifySubstring(prose(cut([-2 10]), cut([-2 15])), 'later extended to -2 to 15 ms');
            tc.verifySubstring(prose(cut([-2 15]), cut([-2 10])), 'later narrowed to -2 to 10 ms');
            tc.verifySubstring(prose(cut([-2 10]), cut([0 10])), 'later narrowed to 0 to 10 ms');
            tc.verifySubstring(prose(cut([-2 15]), cut([-5 12])), 'later changed to -5 to 12 ms');
            b = prose(cut([-2 10]), cut([-2 15]), cut([-2 10]));
            tc.verifySubstring(b, 'later narrowed to -2 to 10 ms');
        end

        function identicalRecutAfterOtherStepsIsKept(tc)
            cut = {'Remove TMS Artifacts (TESA)', struct('cutTimesTMS', [-2 10])};
            b = prose(cut, {'Run TESA ICA', struct()}, cut);
            tc.verifySubstring(b, '(-2 to 10 ms) were removed again');
        end

        function endsWithThePointerToTheSpec(tc)
            b = prose({'Re-Sample', struct('freq', 1000)});
            tc.verifyTrue(endsWith(b, 'pipeline specification saved with the run (spec.mat).'));
        end
    end

    methods (Access = private)
        function assertProse(tc, name, params, keys)
            try
                c = methodsClause(name, params);
            catch err
                tc.verifyFail(sprintf('%s errored: %s', name, err.message));
                return
            end
            tc.verifyTrue(ischar(c) && (isempty(c) || isrow(c)), sprintf('%s: not a char row', name));
            if isempty(c); return; end
            % Lower-case start, unless it opens with an acronym ("TMS-evoked"):
            % pipelineProse capitalises the sentence either way.
            tc.verifyTrue(isstrprop(c(1), 'lower') || ~isempty(regexp(c, '^[A-Z]{2,}', 'once')), ...
                sprintf('%s: "%s" must start lower-case', name, c));
            bad = {'NaN', '[]', '()', '%g', '%s', '%d'};
            for b = bad
                tc.verifyFalse(contains(c, b{1}), sprintf('%s: "%s" contains %s', name, c, b{1}));
            end
            for k = keys
                if numel(k{1}) > 5 && any(isstrprop(k{1}(2:end), 'upper'))
                    tc.verifyFalse(contains(c, k{1}), ...
                        sprintf('%s: "%s" leaks parameter key %s', name, c, k{1}));
                end
            end
        end
    end
end

% ── helpers ─────────────────────────────────────────────────────────────────

function b = prose(varargin)
% pipelineProse over a synthetic run: each argument is {stepName, params}.
    steps = cellfun(@(a) struct('name', a{1}, 'params', a{2}), varargin, 'UniformOutput', false);
    b = pipelineProse(struct('steps', {steps}));
end

function d = stepDefaults(name)
    reg = stepRegistry();
    d = reg(strcmp({reg.name}, name)).defaults;
end
