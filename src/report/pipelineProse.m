
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function body = pipelineProse(report)
% PIPELINEPROSE  Ordered methods-paragraph body for a single report's pipeline.
%   body = PIPELINEPROSE(report) walks report.steps in run order, turns each
%   into a clause via methodsClause, and returns punctuated sentences ending
%   with a pointer to the pipeline specification saved with the run. Returns ''
%   when no step carries parameters (e.g. a legacy report saved before
%   stepRec.params existed), which lets the callers fall back to an
%   outcomes-only sentence.
%
%   REPEATS follow one rule for every step, rather than special cases per
%   step, so a pipeline nobody here has seen still reads correctly:
%     - the same step again with the same settings reads "... again". It used
%       to be dropped as a duplicate, which made a second ICA, the final
%       re-reference and a re-interpolation vanish from the paragraph;
%     - the same step again with different settings, further on, is led by
%       "later,"; straight after itself it needs no lead-in;
%     - an immediate repeat (the same step twice in a row, same settings) is
%       stated once, since it describes nothing new.
%   Manual Command is exempt: each one runs different code, so two of them
%   are two steps, not a repeat.
%   The TMS artifact window keeps its own sentence ("the artifact-removal
%   window was later extended to ..."), with the verb judged from both ends of
%   the window: extended if the new one contains the old, narrowed if it lies
%   inside it, and changed otherwise.
%
%   The text is generated from the report's own parameters - it copies no
%   external source. See methodsClause for the per-step wording.
%
%   See also: methodsClause, methodsSilentSteps, methodsNarrative,
%             methodsParagraphAggregate

    body = '';
    if ~isfield(report, 'steps') || isempty(report.steps)
        return
    end
    clauses  = {};
    lastBy   = containers.Map();   % step name -> its last clause, as written by methodsClause
    prevName = '';                 % step behind the previous sentence, for immediate repeats
    tmsWin   = [];                 % current TMS artifact window

    for i = 1:numel(report.steps)
        s = report.steps{i};
        if ~isstruct(s) || ~isfield(s, 'name') || ~isfield(s, 'params')
            continue            % legacy record without params - skip
        end
        name = canonicalStepName(s.name);
        c    = methodsClause(name, s.params);
        if isempty(c); continue; end

        out = c;
        if isKey(lastBy, name) && ~strcmp(name, 'Manual Command')
            win = tmsWindow(name, s.params);
            if ~isempty(win) && ~isempty(tmsWin) && ~isequal(win, tmsWin)
                out = sprintf('the artifact-removal window was later %s to %g to %g ms', ...
                    windowChange(tmsWin, win), win(1), win(2));
            elseif strcmp(c, lastBy(name))
                if strcmp(prevName, name); continue; end
                out = [c ' again'];
            elseif ~strcmp(prevName, name)
                out = ['later, ' c];
            end
        end
        win = tmsWindow(name, s.params);
        if ~isempty(win); tmsWin = win; end
        lastBy(name) = c;
        prevName     = name;
        clauses{end+1} = out; %#ok<AGROW>
    end

    if isempty(clauses)
        return
    end
    for i = 1:numel(clauses)
        clauses{i} = [upper(clauses{i}(1)), clauses{i}(2:end)];
    end
    % Clauses carry the settings a reader needs; everything else is in the
    % spec saved with every run (runPipelineCore writes batch/spec.mat).
    clauses{end+1} = ['Full step parameters are recorded in the pipeline ' ...
        'specification saved with the run (spec.mat)'];
    body = [strjoin(clauses, '. '), '.'];
end

% ── helpers ───────────────────────────────────────────────────────────────────

function win = tmsWindow(name, params)
    win = [];
    if strcmp(name, 'Remove TMS Artifacts (TESA)') && isfield(params, 'cutTimesTMS')
        v = params.cutTimesTMS;
        if ischar(v) || isstring(v); v = str2num(char(v)); end %#ok<ST2NM>
        v = double(v(:))';
        if numel(v) >= 2; win = v(1:2); end
    end
end

function verb = windowChange(old, new)
% Judged from both ends: a window that grew on one side and shrank on the
% other was neither extended nor narrowed.
    if new(1) <= old(1) && new(2) >= old(2)
        verb = 'extended';
    elseif new(1) >= old(1) && new(2) <= old(2)
        verb = 'narrowed';
    else
        verb = 'changed';
    end
end
