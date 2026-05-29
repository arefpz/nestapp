
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [text, info] = describePipeline(spec, varargin)
% DESCRIBEPIPELINE  Human-readable summary of a pipeline spec.
%
%   Syntax
%     describePipeline(spec)                 % print a numbered step list
%     text = describePipeline(spec)          % also return the Markdown text
%     [text, info] = describePipeline(spec)  % also return a data struct
%     describePipeline(spec, 'Copy', true)   % copy the text to the clipboard
%     describePipeline(spec, 'Quiet', true)  % do not print to the command window
%
%   Inputs
%     spec   struct array - pipeline steps with fields .name and .params
%                           (as produced by makePipelineStep / saved pipelines).
%     Copy   logical (name-value) - copy the description to the clipboard.
%     Quiet  logical (name-value) - suppress printing to the command window.
%
%   Outputs
%     text  char   - a Markdown description: a numbered list of steps, each
%                    with the parameters that DIFFER from the registry default
%                    (so the report shows exactly what was customised).
%     info  struct - .nestappVersion, .nSteps, and .steps (per-step name +
%                    overridden params).
%
%   Pairs with nestappDoctor: that reports the environment, this reports what
%   the user actually ran, so a bug report can be reproduced. Pure and
%   EEGLAB-free; uses stepRegistry only to resolve defaults.
%
%   Side effects   Optionally writes to the clipboard; prints unless Quiet.
%   See also: nestappDoctor, stepRegistry, makePipelineStep, paramsToVarin

p = inputParser;
p.addParameter('Copy',  false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Quiet', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
doCopy  = logical(p.Results.Copy);
beQuiet = logical(p.Results.Quiet);

try
    reg = stepRegistry();
catch
    reg = struct('name', {}, 'defaults', {});
end

info = struct();
info.nestappVersion = nestappVersion();
info.nSteps = numel(spec);
info.steps  = struct('name', {}, 'overrides', {});

L = {};
L{end+1} = '### nestapp pipeline';
L{end+1} = '';
L{end+1} = sprintf('nestapp %s - %d step(s)', info.nestappVersion, numel(spec));
L{end+1} = '';

if isempty(spec)
    L{end+1} = '_(empty pipeline - no steps)_';
else
    for i = 1:numel(spec)
        overrides = overriddenParams(spec(i), reg);
        info.steps(end+1) = struct('name', spec(i).name, 'overrides', {overrides});
        L{end+1} = sprintf('%d. **%s**', i, spec(i).name); %#ok<AGROW>
        if isempty(overrides)
            L{end+1} = '   - _(defaults)_'; %#ok<AGROW>
        else
            keys = fieldnames(overrides);
            for k = 1:numel(keys)
                L{end+1} = sprintf('   - %s = %s', keys{k}, valToStr(overrides.(keys{k}))); %#ok<AGROW>
            end
        end
    end
end

text = strjoin(L, newline);

if ~beQuiet
    fprintf('%s\n', text);
end
if doCopy
    try
        clipboard('copy', text);
        if ~beQuiet; fprintf('(pipeline description copied to clipboard)\n'); end
    catch ME
        warning('describePipeline:clipboard', 'Could not copy to clipboard: %s', ME.message);
    end
end
end

% ── helpers ───────────────────────────────────────────────────────────────────

function ov = overriddenParams(step, reg)
% Return a struct of the step's params that differ from the registry default.
ov = struct();
if ~isfield(step, 'params') || isempty(fieldnames(step.params))
    return
end
idx = find(strcmp({reg.name}, step.name), 1);
if isempty(idx)
    ov = step.params;   % unknown step - show everything
    return
end
defaults = reg(idx).defaults;
keys = fieldnames(step.params);
for k = 1:numel(keys)
    key = keys{k};
    val = step.params.(key);
    if ~isfield(defaults, key) || ~isequaln(val, defaults.(key))
        ov.(key) = val;
    end
end
end

function s = valToStr(v)
% Compact, human-readable rendering of a parameter value.
if ischar(v)
    s = v;
elseif isstring(v)
    s = char(strjoin(v, ', '));
elseif islogical(v)
    s = mat2str(v);
elseif isnumeric(v)
    if isempty(v)
        s = '[]';
    elseif isscalar(v)
        s = num2str(v);
    else
        s = ['[' strtrim(num2str(v(:)')) ']'];
    end
elseif iscell(v)
    parts = cellfun(@valToStr, v, 'UniformOutput', false);
    s = ['{' strjoin(parts, ', ') '}'];
else
    s = sprintf('<%s>', class(v));
end
end
