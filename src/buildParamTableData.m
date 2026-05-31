function data = buildParamTableData(step, registryEntry)
% BUILDPARAMTABLEDATA  Build Nx2 cell array for UITable from a spec step.
%
%   data = BUILDPARAMTABLEDATA(step, registryEntry) returns a cell array
%   with one row per parameter in registryEntry.params:
%     column 1: label string — "friendlyName" or "friendlyName (unit)"
%     column 2: display string for the current value in step.params
%
%   step          - PipelineStep struct with .name and .params fields
%   registryEntry - single element from stepRegistry() output
%
%   See also: refreshParamTable, applyParamEdit, stepRegistry

    params = registryEntry.params;
    n      = numel(params);
    data   = cell(n, 2);

    for r = 1:n
        p = params(r);

        if isempty(p.unit)
            data{r,1} = p.friendlyName;
        else
            data{r,1} = [p.friendlyName ' (' p.unit ')'];
        end

        if isfield(step.params, p.key)
            val = step.params.(p.key);
        else
            val = [];
        end

        data{r,2} = formatParamValue(val, p);
    end
end

% ── local helper ─────────────────────────────────────────────────────────────

function s = formatParamValue(val, paramMeta)
    if isempty(val)
        if ~isempty(paramMeta.placeholder)
            s = paramMeta.placeholder;
        else
            s = '(not set)';
        end
    elseif ischar(val) && ~isrow(val)
        s = [deblank(val(1,:)), ' ...'];   % char matrix — show first row
    elseif ischar(val) || isstring(val)
        s = char(val);
    elseif iscell(val)
        s = strjoin(cellfun(@char, val, 'UniformOutput', false), ', ');
    else
        s = mat2str(val);
    end
end
