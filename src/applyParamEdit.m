function spec = applyParamEdit(spec, stepIdx, row, newData, registryEntry)
% APPLYPARAMEDIT  Apply a UITable cell edit to the pipeline spec.
%
%   spec = APPLYPARAMEDIT(spec, stepIdx, row, newData, registryEntry)
%   parses newData according to the type declared in registryEntry.params(row)
%   and writes the converted value into spec(stepIdx).params. Returns the
%   updated spec unchanged if row exceeds the parameter count.
%
%   spec          - struct array of PipelineStep
%   stepIdx       - 1-based index into spec
%   row           - UITable row that was edited (1-based)
%   newData       - raw value from the UITable cell edit event
%   registryEntry - single element from stepRegistry() with .params
%
%   See also: buildParamTableData, convertParam, stepRegistry

    params = registryEntry.params;
    if row > numel(params)
        return
    end

    paramMeta = params(row);
    val = convertParam(newData, paramMeta.type);
    spec(stepIdx).params.(paramMeta.key) = val;
end
