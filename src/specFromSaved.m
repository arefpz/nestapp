function [spec, warnings] = specFromSaved(data, registry)
% SPECFROMSAVED Load a PipelineSpec from a saved pipeline .mat.
%   [spec, warnings] = SPECFROMSAVED(data, registry)
%
%   data must contain a 'spec' field. Returns a struct array of PipelineStep
%   ({name, params}). Unknown steps are included with a warning.

warnings = {};

if ~isfield(data, 'spec')
    spec = repmat(struct('name','','params',struct()), 0, 1);
    warnings{end+1} = 'File is not in the current format. Re-save the pipeline from the app.';
    return
end

spec = data.spec;
for k = 1:numel(spec)
    if ~any(strcmp({registry.name}, spec(k).name))
        warnings{end+1} = sprintf('Unknown step "%s" (not in registry)', spec(k).name); %#ok<AGROW>
    end
end
end
