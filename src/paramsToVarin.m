function varin = paramsToVarin(params)
% PARAMSTOVARIN Convert a typed params struct to a flat key-value cell.
%   varin = PARAMSTOVARIN(params) returns a 1x2N cell where element 2k-1
%   is the field name and element 2k is the corresponding value, in the
%   same order as fieldnames(params).
%
%   This ordering matches the registry param order because makePipelineStep
%   builds params structs by assigning fields in registry order. Step cases
%   in runPipelineCore that use positional access (vars{2}, vars{2:2:end})
%   depend on this invariant being preserved.

keys  = fieldnames(params);
varin = cell(1, 2*numel(keys));
for k = 1:numel(keys)
    varin{2*k-1} = keys{k};
    varin{2*k}   = params.(keys{k});
end
end
