function entry = buildHistoryEntry(spec, pipelineName)
% BUILDHISTORYENTRY  Build a human-readable provenance string for EEG.history.
%   entry = BUILDHISTORYENTRY(spec, pipelineName)
timestamp = string(datetime('now'), 'yyyy-MM-dd HH:mm:ss');
if isempty(pipelineName)
    pipelineName = '(unsaved)';
end
lines = { ...
    sprintf('%% --- nestapp pipeline  [%s] ---', timestamp), ...
    sprintf('%% Pipeline: %s', pipelineName), ...
    '%  Steps:' ...
};
for si = 1:numel(spec)
    keys   = fieldnames(spec(si).params);
    pairs  = cell(1, numel(keys));
    for ki = 1:numel(keys)
        val = spec(si).params.(keys{ki});
        if isnumeric(val)
            valStr = mat2str(val);
        else
            valStr = char(val);
        end
        pairs{ki} = sprintf('%s=%s', keys{ki}, valStr);
    end
    if isempty(pairs)
        paramStr = '';
    else
        paramStr = ['  [', strjoin(pairs, ', '), ']'];
    end
    lines{end+1} = sprintf('%%  %2d. %s%s', si, spec(si).name, paramStr); %#ok<AGROW>
end
entry = strjoin(lines, newline);
end
