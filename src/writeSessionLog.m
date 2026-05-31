function writeSessionLog(pathName, fileName, stepLog)
% WRITESESSIONLOG  Write a plain-text processing log alongside the data file.
%   writeSessionLog(pathName, fileName, stepLog)
[~, baseName] = fileparts(fileName);
logPath = fullfile(pathName, [baseName, '_nestapp_log.txt']);
fid = fopen(logPath, 'w');
if fid == -1; return; end
fprintf(fid, '=== nestapp session log ===\n');
fprintf(fid, 'File:      %s\n', fileName);
fprintf(fid, 'Processed: %s\n', datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
fprintf(fid, 'MATLAB:    %s\n', version);
fprintf(fid, '\n%-4s  %-35s  %7s  %9s  %9s  %s\n', ...
    '#', 'Step', 'Time(s)', 'Ch before', 'Ch after', 'Note');
fprintf(fid, '%s\n', repmat('-', 1, 80));
for k = 1:numel(stepLog)
    s = stepLog(k);
    note = s.error;
    if s.epochBefore ~= s.epochAfter && isempty(s.error)
        note = sprintf('epochs %d \x2192 %d', s.epochBefore, s.epochAfter);
    end
    fprintf(fid, '%-4d  %-35s  %7.2f  %9d  %9d  %s\n', ...
        k, s.step, s.duration_s, s.chanBefore, s.chanAfter, note);
end
fprintf(fid, '\nTotal steps: %d\n', numel(stepLog));
errSteps = sum(~cellfun(@isempty, {stepLog.error}));
if errSteps > 0
    fprintf(fid, 'Steps with errors: %d\n', errSteps);
end
fclose(fid);
end
