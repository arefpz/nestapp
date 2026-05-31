function v = convertParam(raw, type)
% CONVERTPARAM Convert a raw parameter value to the correct typed form.
%   v = CONVERTPARAM(raw, type) coerces raw (from UITable, TextArea, or a
%   saved .mat) to the type declared in the step registry.
%
%   Recognised type strings:
%     'scalar'     - single double
%     'integer'    - whole-number double
%     'vector'     - numeric row vector
%     'logical'    - 'on'/'off' string (kept as-is)
%     'string'     - plain char
%     'stringlist' - cell array of chars

switch type
    case {'scalar', 'integer'}
        if isnumeric(raw) && isscalar(raw)
            v = raw;
        elseif isnumeric(raw) && ~isempty(raw)
            v = raw(1);
        elseif iscell(raw) && isscalar(raw)
            v = str2double(string(raw{1}));
        else
            v = str2double(string(raw));
        end

    case 'vector'
        if isnumeric(raw)
            v = raw(:)';
        elseif iscell(raw)
            nums = cellfun(@(c) str2double(string(c)), raw(:));
            if all(~isnan(nums))
                v = nums(:)';
            else
                v = raw;
            end
        elseif ischar(raw) || isstring(raw)
            parsed = str2num(char(raw)); %#ok<ST2NM>
            if ~isempty(parsed)
                v = parsed(:)';
            else
                v = char(raw);
            end
        else
            v = raw;
        end

    case 'logical'
        v = raw;   % keep as 'on'/'off' string

    case 'string'
        if ischar(raw) && isrow(raw)
            v = raw;
        elseif ischar(raw)
            v = strjoin(cellstr(raw), newline);   % char matrix -> newline-joined row
        elseif isstring(raw) && isscalar(raw)
            v = char(raw);
        elseif iscell(raw)
            v = strjoin(raw(:)', newline);         % cell of lines -> newline-joined row
        elseif isnumeric(raw) && isempty(raw)
            v = '';
        else
            v = char(string(raw));
        end

    case 'stringlist'
        if iscell(raw)
            v = raw;
        elseif isstring(raw) && ~isscalar(raw)
            % Multi-element string array (e.g. ["Fp1" "Fp2"]) - char() would
            % produce a char matrix that strsplit cannot handle.
            v = cellstr(raw);
        elseif isnumeric(raw) && isempty(raw)
            v = {};
        elseif (ischar(raw) || isstring(raw)) && strcmp(char(raw), '[]')
            v = {};   % old MATLAB "not set" sentinel
        elseif ischar(raw) || isstring(raw)
            % Scalar char or single string, possibly comma-separated.
            parts = strtrim(strsplit(char(raw), ','));
            parts = parts(~cellfun(@isempty, parts));
            v = parts;
        else
            v = {char(string(raw))};
        end

    otherwise
        v = raw;
end
end
