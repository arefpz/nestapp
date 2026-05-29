function s = varinToStruct(varin)
% VARINTOSTRUCT  Convert a {key, val, key, val, ...} cell to a struct.
%   s = VARINTOSTRUCT(varin) is the inverse of paramsToVarin. Used by
%   dispatch cases that pass the unpacked params to a helper expecting an
%   opts struct rather than positional arguments.

if mod(numel(varin), 2) ~= 0
    error('varinToStruct:OddLength', ...
        'varin must have an even length (key-value pairs).');
end

vars = convertContainedStringsToChars(varin);
s = struct();
for k = 1:2:numel(vars)
    s.(vars{k}) = vars{k+1};
end
end
