function s = sanitizeForPath(name)
% SANITIZEFORPATH  Strip non-filename characters from a step name.
%   s = SANITIZEFORPATH(name) returns name with anything outside
%   [A-Za-z0-9_] removed. Used to turn pipeline step names into safe
%   PNG filename fragments.
%
%   Example:
%     sanitizeForPath('Remove ICA Components (TESA)')
%       -> 'RemoveICAComponentsTESA'
s = regexprep(name, '[^A-Za-z0-9_]', '');
end
