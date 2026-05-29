function vars = stripVarinKeys(vars, keysToStrip)
% STRIPVARINKEYS  Remove name-value pairs whose key is in keysToStrip.
%   vars = STRIPVARINKEYS(vars, keysToStrip) drops both the key and its
%   following value for any key (case-insensitive) appearing in
%   keysToStrip. Used by dispatch cases that need to filter out
%   algorithm-specific parameters before forwarding to a generic
%   underlying function.
%
%   Example:
%     vars = {'icatype','runica','approach','symm','g','tanh'};
%     vars = stripVarinKeys(vars, {'approach','g'});
%     % vars is now {'icatype','runica'}

if isempty(vars) || isempty(keysToStrip)
    return
end

keep = true(1, numel(vars));
for k = 1:2:numel(vars)
    if any(strcmpi(vars{k}, keysToStrip))
        keep(k:k+1) = false;
    end
end
vars = vars(keep);
end
