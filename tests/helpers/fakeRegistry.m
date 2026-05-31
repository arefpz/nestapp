function steps = fakeRegistry(extraSteps)
% FAKEREGISTRY  Return a minimal two-step registry for unit tests.
%
%   steps = FAKEREGISTRY() returns a 2-element struct array shaped like the
%   real stepRegistry() output, containing 'Load Data' and 'Save New Set'
%   stubs with one parameter each. Use this in tests that need a registry-
%   shaped object but have no interest in specific step metadata.
%
%   steps = FAKEREGISTRY(extraSteps) appends the additional step structs
%   in extraSteps to the two built-in stubs. Each element of extraSteps
%   must at minimum have a .name field; missing fields are filled with defaults.
%
%   See also: stepRegistry
    arguments
        extraSteps struct = struct.empty(1,0)
    end

    steps = [makeStub('Load Data',    'filepath', 'string'), ...
             makeStub('Save New Set', 'suffix',   'string')];

    for k = 1:numel(extraSteps)
        steps(end+1) = fillDefaults(extraSteps(k)); %#ok<AGROW>
    end
end

% ── helpers ──────────────────────────────────────────────────────────────────

function s = makeStub(name, paramKey, paramType)
    s.name     = name;
    s.defaults = struct(paramKey, '');
    s.info     = '';
    s.requires = {};
    s.params   = makeParam(paramKey, paramType);
end

function p = makeParam(key, type)
    p.key          = key;
    p.friendlyName = key;
    p.unit         = '';
    p.validRange   = [];
    p.description  = '';
    p.placeholder  = '';
    p.type         = type;
end

function s = fillDefaults(s)
    if ~isfield(s, 'defaults'); s.defaults = struct(); end
    if ~isfield(s, 'info');     s.info     = '';       end
    if ~isfield(s, 'requires'); s.requires = {};       end
    if ~isfield(s, 'params');   s.params   = makeParam('value', 'scalar'); end
end
