
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function s = makePipelineStep(name, registry)
% MAKEPIPELINESTEP Construct a typed pipeline step from registry defaults.
%   s = MAKEPIPELINESTEP(name, registry) returns a struct:
%     s.name    - step name (must match registry and switch/case in dispatchStep)
%     s.params  - struct of typed parameter values initialised from defaults
%
%   Throws nestapp:unknownStep if name is not found in the registry.

idx = find(strcmp({registry.name}, name), 1);
if isempty(idx)
    error('nestapp:unknownStep', 'Unknown pipeline step: "%s"', name);
end
s.name   = name;
s.params = registry(idx).defaults;
end
