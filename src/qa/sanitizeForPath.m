
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
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
