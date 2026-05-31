
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function v = nestappVersion()
% NESTAPPVERSION  Return the nestapp application version string.
%
%   Syntax
%     v = nestappVersion()
%
%   Outputs
%     v  char - the semantic version of nestapp, e.g. '1.0.0'.
%
%   This is the single source of truth for the application version. The
%   About dialog, README, CHANGELOG, and the release git tag must all
%   agree with the value returned here (the version-consistency CI check
%   enforces this). Follows Semantic Versioning 2.0.0 (https://semver.org):
%   MAJOR.MINOR.PATCH.
%
%   NOTE: deliberately NOT named version.m - that would shadow MATLAB's
%   built-in version() function on the path.
%
%   See also: version, CHANGELOG.md

v = '1.0.0';
end
