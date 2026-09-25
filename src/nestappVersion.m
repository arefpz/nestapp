
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
%     v  char - the semantic version of nestapp, e.g. '2.1.0'.
%
%   This is the single source of truth for the application version. Everything
%   that runs calls it. The three documents that restate it as text -
%   CHANGELOG.md, CITATION.cff and the README badge - are written from it by
%   tools/setVersion.m, and VersionTest fails if any of them disagree. Change
%   the version with setVersion('major'|'minor'|'patch'), not by editing this
%   line. Follows Semantic Versioning 2.0.0 (https://semver.org):
%   MAJOR.MINOR.PATCH.
%
%   NOTE: deliberately NOT named version.m - that would shadow MATLAB's
%   built-in version() function on the path.
%
%   See also: version, CHANGELOG.md

v = '2.2.0';
end
