
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_versionConsistency
% TEST_VERSIONCONSISTENCY  The app version must agree across its sources.
%
%   nestappVersion() is the single source of truth. This guards that the
%   most recent CHANGELOG.md entry matches it, so a release can't ship with
%   a stale changelog or an un-bumped version. The release git tag must also
%   match (checked by the release process / CI, not here - no git in tests).
%
%   Run: runtests('tests/unit/test_versionConsistency')
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

function test_versionIsSemver(testCase)
v = nestappVersion();
testCase.verifyClass(v, 'char');
testCase.verifyNotEmpty(regexp(v, '^\d+\.\d+\.\d+$', 'once'), ...
    sprintf('nestappVersion() = "%s" is not MAJOR.MINOR.PATCH', v));
end

function test_changelogTopMatchesVersion(testCase)
v   = nestappVersion();
txt = fileread(fullfile(repoRoot(), 'CHANGELOG.md'));
% First "## [x.y.z]" heading (skipping "## [Unreleased]").
tok = regexp(txt, '##\s*\[(\d+\.\d+\.\d+)\]', 'tokens', 'once');
testCase.verifyNotEmpty(tok, 'No versioned section found in CHANGELOG.md');
testCase.verifyEqual(tok{1}, v, sprintf( ...
    ['CHANGELOG.md top release is %s but nestappVersion() is %s. ' ...
     'Bump one to match before releasing.'], tok{1}, v));
end
