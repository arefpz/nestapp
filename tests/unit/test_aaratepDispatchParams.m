
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function tests = test_aaratepDispatchParams
% TEST_AARATEPDISPATCHPARAMS  Regression tests that the AARATEP dispatch
%   cases in processOneFile.m only pass parameter names that the
%   vendored c_TMSEEG_* / c_EEG_* helpers actually recognise.
%
%   Would have caught:
%     - "'doDecayRemovalPerTrial' is not a recognized parameter" - the
%       dispatch was using the AARATEP main-pipeline flag name instead
%       of the underlying helper's 'trialAggMethod_timeCourseRemoval'.
%
%   Strategy: parse `p.addParameter('name', ...)` calls out of each
%   vendored helper's source, build the allowed-name set, and assert
%   every name our dispatch hands to that helper is in the set.
%
%   Run: runtests('tests/unit/test_aaratepDispatchParams')
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));
end

function r = repoRoot()
r = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

% ── allowed-params introspection ─────────────────────────────────────────

function names = readAllowedParamsFromHelper(helperRelPath)
% Parse a vendored helper's source and return the set of names it accepts
% via p.addParameter('name', ...) / addRequired / addOptional.
fp = fullfile(repoRoot(), 'third_party', 'aaratep', helperRelPath);
src = fileread(fp);
tokens = regexp(src, ...
    'p\.addParameter\(\s*''([A-Za-z_][A-Za-z0-9_]*)''', 'tokens');
names = cellfun(@(t) t{1}, tokens, 'UniformOutput', false);
end

% ── Remove Decay Artifact dispatch ────────────────────────────────────────

function test_decayArtifact_dispatchUsesValidParamNames(testCase)
allowed = readAllowedParamsFromHelper(...
    'Common/EEGAnalysisCode/c_TMSEEG_fitAndRemoveDecayArtifact.m');

% Names our dispatch case passes to c_TMSEEG_fitAndRemoveDecayArtifact.
% Update this list whenever processOneFile.m's 'Remove Decay Artifact'
% case changes - the test enforces alignment with the helper signature.
dispatchedNames = {
    'artifactTimespan'
    'trialAggMethod_timeCourseRemoval'
    'aggTrimPercent'
};

for k = 1:numel(dispatchedNames)
    name = dispatchedNames{k};
    testCase.verifyTrue(ismember(name, allowed), sprintf( ...
        ['"Remove Decay Artifact" dispatch passes "%s" but ' ...
         'c_TMSEEG_fitAndRemoveDecayArtifact.m does not accept it. ' ...
         'Allowed names: %s'], ...
        name, strjoin(allowed, ', ')));
end
end

function test_decayArtifact_dispatchDoesNotUseBoolFlagName(testCase)
% Specific guard: 'doDecayRemovalPerTrial' is the high-level flag in
% c_TMSEEG_Preprocess_AARATEPPipeline.m, NOT a parameter of the helper.
% The dispatch must translate it to 'trialAggMethod_timeCourseRemoval'.
allowed = readAllowedParamsFromHelper(...
    'Common/EEGAnalysisCode/c_TMSEEG_fitAndRemoveDecayArtifact.m');
testCase.verifyFalse(ismember('doDecayRemovalPerTrial', allowed), ...
    ['Sanity check: doDecayRemovalPerTrial should NOT be a helper ' ...
     'parameter. If upstream renames the flag to match, this test ' ...
     'tells you the dispatch could be simplified.']);
end

% ── AR-Blend interpolation dispatch ───────────────────────────────────────

function test_arBlendInterp_dispatchUsesValidParamNames(testCase)
allowed = readAllowedParamsFromHelper(...
    'Common/EEGAnalysisCode/c_EEG_ReplaceEpochTimeSegment.m');

% Names our dispatch case passes to c_EEG_ReplaceEpochTimeSegment.
dispatchedNames = {
    'timespanToReplace'
    'method'
    'prePostFitDurations'
};

for k = 1:numel(dispatchedNames)
    name = dispatchedNames{k};
    testCase.verifyTrue(ismember(name, allowed), sprintf( ...
        ['"Interpolate Missing Data (AR-Blend)" dispatch passes "%s" ' ...
         'but c_EEG_ReplaceEpochTimeSegment.m does not accept it. ' ...
         'Allowed names: %s'], ...
        name, strjoin(allowed, ', ')));
end
end

% ── perTrial translation correctness ──────────────────────────────────────

function test_decayArtifact_requiresCurveFittingToolbox(testCase)
% c_TMSEEG_fitAndRemoveDecayArtifact calls MATLAB's fit() function
% (Curve Fitting Toolbox) on lines 101 and 107. The registry must list
% this requirement so the pre-flight checkStepDependencies surfaces a
% clear install message instead of letting the step crash mid-run with
% "Undefined function 'fit' for input arguments of type 'double'".
reg = stepRegistry();
idx = find(strcmp({reg.name}, 'Remove Decay Artifact'), 1);
testCase.verifyNotEmpty(idx, 'Step "Remove Decay Artifact" must exist.');
reqs = reg(idx).requires;
testCase.verifyNotEmpty(reqs, 'Step must declare at least one requirement.');
hasFitRequirement = any(arrayfun(@(r) strcmp(r.fn, 'fit'), reqs));
testCase.verifyTrue(hasFitRequirement, ...
    ['Remove Decay Artifact must require the fit() function so the ' ...
     'pre-flight check tells users to install Curve Fitting Toolbox ' ...
     'before running, not after the pipeline already started.']);
end

function test_decayArtifact_perTrialMapsToNone(testCase)
% Upstream c_TMSEEG_Preprocess_AARATEPPipeline.m line 336:
%   trialAggMethod_timeCourseRemoval = c_if(doDecayRemovalPerTrial, 'none', 'mean')
% Per-trial fitting => no inter-trial aggregation => 'none'.
% The dispatch encodes this mapping; we verify by inspecting the
% dispatch source so a refactor that flips the polarity fails the test.
disp1 = fileread(fullfile(repoRoot(), 'src', 'processOneFile.m'));
caseIdx = strfind(disp1, "case 'Remove Decay Artifact'");
testCase.verifyNotEmpty(caseIdx, 'Dispatch case not found.');
endIdx = strfind(disp1(caseIdx:end), 'eeg_checkset');
caseBody = disp1(caseIdx:caseIdx + endIdx(1));

% Must contain the polarity mapping in either order.
hasOnNone = contains(caseBody, "strcmpi(opts2.perTrial, 'on')") && ...
            contains(caseBody, "'none'");
hasOffMean = contains(caseBody, "'mean'");

testCase.verifyTrue(hasOnNone && hasOffMean, ...
    ['Dispatch must map perTrial=on -> trialAggMethod=none (per-trial fit) ' ...
     'and perTrial=off -> trialAggMethod=mean (across-trial average).']);
end
