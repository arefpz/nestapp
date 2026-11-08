% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef SuiteHygieneTest < NestappTestCase
% SUITEHYGIENETEST  The suite's rules about itself, enforced rather than written down.
%
%   Every rule here corresponds to a way the previous suite decayed. That suite
%   was not written by people who thought duplication was fine; it drifted one
%   reasonable-looking file at a time, and nothing was watching. A convention
%   document would not have stopped it - 50 verbatim copies of one 3-line
%   helper accumulated while the duplication was entirely obvious.
%
%   So the conventions are executable. Each test names the specific decay it
%   prevents, and fails with the list of offending files.
%
%   WHY EACH TEST LOOPS INSTEAD OF USING TestParameter. One test per rule, not
%   one per file: the FACT being asserted is the rule ("no test rolls its own
%   temp dir"), and the files are merely its domain. Parameterising over ~30
%   files x 8 rules would report 240 results for 8 facts, which is precisely
%   the inflation this rewrite exists to undo. The failure message lists every
%   offender, so nothing is lost diagnostically - unlike the old
%   characterization test, which collapsed ten genuinely distinct behaviours
%   into one verifyEmpty and there the granularity really was wrong.
%
%   THE SUITE IS READ ONCE, in TestClassSetup. An earlier version read and
%   comment-stripped every file inside every rule, which came to ELEVEN full
%   passes over the tree per run - a cost that grew with each rule added, not
%   merely with each file. Reading once makes it one.
%
%   EVERY RULE IS ONE SHAPE: verifyNoneOf(noteFcn, why). Three rules used to
%   hand-roll the accumulate-and-verify loop instead of using the helper the
%   other five went through, which is one more way for a ninth rule to be
%   written inconsistently with the rest.
%
%   These rules apply to the NEW suite only. The legacy folders are exempt
%   until they are deleted; they violate nearly all of it, which is the point.

    properties (Constant)
        % The four folders the rules govern. Legacy folders are absent by
        % design - see the class comment.
        SuiteFolders = {'pure', 'eeglab', 'gui', 'eeglab_gui'}

        % Files permitted to read a file as text, each with a reason.
        %
        % RegistryContractTest and MethodsClauseTest read nestapp SOURCE, which
        % is the practice this rewrite is otherwise eliminating. They are
        % allowed because both DERIVE their expected set from the registry and
        % check the code against it - the opposite of the 19 scrapes being
        % deleted, which hardcoded a source fragment and checked the source
        % still contained it, so they failed on rewording and caught nothing.
        %
        % StepGoldenTest reads recorded JSON, not source. That is data the
        % suite owns, so it carries none of the brittleness above - but it
        % still reads a file as text, and the rule cannot tell a .json target
        % from a .m one without parsing the call. Listing it keeps the rule
        % honest rather than loosening it: every text read in the suite is
        % accounted for by name.
        MayReadSource = {'RegistryContractTest', 'MethodsClauseTest', ...
                         'StepGoldenTest', 'AaratepInstallTest', 'VersionTest', ...
                         ... WindowThemeTest: a window created without appWindow or
                         ... outputFigure brings back the dark-mode mix, and there is
                         ... no behavioural way to enumerate every dialog.
                         'WindowThemeTest'}
    end

    properties (Access = private)
        % One entry per suite file: .path .name .folder .code (comments
        % stripped), .effective (.code plus the code of every tests/helpers
        % function it calls - see helperClosure).
        Files
    end

    methods (TestClassSetup)
        function readTheSuiteOnce(tc)
        % Runs after NestappTestCase's own TestClassSetup has fixed the path.
            tc.Files = struct('path', {}, 'name', {}, 'folder', {}, ...
                              'code', {}, 'effective', {});
            root = fullfile(addNestappPath(), 'tests');
            helpers = readHelpers(fullfile(root, 'helpers'));

            % Every suite folder goes on the path for the duration, because
            % everyTestClassInheritsTheBase resolves classes by reflection and
            % meta.class.fromName cannot see a class MATLAB cannot load. Only
            % the folder being RUN is on the path otherwise, so the rule would
            % report every file in the other three as not inheriting the base -
            % which it did, the moment tests/gui gained its first file.
            %
            % A PathFixture rather than addpath: it restores itself, and it is
            % not the bare addpath call this file's own rules ban elsewhere.
            present = {};
            for k = 1:numel(tc.SuiteFolders)
                d = fullfile(root, tc.SuiteFolders{k});
                if isfolder(d); present{end+1} = d; end %#ok<AGROW>
            end
            if ~isempty(present)
                tc.applyFixture(matlab.unittest.fixtures.PathFixture(present));
            end
            for k = 1:numel(tc.SuiteFolders)
                folder = tc.SuiteFolders{k};
                d = fullfile(root, folder);
                if ~isfolder(d); continue; end
                e = dir(fullfile(d, '*.m'));
                for i = 1:numel(e)
                    % This file is exempt from its own scan. It necessarily
                    % CONTAINS every pattern it forbids - the bans are written
                    % here as regex literals - so including it would make each
                    % rule report itself. Found the moment the suite first ran:
                    % noTestRollsItsOwnTempDir failed on its own
                    % '(tempname|tempdir)' pattern. The alternative, parsing
                    % MATLAB well enough to ignore string literals, is far more
                    % machinery than these rules are worth.
                    if strcmp(e(i).name, 'SuiteHygieneTest.m'); continue; end
                    p = fullfile(d, e(i).name);
                    [~, name] = fileparts(e(i).name);
                    code = stripComments(p);
                    tc.Files(end+1) = struct('path', p, 'name', name, ...
                                             'folder', folder, 'code', code, ...
                                             'effective', helperClosure(code, helpers));
                end
            end
        end
    end

    methods (Test)

        function noTestRollsItsOwnPathSetup(tc)
        % 50 verbatim copies of repoRoot() and 107 fileparts^3 chains, which
        % had already drifted: 49 files used a bare addpath(src) and 56 used
        % genpath, so whether src/qa was reachable depended on which file you
        % ran. NestappTestCase does this once.
            tc.verifyNoMatch('(^|\s)function\s+\w+\s*=\s*(repoRoot|srcRoot)\b', ...
                'defines its own repoRoot/srcRoot - inherit NestappTestCase instead');
            tc.verifyNoMatch('(^|[^.\w])addpath\s*\(', ...
                'calls addpath directly - NestappTestCase does that');
            tc.verifyNoMatch('fileparts\s*\(\s*fileparts\s*\(', ...
                'walks up to the repo root by hand - use addNestappPath');
        end

        function noTestRollsItsOwnTempDir(tc)
        % ~30 creation sites in five idioms, of which the old sweep matched
        % four; one of its five hardcoded prefixes matched nothing at all.
            tc.verifyNoMatch('(^|[^.\w])(tempname|tempdir)\s*[\(;,)]', ...
                'builds its own temp dir - use scratchDir(tc)');
        end

        function noTestDecidesItCannotRun(tc)
        % 38 assume* sites, twelve of which skipped on a function shipped in
        % THIS repo - turning a path bug into a green run. The folder declares
        % what a test needs; run_tests checks it once, up front.
            tc.verifyNoMatch('\.assume(Fail|NotEmpty|True|False|Equal|NotEqual)\s*\(', ...
                'skips itself at runtime - the folder already declares its needs');
        end

        function noTestAssertsSomethingThatCannotFail(tc)
        % The old suite contained a literal verifyTrue(true) as a test's only
        % assertion, plus file-size floors standing in for "it rendered".
            tc.verifyNoMatch('verify(True|False)\s*\(\s*(true|false)\s*\)', ...
                'asserts a literal - that cannot fail');
        end

        function everyTestClassInheritsTheBase(tc)
        % Asked of the class hierarchy MATLAB actually resolved, not of the
        % source text: matching '< NestappTestCase' is whitespace-sensitive and
        % blind to a class that inherits through an intermediate base.
            tc.verifyNoneOf(@(f) noteWhen(~inheritsBase(f.name)), ...
                'does not extend NestappTestCase');
        end

        function sourceScrapingIsOptIn(tc)
        % Behavioural over textual, with a named exception list. The old suite
        % had 19 scrapes; test_dispatchCoverage even documents why it stopped
        % ("passes or fails on how the code is WRITTEN rather than what it
        % DOES") while two other files went on scraping the same file.
            allowed = tc.MayReadSource;
            tc.verifyNoneOf(@(f) noteWhen(~ismember(f.name, allowed) && ...
                    ~isempty(regexp(f.code, '(fileread|readlines)\s*\(', 'once'))), ...
                ['reads a file as text; prefer a behavioural test, or add ' ...
                 'it to SuiteHygieneTest.MayReadSource with a reason']);
        end

        function everyTestIsInTheFolderItsDependenciesRequire(tc)
        % The load-bearing rule. The old taxonomy was drawn on "does it call an
        % EEGLAB function" rather than on what a test NEEDS, which left 11
        % files in the documented no-EEGLAB-no-GUI suite needing one or the
        % other, and made its stated contract false.
            tc.verifyNoneOf(@folderNote, 'is in the wrong folder');
        end

        function everySuiteFolderIsTrackable(tc)
        % .gitignore ignores every directory by default and whitelists test
        % folders one at a time, so a folder that is not listed is SILENTLY
        % untracked - git add reports nothing and the commit looks clean. Not
        % hypothetical: tests/pure was created, two test files were written and
        % run green, and both were omitted from the commit without a word.
        %
        % The postmortem comment left in .gitignore asks a human to remember.
        % This is that same rule, enforced - which is this file's whole
        % argument, applied to the one place it had been left as prose.
        %
        % ASKED OF GIT, NOT OF THE FILE. The first version of this rule read
        % .gitignore and looked for the literal '!tests/<folder>/*.m', which
        % asserts how the rule is SPELLED rather than what it DOES - the same
        % source-scraping mistake the old suite made nineteen times, and one
        % that a negation elsewhere in the file, a different but equivalent
        % pattern, or a second .gitignore deeper in the tree would each defeat
        % while the text still matched. check-ignore is the resolver git itself
        % uses when deciding what to add, so it is the only answer that cannot
        % be wrong.
        % EVERY TRACKED FOLDER UNDER tests/, not the four suite folders. The
        % first version checked only those, and the cutover walked straight
        % through the gap it left: tests/golden/ was created by a git mv, so
        % its ten files stayed tracked because they were already in the index -
        % but a NEW golden would have been ignored silently, and this rule was
        % green. Derived from what git actually tracks, so a folder added later
        % is covered without anyone remembering to list it here.
            missing = {};
            folders = trackedTestFolders();
            tc.assertNotEmpty(folders, 'git reports no tracked folders under tests/');
            for k = 1:numel(folders)
                rel = sprintf('%s/ZZHygieneProbe%s', folders{k}, probeExt(folders{k}));
                if tc.gitIgnores(rel)
                    missing{end+1} = rel; %#ok<AGROW>
                end
            end
            tc.verifyEmpty(missing, sprintf( ...
                ['git would not add a new file in %d tracked tests/ folder(s), ' ...
                 'so anything put there is silently untracked. Whitelist in ' ...
                 '.gitignore: %s'], numel(missing), strjoin(missing, ', ')));

        % Positive control, in the same test because it is the same fact: a
        % folder nobody whitelisted MUST come back ignored. Without it a probe
        % that silently stopped resolving - a git that failed to run, a path
        % git could not read - would report every folder safe, which is how
        % three tests in this rewrite passed while asserting nothing.
            tc.verifyTrue(tc.gitIgnores('tests/not_a_suite_folder/ZZProbe.m'), ...
                ['the ignore probe itself is not working: an unlisted folder ' ...
                 'came back tracked, so the check above proves nothing']);
        end

        function everyHelperHasACaller(tc)
        % A shared fixture with no uptake is worse than none - it is a second
        % thing to keep in step, and nobody knows it is dead. Not theoretical:
        % charFixture had ZERO callers in the old unit/ while nine local EEG
        % builders lived beside it, and fakeRegistry had one caller against
        % four private reimplementations.
        %
        % EVERY HELPER, not a list of names. This rule used to name two
        % ('fakeEeg', 'fakeGroupResult'), which is the same list-of-names
        % failure mode that ledger row A3.10 is about - a hardcoded list only
        % covers what someone remembered to add. It cost exactly what that
        % costs: at the cutover, FIVE helpers turned out to have no callers
        % left (assumeDesktop, driveModalDialog, fakeRegistry, hideFromPath,
        % isolateRoiPresets), and this rule was green throughout.
        %
        % Deleting an unused helper is not a loss - git has it, and `git show`
        % brings it back the day a test needs it. What is a loss is a helpers/
        % folder where a reader cannot tell which files are live.
            helpers = readHelpers(fullfile(addNestappPath(), 'tests', 'helpers'));
            % Everything that may count as a caller EXCEPT the helpers
            % themselves - those are added per-helper below, minus the one
            % being judged. A helper's own file necessarily contains
            % "function <name>(", which matches the call pattern, so including
            % it makes every helper its own caller and the rule vacuous. It
            % did exactly that in the first version, and only planting a
            % deliberate orphan showed it.
            fixed = strjoin([{tc.Files.code}, ...
                {stripComments(fullfile(addNestappPath(), 'tests', 'run_tests.m'))}], ...
                newline);

            orphans = {};
            for k = 1:numel(helpers)
                nm = helpers(k).name;
                % A call, or a superclass clause - NestappTestCase is only ever
                % named by "< NestappTestCase", never called.
                % Concatenated, not sprintf'd: sprintf eats \w and \s as format
                % escapes and hands regexp a mangled pattern that matches
                % nothing - which reported all 19 helpers as orphans.
                pattern = ['(^|[^.\w])' nm '\s*[({]|<\s*' nm '\>'];
                % Every OTHER helper, plus the fixed corpus above.
                others = strjoin({helpers([1:k-1, k+1:end]).code}, newline);
                if isempty(regexp([others newline fixed], pattern, 'once'))
                    orphans{end+1} = nm; %#ok<AGROW>
                end
            end
            tc.verifyEmpty(orphans, sprintf( ...
                ['%d helper(s) have no callers: %s. Either delete them (git ' ...
                 'keeps them) or they are the wrong abstraction.'], ...
                numel(orphans), strjoin(orphans, ', ')));
        end
    end

    % ── the one shape every rule uses ────────────────────────────────────────
    methods (Access = private)

        function tf = gitIgnores(tc, relPath)
        % True when git would refuse to add relPath. Exit 0 means ignored, 1
        % means not; anything else is git failing to answer, which must be a
        % failure and not a quiet false either way.
            [status, out] = system(sprintf('git -C "%s" check-ignore -q "%s"', ...
                                           addNestappPath(), relPath));
            tc.assertTrue(ismember(status, [0 1]), sprintf( ...
                'git check-ignore could not answer for %s (exit %d): %s', ...
                relPath, status, strtrim(out)));
            tf = status == 0;
        end

        function verifyNoneOf(tc, noteFcn, why)
        % noteFcn(fileEntry) returns '' when the file is fine, or text to
        % append to its name when it is not.
            bad = {};
            for i = 1:numel(tc.Files)
                note = noteFcn(tc.Files(i));
                if isempty(note); continue; end
                bad{end+1} = [shortPath(tc.Files(i).path) note]; %#ok<AGROW>
            end
            tc.verifyEmpty(bad, listing(bad, why));
        end

        function verifyNoMatch(tc, pattern, why)
            tc.verifyNoneOf(@(f) noteWhen(~isempty(regexp(f.code, pattern, 'once'))), why);
        end
    end
end

% ── file-level helpers ───────────────────────────────────────────────────────

function code = stripComments(p)
% Not cosmetic: an earlier version of this scan matched the licence header
% "% Part of nestapp; see the LICENSE..." in every single file and classified
% the whole suite as needing a display.
lines = strtrim(strsplit(fileread(p), newline));
lines = lines(~startsWith(lines, '%'));
code  = strjoin(regexprep(lines, '%.*$', ''), newline);
end

function tf = inheritsBase(className)
% Walks the real superclass list. Falls back to false when the class cannot be
% resolved at all, which is itself worth failing on - a suite file MATLAB
% cannot load is not a test.
tf = false;
mc = meta.class.fromName(className);
if isempty(mc); return; end
seen = {};
queue = {mc};
while ~isempty(queue)
    c = queue{1}; queue(1) = [];
    if ismember(c.Name, seen); continue; end
    seen{end+1} = c.Name; %#ok<AGROW>
    if strcmp(c.Name, 'NestappTestCase'); tf = true; return; end
    queue = [queue, num2cell(c.SuperclassList')]; %#ok<AGROW>
end
end

function note = noteWhen(isBad)
% The common case: a rule that only needs to say yes or no. Non-empty means
% "offending"; verifyNoneOf appends it to the file name, and a single space
% appends nothing visible.
if isBad; note = ' '; else; note = ''; end
end

function note = folderNote(f)
% A folder may be MORE permissive than the file needs (a pure test parked in
% eeglab/ is wasteful, not wrong), but never less: that is the failure that
% makes a suite's contract a lie.
%
% READS f.effective, NOT f.code, and that distinction is the whole rule. A
% requirement moved behind a shared helper disappears from the calling file's
% text while remaining just as real - which happened here: extracting
% startEeglab and saveFixtureSet took the last 'eeglab(' and 'pop_saveset('
% out of StepGoldenTest, and a text-only version of this rule went from
% correctly identifying it as needing EEGLAB to seeing a file it would happily
% have accepted in pure/. Since using helpers is the thing this suite is built
% to encourage, a rule that goes blind exactly when a helper is used would
% decay to nothing on its own.
needsEeglab = ~isempty(regexp(f.effective, ['(^|[^.\w])(pop_[a-z]\w*|eeg_checkset|' ...
    'eeg_emptyset|eeglab|topoplot|convertlocs|readlocs)\s*\('], 'once'));
needsDisplay = ~isempty(regexp(f.effective, ['(^|[^.\w])(uifigure|figure|uiaxes|axes|' ...
    'subplot|uipanel|exportapp|nestapp)\s*\('], 'once'));

if      needsEeglab &&  needsDisplay; want = 'eeglab_gui';
elseif  needsEeglab;                  want = 'eeglab';
elseif  needsDisplay;                 want = 'gui';
else;                                 want = 'pure';
end

% Explicit membership rather than a rank comparison, because gui and eeglab are
% SIBLINGS, not ordered - a rank would happily accept a display-needing test in
% eeglab/, where there is no display.
switch want
    case 'pure';       ok = true;
    case 'gui';        ok = ismember(f.folder, {'gui', 'eeglab_gui'});
    case 'eeglab';     ok = ismember(f.folder, {'eeglab', 'eeglab_gui'});
    case 'eeglab_gui'; ok = strcmp(f.folder, 'eeglab_gui');
end

note = '';
if ~ok
    note = sprintf(' (in %s/, needs %s/)', f.folder, want);
end
end

function folders = trackedTestFolders()
% Every distinct folder under tests/ that git currently tracks a file in.
[status, out] = system(sprintf('git -C "%s" ls-files tests', addNestappPath()));
if status ~= 0 || isempty(strtrim(out)); folders = {}; return; end
paths   = strtrim(strsplit(strtrim(out), newline));
folders = unique(cellfun(@(p) fileparts(p), paths, 'UniformOutput', false));
folders = folders(~cellfun('isempty', folders));
folders = strrep(folders, filesep, '/');
end

function ext = probeExt(folder)
% Ask about the kind of file that folder holds: the goldens are whitelisted as
% *.json and everything else as *.m, so probing with the wrong extension would
% report a false positive.
if endsWith(folder, 'golden'); ext = '.json'; else; ext = '.m'; end
end

function h = readHelpers(dir_)
% name -> stripped source, for every function in tests/helpers.
h = struct('name', {}, 'code', {});
if ~isfolder(dir_); return; end
e = dir(fullfile(dir_, '*.m'));
for i = 1:numel(e)
    [~, name] = fileparts(e(i).name);
    h(end+1) = struct('name', name, 'code', stripComments(fullfile(dir_, e(i).name))); %#ok<AGROW>
end
end

function txt = helperClosure(code, helpers)
% A file's own code plus the code of every tests/helpers function it names.
%
% ONE LEVEL DEEP, deliberately, and the same choice drawSourceGraph makes for
% the same reason: helpers call each other at most once here (saveFixtureSet
% calls charFixture and scratchDir), and a transitive closure would need cycle
% handling to buy a depth nothing yet uses. If a helper ever hides a
% requirement two levels down, this is where to deepen it - but the honest
% statement is that it is one level, not that it is complete.
parts = {code};
for k = 1:numel(helpers)
    if ~isempty(regexp(code, ['(^|[^.\w])' helpers(k).name '\s*\('], 'once'))
        parts{end+1} = helpers(k).code; %#ok<AGROW>
    end
end
txt = strjoin(parts, newline);
end

function s = shortPath(p)
parts = strsplit(strrep(p, filesep, '/'), '/tests/');
s = parts{end};
end

function msg = listing(items, why)
if isempty(items); msg = ''; return; end
msg = sprintf('%d file(s) %s:\n  %s', numel(items), why, ...
              strjoin(items, sprintf('\n  ')));
end
