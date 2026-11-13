% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [ok, msg] = ensureEeglabReady()
% ENSUREEEGLABREADY  Put EEGLAB and its plugins on the path, once per session.
%   [ok, msg] = ENSUREEEGLABREADY() runs eeglab('nogui') if EEGLAB is not
%   already reachable, and returns whether it is now usable. msg is '' when ok,
%   and a user-facing explanation when not.
%
%   Why this exists: adding the EEGLAB root to the path is NOT enough. EEGLAB
%   ships its functions in subfolders and its plugins in plugins/<name>, and
%   nothing resolves until eeglab() itself runs the path setup and the plugin
%   scan. On a cold MATLAB, with only the root added:
%
%       which('pop_loadset')          -> ''
%       which('pop_tesa_filtbutter')  -> ''
%
%   That matters because availableSteps decides what the picker offers by
%   probing which() for each step's requirements. Asking before EEGLAB has
%   initialised hides every plugin-backed step - 32 of 54 on a stock install -
%   and hides them the same way it hides a genuinely missing plugin, so the
%   user is told nothing is installed when everything is.
%
%   READINESS IS THE SENTINEL, NOT A FLAG. This used to return early whenever
%   `global PLUGINLIST` was non-empty, on the reasoning that PLUGINLIST is what
%   the plugin scan writes and therefore the state we depend on. That is right
%   about whether the scan RAN and wrong about whether its result is still
%   REACHABLE: nothing clears PLUGINLIST, so a restoredefaultpath, a test using
%   hideFromPath, or a pathdef reset left this reporting ready with pop_loadset
%   gone, and the next pop_loadset died with a bare "Undefined function".
%   pathMemo now owns that policy - probe the function we actually need - and
%   it still costs nothing when EEGLAB was initialised by someone else (a prior
%   run, the console, a test), because then the sentinel simply resolves.
%
%   The call is silenced with evalc: eeglab('nogui') prints a banner and a
%   line per plugin, which is noise in a GUI session and drowns the app's own
%   startup logging.
%
%   When EEGLAB is neither on the path nor set in Preferences, findEeglab
%   looks in the usual install locations and the folder it finds is saved as
%   the preference.
%
%   See also: pathMemo, findEeglab, availableSteps, stepAvailability, loadPrefs

% The memo stores the MESSAGE, and ok is derived from it. A cache hit requires
% the sentinel to resolve, so anything read back was computed on a run where
% EEGLAB came up - a stored failure is never returned, because a miss is
% retried. That makes ok exactly isempty(msg) at every point it is read.
msg = pathMemo('pop_loadset', @initialise);
ok  = isempty(msg);
end

% ── helpers ─────────────────────────────────────────────────────────────────

function msg = initialise()
msg = '';

% The Preferences path is where the user points nestapp at their install; it
% is only added if EEGLAB is not already reachable, so an eeglab the user put
% on the path themselves wins over a stale pref.
if isempty(which('eeglab'))
    eeglabPath = '';
    if ispref('nestapp', 'eeglabPath')
        eeglabPath = getpref('nestapp', 'eeglabPath');
    end
    % No usable folder set: look where the installer and a by-hand unzip put
    % EEGLAB, and remember what is found so Preferences shows it.
    if isempty(eeglabPath) || ~isfolder(eeglabPath)
        found = findEeglab();
        if ~isempty(found)
            eeglabPath = found;
            setpref('nestapp', 'eeglabPath', found);
        end
    end
    if ~isempty(eeglabPath) && isfolder(eeglabPath)
        addpath(eeglabPath);
    end
end

if isempty(which('eeglab'))
    msg = ['EEGLAB was not found. Set its folder in Settings > Preferences; ' ...
           'until then, steps that need EEGLAB or its plugins cannot be offered.'];
    return
end

try
    evalc('eeglab nogui');
catch ME
    msg = sprintf(['EEGLAB could not be initialised: %s\nCheck the EEGLAB ' ...
                   'folder in Settings > Preferences.'], ME.message);
end
end
