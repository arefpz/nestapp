
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function logPath = nestDebugLog(action, folder)
% NESTDEBUGLOG  Start/stop teeing nestLog output to a debug log file.
%
%   Syntax
%     logPath = nestDebugLog('start', folder)  % open a log file in folder
%     nestDebugLog('stop')                      % close the log file
%     fid = nestDebugLog('fid')                 % current fid (-1 if none)
%
%   Inputs
%     action  char - 'start', 'stop', or 'fid'.
%     folder  char - directory for the log file (created if needed); 'start' only.
%
%   Outputs
%     logPath  char - the path of the opened log file ('start'), '' otherwise.
%
%   When a log file is open, nestLog writes every line to both the command
%   window and the file, so a run's full step-by-step trace (per-step EEG
%   dimensions, timings, citations, QC notes) is captured for bug reports.
%   Controlled by runPipelineCore via the 'debugLog' preference. The fid is
%   held in a global so the central nestLog chokepoint can see it.
%
%   Side effects   Opens/closes a file; sets the global NESTAPP_DEBUG_FID.
%   See also: nestLog, runPipelineCore

global NESTAPP_DEBUG_FID %#ok<GVMIS>
logPath = '';

switch lower(action)
    case 'start'
        nestDebugLog('stop');   % close any previous handle first
        if nargin < 2 || isempty(folder); folder = pwd; end
        if ~exist(folder, 'dir'); mkdir(folder); end
        stamp   = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        logPath = fullfile(folder, sprintf('nestapp_debug_%s.log', stamp));
        fid = fopen(logPath, 'a');
        if fid < 0
            warning('nestDebugLog:open', 'Could not open debug log at %s', logPath);
            NESTAPP_DEBUG_FID = [];
            logPath = '';
            return
        end
        NESTAPP_DEBUG_FID = fid;
        fprintf(fid, '=== nestapp debug log - %s ===\n', ...
            char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

    case 'stop'
        if ~isempty(NESTAPP_DEBUG_FID)
            try
                if NESTAPP_DEBUG_FID > 2   % never close std streams
                    fclose(NESTAPP_DEBUG_FID);
                end
            catch
            end
        end
        NESTAPP_DEBUG_FID = [];

    case 'fid'
        if isempty(NESTAPP_DEBUG_FID)
            logPath = -1;   % returned via logPath for convenience
        else
            logPath = NESTAPP_DEBUG_FID;
        end

    otherwise
        error('nestDebugLog:badAction', ...
            'Unknown action "%s" (use start, stop, or fid).', action);
end
end
