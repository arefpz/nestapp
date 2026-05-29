
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function nestLog(label, fmt, varargin)
% NESTLOG  Timestamped debug line to the MATLAB command window.
%   nestLog(label, fmt, ...)  writes:
%     [HH:mm:ss.SSS][LABEL] message
%
%   When a debug log file is open (see nestDebugLog), the same line is also
%   appended to that file, capturing a run's full trace for bug reports.
%
%   Safe to call from parallel workers - fprintf in workers is forwarded to
%   the client command window by MATLAB's PCT runtime. The file tee captures
%   client-side logging (the serial path, where step-by-step debugging lives).
global NESTAPP_DEBUG_FID %#ok<GVMIS>
ts   = char(datetime('now', 'Format', 'HH:mm:ss.SSS'));
line = sprintf('[%s][%s] %s', ts, label, sprintf(fmt, varargin{:}));
fprintf('%s\n', line);
if ~isempty(NESTAPP_DEBUG_FID)
    try
        fprintf(NESTAPP_DEBUG_FID, '%s\n', line);
    catch
        % Never let logging-to-file break the run.
    end
end
end
