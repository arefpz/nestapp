
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function nestLog(label, fmt, varargin)
% NESTLOG  Timestamped debug line to the MATLAB command window.
%   nestLog(label, fmt, ...)  writes:
%     [HH:mm:ss.SSS][LABEL] message
%
%   Safe to call from parallel workers - fprintf in workers is forwarded to
%   the client command window by MATLAB's PCT runtime.
ts = char(datetime('now', 'Format', 'HH:mm:ss.SSS'));
fprintf('[%s][%s] %s\n', ts, label, sprintf(fmt, varargin{:}));
end
