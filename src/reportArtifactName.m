
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function name = reportArtifactName(report, ext)
% REPORTARTIFACTNAME  Filename for a per-file report artifact (.mat / .pdf / ...).
%   name = REPORTARTIFACTNAME(report, ext) builds "<base>_report.<ext>" when
%   the 'overwriteReports' pref is true, otherwise "<base>_report_<ts>.<ext>"
%   where <ts> is report.processedAt formatted yyyyMMdd_HHmmss.
%
%   Centralises the naming rule shared by exportReport (.mat) and
%   exportFileReportPDF (.pdf) so they cannot drift apart.

[~, baseName] = fileparts(report.inputFile);
if isempty(baseName), baseName = 'pipeline'; end

ext = char(ext);
if ~isempty(ext) && ext(1) == '.', ext = ext(2:end); end

if getpref('nestapp', 'overwriteReports', false)
    name = sprintf('%s_report.%s', baseName, ext);
else
    ts   = string(report.processedAt, 'yyyyMMdd_HHmmss');
    name = sprintf('%s_report_%s.%s', baseName, ts, ext);
end
end
