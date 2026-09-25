% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function names = methodsSilentSteps()
% METHODSSILENTSTEPS  Registry steps that deliberately add no methods sentence.
%   names = METHODSSILENTSTEPS() lists the steps whose absence from the
%   methods paragraph is a decision, not an omission: data I/O, event
%   detection, labelling without removal, analysis outputs and QC gates. None
%   of them changes the preprocessed data.
%
%   Everything NOT listed here must produce a clause in methodsClause. A step
%   that changes the data and has no case there gets a generic sentence rather
%   than silence, and MethodsProseTest fails for any registry step that falls
%   through to it. The list is the contract that keeps a new step from quietly
%   dropping out of the paragraph, which is how RANSAC, Robust Detrend and
%   Manual Command once went missing.
%
%   Steps that are silent only in some configurations (Load Channel Location
%   without coordCheck=remove, Visualize EEG Data without its reject flag)
%   are NOT listed: their cases in methodsClause decide per parameter set.
%
%   See also: methodsClause, pipelineProse

names = { ...
    'Load Data', ...
    'Save New Set', ...
    'Choose Data Set', ...
    'Find TMS Pulses (TESA)', ...
    'Find TMS Pulses (AARATEP)', ...
    'Label ICA Components', ...        % labels only; removal is its own step
    'Extract TEP (TESA)', ...
    'Find TEP Peaks (TESA)', ...
    'TEP Peak Output', ...
    'Quality Gate'};
end
