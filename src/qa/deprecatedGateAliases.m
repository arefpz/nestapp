
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function aliases = deprecatedGateAliases()
% DEPRECATEDGATEALIASES  Single source of truth for renamed Quality
%   Gate parameter keys. Returns an Nx2 cell array of {oldKey, newKey}
%   pairs. Used by:
%     - qualityGate.applyDefaults  - silent per-call aliasing
%     - runPipelineCore.warnDeprecatedGateParams  - one-time CFG log
%   When you rename a Quality Gate param, add the (old, new) row here
%   and both consumers pick it up automatically.
aliases = { ...
    'maxBadChanPct',        'maxOutlierChanPct'; ...
    'maxBadTrialPct',       'maxOutlierTrialPct'; ...
    'maxBadChanPctWarnAt',  'maxOutlierChanPctWarnAt'; ...
    'maxBadTrialPctWarnAt', 'maxOutlierTrialPctWarnAt'};
end
