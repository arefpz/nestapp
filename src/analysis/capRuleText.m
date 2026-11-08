% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function [status, footer] = capRuleText(montage, mode)
% CAPRULETEXT  What the status line and figure footer say about mixed caps.
%   [status, footer] = CAPRULETEXT(montage, mode) takes res.info.montage from
%   groupCurves and the curve mode ('TEP' | 'GMFP' | 'LMFP').
%
%   status  text to append to the Explore status line, starting with the
%           '  |  ' separator, or '' when every file is on one cap
%   footer  value for the 'caps' entry of an exported figure's footer, or ''
%           when every file is on one cap - so a normal figure is unchanged
%
%   Both rules give something up, so both say what:
%     'exclude'  how many files were left out, and how they fall across groups
%                when that is uneven - lopsided exclusions are how this rule
%                biases a comparison
%     'shared'   how many electrodes remain, and a GMFP warning, because GMFP
%                is a spread ACROSS channels and is not comparable to a
%                full-cap GMFP
%
%   See also: groupCurves, capSummary

status = '';
footer = '';
if ~isequal(fieldOr(montage, 'uniform', true), false); return; end

if strcmp(montage.rule, 'shared')
    n = montage.nChannels;
    status = sprintf('  |  mixed caps: %d shared electrodes', n);
    if strcmpi(mode, 'GMFP')
        status = sprintf('%s  |  GMFP over %d electrodes, not the full cap', status, n);
    end
    dropped = montage.dropped;
    if isempty(dropped)     % the other cap is a superset of the common one
        footer = sprintf('shared electrodes (%d)', n);
    elseif numel(dropped) <= 6
        footer = sprintf('shared electrodes (%d; dropped %s)', n, strjoin(dropped, ' '));
    else
        footer = sprintf('shared electrodes (%d; %d dropped)', n, numel(dropped));
    end
    return
end

nOut = numel(montage.excluded);
status = sprintf('  |  %d file%s excluded (other cap)', nOut, plural(nOut));
byGroup = montage.excludedByGroup;
if numel(unique([byGroup.n])) > 1
    bits = arrayfun(@(b) sprintf('%s %d', b.group, b.n), byGroup, 'UniformOutput', false);
    status = sprintf('%s - uneven: %s', status, strjoin(bits, ', '));
end
status = [status ' - see Caps...'];
footer = sprintf('exclude (%d file%s on another cap left out)', nOut, plural(nOut));
end
