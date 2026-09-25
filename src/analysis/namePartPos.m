% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function pos = namePartPos(idx, K, fromEnd)
% NAMEPARTPOS  Where rule index idx falls in a name of K parts.
%   pos = NAMEPARTPOS(idx, K, fromEnd) is idx itself when counting from the
%   start, and K - idx + 1 when counting from the end, so index 1 from the
%   end is the last part. The mapping is its own inverse: applied to a chip's
%   position in a name it gives the index a rule stores for that chip.
%
%   The one place "from the end" is defined; everything that reads or writes
%   a rule's nameIdx goes through here.
%
%   See also: subjectIdsFromRule, suggestSubjectRule, subjectRuleDialog

pos = idx;
if fromEnd; pos = K - idx + 1; end
end
