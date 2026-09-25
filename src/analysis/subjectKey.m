% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function key = subjectKey(id)
% SUBJECTKEY  The form two subject ids are compared in: case and zero padding ignored.
%   key = SUBJECTKEY(id) lower-cases id and drops leading zeros from every run
%   of digits, so 's01', 'S01' and 'S1' all give 's1'. Takes a char or a
%   cellstr.
%
%   Within one cohort, 'S1' and 'S01' are the same person typed two ways, not
%   two people - nobody numbers participants so that they differ only in
%   padding or case. Treating them as different splits one person into two,
%   which inflates n exactly as silently as a wrong merge deflates it.
%
%   See also: subjectIdsFromRule, suggestSubjectRule

key = regexprep(lower(id), '(?<!\d)0+(?=\d)', '');
end
