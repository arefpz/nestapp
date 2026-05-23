function summary = aggregateGateVerdicts(reports)
% AGGREGATEGATEVERDICTS  Build a files x gates verdict matrix from reports.
%   summary = AGGREGATEGATEVERDICTS(reports)
%
%   reports : cell array of pipeline report structs (as returned by
%             runPipelineCore / exportReport / loaded from .mat). Each
%             report may carry report.quality.gates{*}.
%
%   summary fields:
%     .files     cellstr of file basenames, one per report (rows of the
%                  matrix). Empty inputFile becomes '<unknown>'.
%     .gates     cellstr of unique gate labels (columns, first-seen
%                  order across the batch).
%     .verdicts  numeric matrix (numel(files) x numel(gates)) using
%                  the integer codes:
%                    0 = NotChecked, 1 = Pass, 2 = Marginal,
%                    3 = Fail,       4 = Pending
%     .counts    struct(Pass, Marginal, Fail, Pending) - totals across
%                  the entire matrix (NotChecked is excluded).
%
%   Notes
%   - Gates with the same label but different stepIndex collapse to one
%     column. If a single file ran two gates with the same label the
%     later gate wins (worst-of severity is preserved by the worstOf
%     local helper).
%   - Files without report.quality.gates produce an all-NotChecked row.
%   - Reports that are completely empty or non-struct are skipped.

summary = struct('files', {{}}, 'gates', {{}}, 'verdicts', [], ...
    'counts', struct('Pass', 0, 'Marginal', 0, 'Fail', 0, 'Pending', 0));

if isempty(reports), return, end

% First pass: collect file names and unique gate labels (first-seen order).
files = {};
gates = {};
gateIdx = containers.Map();   % label -> column index
for ri = 1:numel(reports)
    r = reports{ri};
    if ~isstruct(r) || isempty(fieldnames(r)), continue, end
    files{end+1} = reportBasename(r); %#ok<AGROW>
    if ~isfield(r, 'quality') || ~isfield(r.quality, 'gates'), continue, end
    for gi = 1:numel(r.quality.gates)
        label = gateLabel(r.quality.gates{gi});
        if ~isKey(gateIdx, label)
            gateIdx(label) = numel(gates) + 1;
            gates{end+1} = label; %#ok<AGROW>
        end
    end
end

nFiles = numel(files);
nGates = numel(gates);
verdicts = zeros(nFiles, nGates);    % 0 = NotChecked

% Second pass: fill the matrix.
fileRow = 0;
for ri = 1:numel(reports)
    r = reports{ri};
    if ~isstruct(r) || isempty(fieldnames(r)), continue, end
    fileRow = fileRow + 1;
    if ~isfield(r, 'quality') || ~isfield(r.quality, 'gates'), continue, end
    for gi = 1:numel(r.quality.gates)
        g = r.quality.gates{gi};
        col = gateIdx(gateLabel(g));
        v   = verdictCode(g.verdict);
        % Worst-of wins when the same file has duplicate gate labels.
        if v > verdicts(fileRow, col)
            verdicts(fileRow, col) = v;
        end
    end
end

counts = struct('Pass', sum(verdicts(:) == 1), ...
                'Marginal', sum(verdicts(:) == 2), ...
                'Fail',     sum(verdicts(:) == 3), ...
                'Pending',  sum(verdicts(:) == 4));

summary.files    = files;
summary.gates    = gates;
summary.verdicts = verdicts;
summary.counts   = counts;
end

% -- small helpers --------------------------------------------------------

function name = reportBasename(r)
if isfield(r, 'inputFile') && ~isempty(r.inputFile)
    [~, name] = fileparts(r.inputFile);
else
    name = '<unknown>';
end
end

function label = gateLabel(g)
if isfield(g, 'label') && ~isempty(g.label)
    label = char(g.label);
else
    label = 'gate';
end
end

function code = verdictCode(verdict)
switch verdict
    case 'Pass',     code = 1;
    case 'Marginal', code = 2;
    case 'Fail',     code = 3;
    case 'Pending',  code = 4;
    otherwise,       code = 0;   % NotChecked
end
end
