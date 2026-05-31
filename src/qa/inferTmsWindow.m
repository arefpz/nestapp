function win = inferTmsWindow(EEG, fallback, opts)
% INFERTMSWINDOW  Pick a TMS-pulse exclusion window from EEG.event.
%   win = INFERTMSWINDOW(EEG, fallback)
%   win = INFERTMSWINDOW(EEG, fallback, opts)
%
%   Returns [tStart tEnd] in ms (same axis as EEG.times). Falls back
%   to `fallback` when no matching event is found, when EEG has no
%   events, or when EEG.times is unavailable.
%
%   opts.labels    cell array of event-type patterns. Default:
%                  {'TMS','tms','R128','S128','Stimulus'}
%   opts.decayMs   duration of exclusion past pulse onset (default 25)
%
%   For epoched data (size(EEG.data,3) > 1) the first event in trial 1
%   is used (typical convention: stim is at t = 0 ms). For continuous
%   data the earliest event matching opts.labels is used.

if nargin < 2 || isempty(fallback), fallback = [0 25]; end
if nargin < 3 || ~isstruct(opts),  opts = struct();   end
if ~isfield(opts, 'labels')
    opts.labels = {'TMS','tms','R128','S128','Stimulus'};
end
if ~isfield(opts, 'decayMs'), opts.decayMs = 25; end

win = fallback;

if ~isfield(EEG, 'event') || isempty(EEG.event) ...
        || ~isfield(EEG, 'times') || isempty(EEG.times) ...
        || ~isfield(EEG, 'srate') || isempty(EEG.srate)
    return
end

% Pull event types as a cellstr.
rawTypes = {EEG.event.type};
types = cellfun(@toChar, rawTypes, 'UniformOutput', false);

matchMask = false(size(types));
for k = 1:numel(opts.labels)
    matchMask = matchMask | strcmpi(types, opts.labels{k});
end
hits = find(matchMask, 1, 'first');
if isempty(hits), return, end

ev = EEG.event(hits);
if ~isfield(ev, 'latency') || isempty(ev.latency)
    return
end

% latency is in samples relative to the start of the recording (continuous)
% or the start of the epoch (epoched). Convert to ms using EEG.times when
% possible so we share the time axis the rest of the pipeline uses.
latencyMs = latencyToMs(ev.latency, EEG);
if isnan(latencyMs)
    return
end

win = [latencyMs, latencyMs + opts.decayMs];
end

% -- local helpers --------------------------------------------------------

function s = toChar(v)
if isnumeric(v)
    s = num2str(v);
elseif isstring(v)
    s = char(v);
elseif ischar(v)
    s = v;
else
    s = '';
end
end

function ms = latencyToMs(latency, EEG)
% Map a sample-domain latency to the same ms axis EEG.times uses.
latency = double(latency);
if size(EEG.data, 3) > 1
    % Epoched: latency wraps per-trial; the first trial covers samples
    % 1..pnts. Subtract the offset of trial 1 to get within-epoch sample.
    epochSamp = mod(latency - 1, EEG.pnts) + 1;
    if epochSamp < 1 || epochSamp > numel(EEG.times)
        ms = NaN; return
    end
    ms = EEG.times(round(epochSamp));
else
    if latency < 1 || latency > numel(EEG.times)
        ms = NaN; return
    end
    ms = EEG.times(round(latency));
end
end
