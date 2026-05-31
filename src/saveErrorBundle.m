
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function bundleDir = saveErrorBundle(targetDir, ctx)
% SAVEERRORBUNDLE  Write a metadata-only debug bundle when a step fails.
%
%   Syntax
%     bundleDir = saveErrorBundle(targetDir, ctx)
%
%   Inputs
%     targetDir  char   - parent directory for the bundle. A timestamped
%                         subfolder is created inside it. If empty, tempdir.
%     ctx        struct - fields (all optional except err):
%                  .err        MException that was caught (required)
%                  .EEG        the EEG struct at failure (METADATA only is
%                              written - never EEG.data)
%                  .spec       the pipeline spec
%                  .stepName   name of the failing step
%                  .stepIndex  index of the failing step
%                  .fileName   input file being processed
%                  .pipelineName
%
%   Outputs
%     bundleDir  char - the folder written, or '' on failure.
%
%   The bundle is safe to share: it contains the error + stack, the
%   environment (nestappDoctor), the pipeline (describePipeline), and EEG
%   METADATA ONLY (dimensions, srate, channel labels, history) - never the
%   raw recording. Raw data is intentionally excluded for privacy; a user
%   who wants to share data must do so deliberately and separately.
%
%   Side effects   Creates a folder and writes text/Markdown files.
%   See also: nestappDoctor, describePipeline, processOneFile

bundleDir = '';
if nargin < 2 || ~isfield(ctx, 'err') || isempty(ctx.err)
    return
end
ctx = fillDefaults(ctx);
if isempty(targetDir); targetDir = tempdir; end

stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
bundleDir = fullfile(targetDir, sprintf('error_%s', stamp));
try
    if ~exist(bundleDir, 'dir'); mkdir(bundleDir); end
    writeText(fullfile(bundleDir, 'error.txt'),        errorText(ctx));
    writeText(fullfile(bundleDir, 'eeg_metadata.txt'), eegMetaText(ctx.EEG));
    writeText(fullfile(bundleDir, 'environment.md'),   safeCall(@() nestappDoctor('Quiet', true)));
    writeText(fullfile(bundleDir, 'pipeline.md'),      safeCall(@() describePipeline(ctx.spec, 'Quiet', true)));
catch ME
    warning('saveErrorBundle:write', 'Could not write error bundle: %s', ME.message);
    bundleDir = '';
    return
end

try
    nestLog('DEBUG', 'Saved error bundle (metadata only): %s', bundleDir);
catch
end
end

% ── helpers ───────────────────────────────────────────────────────────────────

function ctx = fillDefaults(ctx)
defaults = struct('EEG', struct(), 'spec', struct('name',{},'params',{}), ...
    'stepName', '', 'stepIndex', 0, 'fileName', '', 'pipelineName', '');
f = fieldnames(defaults);
for i = 1:numel(f)
    if ~isfield(ctx, f{i}) || isempty(ctx.(f{i}))
        ctx.(f{i}) = defaults.(f{i});
    end
end
end

function s = errorText(ctx)
L = {};
L{end+1} = 'nestapp error bundle';
L{end+1} = sprintf('Time:      %s', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
L{end+1} = sprintf('Pipeline:  %s', ctx.pipelineName);
L{end+1} = sprintf('File:      %s', ctx.fileName);
L{end+1} = sprintf('Step:      %d (%s)', ctx.stepIndex, ctx.stepName);
L{end+1} = sprintf('Identifier:%s', ctx.err.identifier);
L{end+1} = sprintf('Message:   %s', ctx.err.message);
L{end+1} = '';
L{end+1} = 'Stack:';
try
    L{end+1} = getReport(ctx.err, 'extended', 'hyperlinks', 'off');
catch
    for k = 1:numel(ctx.err.stack)
        L{end+1} = sprintf('  %s (line %d)', ctx.err.stack(k).name, ctx.err.stack(k).line); %#ok<AGROW>
    end
end
s = strjoin(L, newline);
end

function s = eegMetaText(EEG)
% METADATA ONLY - never writes EEG.data.
L = {'EEG metadata (no raw data)'};
if ~isstruct(EEG) || isempty(fieldnames(EEG))
    L{end+1} = '(no EEG available at failure)';
    s = strjoin(L, newline);
    return
end
L{end+1} = sprintf('nbchan:   %s', numOrNA(EEG, 'nbchan'));
L{end+1} = sprintf('trials:   %s', numOrNA(EEG, 'trials'));
L{end+1} = sprintf('pnts:     %s', numOrNA(EEG, 'pnts'));
L{end+1} = sprintf('srate:    %s', numOrNA(EEG, 'srate'));
L{end+1} = sprintf('xmin:     %s', numOrNA(EEG, 'xmin'));
L{end+1} = sprintf('xmax:     %s', numOrNA(EEG, 'xmax'));
if isfield(EEG, 'data')
    L{end+1} = sprintf('data size:%s', mat2str(size(EEG.data)));   % shape only
end
if isfield(EEG, 'ref');     L{end+1} = sprintf('ref:      %s', toStr(EEG.ref)); end
if isfield(EEG, 'setname'); L{end+1} = sprintf('setname:  %s', toStr(EEG.setname)); end
if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs, 'labels')
    L{end+1} = sprintf('channels: %s', strjoin({EEG.chanlocs.labels}, ' '));
end
if isfield(EEG, 'icaweights')
    L{end+1} = sprintf('has ICA:  %d', ~isempty(EEG.icaweights));
end
if isfield(EEG, 'history') && ~isempty(EEG.history)
    L{end+1} = '';
    L{end+1} = 'EEG.history:';
    L{end+1} = char(EEG.history);
end
s = strjoin(L, newline);
end

function out = numOrNA(EEG, f)
if isfield(EEG, f) && ~isempty(EEG.(f)) && isnumeric(EEG.(f))
    out = num2str(EEG.(f));
else
    out = 'n/a';
end
end

function s = toStr(v)
if ischar(v); s = v; elseif isnumeric(v); s = mat2str(v); else; s = class(v); end
end

function out = safeCall(fn)
try
    out = fn();
catch ME
    out = sprintf('(could not collect: %s)', ME.message);
end
end

function writeText(fpath, text)
fid = fopen(fpath, 'w');
if fid < 0; return; end
fwrite(fid, text);
fclose(fid);
end
