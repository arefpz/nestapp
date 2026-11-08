% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef c_FigurePrinter < handle
% C_FIGUREPRINTER  nestapp's replacement for AARATEP's figure-export helper.
%
%   AARATEP writes seven quality-control PNGs during a whole-pipeline run
%   (early channel rejection, eye-ICA components, decay fit and removal,
%   classified components, the final timtopo, and two debug plots). Every one
%   goes through c_FigurePrinter.copyToFile, so replacing this one class fixes
%   all of them at once.
%
%   THIS SHADOWS third_party/aaratep/Common/c_FigurePrinter.m. ensureAaratepOnPath
%   puts this directory on the path AFTER AARATEP's, so it wins. Fixing it here
%   rather than in the vendored copy is deliberate: third_party/aaratep is
%   git-cloned and git-ignored, so any edit there is lost the next time it is
%   re-fetched.
%
%   Three things were wrong with the original, and all three abort a cleaning
%   run that has otherwise succeeded:
%
%   1. A FAILED SCREENSHOT KILLED THE RUN. copyToFile let its exception escape,
%      so a QC image that could not be written took the whole pipeline with it -
%      after the actual cleaning step had already completed. Worse, it threw
%      between the plot and the close() that follows it, orphaning the figure on
%      screen. Here an export failure is a warning: the image is skipped, the
%      pipeline continues, and the caller's close() still runs.
%
%   2. THE REAL ERROR WAS DISCARDED. The original went through export_fig, whose
%      print2array does "catch ex; err = true; end" without ever assigning its
%      output. MATLAB then reports "Output argument A ... not assigned" and the
%      genuine cause is thrown away. exportgraphics is core MATLAB, needs no
%      Ghostscript, and reports its own failures honestly.
%
%   3. initialize() CLEARED GLOBALS MID-RUN. It called CopyFileToClipboard,
%      which calls javaaddpath, which "clears global and persistent variables".
%      Firing that partway through a batch, while EEGLAB globals are live, is
%      not something a screenshot should ever do. initialize() here does
%      nothing, which also stops AARATEP's export_fig from being prepended back
%      onto the path at run time.
%
%   See also: ensureAaratepOnPath, exportgraphics, publicationFigure

    methods (Static)

        function initialize()
        % Deliberately empty. See note 3 above.
        end

        function copyToFile(varargin)
        % COPYTOFILE  Write a figure to file. AARATEP's calling convention.
        %   copyToFile(filepath, 'magnification',N, 'doCrop',tf,
        %              'doTransparent',tf, 'hf',handle)
            p = inputParser();
            p.addRequired('filepath', @ischar);
            p.addParameter('magnification', 1, @isscalar);
            p.addParameter('doCrop', true, @islogical);
            p.addParameter('doTransparent', true, @islogical);
            p.addParameter('hf', []);
            p.parse(varargin{:});
            s = p.Results;

            if isempty(s.hf); s.hf = gcf; end
            if ~isgraphics(s.hf); return; end

            % A quality-control image is never worth failing a run for. Anything
            % that goes wrong below is reported and stepped over, so the caller
            % reaches its close() and the pipeline keeps going.
            try
                writeFigure(s);
            catch ME
                warning('nestapp:qcFigureExportFailed', ...
                    'Could not write the QC figure %s: %s', s.filepath, ME.message);
            end
        end

        function copyToClipboard(varargin) %#ok<VANUS>
        % Not used by the pipeline; a no-op so nothing errors if it is reached.
        end

        function copyMultipleToClipboard(varargin) %#ok<VANUS>
        % Used by c_fig_arrange's interactive shortcuts, not by the pipeline.
        end

        function copyMonitorScreenshotToClipboard(varargin) %#ok<VANUS>
        end

    end
end

% ── helpers ─────────────────────────────────────────────────────────────────

function writeFigure(s)
% Elements AARATEP tags as non-printing are hidden for the shot and restored
% afterwards, matching the original's behaviour.
hidden = gobjects(0);
tagged = findobj(s.hf, 'Tag', 'c_NonPrinting');
if ~isempty(tagged)
    vis    = get(tagged, 'Visible');
    isOn   = ismember(string(vis), "on");
    hidden = tagged(isOn);
    set(hidden, 'Visible', 'off');
end
restore = onCleanup(@() set(hidden, 'Visible', 'on'));

folder = fileparts(s.filepath);
if ~isempty(folder) && ~isfolder(folder)
    mkdir(folder);
end

% AARATEP creates its own figures, outside appWindow/outputFigure, so they
% would follow a dark desktop theme into the file. Pinned here, at the one
% place they are written out.
pinLightTheme(s.hf);

[~, ~, ext] = fileparts(s.filepath);
if strcmpi(ext, '.svg')
    % exportgraphics does not write SVG; saveas does, as the original did.
    saveas(s.hf, s.filepath);
    return
end

% export_fig's -mN means N times the on-screen size, and a MATLAB figure is
% laid out at 96 px per inch, so N maps to a resolution of 96N dpi.
args = {'Resolution', max(round(96 * s.magnification), 1)};
if s.doTransparent
    args = [args, {'BackgroundColor', 'none'}];
else
    args = [args, {'BackgroundColor', 'white'}];
end

exportgraphics(s.hf, s.filepath, args{:});
end
