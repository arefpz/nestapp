
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef test_renderDashboardPanel < matlab.unittest.TestCase
% TEST_RENDERDASHBOARDPANEL  Smoke test for src/qa/renderDashboardPanel.m
%   The renderer is a side-effect function that paints widgets into a
%   uipanel / uifigure. We can't easily assert on the visual output, so
%   the test verifies:
%     - renders without throwing on a representative report set
%     - parent gains expected children (axes / table / labels)
%     - renders cleanly on an empty report set
%     - renders cleanly when ICA-related metrics are NaN
%     - refresh callback fires once when Refresh button is clicked

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            tc.addTeardown(@() rmpath(genpath(fullfile(root, 'src'))));
        end
    end

    methods (Static, Access = private)
        function fig = makeOffscreenFig()
            fig = uifigure('Visible', 'off', 'Position', [100 100 900 600]);
        end

        function r = makeReport(name, gates)
            r = struct( ...
                'inputFile',   [name '.set'], ...
                'processedAt', datetime('now'), ...
                'quality',     struct('figures', {{}}, ...
                                      'gates', {gates}, ...
                                      'worstVerdict', 'Pass'));
        end

        function g = makeGate(label, verdict, thresholds, metrics)
            g = struct( ...
                'label',      label, ...
                'mode',       'absolute', ...
                'verdict',    verdict, ...
                'reasons',    {{'sample reason'}}, ...
                'metrics',    metrics, ...
                'thresholds', thresholds, ...
                'stepIndex',  1);
        end
    end

    methods (Test)
        function renders_three_reports_without_error(tc)
            fig = test_renderDashboardPanel.makeOffscreenFig();
            tc.addTeardown(@() close(fig, 'force'));

            mr = @test_renderDashboardPanel.makeReport;
            mg = @test_renderDashboardPanel.makeGate;
            r1 = mr('a', {mg('raw','Pass', ...
                struct('maxFlatChans', 5), struct('nFlatChans', 1))});
            r2 = mr('b', {mg('raw','Fail', ...
                struct('maxFlatChans', 5), struct('nFlatChans', 8))});
            r3 = mr('c', {mg('raw','Marginal', ...
                struct('maxFlatChans', 5), struct('nFlatChans', 4))});

            renderDashboardPanel(fig, {r1, r2, r3});

            kids = allchild(fig);
            tc.verifyGreaterThan(numel(kids), 0, ...
                'Dashboard should have painted at least one child');
        end

        function renders_empty_report_set_without_error(tc)
            fig = test_renderDashboardPanel.makeOffscreenFig();
            tc.addTeardown(@() close(fig, 'force'));

            renderDashboardPanel(fig, {});
            kids = allchild(fig);
            tc.verifyGreaterThan(numel(kids), 0, ...
                'Header should still render on empty input');
        end

        function refresh_callback_fires(tc)
            fig = test_renderDashboardPanel.makeOffscreenFig();
            tc.addTeardown(@() close(fig, 'force'));

            counter = struct('n', 0);
            cb = @() incCounter();
            renderDashboardPanel(fig, {}, struct('onRefresh', cb));

            % Find and click the Refresh button.
            btn = findobj(fig, 'Type', 'uibutton', 'Text', 'Refresh');
            tc.verifyNotEmpty(btn);
            % uibutton callbacks expect (src, event)
            btn(1).ButtonPushedFcn(btn(1), struct());
            tc.verifyEqual(counter.n, 1);

            function incCounter()
                counter.n = counter.n + 1;
            end
        end

        function NaN_metric_values_handled(tc)
            fig = test_renderDashboardPanel.makeOffscreenFig();
            tc.addTeardown(@() close(fig, 'force'));

            mr = @test_renderDashboardPanel.makeReport;
            mg = @test_renderDashboardPanel.makeGate;
            r = mr('a', {mg('raw', 'Pass', ...
                struct('maxEMGFraction', 0.3), ...
                struct('emgFraction', NaN))});

            % NaN metrics should be silently filtered; no crash.
            renderDashboardPanel(fig, {r});
            tc.verifyTrue(true);
        end

        function rerender_does_not_double_children(tc)
            fig = test_renderDashboardPanel.makeOffscreenFig();
            tc.addTeardown(@() close(fig, 'force'));

            renderDashboardPanel(fig, {});
            n1 = numel(allchild(fig));
            renderDashboardPanel(fig, {});
            n2 = numel(allchild(fig));
            % Same logical render should produce a similar child count -
            % allow slack for MATLAB's invisible internal handles, but
            % refuse to silently double the visible widget tree.
            tc.verifyLessThan(n2, 2 * max(n1, 1));
        end
    end
end
