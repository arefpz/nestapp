% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef SelectDataTreeTest < NestappTestCase
% SELECTDATATREETEST  The data browser opens once, and one window at a time.
%
%   Reported on Linux: picking data opened a confusing number of windows.
%   The browser blocks only its caller, so a second click on Browse... (or
%   File > Select Data, a recent-data item, Add group...) nested another
%   browser - and with no remembered folder each one also opened the system
%   folder picker at once.
%
%   The browser is modal-by-waitfor, so every test here arms a FAILSAFE timer
%   that deletes any open browser after a few seconds: a regression fails
%   the test instead of leaving MATLAB waiting for a human. uigetdir is
%   shadowed for the same reason - a real system picker cannot be closed
%   from code.
%
%   In tests/gui: it creates real windows, but needs no EEGLAB and no app.

    properties (Constant)
        NAME = 'Select Data Files'
        KEY  = 'nestappOpenDataBrowser'
    end

    methods (TestMethodSetup)
        function safety(tc)
            % The real picker must never open; if it were called it would
            % answer with this folder, which the tests can then detect.
            shadowFunction(tc, 'uigetdir', {tempdir});
            t = timer('StartDelay', 8, 'TimerFcn', @(~,~) closeBrowsers(tc.NAME));
            start(t);
            tc.addTeardown(@() stopAndDelete(t));
            tc.addTeardown(@() closeBrowsers(tc.NAME));
            tc.addTeardown(@() clearRegistration(tc.KEY));
        end
    end

    methods (Test)

        function aSecondOpenBringsTheFirstForwardInsteadOfNesting(tc)
            first = appWindow('Name', tc.NAME, 'Visible', 'off');
            tc.addTeardown(@() delete(first));
            setappdata(groot, tc.KEY, first);   % as if a browser were open

            t0 = tic;
            got = selectDataTree(scratchDir(tc));
            tc.verifyLessThan(toc(t0), 3, 'returned at once, without blocking');
            tc.verifyEmpty(got, 'a refused open reads as cancelled');
            tc.verifyNumElements(findall(groot, 'Type', 'figure', 'Name', tc.NAME), 1, ...
                'still exactly one browser');
        end

        function withNoRememberedFolderItWaitsForChoose(tc)
        % It used to open the system folder picker at the same moment as
        % itself. Now it opens empty and points at Choose...
            report = struct();
            t = timer('StartDelay', 1.5, 'TimerFcn', @(~,~) inspectThenClose());
            start(t);
            tc.addTeardown(@() stopAndDelete(t));
            got = selectDataTree('');
            tc.verifyEmpty(got);
            tc.verifyTrue(report.seen, 'the browser opened');
            tc.verifySubstring(report.parent, 'Choose...', ...
                'it asks for a folder instead of choosing one for the user');
            tc.verifyEqual(char(report.weight), 'bold', 'Choose... is highlighted');

            function inspectThenClose()
                f = findall(groot, 'Type', 'figure', 'Name', tc.NAME);
                report.seen = ~isempty(f);
                if ~report.seen; return; end
                fld = findall(f, 'Type', 'uieditfield');
                report.parent = fld(~logical([fld.Editable])).Value;
                report.weight = findall(f, 'Type', 'uibutton', 'Text', 'Choose...').FontWeight;
                delete(f);
            end
        end

        function closingTheBrowserFreesTheNextOpen(tc)
        % The registration must not outlive the window, or the browser could
        % never be opened again after the first time.
            t = timer('StartDelay', 1.5, 'TimerFcn', @(~,~) closeBrowsers(tc.NAME));
            start(t);
            tc.addTeardown(@() stopAndDelete(t));
            selectDataTree(scratchDir(tc));
            tc.verifyFalse(isappdata(groot, tc.KEY), 'unregistered on close');
        end
    end
end

function closeBrowsers(name)
delete(findall(groot, 'Type', 'figure', 'Name', name));
end

function clearRegistration(key)
if isappdata(groot, key); rmappdata(groot, key); end
end

function stopAndDelete(t)
if isvalid(t); stop(t); delete(t); end
end
