% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef TestProgressLogger < matlab.unittest.plugins.TestRunnerPlugin
% TESTPROGRESSLOGGER  Print each test's start/end to the console and a log file.
%
%   Used by run_tests to pinpoint a hanging or slow test: every test logs a
%   timestamped START before it runs and an END (+ elapsed) after. The log is
%   teed to a file, opened/closed per line so the last line survives even if a
%   test hangs and the session has to be killed - the final START with no
%   matching END is the culprit.
%
%   runner = matlab.unittest.TestRunner.withTextOutput;
%   runner.addPlugin(TestProgressLogger('C:\path\progress.log'));
%
%   See also: run_tests, matlab.unittest.plugins.TestRunnerPlugin

    properties
        LogFile = '';
    end

    methods
        function plugin = TestProgressLogger(logFile)
            if nargin >= 1
                plugin.LogFile = logFile;
            end
        end
    end

    methods (Access = protected)
        function runTest(plugin, pluginData)
            plugin.emit(sprintf('START %s', pluginData.Name));
            tStart = tic;
            runTest@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
            plugin.emit(sprintf('END   %s  (%.2fs)', pluginData.Name, toc(tStart)));
        end
    end

    methods (Access = private)
        function emit(plugin, msg)
            line = sprintf('[%s] %s', ...
                char(datetime('now', 'Format', 'HH:mm:ss.SSS')), msg);
            fprintf('%s\n', line);
            if ~isempty(plugin.LogFile)
                fid = fopen(plugin.LogFile, 'a');
                if fid > 0
                    fprintf(fid, '%s\n', line);
                    fclose(fid);   % flush per line so a hang can't lose it
                end
            end
        end
    end
end
