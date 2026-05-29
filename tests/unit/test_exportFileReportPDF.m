
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
classdef test_exportFileReportPDF < matlab.unittest.TestCase
% TEST_EXPORTFILEREPORTPDF  Smoke test for src/qa/exportFileReportPDF.m
%   PDF rendering is hard to assert on visually, so the test verifies:
%     - a PDF appears at the expected path
%     - file size grows when QC PNG pages are added (proxy for "multi-page")
%     - missing-PNG case still produces a PDF (placeholder page)
%     - a re-run on the same file overwrites cleanly (no stale append)

    methods (TestClassSetup)
        function addSrcPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            tc.addTeardown(@() rmpath(genpath(fullfile(root, 'src'))));
        end
    end

    methods (Static, Access = private)
        function r = makeReport(inputFile, figures)
            r = initPipelineReport(inputFile);
            r.processedAt = datetime(2026, 5, 23, 12, 0, 0);
            r.channels.original = 16;
            r.channels.final    = 16;
            r.quality.figures   = figures;
        end

        function p = writePNG(dir, name)
            p = fullfile(dir, name);
            imwrite(uint8(255 * rand(64, 64, 3)), p);
        end
    end

    methods (Test)
        function text_only_PDF_is_produced(tc)
            tmpDir = fullfile(tempdir, ['pdf_text_', char(matlab.lang.internal.uuid())]);
            mkdir(tmpDir);
            tc.addTeardown(@() rmdir(tmpDir, 's'));

            r = test_exportFileReportPDF.makeReport( ...
                fullfile(tmpDir, 'sub.set'), {});
            pdfPath = exportFileReportPDF(r, tmpDir);

            tc.verifyTrue(exist(pdfPath, 'file') == 2);
            d = dir(pdfPath);
            tc.verifyGreaterThan(d.bytes, 1000, ...
                'Text-only PDF should be > 1 KB');
        end

        function PDF_grows_when_PNG_pages_added(tc)
            tmpDir = fullfile(tempdir, ['pdf_grow_', char(matlab.lang.internal.uuid())]);
            mkdir(tmpDir);
            tc.addTeardown(@() rmdir(tmpDir, 's'));

            png1 = test_exportFileReportPDF.writePNG(tmpDir, 'p1.png');
            png2 = test_exportFileReportPDF.writePNG(tmpDir, 'p2.png');

            rTextOnly = test_exportFileReportPDF.makeReport( ...
                fullfile(tmpDir, 'sub.set'), {});
            pTextOnly = exportFileReportPDF(rTextOnly, tmpDir);
            sizeTextOnly = dir(pTextOnly).bytes;
            delete(pTextOnly);

            rWithImages = test_exportFileReportPDF.makeReport( ...
                fullfile(tmpDir, 'sub.set'), {png1, png2});
            pWithImages = exportFileReportPDF(rWithImages, tmpDir);
            sizeWithImages = dir(pWithImages).bytes;

            tc.verifyGreaterThan(sizeWithImages, sizeTextOnly, ...
                'Adding PNG pages should grow the PDF.');
        end

        function missing_PNG_still_produces_PDF(tc)
            tmpDir = fullfile(tempdir, ['pdf_miss_', char(matlab.lang.internal.uuid())]);
            mkdir(tmpDir);
            tc.addTeardown(@() rmdir(tmpDir, 's'));

            % Reference a PNG that doesn't exist; renderImagePage falls
            % back to a placeholder page rather than throwing.
            phantom = fullfile(tmpDir, 'does_not_exist.png');
            r = test_exportFileReportPDF.makeReport( ...
                fullfile(tmpDir, 'sub.set'), {phantom});
            pdfPath = exportFileReportPDF(r, tmpDir);

            tc.verifyTrue(exist(pdfPath, 'file') == 2);
        end

        function rerun_overwrites_existing_pdf(tc)
            tmpDir = fullfile(tempdir, ['pdf_over_', char(matlab.lang.internal.uuid())]);
            mkdir(tmpDir);
            tc.addTeardown(@() rmdir(tmpDir, 's'));

            png1 = test_exportFileReportPDF.writePNG(tmpDir, 'p1.png');
            r = test_exportFileReportPDF.makeReport( ...
                fullfile(tmpDir, 'sub.set'), {png1});

            pdf1 = exportFileReportPDF(r, tmpDir);
            size1 = dir(pdf1).bytes;
            pdf2 = exportFileReportPDF(r, tmpDir);
            size2 = dir(pdf2).bytes;

            tc.verifyEqual(pdf1, pdf2, ...
                'Second run should target the same filename (timestamp matches report.processedAt)');
            % Sizes should be similar - certainly not 2x (which would
            % indicate a stale-append).
            tc.verifyLessThan(size2, 1.5 * size1, ...
                'Re-run should overwrite, not double the PDF size.');
        end
    end
end
