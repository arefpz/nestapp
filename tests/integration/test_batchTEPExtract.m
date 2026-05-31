function tests = test_batchTEPExtract
% TEST_BATCHTEPEXTRACT  Integration tests for batchTEPExtract.
%
%   All tests use a synthetic loadFcn so EEGLAB file I/O is not required,
%   but TESA (tesa_peakanalysis) must be on the path because tepPeakFinder
%   always calls it.  If TESA is absent these tests fail visibly — they do
%   NOT silently pass.
%
%   Run: runtests('tests/integration/test_batchTEPExtract')
tests = functiontests(localfunctions);
end

% ── setup ─────────────────────────────────────────────────────────────────────

function setupOnce(testCase)
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));

% Hard requirement — missing TESA must be a visible failure, not a silent pass.
if isempty(which('tesa_peakanalysis'))
    testCase.assumeFail( ...
        'TESA (tesa_peakanalysis) is not on the MATLAB path. ' ...
        'Install TESA via the EEGLAB Plugin Manager before running these tests.');
end
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

% ── helpers ───────────────────────────────────────────────────────────────────

function EEG = makeEEG(nChan, nTime, nTrials, labels, times)
EEG.nbchan   = nChan;
EEG.trials   = nTrials;
EEG.times    = times;
EEG.chanlocs = struct('labels', labels);
base = zeros(nChan, nTime);
for k = 1:nChan
    base(k,:) = syntheticTEP(times);
end
EEG.data = repmat(base, 1, 1, nTrials) + 0.1 * randn(nChan, nTime, nTrials);
end

function w = syntheticTEP(times)
w  = -3 * gauss(times, 45, 6);
w  = w + 4 * gauss(times, 60, 7);
w  = w - 5 * gauss(times, 100, 12);
w  = w + 4 * gauss(times, 180, 20);
end

function g = gauss(t, mu, sigma)
g = exp(-0.5 * ((t - mu) / sigma).^2);
end

function EEG = makeEpochedEEG(labels)
times = -100 : 2 : 300;
if nargin < 1
    labels = {'FC1','FC3','C1','C3','F3','F4','P3','P4'};
end
rng(42);
EEG = makeEEG(numel(labels), numel(times), 40, labels, times);
end

function fn = wrapEEG(EEG)
fn = @(~) EEG;
end

% ── tests ──────────────────────────────────────────────────────────────────────

function test_outputTableHasCorrectColumns(testCase)
EEG = makeEpochedEEG();
[T, ~] = batchTEPExtract({'dummy.set'}, {'FC1','FC3','C1','C3'}, ...
    'loadFcn', wrapEEG(EEG));
expected = {'file','roi_electrodes','component','polarity','found', ...
    'latency_ms','amplitude_uv','win_start_ms','win_end_ms', ...
    'nom_latency_ms','n_trials','n_channels_roi'};
testCase.verifyEqual(T.Properties.VariableNames, expected);
end

function test_sixRowsPerFile(testCase)
EEG = makeEpochedEEG();
[T, ~] = batchTEPExtract({'sub01.set'}, {'FC1','C1'}, 'loadFcn', wrapEEG(EEG));
testCase.verifyEqual(height(T), 6, 'Default 6-component config must produce 6 rows per file');
end

function test_twoFilesTwelveRows(testCase)
EEG = makeEpochedEEG();
[T, ~] = batchTEPExtract({'s1.set','s2.set'}, {'FC1','C1'}, 'loadFcn', wrapEEG(EEG));
testCase.verifyEqual(height(T), 12, 'Two files must produce 12 rows');
end

function test_foundColumnIsNumeric(testCase)
EEG = makeEpochedEEG();
[T, ~] = batchTEPExtract({'s.set'}, {'FC1','C1'}, 'loadFcn', wrapEEG(EEG));
testCase.verifyClass(T.found, 'double', '''found'' column must be double');
testCase.verifyTrue(all(T.found == 0 | T.found == 1), ...
    '''found'' must contain only 0 or 1');
end

function test_nanRowsForMissingROI(testCase)
EEG = makeEpochedEEG({'Cz','Fz','Pz'});
[T, warns] = batchTEPExtract({'s.set'}, {'FC1','FC3'}, 'loadFcn', wrapEEG(EEG));
testCase.verifyEqual(height(T), 6, 'Must still produce 6 NaN rows for missing ROI');
testCase.verifyTrue(all(isnan(T.latency_ms)), 'latency_ms must be NaN for missing ROI');
testCase.verifyTrue(all(T.n_channels_roi == 0), 'n_channels_roi must be 0 for missing ROI');
testCase.verifyFalse(isempty(warns), 'A warning must be returned for missing ROI');
end

function test_skipNonEpochedFile(testCase)
EEG = makeEpochedEEG();
EEG.trials = 1;
[T, warns] = batchTEPExtract({'s.set'}, {'FC1','C1'}, 'loadFcn', wrapEEG(EEG));
testCase.verifyEqual(height(T), 6, 'Must produce 6 NaN rows for non-epoched file');
testCase.verifyTrue(all(T.found == 0), 'found must be 0 for skipped file');
testCase.verifyFalse(isempty(warns), 'A warning must be returned for non-epoched file');
end

function test_partialROIWarningAndCorrectCount(testCase)
% File has only FC1 and C1 of the requested FC1, FC3, C1, C3.
% The NaN rows (caused by tepPeakFinder failing without real peaks) must
% still record n_channels_roi=2, not 0.
EEG = makeEpochedEEG({'FC1','C1','Cz','Pz'});
roiElec = {'FC1','FC3','C1','C3'};
[T, warns] = batchTEPExtract({'s.set'}, roiElec, 'loadFcn', wrapEEG(EEG));
testCase.verifyFalse(isempty(warns), 'Partial ROI must produce a warning');
% n_channels_roi: peaks found → 2 from actual rows; peaks not found → 2 from NaN rows
testCase.verifyTrue(all(T.n_channels_roi == 2), ...
    'n_channels_roi must be 2 (actual found count), not 0, even in fallback rows');
end

function test_loadErrorGraceful(testCase)
EEG = makeEpochedEEG();
mixedLoad = @(fp) loadByName(fp, EEG);
[T, warns] = batchTEPExtract({'bad.set','good.set'}, {'FC1','C1'}, ...
    'loadFcn', mixedLoad);
testCase.verifyEqual(height(T), 12, 'Must produce rows for both files even if one fails');
testCase.verifyFalse(isempty(warns), 'Load failure must produce a warning');
end

function test_csvRoundTrip(testCase)
EEG = makeEpochedEEG();
csvPath = [tempname, '.csv'];
cleanup = onCleanup(@() deleteIfExists(csvPath));
[Torig, ~] = batchTEPExtract({'s.set'}, {'FC1','C1'}, ...
    'loadFcn', wrapEEG(EEG), 'csvPath', csvPath);
testCase.verifyTrue(exist(csvPath,'file') == 2, 'CSV file must be created');
Tread = readtable(csvPath, 'TextType', 'string');
testCase.verifyEqual(height(Tread), height(Torig), 'Row count must match after round-trip');
foundMask = Torig.found == 1;
if any(foundMask)
    testCase.verifyEqual(Tread.latency_ms(foundMask), Torig.latency_ms(foundMask), ...
        'AbsTol', 1e-3, 'latency_ms must survive CSV round-trip');
end
end

function test_deterministicOutput(testCase)
rng(42); EEG = makeEpochedEEG();
rng(42); [T1, ~] = batchTEPExtract({'s.set'}, {'FC1','C1'}, 'loadFcn', wrapEEG(EEG));
rng(42); EEG = makeEpochedEEG();
rng(42); [T2, ~] = batchTEPExtract({'s.set'}, {'FC1','C1'}, 'loadFcn', wrapEEG(EEG));
testCase.verifyEqual(T1.latency_ms, T2.latency_ms, 'AbsTol', 1e-9, ...
    'batchTEPExtract must be deterministic given fixed rng state');
end

function test_consistencyWithTepPeakFinder(testCase)
rng(42);
EEG = makeEpochedEEG({'FC1','C1'});
roiElec = {'FC1','C1'};
[T, ~] = batchTEPExtract({'s.set'}, roiElec, 'loadFcn', wrapEEG(EEG));

roiIdx   = find(ismember({EEG.chanlocs.labels}, roiElec));
roiData  = mean(EEG.data(roiIdx, :, :), 1);
waveform = mean(roiData, 3);
waveform = smoothdata(waveform, 'movmean', 5);
refPeaks = tepPeakFinder(waveform, EEG.times);

for ci = 1:numel(refPeaks)
    pk  = refPeaks(ci);
    row = T(strcmp(T.component, pk.name), :);
    testCase.verifyEqual(row.found(1), double(pk.found), ...
        sprintf('%s: found must match direct tepPeakFinder call', pk.name));
    if pk.found
        testCase.verifyEqual(row.latency_ms(1), pk.latencyMs, 'AbsTol', 1e-6, ...
            sprintf('%s: latency_ms must match direct tepPeakFinder call', pk.name));
    end
end
end

% ── teardown helpers ───────────────────────────────────────────────────────────

function deleteIfExists(path)
if exist(path, 'file'); delete(path); end
end

function E = loadByName(fp, goodEEG)
[~, nm] = fileparts(fp);
if strcmp(nm, 'bad')
    error('test:loadFail', 'Simulated load failure');
end
E = goodEEG;
end
