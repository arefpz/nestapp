function tests = test_nestappCallbacks
% TEST_NESTAPPCALLBACKS  Unit tests for GUI callback helpers extracted from @nestapp.
%
%   Tests buildParamTableData and applyParamEdit — the pure-function cores
%   of the refreshParamTable and UITableCellEdit callbacks — without
%   instantiating the App Designer app.
%
%   Run: runtests('tests/unit/test_nestappCallbacks')
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
r = repoRoot();
addpath(r);
addpath(fullfile(r, 'src'));
end

function r = repoRoot()
r = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

% ══════════════════════════════════════════════════════════════════════════
%% buildParamTableData — column structure
% ══════════════════════════════════════════════════════════════════════════

function test_buildParamTableData_rowCountMatchesParamCount(testCase)
reg   = fakeRegistry();
step  = makeStep('Load Data', struct('filepath', ''));
data  = buildParamTableData(step, reg(1));
testCase.verifyEqual(size(data, 1), numel(reg(1).params));
testCase.verifyEqual(size(data, 2), 2);
end

function test_buildParamTableData_labelNoUnit(testCase)
reg  = fakeRegistry();
step = makeStep('Load Data', struct('filepath', 'test.set'));
data = buildParamTableData(step, reg(1));
testCase.verifyEqual(data{1,1}, 'filepath');
end

function test_buildParamTableData_labelWithUnit(testCase)
reg = fakeRegistry(makeRegistryStub('Resample', 'freq', 'Hz'));
step = makeStep('Resample', struct('freq', 1000));
data = buildParamTableData(step, reg(3));
testCase.verifyEqual(data{1,1}, 'freq (Hz)');
end

function test_buildParamTableData_valueSetShowsFormatted(testCase)
reg  = fakeRegistry();
step = makeStep('Load Data', struct('filepath', 'sub01.set'));
data = buildParamTableData(step, reg(1));
testCase.verifyEqual(data{1,2}, 'sub01.set');
end

function test_buildParamTableData_missingParamShowsNotSet(testCase)
reg  = fakeRegistry();
step = makeStep('Load Data', struct());   % filepath absent from params
data = buildParamTableData(step, reg(1));
testCase.verifyEqual(data{1,2}, '(not set)');
end

function test_buildParamTableData_placeholderShownWhenEmpty(testCase)
stub = makeRegistryStub('MyStep', 'thr', '');
stub.params.placeholder = 'e.g. 0.8';
reg  = fakeRegistry(stub);
step = makeStep('MyStep', struct());
data = buildParamTableData(step, reg(3));
testCase.verifyEqual(data{1,2}, 'e.g. 0.8');
end

function test_buildParamTableData_numericVectorFormatted(testCase)
reg  = fakeRegistry();
step = makeStep('Load Data', struct('filepath', [1 2 3]));
data = buildParamTableData(step, reg(1));
testCase.verifyEqual(data{1,2}, '[1 2 3]');
end

function test_buildParamTableData_cellArrayJoined(testCase)
reg  = fakeRegistry();
step = makeStep('Load Data', struct('filepath', {{'Cz','Pz','Fz'}}));
data = buildParamTableData(step, reg(1));
testCase.verifyEqual(data{1,2}, 'Cz, Pz, Fz');
end

% ══════════════════════════════════════════════════════════════════════════
%% applyParamEdit — spec mutation
% ══════════════════════════════════════════════════════════════════════════

function test_applyParamEdit_updatesCorrectKey(testCase)
reg    = fakeRegistry();
spec   = makeStep('Load Data', struct('filepath', ''));
spec   = applyParamEdit(spec, 1, 1, '/data/sub01.set', reg(1));
testCase.verifyEqual(spec(1).params.filepath, '/data/sub01.set');
end

function test_applyParamEdit_parsesNumericStringToScalar(testCase)
stub = makeRegistryStub('Resample', 'freq', 'Hz');
stub.params.type = 'scalar';          % unit is 'Hz'; type must be set explicitly
reg  = fakeRegistry(stub);
spec(1).name   = 'Resample';
spec(1).params = struct('freq', 0);
spec = applyParamEdit(spec, 1, 1, '1000', reg(3));
testCase.verifyEqual(spec(1).params.freq, 1000);
testCase.verifyTrue(isnumeric(spec(1).params.freq));
end

function test_applyParamEdit_passesthroughStringValue(testCase)
reg  = fakeRegistry();
spec = makeStep('Save New Set', struct('suffix', ''));
spec = applyParamEdit(spec, 1, 1, '_cleaned', reg(2));
testCase.verifyEqual(spec(1).params.suffix, '_cleaned');
end

function test_applyParamEdit_noopWhenRowExceedsParams(testCase)
reg    = fakeRegistry();
spec   = makeStep('Load Data', struct('filepath', 'original.set'));
before = spec;
spec   = applyParamEdit(spec, 1, 99, 'anything', reg(1));   % row 99 > 1 param
testCase.verifyEqual(spec(1).params.filepath, before(1).params.filepath);
end

function test_applyParamEdit_doesNotMutateOtherSteps(testCase)
reg     = fakeRegistry();
spec(1) = makeStep('Load Data',    struct('filepath', 'a.set'));
spec(2) = makeStep('Save New Set', struct('suffix',   '_orig'));
spec    = applyParamEdit(spec, 1, 1, 'b.set', reg(1));
testCase.verifyEqual(spec(2).params.suffix, '_orig');
end

% ══════════════════════════════════════════════════════════════════════════
%% helpers
% ══════════════════════════════════════════════════════════════════════════

function step = makeStep(name, params)
step.name   = name;
step.params = params;
end

function stub = makeRegistryStub(name, paramKey, unitOrType)
% makeRegistryStub  Build a single registry entry for test data.
%   unitOrType is used as both unit and type for simplicity; caller can
%   override individual fields after the fact.
stub.name     = name;
stub.defaults = struct(paramKey, []);
stub.info     = '';
stub.requires = {};
p.key          = paramKey;
p.friendlyName = paramKey;
p.unit         = unitOrType;
p.validRange   = [];
p.description  = '';
p.placeholder  = '';
p.type         = 'string';
stub.params    = p;
end
