function tests = test_solver
%TEST_SOLVER  Interface and bookkeeping of solve_burgers.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
oldPath = path;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
testCase.addTeardown(@path, oldPath);
end

function testGrid(testCase)
[x, u] = solve_burgers(@sin, 'L', 4, 'N', 40, 'tFinal', 0.1);
verifyEqual(testCase, size(x), [40 1]);
verifyEqual(testCase, size(u), [40 1]);
verifyEqual(testCase, x(1), 0);
verifyEqual(testCase, x(end), 4 - 4/40, 'AbsTol', 1e-14);
end

function testZeroFinalTimeReturnsInitialCondition(testCase)
[x, u, info] = solve_burgers(@sin, 'N', 50, 'tFinal', 0);
verifyEqual(testCase, u, sin(x));
verifyEqual(testCase, info.steps, 0);
verifyEqual(testCase, info.dtMin, 0);
end

function testSnapshotsHitRequestedTimes(testCase)
[~, u, info] = solve_burgers(@sin, 'N', 80, 'tFinal', 0.6, 'tSave', [0 0.2 0.45], ...
    'scheme', 'godunov');
verifyEqual(testCase, info.t, [0; 0.2; 0.45; 0.6], 'AbsTol', 1e-14);
verifyEqual(testCase, size(info.U), [80 4]);
verifyEqual(testCase, info.U(:, end), u);
verifyEqual(testCase, info.U(:, 1), sin((0:79)' * 2*pi/80), 'AbsTol', 1e-14);
end

function testInfoBookkeeping(testCase)
[~, ~, info] = solve_burgers(@sin, 'N', 100, 'tFinal', 0.3, 'scheme', 'maccormack');
verifyEqual(testCase, info.scheme, 'maccormack');
verifyGreaterThan(testCase, info.steps, 0);
verifyLessThanOrEqual(testCase, info.dtMin, info.dtMax);
verifyEqual(testCase, info.newtonIters, 0);
end

function testMaccormackMatchesCharacteristicsBeforeShock(testCase)
[x, u] = solve_burgers(@sin, 'N', 400, 'tFinal', 0.5, 'scheme', 'maccormack');
e = error_norms(u, burgers_characteristics(x, 0.5, @sin, 1), 2*pi/400);
verifyLessThan(testCase, e.L1, 5e-4);
end

function testErrorNormsOfKnownDifference(testCase)
e = error_norms([1; 2; 3; 4], [1; 1; 1; 1], 0.5);
verifyEqual(testCase, e.L1, 0.5 * (0 + 1 + 2 + 3), 'AbsTol', 1e-14);
verifyEqual(testCase, e.L2, sqrt(0.5 * (0 + 1 + 4 + 9)), 'AbsTol', 1e-14);
verifyEqual(testCase, e.Linf, 3);
end

function testNegativeViscosityRejected(testCase)
verifyError(testCase, @() solve_burgers(@sin, 'nu', -1), ?MException);
end

function testTooFewPointsRejected(testCase)
verifyError(testCase, @() solve_burgers(@sin, 'N', 3), ?MException);
end

function testSnapshotTimeBeyondFinalTimeRejected(testCase)
verifyError(testCase, @() solve_burgers(@sin, 'tFinal', 1, 'tSave', 2), 'burgers:invalidInput');
end

function testUnknownSchemeRejected(testCase)
verifyError(testCase, @() solve_burgers(@sin, 'scheme', 'rk99'), 'burgers:unknownScheme');
end

function testBadInitialConditionRejected(testCase)
verifyError(testCase, @() solve_burgers(@(x) 1, 'N', 20), 'burgers:invalidInitialCondition');
verifyError(testCase, @() solve_burgers(@(x) NaN(size(x)), 'N', 20), ...
    'burgers:invalidInitialCondition');
end
