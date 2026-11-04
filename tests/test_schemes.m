function tests = test_schemes
%TEST_SCHEMES  Properties of the single-step update functions and of the solver runs.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
oldPath = path;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
testCase.addTeardown(@path, oldPath);
end

function testGetSchemeListsAllSchemes(testCase)
names = get_scheme();
verifyEqual(testCase, numel(names), 7);
for k = 1:numel(names)
    step = get_scheme(names{k});
    verifyClass(testCase, step, 'function_handle');
end
end

function testUnknownSchemeErrors(testCase)
verifyError(testCase, @() get_scheme('nonsense'), 'burgers:unknownScheme');
end

function testUniformStateIsPreserved(testCase)
N = 32;
dx = 2*pi / N;
u = 3 * ones(N, 1);
names = get_scheme();
for k = 1:numel(names)
    nu = 0.01 * ~strcmp(names{k}, 'lax_friedrichs');   % LF is inviscid-only
    step = get_scheme(names{k});
    unew = step(u, 0.01, dx, nu);
    verifyEqual(testCase, unew, u, 'AbsTol', 1e-12);
end
end

function testConservativeStepsPreserveMass(testCase)
N = 64;
dx = 2*pi / N;
x = (0:N-1)' * dx;
u = 1 + 0.5*sin(x) + 0.2*cos(3*x);
names = get_scheme();
for k = 1:numel(names)
    if strcmp(names{k}, 'upwind_nc'), continue, end
    nu = 0.02 * ~strcmp(names{k}, 'lax_friedrichs');
    step = get_scheme(names{k});
    unew = step(u, 0.2*dx, dx, nu);
    verifyEqual(testCase, dx*sum(unew), dx*sum(u), 'AbsTol', 1e-10);
end
end

function testNonConservativeUpwindDriftsAcrossShock(testCase)
% Plateau 1 -> 2 -> 1: the conservative schemes keep the mass, upwind_nc loses it,
% and its shock lags the exact position 1.4.
f = ic_plateau(1, 2, 0.2, 0.8);
[x, uNc, infoNc] = solve_burgers(f, 'L', 2, 'N', 400, 'tFinal', 0.4, 'scheme', 'upwind_nc');
[~, uG, infoG] = solve_burgers(f, 'L', 2, 'N', 400, 'tFinal', 0.4, 'scheme', 'godunov');
verifyGreaterThan(testCase, abs(infoNc.massEnd - infoNc.mass0), 1e-3);
verifyLessThan(testCase, abs(infoG.massEnd - infoG.mass0), 1e-10);
verifyGreaterThan(testCase, abs(shock_position(x, uNc) - 1.4), 0.02);
verifyLessThan(testCase, abs(shock_position(x, uG) - 1.4), 0.01);
end

function testShockSpeedOfConservativeSchemes(testCase)
f = ic_plateau(1, 2, 0.2, 0.8);
for name = {'godunov', 'lax_friedrichs', 'lax_wendroff', 'maccormack'}
    [x, u] = solve_burgers(f, 'L', 2, 'N', 400, 'tFinal', 0.4, 'scheme', name{1});
    verifyLessThan(testCase, abs(shock_position(x, u) - 1.4), 0.015, name{1});
end
end

function testMonotoneSchemesRespectMaximumPrinciple(testCase)
f = @(x) 1.5 + 0.5*sin(x);
for name = {'godunov', 'lax_friedrichs'}
    [~, u] = solve_burgers(f, 'N', 100, 'tFinal', 1, 'scheme', name{1});
    verifyGreaterThanOrEqual(testCase, min(u), 1 - 1e-12, name{1});
    verifyLessThanOrEqual(testCase, max(u), 2 + 1e-12, name{1});
end
end

function testGodunovPreservesSymmetry(testCase)
% sin(x) is odd about x = pi, and so must the solution stay.
N = 100;
[~, u] = solve_burgers(@sin, 'N', N, 'tFinal', 2, 'scheme', 'godunov');
mirror = mod(N - (0:N-1), N) + 1;
verifyEqual(testCase, u + u(mirror), zeros(N, 1), 'AbsTol', 1e-10);
end

function testFtcsIsUnstableWithoutViscosity(testCase)
verifyError(testCase, ...
    @() solve_burgers(@sin, 'N', 200, 'nu', 0, 'tFinal', 10, 'scheme', 'ftcs'), ...
    'burgers:unstable');
end

function testExplicitSchemeBlowsUpAboveCflLimit(testCase)
verifyError(testCase, ...
    @() solve_burgers(@sin, 'N', 100, 'tFinal', 5, 'scheme', 'maccormack', 'cfl', 3), ...
    'burgers:unstable');
end

function testLaxFriedrichsRejectsViscosity(testCase)
verifyError(testCase, ...
    @() solve_burgers(@sin, 'nu', 0.01, 'scheme', 'lax_friedrichs'), ...
    'burgers:schemeInviscidOnly');
end

function testCrankNicolsonNewtonConvergesQuickly(testCase)
[~, ~, info] = solve_burgers(ic_cole_hopf(0.1), 'N', 200, 'nu', 0.1, 'tFinal', 0.5, ...
    'scheme', 'crank_nicolson');
verifyGreaterThan(testCase, info.newtonIters, 0);
verifyLessThan(testCase, info.newtonIters / info.steps, 6);
end

function pos = shock_position(x, u)
%SHOCK_POSITION  Midpoint of the cell pair with the steepest decrease.
[~, j] = min(diff(u));
pos = 0.5 * (x(j) + x(j+1));
end
