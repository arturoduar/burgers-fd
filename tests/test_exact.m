function tests = test_exact
%TEST_EXACT  Tests for the exact solutions used as references.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
oldPath = path;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
testCase.addTeardown(@path, oldPath);
end

function testColeHopfIsPeriodic(testCase)
x = linspace(0, 2*pi, 50)';
verifyEqual(testCase, burgers_cole_hopf(x + 2*pi, 0.3, 0.1), ...
    burgers_cole_hopf(x, 0.3, 0.1), 'AbsTol', 1e-10);
end

function testColeHopfMeanIsFour(testCase)
% The Galilean shift makes the period mean exactly 4 at every time.
x = (0:1999)' * 2*pi / 2000;
verifyEqual(testCase, mean(burgers_cole_hopf(x, 0.5, 0.1)), 4, 'AbsTol', 1e-9);
end

function testColeHopfSatisfiesPde(testCase)
% Residual of u_t + u*u_x - nu*u_xx with fourth-order central differences.
nu = 0.2;
t = 0.5;
h = 1e-3;
x = linspace(0, 2*pi, 31)';
f = @(xx, tt) burgers_cole_hopf(xx, tt, nu);
d1 = @(fm2, fm1, fp1, fp2) (-fp2 + 8*fp1 - 8*fm1 + fm2) / (12*h);
ut = d1(f(x, t-2*h), f(x, t-h), f(x, t+h), f(x, t+2*h));
ux = d1(f(x-2*h, t), f(x-h, t), f(x+h, t), f(x+2*h, t));
uxx = (-f(x+2*h, t) + 16*f(x+h, t) - 30*f(x, t) + 16*f(x-h, t) - f(x-2*h, t)) / (12*h^2);
residual = ut + f(x, t) .* ux - nu * uxx;
verifyLessThan(testCase, max(abs(residual)), 1e-4);
end

function testColeHopfRejectsZeroViscosity(testCase)
verifyError(testCase, @() burgers_cole_hopf(0, 0.1, 0), 'burgers:invalidInput');
end

function testCharacteristicsSolveImplicitEquation(testCase)
x = linspace(0, 2*pi, 101)';
t = 0.5;
u = burgers_characteristics(x, t, @sin, 1);
verifyEqual(testCase, u, sin(x - t*u), 'AbsTol', 1e-12);
end

function testCharacteristicsAtTimeZero(testCase)
x = linspace(0, 2*pi, 51)';
verifyEqual(testCase, burgers_characteristics(x, 0, @sin, 1), sin(x), 'AbsTol', 1e-14);
end

function testBreakingTimeOfSine(testCase)
verifyEqual(testCase, burgers_breaking_time(@sin, 2*pi), 1, 'RelTol', 1e-6);
end

function testBreakingTimeScalesWithAmplitude(testCase)
verifyEqual(testCase, burgers_breaking_time(ic_sine(2*pi, 2), 2*pi), 0.5, 'RelTol', 1e-6);
end

function testBreakingTimeInfiniteWithoutCompression(testCase)
verifyEqual(testCase, burgers_breaking_time(@(x) 2 + 0*x, 2*pi), Inf);
end

function testPlateauSolutionValues(testCase)
% ub = 1, uh = 2 on (0.2, 0.8); at t = 0.4 the fan spans [0.6, 1.0] and the shock sits at 1.4.
x = [0.1; 0.8; 1.2; 1.5];
u = burgers_exact_plateau(x, 0.4, 1, 2, 0.2, 0.8);
verifyEqual(testCase, u, [1; 1.5; 2; 1], 'AbsTol', 1e-14);
end

function testPlateauAtTimeZeroIsInitialCondition(testCase)
x = linspace(0, 2, 41)';
u0 = ic_plateau(1, 2, 0.2, 0.8);
verifyEqual(testCase, burgers_exact_plateau(x, 0, 1, 2, 0.2, 0.8), u0(x));
end

function testPlateauConservesMass(testCase)
x = linspace(0, 2, 200001)';
u = burgers_exact_plateau(x, 0.4, 1, 2, 0.2, 0.8);
mass = trapz(x, u);
verifyEqual(testCase, mass, 1*2 + (2-1)*(0.8-0.2), 'AbsTol', 1e-4);
end

function testPlateauInvalidAfterFanCatchesShock(testCase)
% t* = 2*(0.8-0.2)/(2-1) = 1.2
verifyError(testCase, @() burgers_exact_plateau(1, 1.3, 1, 2, 0.2, 0.8), ...
    'burgers:exactPlateauInvalid');
end
