function tests = test_convergence
%TEST_CONVERGENCE  Observed orders of accuracy against exact solutions.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
oldPath = path;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
testCase.addTeardown(@path, oldPath);
end

function testViscousFirstOrderScheme(testCase)
r = convergence_study('cole_hopf', 'godunov', [100 200 400]);
verifyGreaterThan(testCase, r.order(end), 0.8);
end

function testViscousSecondOrderSchemes(testCase)
for name = {'lax_wendroff', 'maccormack', 'crank_nicolson'}
    r = convergence_study('cole_hopf', name{1}, [100 200 400]);
    verifyGreaterThan(testCase, r.order(end), 1.8, name{1});
end
end

function testSmoothInviscidFirstOrderSchemes(testCase)
for name = {'upwind_nc', 'godunov', 'lax_friedrichs'}
    r = convergence_study('smooth', name{1}, [100 200 400]);
    verifyGreaterThan(testCase, r.order(end), 0.85, name{1});
    verifyLessThan(testCase, r.order(end), 1.3, name{1});
end
end

function testSmoothInviscidSecondOrderSchemes(testCase)
for name = {'lax_wendroff', 'maccormack'}
    r = convergence_study('smooth', name{1}, [100 200 400]);
    verifyGreaterThan(testCase, r.order(end), 1.8, name{1});
end
end

function testErrorsDecreaseWithResolution(testCase)
r = convergence_study('smooth', 'maccormack', [50 100 200]);
verifyTrue(testCase, all(diff(r.L1) < 0));
verifyEqual(testCase, r.N, [50; 100; 200]);
end

function testUnknownProblemErrors(testCase)
verifyError(testCase, @() convergence_study('nope', 'godunov', 50), 'burgers:invalidInput');
end
