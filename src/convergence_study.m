function result = convergence_study(problem, scheme, Ns)
%CONVERGENCE_STUDY  Measure a scheme's error and observed order of accuracy.
%
%   result = convergence_study(problem, scheme, Ns)
%
%   problem  'cole_hopf'  viscous Burgers, nu = 0.1, exact Cole-Hopf solution, t = 0.5
%            'smooth'     inviscid Burgers, u0 = sin(x), exact solution by
%                         characteristics at t = 0.5 (before the shock forms at t = 1)
%   scheme   any name from get_scheme (lax_friedrichs only on 'smooth')
%   Ns       vector of grid sizes, e.g. [50 100 200 400]
%
%   result has the columns N, dx, L1, L2, Linf and order, where order(k) is the
%   observed rate log(L1(k-1)/L1(k)) / log(N(k)/N(k-1)); order(1) is NaN.
%   Time steps scale with the grid automatically (see solve_burgers).

Ns = Ns(:);
n = numel(Ns);
L1 = zeros(n, 1);
L2 = zeros(n, 1);
Linf = zeros(n, 1);
T = 0.5;
nu = 0.1;

for k = 1:n
    N = Ns(k);
    switch problem
        case 'cole_hopf'
            [x, u] = solve_burgers(ic_cole_hopf(nu), 'N', N, 'nu', nu, ...
                'tFinal', T, 'scheme', scheme);
            uExact = burgers_cole_hopf(x, T, nu);
        case 'smooth'
            [x, u] = solve_burgers(@sin, 'N', N, 'nu', 0, 'tFinal', T, 'scheme', scheme);
            uExact = burgers_characteristics(x, T, @sin, 1);
        otherwise
            error('burgers:invalidInput', ...
                'Unknown problem "%s"; use ''cole_hopf'' or ''smooth''.', problem);
    end
    e = error_norms(u, uExact, 2*pi/N);
    L1(k) = e.L1;
    L2(k) = e.L2;
    Linf(k) = e.Linf;
end

order = [NaN; log(L1(1:end-1) ./ L1(2:end)) ./ log(Ns(2:end) ./ Ns(1:end-1))];
result = struct('problem', problem, 'scheme', scheme, 'N', Ns, 'dx', 2*pi ./ Ns, ...
    'L1', L1, 'L2', L2, 'Linf', Linf, 'order', order);
end
