function [x, u, info] = solve_burgers(u0fun, varargin)
%SOLVE_BURGERS  Finite-difference solution of the 1-D Burgers equation on a periodic domain.
%
%   u_t + u*u_x = nu*u_xx   on  0 <= x < L,  periodic.
%
%   [x, u, info] = solve_burgers(u0fun, Name, Value, ...)
%
%   u0fun is a function handle returning the initial condition for a column
%   vector of grid points, e.g. @(x) sin(x).  The grid is x_j = (j-1)*L/N for
%   j = 1..N (the right end point is the periodic image of the left one).
%
%   Name-value options
%     'L'       domain length                              (default 2*pi)
%     'N'       number of grid points, N >= 4              (default 200)
%     'nu'      viscosity, nu >= 0                         (default 0)
%     'tFinal'  final time                                 (default 1)
%     'scheme'  see get_scheme                             (default 'maccormack')
%     'cfl'     safety factor for the time step            (default 0.5)
%     'tSave'   extra times at which to record the solution (default [])
%
%   The time step is chosen automatically from the largest current |u| and nu
%   (see choose_dt below), and shortened so every requested time is hit exactly.
%
%   Outputs
%     x     grid (N-by-1)
%     u     solution at tFinal (N-by-1)
%     info  struct with fields scheme, nu, N, dx, steps, dtMin, dtMax,
%           t (recorded times, column), U (N-by-numel(t) solution snapshots,
%           the last column equals u), mass0, massEnd (integral of u, which
%           conservative schemes keep constant) and newtonIters (total Newton
%           iterations, crank_nicolson only).
%
%   Errors
%     burgers:unknownScheme, burgers:schemeInviscidOnly, burgers:invalidInput,
%     burgers:invalidInitialCondition, burgers:unstable (the solution blew up).
%
%   Example
%     [x, u] = solve_burgers(@sin, 'N', 200, 'tFinal', 2, 'scheme', 'godunov');

p = inputParser;
p.FunctionName = 'solve_burgers';
addParameter(p, 'L', 2*pi, @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'N', 200, @(v) isnumeric(v) && isscalar(v) && v == round(v) && v >= 4);
addParameter(p, 'nu', 0, @(v) isnumeric(v) && isscalar(v) && v >= 0);
addParameter(p, 'tFinal', 1, @(v) isnumeric(v) && isscalar(v) && v >= 0);
addParameter(p, 'scheme', 'maccormack', @(v) ischar(v) || isstring(v));
addParameter(p, 'cfl', 0.5, @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'tSave', [], @(v) isnumeric(v));
parse(p, varargin{:});
opt = p.Results;

L = opt.L;
N = opt.N;
nu = opt.nu;
cfl = opt.cfl;
scheme = char(opt.scheme);
step = get_scheme(scheme);

if strcmp(scheme, 'lax_friedrichs') && nu > 0
    error('burgers:schemeInviscidOnly', ...
        ['lax_friedrichs cannot be combined with an explicit viscous term: its amplification ' ...
         'factor at the Nyquist wavenumber is -1-4*nu*dt/dx^2, so it is unstable for any nu > 0. ' ...
         'Use nu = 0, or choose godunov, lax_wendroff, maccormack or crank_nicolson.']);
end
if any(opt.tSave(:) < 0 | opt.tSave(:) > opt.tFinal)
    error('burgers:invalidInput', 'Every time in tSave must lie between 0 and tFinal.');
end

dx = L / N;
x = (0:N-1)' * dx;
u = u0fun(x);
u = u(:);
if numel(u) ~= N || any(~isfinite(u))
    error('burgers:invalidInitialCondition', ...
        'u0fun must return N finite values for the N-by-1 grid it is given.');
end

targets = unique([opt.tSave(:); opt.tFinal]);   % sorted, increasing
blowup = 1e3 * max(1, max(abs(u)));
mass0 = dx * sum(u);
U = zeros(N, numel(targets));
t = 0;
steps = 0;
dtMin = Inf;
dtMax = 0;
newtonIters = 0;

for k = 1:numel(targets)
    target = targets(k);
    while t < target - 1e-12 * max(1, target)
        dt = min(choose_dt(scheme, u, dx, nu, cfl), target - t);
        if strcmp(scheme, 'crank_nicolson')
            [u, its] = step(u, dt, dx, nu);
            newtonIters = newtonIters + its;
        else
            u = step(u, dt, dx, nu);
        end
        t = t + dt;
        steps = steps + 1;
        dtMin = min(dtMin, dt);
        dtMax = max(dtMax, dt);
        if ~all(isfinite(u)) || max(abs(u)) > blowup
            error('burgers:unstable', ...
                ['Scheme "%s" became unstable at t = %.4g (step %d). Try a smaller cfl, add ' ...
                 'viscosity, or choose a different scheme.'], scheme, t, steps);
        end
    end
    t = target;
    U(:, k) = u;
end

if steps == 0
    dtMin = 0;
end
info = struct('scheme', scheme, 'nu', nu, 'N', N, 'dx', dx, 'steps', steps, ...
    'dtMin', dtMin, 'dtMax', dtMax, 't', targets, 'U', U, ...
    'mass0', mass0, 'massEnd', dx * sum(u), 'newtonIters', newtonIters);
end

function dt = choose_dt(scheme, u, dx, nu, cfl)
%CHOOSE_DT  Largest stable-looking time step for the current state.
%   Explicit schemes combine the advective and diffusive limits,
%       dt = cfl / (max|u|/dx + 2*nu/dx^2),
%   so that cfl = 1 sits on the stability boundary of first-order upwinding.
%   FTCS also needs dt <= 2*nu/max|u|^2.  Crank-Nicolson is implicit, so only
%   accuracy limits its step: dt = cfl*dx/max|u|.
umax = max(abs(u));
if strcmp(scheme, 'crank_nicolson')
    dt = cfl * dx / max(umax, eps);
    return
end
dt = cfl / (umax/dx + 2*nu/dx^2);
if strcmp(scheme, 'ftcs') && nu > 0
    dt = min(dt, 0.9 * 2 * nu / max(umax^2, eps));
end
end
