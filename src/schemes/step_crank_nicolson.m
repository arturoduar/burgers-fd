function [unew, iters] = step_crank_nicolson(u, dt, dx, nu)
%STEP_CRANK_NICOLSON  One implicit Crank-Nicolson step, solved with Newton's method.
%
%   Solves for w = u^{n+1}:
%       w - u + dt/2*D(w.^2/2) + dt/2*D(u.^2/2) - nu*dt/2*(L*w + L*u) = 0
%   where D is the periodic central first-derivative matrix and L the periodic
%   central second-derivative matrix.  Newton's method uses the sparse Jacobian
%       J = I + dt/2*D*diag(w) - nu*dt/2*L
%   and starts from w = u.  Second-order in time and space, conservative, and
%   not limited by a CFL condition (solve_burgers uses dt = cfl*dx/max|u|
%   purely for accuracy).  Central differencing makes it oscillatory for
%   under-resolved shocks.
%
%   [unew, iters] = step_crank_nicolson(...) also returns the Newton iteration count.
%   Raises burgers:newtonFailed if Newton does not converge in 25 iterations.

N = numel(u);
idx = (1:N)';
ip = [2:N, 1]';
im = [N, 1:N-1]';
D = sparse(idx, ip, 1/(2*dx), N, N) + sparse(idx, im, -1/(2*dx), N, N);
Lap = (sparse(idx, ip, 1, N, N) + sparse(idx, im, 1, N, N) - 2*speye(N)) / dx^2;
Id = speye(N);

fixed = dt/2 * (D * (0.5 * u.^2)) - nu*dt/2 * (Lap * u) - u;   % terms that do not involve w
w = u;
maxIter = 25;
tol = 1e-12;
for iters = 1:maxIter
    R = w + dt/2 * (D * (0.5 * w.^2)) - nu*dt/2 * (Lap * w) + fixed;
    J = Id + dt/2 * D * spdiags(w, 0, N, N) - nu*dt/2 * Lap;
    dw = J \ R;
    w = w - dw;
    if max(abs(dw)) < tol * max(1, max(abs(w)))
        unew = w;
        return
    end
end
error('burgers:newtonFailed', ...
    'Newton iteration did not converge in %d iterations (dt = %g).', maxIter, dt);
end
