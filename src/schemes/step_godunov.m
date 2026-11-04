function unew = step_godunov(u, dt, dx, nu)
%STEP_GODUNOV  One step of Godunov's conservative first-order method.
%
%   unew_j = u_j - dt/dx*(F_{j+1/2} - F_{j-1/2}) + nu*dt/dx^2*(u_{j+1} - 2u_j + u_{j-1})
%
%   with the exact Riemann-problem flux of f(u) = u^2/2 between u_j and u_{j+1}.
%   Monotone (no new extrema for inviscid flow under dt*max|u|/dx <= 1),
%   conservative, and it picks the entropy-satisfying solution, so shocks have
%   the correct speed and rarefactions are not turned into shocks.

up = circshift(u, -1);
F = godunov_flux(u, up);            % F_{j+1/2}
unew = u - dt/dx * (F - circshift(F, 1)) + nu*dt/dx^2 * (up - 2*u + circshift(u, 1));
end

function F = godunov_flux(a, b)
%GODUNOV_FLUX  Exact Godunov flux for f(u) = u^2/2 with left state a, right state b.
%   a > b (shock):        max(f(a), f(b))
%   a <= b (rarefaction): min(f(a), f(b)), or 0 when the fan contains u = 0
fa = 0.5 * a.^2;
fb = 0.5 * b.^2;
F = max(fa, fb);
rare = a <= b;
F(rare) = min(fa(rare), fb(rare));
F(a < 0 & b > 0) = 0;
end
