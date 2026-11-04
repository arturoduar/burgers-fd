function unew = step_ftcs(u, dt, dx, nu)
%STEP_FTCS  One step of forward-time, central-space differencing.
%
%   unew_j = u_j - dt/(2*dx)*(f_{j+1} - f_{j-1}) + nu*dt/dx^2*(u_{j+1} - 2u_j + u_{j-1})
%
%   with f = u^2/2.  Only usable for viscous flow: with nu = 0 the scheme is
%   unconditionally unstable (solve_burgers then raises burgers:unstable), and
%   with nu > 0 it needs dt <= 2*nu/max|u|^2 as well as the diffusive limit.

f = 0.5 * u.^2;
up = circshift(u, -1);
um = circshift(u, 1);
unew = u - dt/(2*dx) * (circshift(f, -1) - circshift(f, 1)) + nu*dt/dx^2 * (up - 2*u + um);
end
