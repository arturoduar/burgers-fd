function unew = step_lax_wendroff(u, dt, dx, nu)
%STEP_LAX_WENDROFF  One step of the two-step Richtmyer form of Lax-Wendroff.
%
%   Half step (cell interfaces):
%       u*_{j+1/2} = (u_j + u_{j+1})/2 - dt/(2*dx)*(f_{j+1} - f_j)
%   Full step (conservative):
%       unew_j = u_j - dt/dx*(f(u*_{j+1/2}) - f(u*_{j-1/2})) + nu*dt/dx^2*(u_{j+1} - 2u_j + u_{j-1})
%
%   Second-order accurate for smooth flow, but oscillates near shocks.

f = 0.5 * u.^2;
up = circshift(u, -1);                        % u_{j+1}
um = circshift(u, 1);                         % u_{j-1}
uh = 0.5 * (u + up) - dt/(2*dx) * (circshift(f, -1) - f);   % u*_{j+1/2}
Fh = 0.5 * uh.^2;
unew = u - dt/dx * (Fh - circshift(Fh, 1)) + nu*dt/dx^2 * (up - 2*u + um);
end
