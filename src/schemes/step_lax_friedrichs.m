function unew = step_lax_friedrichs(u, dt, dx, nu)
%STEP_LAX_FRIEDRICHS  One step of the Lax-Friedrichs scheme (inviscid only).
%
%   unew_j = (u_{j+1} + u_{j-1})/2 - dt/(2*dx)*(f_{j+1} - f_{j-1}),   f = u^2/2
%
%   Conservative and monotone but very diffusive.  Its numerical viscosity
%   cannot be combined with an explicit physical viscosity (the Nyquist mode
%   then has amplification -1 - 4*nu*dt/dx^2), so nu must be 0.

if nu > 0
    error('burgers:schemeInviscidOnly', ...
        'lax_friedrichs supports nu = 0 only (it is unstable with an explicit viscous term).');
end
up = circshift(u, -1);
um = circshift(u, 1);
unew = 0.5 * (up + um) - dt/(4*dx) * (up.^2 - um.^2);
end
