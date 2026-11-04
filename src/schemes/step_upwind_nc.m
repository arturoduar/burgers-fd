function unew = step_upwind_nc(u, dt, dx, nu)
%STEP_UPWIND_NC  One step of first-order upwinding of the NON-conservative form.
%
%   Discretizes u_t + u*u_x = nu*u_xx with u*u_x written as u times a one-sided
%   difference taken on the side the wave comes from:
%       u_x ~ (u_j - u_{j-1})/dx  if u_j >= 0,   (u_{j+1} - u_j)/dx  if u_j < 0
%   Viscosity uses the central second difference.
%
%   This scheme is deliberately NOT conservative: sum(u) drifts, and across a
%   shock it converges to a discontinuity that moves at the wrong speed.  It is
%   included to show why conservation form matters (see example_plateau_shock).

up = circshift(u, -1);
um = circshift(u, 1);
ux = u - um;                 % backward difference (numerator)
neg = u < 0;
ux(neg) = up(neg) - u(neg);  % forward difference where the flow goes left
unew = u - dt/dx * u .* ux + nu*dt/dx^2 * (up - 2*u + um);
end
