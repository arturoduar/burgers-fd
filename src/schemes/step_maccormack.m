function unew = step_maccormack(u, dt, dx, nu)
%STEP_MACCORMACK  One step of the MacCormack predictor-corrector scheme.
%
%   Predictor (forward difference):
%       u*_j = u_j - dt/dx*(f_{j+1} - f_j) + nu*dt/dx^2*lap(u)_j
%   Corrector (backward difference on the predicted values):
%       unew_j = (u_j + u*_j - dt/dx*(f*_j - f*_{j-1}) + nu*dt/dx^2*lap(u*)_j)/2
%
%   with f = u^2/2 and lap(v)_j = v_{j+1} - 2v_j + v_{j-1}.  Second-order
%   accurate and conservative; like Lax-Wendroff it oscillates at shocks.

f = 0.5 * u.^2;
ustar = u - dt/dx * (circshift(f, -1) - f) ...
    + nu*dt/dx^2 * (circshift(u, -1) - 2*u + circshift(u, 1));
fs = 0.5 * ustar.^2;
unew = 0.5 * (u + ustar - dt/dx * (fs - circshift(fs, 1)) ...
    + nu*dt/dx^2 * (circshift(ustar, -1) - 2*ustar + circshift(ustar, 1)));
end
