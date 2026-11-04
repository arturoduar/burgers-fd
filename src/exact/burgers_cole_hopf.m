function u = burgers_cole_hopf(x, t, nu)
%BURGERS_COLE_HOPF  Exact periodic solution of the viscous Burgers equation.
%
%   u = burgers_cole_hopf(x, t, nu) evaluates, on the period 2*pi, the solution of
%       u_t + u*u_x = nu*u_xx
%   obtained from the Cole-Hopf transform u = -2*nu*phi_x/phi with the periodized
%   heat-kernel solution phi.  With z_k = x - 4t - 2*pi*k,
%       u(x,t) = 4 + sum_k w_k*z_k/(t+1),    w_k ~ exp(-z_k^2/(4*nu*(t+1))),
%   the weights normalized to sum to one (computed with the maximum subtracted so
%   the exponentials cannot underflow).  The "4" is a Galilean shift; the mean of
%   u over a period is exactly 4 for all t.  Image sums k = -4..4 are far more
%   than enough for nu >= 0.001.
%
%   x may have any shape and any value (the solution is 2*pi-periodic); t >= 0
%   and nu > 0 are scalars.  At t = 0 this gives the initial condition used by
%   ic_cole_hopf.

if ~(isscalar(nu) && nu > 0)
    error('burgers:invalidInput', 'burgers_cole_hopf needs a scalar nu > 0.');
end
sz = size(x);
xr = mod(x(:) - 4*t, 2*pi);          % reduce to one period
k = -4:4;
z = xr - 2*pi*k;                     % numel(x)-by-9 (implicit expansion)
e = -z.^2 / (4*nu*(t + 1));
e = e - max(e, [], 2);
w = exp(e);
u = 4 + (sum(w .* z, 2) ./ sum(w, 2)) / (t + 1);
u = reshape(u, sz);
end
