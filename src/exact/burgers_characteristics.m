function u = burgers_characteristics(x, t, u0fun, umax)
%BURGERS_CHARACTERISTICS  Exact inviscid Burgers solution by the method of characteristics.
%
%   u = burgers_characteristics(x, t, u0fun, umax)
%
%   Before shocks form, u is constant along the lines x = xi + t*u0(xi), so
%   u(x,t) = u0(xi) where xi solves xi + t*u0(xi) = x.  The left side increases
%   with xi when t < breaking time (see burgers_breaking_time), so xi is found
%   by bisection on [x - t*umax, x + t*umax], vectorized over x.
%
%   u0fun  vectorized, periodic initial condition handle
%   umax   any bound with umax >= max|u0|
%
%   Valid only for 0 <= t <= breaking time; no check is made.

if t < 0
    error('burgers:invalidInput', 't must be nonnegative.');
end
sz = size(x);
xv = x(:);
lo = xv - t*umax;
hi = xv + t*umax;
for iter = 1:80
    mid = 0.5 * (lo + hi);
    g = mid + t * u0fun(mid) - xv;
    pos = g > 0;
    hi(pos) = mid(pos);
    lo(~pos) = mid(~pos);
end
u = reshape(u0fun(0.5 * (lo + hi)), sz);
end
