function tb = burgers_breaking_time(u0fun, L)
%BURGERS_BREAKING_TIME  Time at which the inviscid solution first forms a shock.
%
%   tb = burgers_breaking_time(u0fun, L) returns -1/min(u0'), estimated with
%   central differences on a fine grid over [0, L], or Inf if u0 is
%   nondecreasing (no shock ever forms).  For u0 = sin(x) the result is 1.

xs = linspace(0, L, 20001)';
h = 1e-6 * L;
d = (u0fun(xs + h) - u0fun(xs - h)) / (2*h);
dmin = min(d);
if dmin >= 0
    tb = Inf;
else
    tb = -1 / dmin;
end
end
