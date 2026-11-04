function f = ic_plateau(ub, uh, xa, xb)
%IC_PLATEAU  A raised plateau on a constant background, as a function handle.
%
%   u0(x) = uh for xa < x < xb, and ub elsewhere.
%
%   With ub < uh the left edge becomes a rarefaction fan and the right edge a
%   shock; burgers_exact_plateau gives the exact solution until they meet.
%   ub = 0 puts the fan's tail at the sonic point u = 0, which is a hard case
%   for schemes that are not entropy-satisfying.

f = @(x) ub + (uh - ub) * (x > xa & x < xb);
end
