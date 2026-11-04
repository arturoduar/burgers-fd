function u = burgers_exact_plateau(x, t, ub, uh, xa, xb)
%BURGERS_EXACT_PLATEAU  Exact inviscid solution for a raised plateau on a background.
%
%   u = burgers_exact_plateau(x, t, ub, uh, xa, xb)
%
%   Initial data: u0 = uh on (xa, xb) and ub elsewhere, with ub < uh.  Then
%     * the left edge (ub | uh) is a rarefaction fan,
%           u = ub               for x <= xa + ub*t
%           u = (x - xa)/t       for xa + ub*t < x <= xa + uh*t
%     * the right edge (uh | ub) is a shock moving at the Rankine-Hugoniot
%       speed s = (uh + ub)/2, from xb to xb + s*t,
%     * between them u = uh.
%   The fan head xa + uh*t reaches the shock at t* = 2*(xb - xa)/(uh - ub); after
%   that the plateau is eaten away and this formula no longer applies
%   (burgers:exactPlateauInvalid).  Choose the domain wide enough that no wave
%   reaches the periodic boundary before t.

if ~(ub < uh)
    error('burgers:invalidInput', 'burgers_exact_plateau needs ub < uh.');
end
tStar = 2 * (xb - xa) / (uh - ub);
if t >= tStar
    error('burgers:exactPlateauInvalid', ...
        'The fan has caught the shock at t* = %g; the exact formula holds only for t < t*.', tStar);
end
u = ub + zeros(size(x));
if t == 0
    u(x > xa & x < xb) = uh;
    return
end
s = 0.5 * (uh + ub);
fan = x > xa + ub*t & x <= xa + uh*t;
u(fan) = (x(fan) - xa) / t;
plateau = x > xa + uh*t & x < xb + s*t;
u(plateau) = uh;
end
