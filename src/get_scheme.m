function [step, names] = get_scheme(name)
%GET_SCHEME  Look up a finite-difference time-stepping function by name.
%
%   step = get_scheme(name) returns a handle to step_<name>, a function
%       unew = step(u, dt, dx, nu)
%   advancing the periodic solution u (column vector) by one step dt.
%
%   names = the list of available schemes (also returned as the second output;
%   call with no input to get just the list).
%
%   Scheme            Type                 Order (space/time)   Use it for
%   upwind_nc         explicit, NOT cons.  1 / 1                demonstrating why conservation matters
%   godunov           explicit, cons.      1 / 1                shocks (monotone, entropy-satisfying)
%   lax_friedrichs    explicit, cons.      1 / 1                shocks, inviscid only (very diffusive)
%   lax_wendroff      explicit, cons.      2 / 2 (inviscid)     smooth flow; oscillates at shocks
%   maccormack        explicit, cons.      2 / 2                smooth or viscous flow; oscillates at shocks
%   ftcs              explicit, cons.      2 / 1                viscous flow only (unstable when nu = 0)
%   crank_nicolson    implicit, cons.      2 / 2                smooth or viscous flow, larger time steps

names = {'upwind_nc', 'godunov', 'lax_friedrichs', 'lax_wendroff', ...
         'maccormack', 'ftcs', 'crank_nicolson'};
if nargin == 0
    step = names;
    return
end
name = char(name);
if ~any(strcmp(name, names))
    error('burgers:unknownScheme', 'Unknown scheme "%s". Available schemes: %s.', ...
        name, strjoin(names, ', '));
end
step = str2func(['step_' name]);
end
