function f = ic_cole_hopf(nu)
%IC_COLE_HOPF  Initial condition of the Cole-Hopf test problem, as a function handle.
%
%   The t = 0 state of burgers_cole_hopf(x, t, nu), a smooth periodic profile
%   on [0, 2*pi) that steepens into a viscous shock layer.

f = @(x) burgers_cole_hopf(x, 0, nu);
end
