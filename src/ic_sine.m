function f = ic_sine(L, amplitude)
%IC_SINE  Initial condition u0(x) = amplitude*sin(2*pi*x/L), as a function handle.
%
%   With L = 2*pi and amplitude = 1 this is the classic sin(x) problem whose
%   characteristics cross, forming a shock at t = 1.

if nargin < 1, L = 2*pi; end
if nargin < 2, amplitude = 1; end
f = @(x) amplitude * sin(2*pi*x / L);
end
