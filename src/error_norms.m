function e = error_norms(u, uExact, dx)
%ERROR_NORMS  Grid-weighted L1, L2 and maximum norms of the error u - uExact.
%
%   e = error_norms(u, uExact, dx) returns a struct with fields
%       L1   = dx * sum(|u - uExact|)
%       L2   = sqrt(dx * sum((u - uExact).^2))
%       Linf = max(|u - uExact|)

d = u(:) - uExact(:);
e.L1 = dx * sum(abs(d));
e.L2 = sqrt(dx * sum(d.^2));
e.Linf = max(abs(d));
end
