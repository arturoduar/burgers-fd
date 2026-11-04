%EXAMPLE_PLATEAU_SHOCK  Rarefaction plus shock: why conservation form and entropy matter.
%   Plateau u = 2 on (0.2, 0.8) over background u = 1; exact solution known until t = 1.2.
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root, 'burgers_setup.m'));

L = 2; N = 400; T = 0.4;
ub = 1; uh = 2; xa = 0.2; xb = 0.8;
u0 = ic_plateau(ub, uh, xa, xb);
schemes = {'upwind_nc', 'godunov', 'lax_friedrichs', 'lax_wendroff', 'maccormack', 'crank_nicolson'};

fprintf('Exact shock position at t = %g: %.4f\n', T, xb + 0.5*(uh + ub)*T);
fprintf('%-16s %-10s %-10s %-11s %-8s %-8s\n', 'scheme', 'L1 error', 'shock x', 'mass drift', 'min u', 'max u');

fig = figure('Visible', 'off', 'Position', [100 100 1000 700]);
for s = 1:numel(schemes)
    [x, u, info] = solve_burgers(u0, 'L', L, 'N', N, 'tFinal', T, 'scheme', schemes{s});
    e = error_norms(u, burgers_exact_plateau(x, T, ub, uh, xa, xb), info.dx);
    [~, j] = min(diff(u));
    fprintf('%-16s %-10.3e %-10.4f %-11.2e %-8.3f %-8.3f\n', schemes{s}, e.L1, ...
        0.5*(x(j) + x(j+1)), info.massEnd - info.mass0, min(u), max(u));
    subplot(2, 3, s);
    xf = linspace(0, L, 2000)';
    plot(xf, burgers_exact_plateau(xf, T, ub, uh, xa, xb), 'k-'); hold on;
    plot(x, u, 'r.', 'MarkerSize', 6);
    title(strrep(schemes{s}, '_', '\_'));
    xlim([0 L]); ylim([0.7 2.4]); grid on;
end
file = fullfile(results_dir(), 'plateau_shock.png');
saveas(fig, file);
close(fig);
fprintf('Saved %s\n', file);
