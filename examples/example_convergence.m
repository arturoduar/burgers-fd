%EXAMPLE_CONVERGENCE  Observed order of accuracy of each scheme on two smooth problems.
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root, 'burgers_setup.m'));

Ns = [50 100 200 400 800];
cases = {'cole_hopf', {'godunov', 'lax_wendroff', 'maccormack', 'crank_nicolson'}; ...
         'smooth',    {'upwind_nc', 'godunov', 'lax_friedrichs', 'lax_wendroff', 'maccormack'}};

fig = figure('Visible', 'off', 'Position', [100 100 1000 450]);
for c = 1:size(cases, 1)
    problem = cases{c, 1};
    names = cases{c, 2};
    fprintf('\nProblem "%s": L1 error and observed order\n', problem);
    subplot(1, 2, c);
    for s = 1:numel(names)
        r = convergence_study(problem, names{s}, Ns);
        fprintf('%-16s', names{s});
        fprintf('%9.2e(%4.2f) ', [r.L1, r.order]');
        fprintf('\n');
        loglog(r.dx, r.L1, 'o-', 'DisplayName', strrep(names{s}, '_', '\_'));
        hold on;
    end
    dx = 2*pi ./ Ns;
    loglog(dx, 0.5*dx, 'k--', 'DisplayName', 'O(\Deltax)');
    loglog(dx, 0.5*dx.^2, 'k:', 'DisplayName', 'O(\Deltax^2)');
    xlabel('\Deltax'); ylabel('L1 error');
    title(strrep(problem, '_', '\_'));
    legend('Location', 'southeast');
    grid on;
end
file = fullfile(results_dir(), 'convergence.png');
saveas(fig, file);
close(fig);
fprintf('\nSaved %s\n', file);
