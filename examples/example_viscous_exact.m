%EXAMPLE_VISCOUS_EXACT  Compare three schemes with the exact Cole-Hopf solution (nu = 0.1).
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root, 'burgers_setup.m'));

nu = 0.1;
N = 200;
times = [0.25 0.5 1];
schemes = {'godunov', 'maccormack', 'crank_nicolson'};
u0 = ic_cole_hopf(nu);

fprintf('Viscous Burgers, nu = %g, N = %d: L1 error against the exact solution\n', nu, N);
fprintf('%-16s', 'scheme');
fprintf('t = %-8g', times);
fprintf('\n');

fig = figure('Visible', 'off', 'Position', [100 100 900 450]);
hold on;
xf = linspace(0, 2*pi, 1000)';
for k = 1:numel(times)
    h = plot(xf, burgers_cole_hopf(xf, times(k), nu), 'k-', 'LineWidth', 1);
    if k == 1
        set(h, 'DisplayName', 'exact');
    else
        set(h, 'HandleVisibility', 'off');
    end
end
colors = lines(numel(schemes));
for s = 1:numel(schemes)
    [x, ~, info] = solve_burgers(u0, 'N', N, 'nu', nu, 'tFinal', times(end), ...
        'tSave', times(1:end-1), 'scheme', schemes{s});
    fprintf('%-16s', schemes{s});
    for k = 1:numel(times)
        e = error_norms(info.U(:, k), burgers_cole_hopf(x, times(k), nu), info.dx);
        fprintf('%-12.3e', e.L1);
    end
    fprintf('\n');
    plot(x, info.U(:, end), '.', 'Color', colors(s, :), 'MarkerSize', 10, ...
        'DisplayName', strrep(schemes{s}, '_', '\_'));
end
xlabel('x'); ylabel('u');
title(sprintf('Viscous Burgers, \\nu = %g: solutions at t = %s', nu, mat2str(times)));
legend('Location', 'best');
grid on;
file = fullfile(results_dir(), 'viscous_exact.png');
saveas(fig, file);
close(fig);
fprintf('Saved %s\n', file);
