%EXAMPLE_SHOCK_FORMATION  Inviscid sin(x): smooth steepening, shock at t = 1, and what follows.
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root, 'burgers_setup.m'));

N = 200;
times = [0.5 1 2];
schemes = {'godunov', 'maccormack'};
fprintf('Breaking time of sin(x): %.4f\n', burgers_breaking_time(@sin, 2*pi));

fig = figure('Visible', 'off', 'Position', [100 100 1000 450]);
for s = 1:numel(schemes)
    [x, ~, info] = solve_burgers(@sin, 'N', N, 'tFinal', times(end), ...
        'tSave', times(1:end-1), 'scheme', schemes{s});
    subplot(1, 2, s);
    hold on;
    plot(x, sin(x), 'k:', 'DisplayName', 't = 0');
    for k = 1:numel(times)
        plot(x, info.U(:, k), 'LineWidth', 1.2, 'DisplayName', sprintf('t = %g', times(k)));
    end
    xlabel('x'); ylabel('u');
    title(strrep(schemes{s}, '_', '\_'));
    legend('Location', 'best');
    grid on;
    fprintf('%-12s at t = %g: min u = %+.4f, max u = %+.4f, mass drift = %.1e\n', ...
        schemes{s}, times(end), min(info.U(:, end)), max(info.U(:, end)), ...
        info.massEnd - info.mass0);
end
file = fullfile(results_dir(), 'shock_formation.png');
saveas(fig, file);
close(fig);
fprintf('Saved %s\n', file);
