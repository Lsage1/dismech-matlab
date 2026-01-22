clc
clear all
close all
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesColor','w');

%% Experimental Data
weights = [290, 436, 186.5, 480, 520, 340, 490, 550, 610, 571, 610];
displacements = [.175, .275, .125, .375, .475, .225, .425, .525, .675, .625, .725];
force_N = weights * 9.81 / 1000;
displacement_mm = displacements * 25.4;

% Fit experimental data
p = polyfit(displacement_mm, force_N, 2);
x_fit = linspace(min(displacement_mm), max(displacement_mm), 100);
fit_line = polyval(p, x_fit);

%% Load Simulation Results
load('force_displacement_results.mat', 'results');
sim_force = results.force;
sim_displacement_mm = -1 * results.displacement * 1000; % Convert to mm and make positive

%% Overlay Plot
scatter(displacement_mm, force_N, 'MarkerEdgeColor', [0.7 0.7 0.7], 'MarkerFaceColor', [0.7 0.7 0.7]);
hold on;
plot(x_fit, fit_line, 'Color', [0.5 0 0.5], 'LineWidth', 3);
plot(sim_displacement_mm, abs(sim_force), 'Color', [0 0.5 0.8], 'LineWidth', 3);
scatter(sim_displacement_mm, abs(sim_force), 'MarkerEdgeColor', [0 0.3 0.6], 'MarkerFaceColor', [0.3 0.6 0.9]);

xlabel('Displacement (mm)', 'FontSize', 12);
ylabel('Force (N)', 'FontSize', 12);
set(gca, 'Color', [0.95 0.95 0.95], 'XColor', 'k', 'YColor', 'k', 'TickDir', 'out', 'TickLength', [0 0], 'FontSize', 11);
grid on;
set(gca, 'GridColor', 'w', 'GridAlpha', 1);
legend('Experimental Data', 'Experimental Fit', 'Simulation Curve', 'Simulation Points', 'Location', 'northwest');
hold off;