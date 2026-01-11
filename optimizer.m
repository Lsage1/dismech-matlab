clear all
close all
clc

addpath experiments/hapticDevice
addpath springs/
addpath util_functions/
addpath contact_functions/
addpath rod_dynamics/
addpath shell_dynamics/
addpath external_forces/
addpath adaptive_stepping/
addpath logging/
addpath(genpath('experiments')); 

%% First, run the robot description to set up all variables
robotDescriptionHapticDevice;

%% Define simulation parameters
edges_to_actuate = [1, 2, 3];  % Which edges to actuate
actuation_factors = [1, 0.5, 1];  % Multiply initial rest length by these factors
top_verts_ind = [5, 35, 65];  % Nodes to track (the top vertices)

%% Run single simulation
%[tilt_x, tilt_y] = simulate_geometry(edges_to_actuate, actuation_factors, top_verts_ind, ...
%    material, geom, env, sim_params, nodes, edges, face_nodes, ...
%    fixed_node_indices, fixed_edge_indices);

%% Display results
%fprintf('\n=== Simulation Results ===\n');
%fprintf('Tilt about X-axis (pitch): %.2f degrees\n', tilt_x);
%fprintf('Tilt about Y-axis (roll):  %.2f degrees\n', tilt_y);
%fprintf('Total tilt magnitude:      %.2f degrees\n', sqrt(tilt_x^2 + tilt_y^2));


%% Example: Test multiple actuation scenarios
fprintf('\n=== Testing Multiple Actuation Scenarios ===\n');

test_scenarios = [
    0.5867, 0.8330, 0.5867;   % Scenario 1: single tendon shortened
    0.3333, 0.9003, 0.9003;   % Scenario 2: double tendon shortened
    0.6946, 0.6946, 0.6946;   % Scenario 3: equal reduction for 3 tendons
];

results_x = zeros(size(test_scenarios, 1), 1);
results_y = zeros(size(test_scenarios, 1), 1);
results_mag = zeros(size(test_scenarios, 1), 1);

for i = 1:size(test_scenarios, 1)
    fprintf('Running scenario %d: [%.2f, %.2f, %.2f]... ', i, test_scenarios(i,:));
    [tx, ty] = simulate_geometry(edges_to_actuate, test_scenarios(i,:), top_verts_ind, ...
        material, geom, env, sim_params, nodes, edges, face_nodes, ...
        fixed_node_indices, fixed_edge_indices);
    results_x(i) = tx;
    results_y(i) = ty;
    results_mag(i) = sqrt(tx^2 + ty^2);
    fprintf('Tilt = [%.2f°, %.2f°], mag = %.2f°\n', tx, ty, results_mag(i));
end

%% Find best scenario (minimize total tilt magnitude)
[min_tilt, best_idx] = min(results_mag);
fprintf('\nBest scenario: #%d with total tilt = %.2f°\n', best_idx, min_tilt);
fprintf('Best actuation factors: [%.2f, %.2f, %.2f]\n', test_scenarios(best_idx,:));
fprintf('Tilt components: X=%.2f°, Y=%.2f°\n', results_x(best_idx), results_y(best_idx));

%% Plot results
figure('Name', 'Tilt Comparison');

subplot(2,2,1);
bar(results_x);
xlabel('Scenario');
ylabel('X-Tilt (degrees)');
title('Tilt about X-axis (Pitch)');
grid on;

subplot(2,2,2);
bar(results_y);
xlabel('Scenario');
ylabel('Y-Tilt (degrees)');
title('Tilt about Y-axis (Roll)');
grid on;

subplot(2,2,3);
bar(results_mag);
xlabel('Scenario');
ylabel('Total Tilt (degrees)');
title('Total Tilt Magnitude');
grid on;
hold on;
plot(best_idx, min_tilt, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
legend('Tilt', 'Best');
hold off;

subplot(2,2,4);
plot(results_x, results_y, 'bo-', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
plot(results_x(best_idx), results_y(best_idx), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
plot(0, 0, 'gx', 'MarkerSize', 15, 'LineWidth', 3);
grid on;
xlabel('X-Tilt (degrees)');
ylabel('Y-Tilt (degrees)');
title('Tilt Components (2D view)');
legend('Scenarios', 'Best', 'Target (0,0)');
axis equal;
hold off;
