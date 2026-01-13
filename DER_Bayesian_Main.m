% Bayesian Optimization for Discrete Elastic Rod Geometry
% Optimizes control points (x,y,z) to minimize average tilt across multiple scenarios
% First point is fixed, last point Z is fixed (20 dimensions total)
clear; close all; clc;

%% Add paths
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

%% ==================== CONFIGURATION ====================

% Initial control points (baseline geometry)
control_points_base = [
     0,      0,       .015;
    .007,    .005,       .015;
    .019,    .0040,       .007;
    .007,    .0011,       .005;
    .004,    .0007,       .003;
    .004,    .0005,   .002;
    .009,    .0005,   .002;
    .009,    0,     0;
];

% Define bounds for optimization
% X and Y: [0, 0.018]
% Z: [0, 0.015]
x_y_min = 0;
x_y_max = 0.018;
z_min = 0;
z_max = 0.015;

% Build parameter configuration
% Skip point 1 (fixed) and point 8 Z-coordinate (fixed)
param_config = [];
param_names = {};
x_min_vec = [];
x_max_vec = [];

for pt = 1:size(control_points_base, 1)
    for dim = 1:3  % x, y, z
        % Skip first point entirely
        if pt == 1
            continue;
        end
        
        % Skip last point's Z coordinate
        if pt == 8 && dim == 3
            continue;
        end
        
        % Add this parameter to optimization
        param_config = [param_config; pt, dim];
        param_names{end+1} = sprintf('P%d.%s', pt, char('X'+dim-1));
        
        % Set bounds based on dimension
        if dim == 3  % Z coordinate
            x_min_vec = [x_min_vec, z_min];
            x_max_vec = [x_max_vec, z_max];
        else  % X or Y coordinate
            x_min_vec = [x_min_vec, x_y_min];
            x_max_vec = [x_max_vec, x_y_max];
        end
    end
end

% Test scenarios (actuation factors for 3 tendons)
test_scenarios = [
    0.5867, 0.8330, 0.5867;   % Scenario 1: single tendon shortened
    0.3333, 0.9003, 0.9003;   % Scenario 2: double tendon shortened
    0.6946, 0.6946, 0.6946;   % Scenario 3: equal reduction for 3 tendons
];

% Simulation settings
edges_to_actuate = [1, 2, 3];
top_verts_ind = [5, 35, 65];

% Bayesian optimization settings
n_iterations = 1000;
kappa = 2.0;           % UCB exploration parameter (higher = more exploration)
epsilon = 0.2;         % Random exploration probability
alpha = 0.3;           % Constraint weight (0=ignore failures, 1=only feasibility)

% Create output directory for plots
plot_dir = 'optimization_plots';
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

%% ==================== INITIALIZATION ====================

n_params = size(param_config, 1);
x_min = x_min_vec;
x_max = x_max_vec;

fprintf('=== BAYESIAN OPTIMIZATION FOR ROD GEOMETRY ===\n');
fprintf('Optimizing %d control point parameters\n', n_params);
fprintf('  - Point 1: FIXED (all dimensions)\n');
fprintf('  - Points 2-7: X,Y,Z optimized\n');
fprintf('  - Point 8: X,Y optimized, Z FIXED\n');
fprintf('Bounds: X,Y ∈ [%.4f, %.4f], Z ∈ [%.4f, %.4f]\n', x_y_min, x_y_max, z_min, z_max);
fprintf('Testing %d actuation scenarios and averaging tilt\n', size(test_scenarios, 1));
fprintf('Running %d optimization iterations\n', n_iterations);
fprintf('Plots will be saved to: %s/\n\n', plot_dir);

% Start at baseline values for optimized parameters
x_initial = zeros(1, n_params);
for i = 1:n_params
    pt_idx = param_config(i, 1);
    dim_idx = param_config(i, 2);
    x_initial(i) = control_points_base(pt_idx, dim_idx);
end

% Evaluate initial geometry
fprintf('=== Initial Evaluation ===\n');
[avg_tilt_init, success_init, scenario_tilts_init] = evaluate_geometry(x_initial, ...
    control_points_base, param_config, edges_to_actuate, test_scenarios, top_verts_ind);

fprintf('Initial average tilt: %.4f°\n', avg_tilt_init);
fprintf('Scenario tilts: [%.4f°, %.4f°, %.4f°]\n', scenario_tilts_init);

% Initialize storage
X_observed = x_initial;
Y_observed = avg_tilt_init;
Success_flags = success_init;
Scenario_tilts = scenario_tilts_init;  % Store individual scenario results

%% ==================== OPTIMIZATION LOOP ====================

% Create invisible figure for plotting (won't display on screen)
fig_progress = figure('Visible', 'off', 'Position', [100, 100, 1600, 900]);

for iter = 1:n_iterations
    fprintf('\n========== Iteration %d/%d ==========\n', iter, n_iterations);
    
    % Get successful observations
    idx_success = Success_flags & ~isnan(Y_observed);
    n_success = sum(idx_success);
    
    fprintf('Success rate: %.1f%% (%d/%d)\n', ...
        100*n_success/length(Success_flags), n_success, length(Success_flags));
    
    if n_success > 0
        fprintf('Current best average tilt: %.4f°\n', min(Y_observed(idx_success)));
    end
    
    %% Select next point to evaluate
    if n_success < 3
        % Not enough data - random exploration
        fprintf('Strategy: RANDOM (insufficient data)\n');
        x_next = x_min + rand(1, n_params) .* (x_max - x_min);
        
    else
        % Fit Gaussian Processes
        X_success = X_observed(idx_success, :);
        Y_success = Y_observed(idx_success);
        
        % Objective GP (minimize average tilt)
        gprMdl_objective = fitrgp(X_success, Y_success, ...
            'KernelFunction', 'squaredexponential', 'Sigma', 0.1);
        
        % Constraint GP (predict success probability)
        Y_constraint = double(Success_flags);
        gprMdl_constraint = fitrgp(X_observed, Y_constraint, ...
            'KernelFunction', 'squaredexponential', 'Sigma', 0.05, ...
            'Beta', 1.0, 'FitMethod', 'exact');
        
        % Generate candidate points
        n_candidates = 2000;
        X_candidates = repmat(x_min, n_candidates, 1) + ...
            rand(n_candidates, n_params) .* repmat(x_max - x_min, n_candidates, 1);
        
        % Predict objective
        [y_pred, y_std] = predict(gprMdl_objective, X_candidates);
        
        % Predict constraint with optimistic prior for unexplored regions
        [p_success_raw, ~] = predict(gprMdl_constraint, X_candidates);
        
        % Calculate minimum distance to any observed point
        all_distances = pdist2(X_candidates, X_observed, 'euclidean');
        min_dist = min(all_distances, [], 2);
        
        % Ensure column vectors
        p_success_raw = p_success_raw(:);
        min_dist = min_dist(:);
        
        max_dist = norm(x_max - x_min);
        optimism = min(1, (min_dist / max_dist) * 3);
        p_success = (1 - optimism) .* p_success_raw + optimism;
        p_success = max(0, min(1, p_success));
        
        % Epsilon-greedy strategy
        if rand() < epsilon
            % Random exploration
            fprintf('Strategy: RANDOM EXPLORATION\n');
            prob_sample = p_success / sum(p_success);
            idx_next = randsample(n_candidates, 1, true, prob_sample);
            x_next = X_candidates(idx_next, :);
            fprintf('  P(success) = %.2f\n', p_success(idx_next));
        else
            % UCB acquisition
            fprintf('Strategy: UCB ACQUISITION\n');
            acquisition = y_pred - kappa * y_std;
            
            % Combine with constraint
            acq_norm = (acquisition - min(acquisition)) / (max(acquisition) - min(acquisition));
            constraint_term = 1 - p_success;
            combined = (1 - alpha) * acq_norm + alpha * constraint_term;
            
            [~, idx_next] = min(combined);
            x_next = X_candidates(idx_next, :);
            fprintf('  Predicted avg tilt: %.4f ± %.4f°\n', y_pred(idx_next), y_std(idx_next));
            fprintf('  P(success): %.2f\n', p_success(idx_next));
        end
    end
    
    %% Evaluate new geometry
    [avg_tilt_next, success_next, scenario_tilts_next] = evaluate_geometry(x_next, ...
        control_points_base, param_config, edges_to_actuate, test_scenarios, top_verts_ind);
    
    % Store results
    X_observed = [X_observed; x_next];
    Y_observed = [Y_observed; avg_tilt_next];
    Success_flags = [Success_flags; success_next];
    Scenario_tilts = [Scenario_tilts; scenario_tilts_next];
    
    %% Save plot periodically
    if mod(iter, 50) == 0 && sum(idx_success) >= 3
        fprintf('Saving progress plot...\n');
        plot_progress(fig_progress, X_observed, Y_observed, Success_flags, Scenario_tilts, ...
            iter, n_iterations, param_config);
        saveas(fig_progress, fullfile(plot_dir, sprintf('progress_iter_%04d.png', iter)));
        fprintf('  Saved to %s/progress_iter_%04d.png\n', plot_dir, iter);
    end
    
    %% Save checkpoint periodically
    if mod(iter, 100) == 0
        checkpoint.X_observed = X_observed;
        checkpoint.Y_observed = Y_observed;
        checkpoint.Success_flags = Success_flags;
        checkpoint.Scenario_tilts = Scenario_tilts;
        checkpoint.iteration = iter;
        checkpoint.timestamp = datetime('now');
        save(fullfile(plot_dir, sprintf('checkpoint_iter_%04d.mat', iter)), 'checkpoint');
        fprintf('Checkpoint saved at iteration %d\n', iter);
    end
end

close(fig_progress);  % Close the invisible figure

%% ==================== FINAL RESULTS ====================

fprintf('\n\n========================================\n');
fprintf('===   OPTIMIZATION COMPLETE   ===\n');
fprintf('========================================\n');

idx_success = Success_flags & ~isnan(Y_observed);
fprintf('Final success rate: %.1f%% (%d/%d)\n', ...
    100*sum(idx_success)/length(Success_flags), sum(idx_success), length(Success_flags));

if sum(idx_success) > 0
    [best_avg_tilt, best_idx_rel] = min(Y_observed(idx_success));
    X_success = X_observed(idx_success, :);
    best_params = X_success(best_idx_rel, :);
    
    % Find which iteration this was
    success_indices = find(idx_success);
    best_iter_abs = success_indices(best_idx_rel);
    best_scenario_tilts = Scenario_tilts(best_iter_abs, :);
    
    fprintf('\n*** BEST SOLUTION ***\n');
    fprintf('Minimum average tilt: %.4f°\n', best_avg_tilt);
    fprintf('Individual scenario tilts:\n');
    fprintf('  Scenario 1: %.4f°\n', best_scenario_tilts(1));
    fprintf('  Scenario 2: %.4f°\n', best_scenario_tilts(2));
    fprintf('  Scenario 3: %.4f°\n', best_scenario_tilts(3));
    
    % Reconstruct best control points (preserve fixed values)
    best_control_points = control_points_base;
    for i = 1:n_params
        pt_idx = param_config(i, 1);
        dim_idx = param_config(i, 2);
        best_control_points(pt_idx, dim_idx) = best_params(i);
    end
    
    fprintf('\nOptimized control points:\n');
    disp(array2table(best_control_points, 'VariableNames', {'X', 'Y', 'Z'}));
    
    fprintf('\nChange from baseline:\n');
    change = best_control_points - control_points_base;
    disp(array2table(change, 'VariableNames', {'dX', 'dY', 'dZ'}));
    
    % Highlight fixed values
    fprintf('\nFixed values (not optimized):\n');
    fprintf('  Point 1: [%.4f, %.4f, %.4f] (entire point fixed)\n', ...
        control_points_base(1,1), control_points_base(1,2), control_points_base(1,3));
    fprintf('  Point 8 Z: %.4f (z-coordinate fixed)\n', control_points_base(8,3));
    
    % Save results
    results.best_avg_tilt = best_avg_tilt;
    results.best_scenario_tilts = best_scenario_tilts;
    results.best_params = best_params;
    results.best_control_points = best_control_points;
    results.baseline_control_points = control_points_base;
    results.X_observed = X_observed;
    results.Y_observed = Y_observed;
    results.Success_flags = Success_flags;
    results.Scenario_tilts = Scenario_tilts;
    results.param_config = param_config;
    results.test_scenarios = test_scenarios;
    results.bounds.x_y_min = x_y_min;
    results.bounds.x_y_max = x_y_max;
    results.bounds.z_min = z_min;
    results.bounds.z_max = z_max;
    
    save('optimization_results.mat', 'results');
    fprintf('\nResults saved to optimization_results.mat\n');
    
    % Create and save final visualizations
    fprintf('\nGenerating final plots...\n');
    
    % Final progress plot
    fig_final_progress = figure('Visible', 'on', 'Position', [100, 100, 1600, 900]);
    plot_progress(fig_final_progress, X_observed, Y_observed, Success_flags, ...
        Scenario_tilts, n_iterations, n_iterations, param_config);
    saveas(fig_final_progress, fullfile(plot_dir, 'final_progress.png'));
    saveas(fig_final_progress, fullfile(plot_dir, 'final_progress.fig'));
    fprintf('  Saved final progress plot\n');
    
    % Final results summary
    fig_final_results = create_final_plot(results);
    saveas(fig_final_results, fullfile(plot_dir, 'final_results.png'));
    saveas(fig_final_results, fullfile(plot_dir, 'final_results.fig'));
    fprintf('  Saved final results plot\n');
    
    % Geometry comparison (detailed 3D view)
    fig_geometry = create_geometry_plot(results);
    saveas(fig_geometry, fullfile(plot_dir, 'geometry_comparison.png'));
    saveas(fig_geometry, fullfile(plot_dir, 'geometry_comparison.fig'));
    fprintf('  Saved geometry comparison plot\n');
    
    fprintf('\nAll plots saved to: %s/\n', plot_dir);
    fprintf('  - progress_iter_XXXX.png (periodic snapshots)\n');
    fprintf('  - final_progress.png/fig (complete optimization history)\n');
    fprintf('  - final_results.png/fig (summary of best solution)\n');
    fprintf('  - geometry_comparison.png/fig (3D geometry visualization)\n');
end

%% ==================== HELPER FUNCTIONS ====================

function [avg_tilt_magnitude, success, scenario_tilts] = evaluate_geometry(params, ...
    control_points_base, param_config, edges_to_actuate, test_scenarios, top_verts_ind)
    
    % Reconstruct control points from parameter vector
    control_points = control_points_base;
    
    % Update only the optimized parameters
    for i = 1:length(params)
        pt_idx = param_config(i, 1);
        dim_idx = param_config(i, 2);
        control_points(pt_idx, dim_idx) = params(i);
    end
    
    % Generate geometry
    try
        create_haptic_node_fun(control_points);
    catch ME
        fprintf('  ❌ FAILED (geometry): %s\n', ME.message);
        avg_tilt_magnitude = NaN;
        success = false;
        scenario_tilts = [NaN, NaN, NaN];
        return;
    end
    
    % Load robot description
    try
        evalin('base', 'robotDescriptionHapticDevice;');
    catch ME
        fprintf('  ❌ FAILED (robot description): %s\n', ME.message);
        avg_tilt_magnitude = NaN;
        success = false;
        scenario_tilts = [NaN, NaN, NaN];
        return;
    end
    
    % Get variables from base workspace
    try
        material = evalin('base', 'material');
        geom = evalin('base', 'geom');
        env = evalin('base', 'env');
        sim_params = evalin('base', 'sim_params');
        nodes = evalin('base', 'nodes');
        edges = evalin('base', 'edges');
        face_nodes = evalin('base', 'face_nodes');
        fixed_node_indices = evalin('base', 'fixed_node_indices');
        fixed_edge_indices = evalin('base', 'fixed_edge_indices');
    catch ME
        fprintf('  ❌ FAILED (getting variables from workspace): %s\n', ME.message);
        avg_tilt_magnitude = NaN;
        success = false;
        scenario_tilts = [NaN, NaN, NaN];
        return;
    end
    
    % Run simulation for each scenario
    n_scenarios = size(test_scenarios, 1);
    tilt_magnitudes = zeros(n_scenarios, 1);
    
    fprintf('  Testing %d scenarios:\n', n_scenarios);
    
    for s = 1:n_scenarios
        try
            actuation_factors = test_scenarios(s, :);
            
            % Enable warning tracking
            lastwarn('');
            warning('on', 'all');
            
            [tilt_x, tilt_y] = simulate_geometry(edges_to_actuate, actuation_factors, ...
                top_verts_ind, material, geom, env, sim_params, nodes, edges, ...
                face_nodes, fixed_node_indices, fixed_edge_indices);
            
            % Check for warnings indicating numerical problems
            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg)
                if contains(warnMsg, 'singular') || contains(warnMsg, 'badly scaled') || ...
                   contains(warnMsg, 'not converge') || contains(warnId, 'MATLAB:singularMatrix')
                    fprintf('  ❌ FAILED (scenario %d): Numerical instability detected - %s\n', s, warnMsg);
                    avg_tilt_magnitude = NaN;
                    success = false;
                    scenario_tilts = [NaN, NaN, NaN];
                    return;
                end
            end
            
            % Check for NaN or Inf in results
            if isnan(tilt_x) || isnan(tilt_y) || isinf(tilt_x) || isinf(tilt_y)
                fprintf('  ❌ FAILED (scenario %d): Invalid tilt values (NaN or Inf)\n', s);
                avg_tilt_magnitude = NaN;
                success = false;
                scenario_tilts = [NaN, NaN, NaN];
                return;
            end
            
            % Combine tilt components into magnitude
            tilt_magnitudes(s) = sqrt(tilt_x^2 + tilt_y^2);
            
            fprintf('    Scenario %d [%.2f, %.2f, %.2f]: Tilt = %.4f° (x=%.2f°, y=%.2f°)\n', ...
                s, actuation_factors, tilt_magnitudes(s), tilt_x, tilt_y);
            
        catch ME
            fprintf('  ❌ FAILED (scenario %d simulation): %s\n', s, ME.message);
            avg_tilt_magnitude = NaN;
            success = false;
            scenario_tilts = [NaN, NaN, NaN];
            return;
        end
    end
    
    % Calculate average tilt across scenarios
    avg_tilt_magnitude = mean(tilt_magnitudes);
    scenario_tilts = tilt_magnitudes';
    success = true;
    
    fprintf('  ✓ SUCCESS: Average tilt = %.4f°\n', avg_tilt_magnitude);
end

function plot_progress(fig_handle, X_obs, Y_obs, Success, Scenario_tilts, iter, n_iter, param_config)
    figure(fig_handle);
    clf;
    
    idx_ok = Success & ~isnan(Y_obs);
    Y_ok = Y_obs(idx_ok);
    Scenario_ok = Scenario_tilts(idx_ok, :);
    n_success = sum(idx_ok);
    
    % Convergence plot
    subplot(2, 3, 1);
    if n_success > 0
        Y_cummin = cummin(Y_ok);
        plot(1:length(Y_cummin), Y_cummin, 'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'g');
        xlabel('Successful Iteration');
        ylabel('Best Avg Tilt (°)');
        title(sprintf('Convergence (Iter %d/%d)', iter, n_iter));
        grid on;
    else
        text(0.5, 0.5, 'No successful runs yet', 'HorizontalAlignment', 'center');
        title('Convergence');
    end
    
    % Success rate
    subplot(2, 3, 2);
    window = min(5, length(Success));
    if length(Success) >= window
        success_rate = movmean(double(Success), window) * 100;
        plot(success_rate, 'b-', 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Success Rate (%)');
        title('Rolling Success Rate');
        ylim([0, 105]);
        grid on;
    else
        text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center');
        title('Rolling Success Rate');
    end
    
    % Tilt histogram
    subplot(2, 3, 3);
    if n_success > 0
        histogram(Y_ok, 10, 'FaceColor', [0.2 0.8 0.2]);
        xlabel('Average Tilt (°)');
        ylabel('Count');
        title('Avg Tilt Distribution');
        grid on;
    else
        text(0.5, 0.5, 'No successful runs yet', 'HorizontalAlignment', 'center');
        title('Avg Tilt Distribution');
    end
    
    % All observations timeline
    subplot(2, 3, 4);
    hold on;
    if any(~idx_ok)
        plot(find(~idx_ok), Y_obs(~idx_ok), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    end
    if any(idx_ok)
        plot(find(idx_ok), Y_ok, 'go-', 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
    end
    xlabel('Iteration');
    ylabel('Avg Tilt (°)');
    title('All Observations');
    if any(~idx_ok) && any(idx_ok)
        legend('Failed', 'Success', 'Location', 'best');
    end
    grid on;
    hold off;
    
    % Scenario comparison
    subplot(2, 3, 5);
    if n_success > 0
        boxplot(Scenario_ok, 'Labels', {'Scenario 1', 'Scenario 2', 'Scenario 3'});
        ylabel('Tilt (°)');
        title('Tilt Distribution by Scenario');
        grid on;
    else
        text(0.5, 0.5, 'No successful runs yet', 'HorizontalAlignment', 'center');
        title('Tilt Distribution by Scenario');
    end
    
    % Parameter changes
    subplot(2, 3, 6);
    if n_success > 3
        X_ok = X_obs(idx_ok, :);
        param_std = std(X_ok);
        [~, sorted_idx] = sort(param_std, 'descend');
        
        n_show = min(6, length(sorted_idx));
        boxplot(X_ok(:, sorted_idx(1:n_show)));
        
        labels = cell(n_show, 1);
        for i = 1:n_show
            idx = sorted_idx(i);
            labels{i} = sprintf('P%d.%d', param_config(idx, 1), param_config(idx, 2));
        end
        set(gca, 'XTickLabel', labels);
        xtickangle(45);
        ylabel('Parameter Value');
        title('Most Variable Parameters');
        grid on;
    else
        text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center');
        title('Most Variable Parameters');
    end
    
    drawnow;
end

function fig = create_final_plot(results)
    fig = figure('Position', [150, 150, 1600, 600], 'Name', 'Final Results');
    
    idx_ok = results.Success_flags & ~isnan(results.Y_observed);
    Y_ok = results.Y_observed(idx_ok);
    
    % Convergence
    subplot(1, 4, 1);
    Y_cummin = cummin(Y_ok);
    plot(1:length(Y_cummin), Y_cummin, 'g-o', 'LineWidth', 2.5, ...
        'MarkerFaceColor', 'g', 'MarkerSize', 7);
    xlabel('Successful Iteration', 'FontSize', 11);
    ylabel('Best Avg Tilt (°)', 'FontSize', 11);
    title('Optimization Progress', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Control points comparison (3D)
    subplot(1, 4, 2);
    plot3(results.baseline_control_points(:,1), results.baseline_control_points(:,2), ...
        results.baseline_control_points(:,3), 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8);
    hold on;
    plot3(results.best_control_points(:,1), results.best_control_points(:,2), ...
        results.best_control_points(:,3), 'r-^', 'LineWidth', 2.5, 'MarkerSize', 8);
    
    % Highlight fixed point 1
    plot3(results.baseline_control_points(1,1), results.baseline_control_points(1,2), ...
        results.baseline_control_points(1,3), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k');
    
    xlabel('X', 'FontSize', 11);
    ylabel('Y', 'FontSize', 11);
    zlabel('Z', 'FontSize', 11);
    title('Control Points', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Baseline', 'Optimized', 'Fixed Pt 1', 'Location', 'best');
    grid on;
    view(45, 30);
    hold off;
    
    % Scenario comparison
    subplot(1, 4, 3);
    baseline_scenarios = results.Scenario_tilts(1, :);
    best_scenarios = results.best_scenario_tilts;
    
    x_pos = 1:3;
    bar_data = [baseline_scenarios; best_scenarios]';
    b = bar(x_pos, bar_data);
    b(1).FaceColor = [0.3 0.3 0.8];
    b(2).FaceColor = [0.8 0.2 0.2];
    
    set(gca, 'XTickLabel', {'Scenario 1', 'Scenario 2', 'Scenario 3'});
    ylabel('Tilt (°)', 'FontSize', 11);
    title('Scenario Comparison', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Baseline', 'Optimized', 'Location', 'best');
    grid on;
    xtickangle(45);
    
    % Largest parameter changes
    subplot(1, 4, 4);
    all_changes = zeros(size(results.baseline_control_points));
    all_changes(:) = results.best_control_points(:) - results.baseline_control_points(:);
    
    param_changes = [];
    param_labels_all = {};
    for i = 1:size(results.param_config, 1)
        pt = results.param_config(i, 1);
        dim = results.param_config(i, 2);
        param_changes(end+1) = abs(all_changes(pt, dim));
        param_labels_all{end+1} = sprintf('P%d.%s', pt, char('X'+dim-1));
    end
    
    [sorted_changes, sorted_idx] = sort(param_changes, 'descend');
    
    n_show = min(8, length(sorted_changes));
    bar(sorted_changes(1:n_show));
    
    labels = param_labels_all(sorted_idx(1:n_show));
    set(gca, 'XTickLabel', labels);
    ylabel('Absolute Change', 'FontSize', 11);
    title('Largest Parameter Changes', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    xtickangle(45);
end

function fig = create_geometry_plot(results)
    fig = figure('Position', [200, 200, 1200, 800], 'Name', 'Geometry Comparison');
    
    % Multiple views of the geometry
    views = [45, 30; 0, 0; 90, 0; 0, 90];
    view_names = {'Isometric', 'Front', 'Side', 'Top'};
    
    for v = 1:4
        subplot(2, 2, v);
        
        % Plot baseline
        plot3(results.baseline_control_points(:,1), ...
              results.baseline_control_points(:,2), ...
              results.baseline_control_points(:,3), ...
              'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        hold on;
        
        % Plot optimized
        plot3(results.best_control_points(:,1), ...
              results.best_control_points(:,2), ...
              results.best_control_points(:,3), ...
              'r-^', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        
        % Highlight fixed point 1
        plot3(results.baseline_control_points(1,1), ...
              results.baseline_control_points(1,2), ...
              results.baseline_control_points(1,3), ...
              'ko', 'MarkerSize', 14, 'MarkerFaceColor', 'k', 'LineWidth', 2);
        
        % Draw vectors showing movement
        for i = 2:size(results.baseline_control_points, 1)
            quiver3(results.baseline_control_points(i,1), ...
                   results.baseline_control_points(i,2), ...
                   results.baseline_control_points(i,3), ...
                   results.best_control_points(i,1) - results.baseline_control_points(i,1), ...
                   results.best_control_points(i,2) - results.baseline_control_points(i,2), ...
                   results.best_control_points(i,3) - results.baseline_control_points(i,3), ...
                   0, 'g', 'LineWidth', 1.5, 'MaxHeadSize', 2);
        end
        
        % Label points
        for i = 1:size(results.baseline_control_points, 1)
            text(results.baseline_control_points(i,1), ...
                 results.baseline_control_points(i,2), ...
                 results.baseline_control_points(i,3), ...
                 sprintf(' %d', i), 'FontSize', 9, 'Color', 'b', 'FontWeight', 'bold');
        end
        
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        title(sprintf('%s View', view_names{v}), 'FontWeight', 'bold');
        legend('Baseline', 'Optimized', 'Fixed', 'Movement', 'Location', 'best');
        grid on;
        axis equal;
        view(views(v, :));
        hold off;
    end
    
    % Add overall title
    sgtitle(sprintf('Geometry Comparison: %.4f° → %.4f° (%.1f%% improvement)', ...
        results.Y_observed(1), results.best_avg_tilt, ...
        100*(results.Y_observed(1) - results.best_avg_tilt)/results.Y_observed(1)), ...
        'FontSize', 14, 'FontWeight', 'bold');
end