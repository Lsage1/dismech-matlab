% 2D Bayesian Optimization with Failure Handling
% Objective: Minimize f(x,y) = 0.5*x*sin(x) + sin(y) on [0,8] x [0,8]
% Constraint: Simulations fail when (x-4)^2 + (y-4)^2 > 10
clear; close all; clc;

%% Define the true objective function
f = @(x, y) 0.5 * x .* sin(x) + sin(y);

%% Define failure region
% Returns true if simulation will fail
will_fail = @(x, y) (x - 4).^2 + (y - 4).^2 > 10;

%% Simulation function with failure handling
function [value, success] = simulate(x, y, f, will_fail)
    if will_fail(x, y)
        value = NaN;
        success = false;
        fprintf('  FAILED at [%.3f, %.3f]\n', x, y);
    else
        value = f(x, y);
        success = true;
        fprintf('  SUCCESS at [%.3f, %.3f], f = %.4f\n', x, y, value);
    end
end

%% Settings
x_min = 0;
x_max = 8;
y_min = 0;
y_max = 8;
n_iterations = 40;
initial_guess = [2, 2];  % Starting point [x, y]

% EXPLORATION PARAMETERS
kappa = 2.0;          % UCB exploration parameter
epsilon = 0.15;       % Probability of pure random exploration

% FAILURE HANDLING STRATEGY
USE_CONSTRAINT_MODEL = true;  % Model probability of success
alpha = 0.1;                  % Weight for constraint vs objective (0=ignore constraints, 1=only feasibility)

%% Initialize with the starting point
X_observed = initial_guess;
[z_init, success_init] = simulate(initial_guess(1), initial_guess(2), f, will_fail);
Y_observed = z_init;
Success_flags = success_init;

%% Create fine grid for plotting
n_grid = 50;
x_plot = linspace(x_min, x_max, n_grid);
y_plot = linspace(y_min, y_max, n_grid);
[X_grid, Y_grid] = meshgrid(x_plot, y_plot);

% Evaluate true function and failure region on grid
Z_true = f(X_grid, Y_grid);
Failure_region = will_fail(X_grid, Y_grid);

%% Bayesian Optimization Loop
figure('Position', [50, 50, 1600, 900]);

for iter = 1:n_iterations
    fprintf('\n--- Iteration %d ---\n', iter);
    
    % Separate successful and failed observations
    idx_success = Success_flags;
    X_success = X_observed(idx_success, :);
    Y_success = Y_observed(idx_success);
    X_failed = X_observed(~idx_success, :);
    
    fprintf('Success rate: %.1f%% (%d/%d)\n', ...
            100*sum(Success_flags)/length(Success_flags), sum(Success_flags), length(Success_flags));
    
    % Check if we have enough successful points
    if sum(idx_success) < 3
        fprintf('Not enough successful points. Random exploration.\n');
        % Random exploration near successful points
        if sum(idx_success) > 0
            random_idx = randi(sum(idx_success));
            x_next = X_success(random_idx, 1) + randn() * 0.5;
            y_next = X_success(random_idx, 2) + randn() * 0.5;
        else
            x_next = x_min + rand() * (x_max - x_min);
            y_next = y_min + rand() * (y_max - y_min);
        end
        x_next = max(x_min, min(x_max, x_next));
        y_next = max(y_min, min(y_max, y_next));
    else
        % Fit GP to SUCCESSFUL observations only
        gprMdl_objective = fitrgp(X_success, Y_success, ...
            'KernelFunction', 'squaredexponential', ...
            'Sigma', 0.1);
        
        % Fit CONSTRAINT MODEL: GP to predict probability of success
        % Use 1 for success, 0 for failure
        Y_constraint = double(Success_flags);
        gprMdl_constraint = fitrgp(X_observed, Y_constraint, ...
            'KernelFunction', 'squaredexponential', ...
            'Sigma', 0.05, ...
            'Beta', 1.0, ...  % Prior mean = 1 (optimistic: assume safe until proven otherwise)
            'FitMethod', 'exact');
        
        % Predict over grid
        X_test = [X_grid(:), Y_grid(:)];
        
        % Objective predictions (from successful points only)
        [z_pred_vec, z_std_vec] = predict(gprMdl_objective, X_test);
        Z_pred = reshape(z_pred_vec, size(X_grid));
        Z_std = reshape(z_std_vec, size(X_grid));
        
        % Constraint predictions (probability of success)
        % Optimistic prior: starts at 1, decreases near observed failures
        [p_success_vec, p_std_vec] = predict(gprMdl_constraint, X_test);
        
        % Apply optimistic adjustment: in unexplored regions (high uncertainty), assume success
        % Where uncertainty is high and no nearby observations, push toward 1
        p_success_vec_adjusted = p_success_vec;
        
        % Calculate distance to nearest observation for each test point
        min_dist_to_obs = zeros(size(p_success_vec));
        for i = 1:length(p_success_vec)
            distances = sqrt(sum((X_observed - X_test(i, :)).^2, 2));
            min_dist_to_obs(i) = min(distances);
        end
        
        % Normalize distances
        max_dist = sqrt((x_max - x_min)^2 + (y_max - y_min)^2);
        norm_dist = min_dist_to_obs / max_dist;
        
        % In regions far from observations, blend toward optimistic prior (1.0)
        % weight: 0 = use GP prediction, 1 = use prior
        optimism_weight = min(1, norm_dist * 3);  % Increase optimism with distance
        p_success_vec_adjusted = (1 - optimism_weight) .* p_success_vec + optimism_weight * 1.0;
        
        P_success = reshape(p_success_vec_adjusted, size(X_grid));
        P_success = max(0, min(1, P_success));  % Clamp to [0, 1]
        
        % Epsilon-greedy
        use_random = rand() < epsilon;
        
        if use_random
            % Random exploration, but biased away from failure regions
            % Sample with probability proportional to P_success
            prob_sample = P_success(:) / sum(P_success(:));
            idx_next = randsample(length(prob_sample), 1, true, prob_sample);
            x_next = X_test(idx_next, 1);
            y_next = X_test(idx_next, 2);
            fprintf('RANDOM EXPLORATION (weighted by success probability)\n');
        else
            % UCB Acquisition with constraint handling
            % For MINIMIZATION: Lower Confidence Bound (LCB)
            Z_acquisition = Z_pred - kappa * Z_std;  % Lower values are better
            
            % NORMALIZE both terms to [0, 1] so they're on the same scale
            % Normalize acquisition function
            Z_acq_min = min(Z_acquisition(:));
            Z_acq_max = max(Z_acquisition(:));
            Z_acquisition_norm = (Z_acquisition - Z_acq_min) / (Z_acq_max - Z_acq_min);
            
            % Constraint term (already in [0,1] but invert: 1-P_success)
            constraint_term = 1 - P_success;
            
            % METHOD 1: Weighted combination (both normalized to [0,1])
            % Now alpha is interpretable: how much we care about feasibility vs optimality
            alpha = 0.5;  % Weight for constraint (0 = ignore constraints, 1 = only care about feasibility)
            Z_acquisition_constrained = (1 - alpha) * Z_acquisition_norm + alpha * constraint_term;
            
            % METHOD 2: Hard constraint (only consider points with P_success > threshold)
            % Uncomment to use this instead:
            % threshold = 0.5;
            % Z_acquisition_constrained = Z_acquisition;
            % Z_acquisition_constrained(P_success < threshold) = inf;
            
            % Find minimum
            [~, idx_next] = min(Z_acquisition_constrained(:));
            x_next = X_test(idx_next, 1);
            y_next = X_test(idx_next, 2);
            fprintf('UCB EXPLOITATION (P_success = %.2f at selected point)\n', P_success(idx_next));
        end
    end
    
    % Evaluate
    [z_next, success_next] = simulate(x_next, y_next, f, will_fail);
    
    % Add observation
    X_observed = [X_observed; x_next, y_next];
    Y_observed = [Y_observed; z_next];
    Success_flags = [Success_flags; success_next];
    
    % Track best so far
    if success_next
        current_best = min(Y_observed(Success_flags));
        fprintf('  Current best objective: %.4f\n', current_best);
    end
    
    %% Plotting
    if sum(idx_success) >= 3
        clf;
        
        % Update indices
        idx_success = Success_flags;
        X_success = X_observed(idx_success, :);
        Y_success = Y_observed(idx_success);
        X_failed = X_observed(~idx_success, :);
        
        % Subplot 1: True function with failure region
        subplot(2, 3, 1);
        surf(X_grid, Y_grid, Z_true, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on;
        % Show failure region
        Z_failure_viz = Z_true;
        Z_failure_viz(~Failure_region) = NaN;
        surf(X_grid, Y_grid, Z_failure_viz, 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot3(X_success(:, 1), X_success(:, 2), Y_success, ...
              'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
        if ~isempty(X_failed)
            plot3(X_failed(:, 1), X_failed(:, 2), zeros(size(X_failed, 1), 1), ...
                  'rx', 'MarkerSize', 10, 'LineWidth', 2);
        end
        xlabel('x'); ylabel('y'); zlabel('f(x,y)');
        title('True Function (Red = Failure Region)');
        view(45, 30); grid on;
        
        % Subplot 2: GP Mean
        subplot(2, 3, 2);
        surf(X_grid, Y_grid, Z_pred, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        hold on;
        plot3(X_success(:, 1), X_success(:, 2), Y_success, ...
              'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
        xlabel('x'); ylabel('y'); zlabel('GP Mean');
        title(sprintf('Iter %d: Objective Model', iter));
        view(45, 30); grid on;
        
        % Subplot 3: Probability of Success
        subplot(2, 3, 3);
        surf(X_grid, Y_grid, P_success, 'EdgeColor', 'none');
        hold on;
        plot3(X_success(:, 1), X_success(:, 2), ones(size(X_success, 1), 1), ...
              'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
        if ~isempty(X_failed)
            plot3(X_failed(:, 1), X_failed(:, 2), zeros(size(X_failed, 1), 1), ...
                  'rx', 'MarkerSize', 10, 'LineWidth', 2);
        end
        xlabel('x'); ylabel('y'); zlabel('P(success)');
        title('Constraint Model (Success Probability)');
        colorbar;
        view(45, 30); grid on;
        
        % Subplot 4: Acquisition function
        subplot(2, 3, 4);
        if ~use_random
            surf(X_grid, Y_grid, Z_acquisition_constrained, 'EdgeColor', 'none');
            hold on;
            plot3(x_next, y_next, Z_acquisition_constrained(idx_next), ...
                  'g^', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
        end
        xlabel('x'); ylabel('y'); zlabel('Acquisition');
        title('Constrained Acquisition (Lower is Better)');
        colorbar;
        view(45, 30); grid on;
        
        % Subplot 5: Top view of probability of success
        subplot(2, 3, 5);
        contourf(X_grid, Y_grid, P_success, 20);
        hold on;
        plot(X_success(:, 1), X_success(:, 2), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
        if ~isempty(X_failed)
            plot(X_failed(:, 1), X_failed(:, 2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
        end
        plot(x_next, y_next, 'w^', 'MarkerSize', 12, 'MarkerFaceColor', 'w');
        xlabel('x'); ylabel('y');
        title('Success Probability (Top View)');
        colorbar;
        grid on;
        
        % Subplot 6: Convergence
        subplot(2, 3, 6);
        if sum(idx_success) > 0
            Y_success_cummin = zeros(sum(idx_success), 1);
            Y_success_cummin(1) = Y_success(1);
            for i = 2:sum(idx_success)
                Y_success_cummin(i) = min(Y_success_cummin(i-1), Y_success(i));
            end
            plot(1:sum(idx_success), Y_success_cummin, 'g-o', ...
                 'LineWidth', 2, 'MarkerFaceColor', 'g');
        end
        xlabel('Successful Iteration');
        ylabel('Best f(x,y) Found');
        title(sprintf('Convergence (%.1f%% success)', 100*sum(idx_success)/length(Success_flags)));
        grid on;
        
        drawnow;
        pause(0.3);
    end
end

%% Final Results
fprintf('\n\n=== Final Results ===\n');
fprintf('Total iterations: %d\n', n_iterations);
fprintf('Successful: %d (%.1f%%)\n', sum(Success_flags), 100*sum(Success_flags)/length(Success_flags));
fprintf('Failed: %d (%.1f%%)\n', sum(~Success_flags), 100*sum(~Success_flags)/length(Success_flags));

if sum(Success_flags) > 0
    idx_success = Success_flags;
    Y_success = Y_observed(idx_success);
    X_success = X_observed(idx_success, :);
    
    [y_min, idx_min] = min(Y_success);
    x_min_found = X_success(idx_min, :);
    
    fprintf('\nBest solution found:\n');
    fprintf('  [x, y] = [%.4f, %.4f]\n', x_min_found(1), x_min_found(2));
    fprintf('  f(x,y) = %.4f\n', y_min);
    
    % Check if it's in safe region
    if ~will_fail(x_min_found(1), x_min_found(2))
        fprintf('  Status: In SAFE region ✓\n');
    else
        fprintf('  Status: In FAILURE region (should not happen!)\n');
    end
end

%% Final visualization
figure('Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
contourf(X_grid, Y_grid, Z_true, 20);
hold on;
% Draw failure boundary
theta = linspace(0, 2*pi, 100);
r = sqrt(10);
x_circle = 4 + r * cos(theta);
y_circle = 4 + r * sin(theta);
plot(x_circle, y_circle, 'r-', 'LineWidth', 3);

idx_success = Success_flags;
X_success = X_observed(idx_success, :);
X_failed = X_observed(~idx_success, :);

plot(X_success(:, 1), X_success(:, 2), 'go-', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'LineWidth', 1);
if ~isempty(X_failed)
    plot(X_failed(:, 1), X_failed(:, 2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end
if sum(Success_flags) > 0
    [~, idx_min] = min(Y_observed(Success_flags));
    best_pt = X_success(idx_min, :);
    plot(best_pt(1), best_pt(2), 'r*', 'MarkerSize', 20, 'LineWidth', 3);
end
xlabel('x'); ylabel('y');
title('Optimization Path (Red Circle = Failure Boundary)');
colorbar;
legend('Success Path', 'Failures', 'Best Found', 'Location', 'best');
grid on;

subplot(1, 2, 2);
if sum(Success_flags) >= 3
    % Show final constraint model
    contourf(X_grid, Y_grid, P_success, 20);
    hold on;
    plot(x_circle, y_circle, 'r-', 'LineWidth', 3);
    plot(X_success(:, 1), X_success(:, 2), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    if ~isempty(X_failed)
        plot(X_failed(:, 1), X_failed(:, 2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    end
    xlabel('x'); ylabel('y');
    title('Learned Constraint Model (P(success))');
    colorbar;
    grid on;
end