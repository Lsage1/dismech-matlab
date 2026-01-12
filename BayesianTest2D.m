% 2D Bayesian Optimization Demo
% Objective: Minimize f(x,y) = 0.5*x*sin(x) + sin(y) on [0,8] x [0,8]
clear; close all; clc;

%% Define the true objective function
f = @(x, y) 0.5 * x .* sin(x) + sin(y);

%% Settings
x_min = 0;
x_max = 8;
y_min = 0;
y_max = 8;
n_iterations = 30;  % Number of optimization iterations
initial_guess = [2, 2];  % Starting point [x, y]

% EXPLORATION PARAMETERS
kappa = 4.0;          % UCB exploration parameter
epsilon = 0.2;        % Probability of pure random exploration

%% Initialize with the starting point
X_observed = initial_guess;  % [x, y] pair
Y_observed = f(initial_guess(1), initial_guess(2));  % Function value

%% Create fine grid for plotting
n_grid = 50;  % Grid resolution
x_plot = linspace(x_min, x_max, n_grid);
y_plot = linspace(y_min, y_max, n_grid);
[X_grid, Y_grid] = meshgrid(x_plot, y_plot);

% Evaluate true function on grid
Z_true = f(X_grid, Y_grid);

%% Bayesian Optimization Loop
figure('Position', [100, 100, 1400, 600]);

for iter = 1:n_iterations
    fprintf('\n--- Iteration %d ---\n', iter);
    
    % Fit Gaussian Process to observed data
    % For 2D input, we concatenate [x, y] as rows
    gprMdl = fitrgp(X_observed, Y_observed, ...
        'KernelFunction', 'squaredexponential', ...
        'Sigma', 0.1);  % Noise parameter
    
    % Predict over the entire grid
    % Reshape grid to Nx2 matrix for GP prediction
    X_test = [X_grid(:), Y_grid(:)];  % Each row is [x, y]
    [z_pred_vec, z_std_vec] = predict(gprMdl, X_test);
    
    % Reshape predictions back to grid
    Z_pred = reshape(z_pred_vec, size(X_grid));
    Z_std = reshape(z_std_vec, size(X_grid));
    
    % Epsilon-greedy: random exploration with probability epsilon
    use_random = rand() < epsilon;
    
    if use_random
        % Pure random exploration
        x_next = x_min + rand() * (x_max - x_min);
        y_next = y_min + rand() * (y_max - y_min);
        fprintf('RANDOM EXPLORATION at [%.3f, %.3f]\n', x_next, y_next);
    else
        % Acquisition function: Upper Confidence Bound (UCB)
        Z_acquisition = Z_pred - kappa * Z_std;
        
        % Find minimum acquisition value
        [~, idx_min] = min(Z_acquisition(:));
        x_next = X_test(idx_min, 1);
        y_next = X_test(idx_min, 2);
        fprintf('UCB EXPLOITATION at [%.3f, %.3f]\n', x_next, y_next);
    end
    
    % Evaluate the true function
    z_next = f(x_next, y_next);
    fprintf('Function value: %.4f\n', z_next);
    
    % Add new observation
    X_observed = [X_observed; x_next, y_next];
    Y_observed = [Y_observed; z_next];
    
    %% Plotting
    clf;
    
    % Subplot 1: True function with observations
    subplot(1, 3, 1);
    surf(X_grid, Y_grid, Z_true, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on;
    % Plot all observed points
    plot3(X_observed(1:end-1, 1), X_observed(1:end-1, 2), Y_observed(1:end-1), ...
          'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
    % Plot new point
    if use_random
        plot3(x_next, y_next, z_next, 'm^', 'MarkerSize', 12, 'MarkerFaceColor', 'm', 'LineWidth', 1.5);
    else
        plot3(x_next, y_next, z_next, 'g^', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
    end
    xlabel('x');
    ylabel('y');
    zlabel('f(x,y)');
    title('True Function');
    colorbar;
    view(45, 30);
    grid on;
    
    % Subplot 2: GP Mean Prediction
    subplot(1, 3, 2);
    surf(X_grid, Y_grid, Z_pred, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on;
    plot3(X_observed(1:end-1, 1), X_observed(1:end-1, 2), Y_observed(1:end-1), ...
          'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
    if use_random
        plot3(x_next, y_next, z_next, 'm^', 'MarkerSize', 12, 'MarkerFaceColor', 'm', 'LineWidth', 1.5);
    else
        plot3(x_next, y_next, z_next, 'g^', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
    end
    xlabel('x');
    ylabel('y');
    zlabel('GP Mean');
    title(sprintf('Iteration %d: GP Prediction', iter));
    colorbar;
    view(45, 30);
    grid on;
    
    % Subplot 3: Acquisition Function
    subplot(1, 3, 3);
    if use_random
        % Show uniform distribution during random exploration
        Z_acquisition = ones(size(Z_pred));
        surf(X_grid, Y_grid, Z_acquisition, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        title('RANDOM EXPLORATION');
    else
        Z_acquisition = Z_pred - kappa * Z_std;
        surf(X_grid, Y_grid, Z_acquisition, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        hold on;
        plot3(x_next, y_next, Z_acquisition(idx_min), 'g^', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
        title('UCB Acquisition (Lower is Better)');
    end
    xlabel('x');
    ylabel('y');
    zlabel('Acquisition');
    colorbar;
    view(45, 30);
    grid on;
    
    drawnow;
    pause(0.5);
end

%% Final results
[y_min, idx_min] = min(Y_observed);
x_min_found = X_observed(idx_min, :);

fprintf('\n=== Bayesian Optimization Results ===\n');
fprintf('Number of iterations: %d\n', n_iterations);
fprintf('Best point found: [%.4f, %.4f]\n', x_min_found(1), x_min_found(2));
fprintf('Best f(x,y) found: %.4f\n', y_min);

% Find true minimum from grid
[z_true_min, idx_true] = min(Z_true(:));
x_true_min = X_grid(idx_true);
y_true_min = Y_grid(idx_true);
fprintf('\nTrue minimum (from grid):\n');
fprintf('Point: [%.4f, %.4f]\n', x_true_min, y_true_min);
fprintf('Value: %.4f\n', z_true_min);

%% Plot convergence
figure('Position', [100, 100, 800, 400]);
plot(1:length(Y_observed), cummin(Y_observed), 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;
yline(z_true_min, 'r--', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Best f(x,y) Found');
title('Convergence of 2D Bayesian Optimization');
legend('Best Found', 'True Minimum', 'Location', 'best');
grid on;

%% Final 3D visualization
figure('Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
surf(X_grid, Y_grid, Z_true, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on;
plot3(X_observed(:, 1), X_observed(:, 2), Y_observed, ...
      'ko-', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'LineWidth', 1);
plot3(x_min_found(1), x_min_found(2), y_min, ...
      'r*', 'MarkerSize', 20, 'LineWidth', 3);
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Optimization Path on True Function');
colorbar;
view(45, 30);
grid on;
legend('True Function', 'Sample Path', 'Best Found', 'Location', 'best');

subplot(1, 2, 2);
contourf(X_grid, Y_grid, Z_true, 20);
hold on;
plot(X_observed(:, 1), X_observed(:, 2), 'ko-', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'LineWidth', 1);
plot(x_min_found(1), x_min_found(2), 'r*', 'MarkerSize', 20, 'LineWidth', 3);
xlabel('x');
ylabel('y');
title('Optimization Path (Top View)');
colorbar;
grid on;
legend('Sample Path', 'Best Found', 'Location', 'best');