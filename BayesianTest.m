% Simple Bayesian Optimization Demo with Epsilon-Greedy Exploration
% Objective: Minimize f(x) = x * sin(x) on interval [0, 10]
clear; close all; clc;
%% Define the true objective function
f = @(x) x .* sin(x);
%% Settings
x_min = 0;
x_max = 8;
n_iterations = 20;  % Number of optimization iterations
initial_guess = 2;

% EXPLORATION PARAMETERS
kappa = 2.0;          % UCB exploration parameter
epsilon = 0.2;        % Probability of pure random exploration (0.2 = 20% chance)
add_distance_bonus = true;  % Add bonus for unexplored regions

%% Initialize with the starting point
X_observed = initial_guess;
Y_observed = f(initial_guess);
%% Create fine grid for plotting
x_plot = linspace(x_min, x_max, 200)';
y_true = f(x_plot);
%% Bayesian Optimization Loop
figure('Position', [100, 100, 1200, 800]);
for iter = 1:n_iterations
    % Fit Gaussian Process to observed data
    gprMdl = fitrgp(X_observed, Y_observed, ...
        'KernelFunction', 'squaredexponential', ...
        'Sigma', 0.1);  % Moderate noise to maintain uncertainty
    
    % Predict mean and standard deviation over the domain
    [y_pred, y_std] = predict(gprMdl, x_plot);
    
    % Epsilon-greedy: random exploration with probability epsilon
    use_random = rand() < epsilon;
    
    if use_random
        % Pure random exploration
        idx_next = randi(length(x_plot));
        x_next = x_plot(idx_next);
        fprintf('Iteration %d: RANDOM EXPLORATION at x = %.3f\n', iter, x_next);
    else
        % Acquisition function: Upper Confidence Bound (UCB)
        acquisition = y_pred - kappa * y_std;
        
        % OPTIONAL: Add distance bonus to favor unexplored regions
        if add_distance_bonus
            % Calculate minimum distance to any observed point
            min_distances = zeros(size(x_plot));
            for i = 1:length(x_plot)
                min_distances(i) = min(abs(x_plot(i) - X_observed));
            end
            % Normalize to [0, 1]
            distance_bonus = min_distances / max(min_distances);
            
            % Weight to control distance bonus influence (negative = favor distant points)
            distance_weight = -0.5;
            acquisition = acquisition + distance_weight * distance_bonus;
        end
        
        % Find next point to sample (minimize acquisition)
        [~, idx_next] = min(acquisition);
        x_next = x_plot(idx_next);
        fprintf('Iteration %d: UCB EXPLOITATION at x = %.3f\n', iter, x_next);
    end
    
    y_next = f(x_next);
    
    % Add new observation
    X_observed = [X_observed; x_next];
    Y_observed = [Y_observed; y_next];
    
    % Plotting
    clf;
    
    % Subplot 1: GP mean and uncertainty
    subplot(2,1,1);
    hold on;
    plot(x_plot, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True Function');
    plot(x_plot, y_pred, 'b-', 'LineWidth', 1.5, 'DisplayName', 'GP Mean');
    fill([x_plot; flipud(x_plot)], ...
         [y_pred + 2*y_std; flipud(y_pred - 2*y_std)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% CI');
    plot(X_observed(1:end-1), Y_observed(1:end-1), 'ro', ...
         'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Observed Points');
    if use_random
        plot(x_next, y_next, 'm^', 'MarkerSize', 12, ...
             'MarkerFaceColor', 'm', 'DisplayName', 'Random Exploration');
    else
        plot(x_next, y_next, 'g^', 'MarkerSize', 12, ...
             'MarkerFaceColor', 'g', 'DisplayName', 'UCB Selection');
    end
    xlabel('x');
    ylabel('f(x) = x sin(x)');
    title(sprintf('Iteration %d: GP Model (kappa=%.1f, epsilon=%.2f)', iter, kappa, epsilon));
    legend('Location', 'best');
    grid on;
    
    % Subplot 2: Acquisition function (UCB)
    subplot(2,1,2);
    if ~use_random
        plot(x_plot, acquisition, 'r-', 'LineWidth', 1.5);
        hold on;
        plot(x_next, acquisition(idx_next), 'g^', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
        ylabel('UCB Acquisition');
        title('Acquisition Function (Lower is Better)');
    else
        % Show all points have equal probability during random exploration
        plot(x_plot, ones(size(x_plot)), 'r-', 'LineWidth', 1.5);
        hold on;
        plot(x_next, 1, 'm^', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
        ylabel('Random Selection Probability');
        title('RANDOM EXPLORATION MODE');
    end
    xlabel('x');
    grid on;
    
    drawnow;
    pause(0.5);  % Pause to see the progression
end
%% Final results
[y_min, idx_min] = min(Y_observed);
x_min_found = X_observed(idx_min);
fprintf('\n=== Bayesian Optimization Results ===\n');
fprintf('Number of iterations: %d\n', n_iterations);
fprintf('Best x found: %.4f\n', x_min_found);
fprintf('Best f(x) found: %.4f\n', y_min);
% Find true minimum for comparison
[y_true_min, idx_true] = min(y_true);
x_true_min = x_plot(idx_true);
fprintf('\nTrue minimum (from fine grid):\n');
fprintf('x = %.4f, f(x) = %.4f\n', x_true_min, y_true_min);
%% Plot convergence
figure('Position', [100, 100, 800, 400]);
plot(1:length(Y_observed), cummin(Y_observed), 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;
yline(y_true_min, 'r--', 'LineWidth', 1.5, 'Label', 'True Minimum');
xlabel('Iteration');
ylabel('Best f(x) Found');
title('Convergence of Bayesian Optimization');
grid on;