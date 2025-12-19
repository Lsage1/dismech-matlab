% Nurbs Curve test

% Generate NURBS curve with evenly spaced points and plot in 3D

%% 1. Define control points (example: your cone branch)
% You can replace these with your own points
control_points = [
     0,      0,       .015;
    .007,    .005,       .015;
    .019,    .0040,       .007;
    .007,    .0011,       .005;
    .004,    .0009,       .004;
    .004,    .0007,       .003;
    .004,    .0005,   .002;
    .009,    .0005,   .002;
    .009,    0.00025, .001;
    .009,    0,     0;
    
    
];

% Or use your quarter circle points
% control_points = [...your points here...];

%% 2. Create NURBS curve
% Number of control points
n_ctrl = size(control_points, 1);

% Degree of the curve (typically 3 for cubic NURBS)
degree = min(3, n_ctrl - 1);

% Create CLAMPED knot vector (ensures curve passes through endpoints)
n_internal = n_ctrl - degree - 1;
if n_internal > 0
    internal_knots = linspace(0, 1, n_internal + 2);
    internal_knots = internal_knots(2:end-1); % Remove endpoints
else
    internal_knots = [];
end
knots = [zeros(1, degree+1), internal_knots, ones(1, degree+1)];

% Weights (all equal to 1 for non-rational B-spline, i.e., standard spline)
weights = ones(n_ctrl, 1);

%% 3. Evaluate NURBS curve at many points for smooth visualization
n_eval = 1000;  % Number of evaluation points for smooth curve
u_eval = linspace(0, 1, n_eval);

% Evaluate curve using De Boor's algorithm or basis functions
curve_points = evaluateNURBS(control_points, weights, knots, degree, u_eval);

%% 4. Place evenly spaced points along the curve
n_spaced = 25;  % Number of evenly spaced points desired

% Calculate arc length parameterization
arc_lengths = zeros(n_eval, 1);
for i = 2:n_eval
    arc_lengths(i) = arc_lengths(i-1) + norm(curve_points(i,:) - curve_points(i-1,:));
end
total_length = arc_lengths(end);

% Find parameters for evenly spaced points
target_lengths = linspace(0, total_length, n_spaced);
u_spaced = zeros(n_spaced, 1);

for i = 1:n_spaced
    % Find parameter u that gives the target arc length
    idx = find(arc_lengths >= target_lengths(i), 1, 'first');
    if isempty(idx)
        u_spaced(i) = 1;
    else
        u_spaced(i) = u_eval(idx);
    end
end

% Evaluate curve at evenly spaced parameters
spaced_points = evaluateNURBS(control_points, weights, knots, degree, u_spaced);

% FORCE first and last points to be exactly on control points
spaced_points(1, :) = control_points(1, :);
spaced_points(end, :) = control_points(end, :);

%% 5. Plot in 3D
figure('Color', 'w');
hold on; grid on; axis equal;

% Define plot limits
plot_x = [-0.02, 0.02];
plot_y = [-0.02, 0.02];
plot_z = [-0.02, 0.02];

% Plot axis lines
plot3([plot_x(1), plot_x(2)], [0, 0], [0, 0], 'k-', 'LineWidth', 1.5, 'DisplayName', 'X-axis');
plot3([0, 0], [plot_y(1), plot_y(2)], [0, 0], 'k-', 'LineWidth', 1.5, 'DisplayName', 'Y-axis');
plot3([0, 0], [0, 0], [plot_z(1), plot_z(2)], 'k-', 'LineWidth', 1.5, 'DisplayName', 'Z-axis');

% Plot control points
plot3(control_points(:,1), control_points(:,2), control_points(:,3), ...
    'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Control Points');

% Plot control polygon
plot3(control_points(:,1), control_points(:,2), control_points(:,3), ...
    'r--', 'LineWidth', 1, 'DisplayName', 'Control Polygon');

% Plot smooth NURBS curve
plot3(curve_points(:,1), curve_points(:,2), curve_points(:,3), ...
    'b-', 'LineWidth', 2, 'DisplayName', 'NURBS Curve');

% Plot evenly spaced points
plot3(spaced_points(:,1), spaced_points(:,2), spaced_points(:,3), ...
    'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Evenly Spaced Points');

% Highlight first and last points with larger markers
plot3(spaced_points(1,1), spaced_points(1,2), spaced_points(1,3), ...
    'ko', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Start/End Points');
plot3(spaced_points(end,1), spaced_points(end,2), spaced_points(end,3), ...
    'ko', 'MarkerSize', 12, 'LineWidth', 2, 'HandleVisibility', 'off');

% Labels and legend
xlabel('X'); ylabel('Y'); zlabel('Z');
title('NURBS Curve with Evenly Spaced Points');
legend('Location', 'best');
view(3);
rotate3d on;

fprintf('Total curve length: %.4f\n', total_length);
fprintf('Spacing between points: %.4f\n', total_length / (n_spaced - 1));

% Verify endpoints match
endpoint_error_start = norm(spaced_points(1,:) - control_points(1,:));
endpoint_error_end = norm(spaced_points(end,:) - control_points(end,:));
fprintf('Endpoint error (start): %.6e\n', endpoint_error_start);
fprintf('Endpoint error (end): %.6e\n', endpoint_error_end);

%% Helper function: Evaluate NURBS curve
function points = evaluateNURBS(ctrl_pts, weights, knots, degree, u)
    n_eval = length(u);
    n_ctrl = size(ctrl_pts, 1);
    points = zeros(n_eval, 3);
    
    for i = 1:n_eval
        % Find the knot span
        span = findSpan(n_ctrl - 1, degree, u(i), knots);
        
        % Compute basis functions
        N = basisFunctions(span, u(i), degree, knots);
        
        % Compute point on curve
        numerator = [0, 0, 0];
        denominator = 0;
        
        for j = 0:degree
            idx = span - degree + j + 1;
            if idx >= 1 && idx <= n_ctrl
                numerator = numerator + N(j+1) * weights(idx) * ctrl_pts(idx, :);
                denominator = denominator + N(j+1) * weights(idx);
            end
        end
        
        if denominator > 0
            points(i, :) = numerator / denominator;
        end
    end
end

%% Helper function: Find knot span
function span = findSpan(n, p, u, U)
    if u >= U(n+2)
        span = n;
        return;
    end
    
    low = p;
    high = n + 1;
    mid = floor((low + high) / 2);
    
    while u < U(mid+1) || u >= U(mid+2)
        if u < U(mid+1)
            high = mid;
        else
            low = mid;
        end
        mid = floor((low + high) / 2);
    end
    
    span = mid;
end

%% Helper function: Compute basis functions
function N = basisFunctions(i, u, p, U)
    N = zeros(p+1, 1);
    left = zeros(p+1, 1);
    right = zeros(p+1, 1);
    
    N(1) = 1.0;
    
    for j = 1:p
        left(j+1) = u - U(i+2-j);
        right(j+1) = U(i+j+1) - u;
        saved = 0.0;
        
        for r = 0:j-1
            temp = N(r+1) / (right(r+2) + left(j-r+1));
            N(r+1) = saved + right(r+2) * temp;
            saved = left(j-r+1) * temp;
        end
        
        N(j+1) = saved;
    end
end