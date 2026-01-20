% Define pyramid parameters
tip_position = [0, 0, 0.0195];
base_center = [0, 0, 0];
base_radius = 0.0235;

% Hardcoded test point
test_point = [0, 0, 0.005];

% Define the 4 points to plot
% Point 1: Tip position
point1 = tip_position;

% Points 2-4: Three equally spaced points on the base circle (120 degrees apart)
% Rotate by 90 degrees (pi/2) so one point is at (0, 0.011, 0)
angles = [pi/2, pi/2 + 2*pi/3, pi/2 + 4*pi/3];
base_points = zeros(3, 3);
for i = 1:3
    angle = angles(i);
    x = base_radius * cos(angle);
    y = base_radius * sin(angle);
    z = 0;
    base_points(i, :) = [x, y, z];
end

point2 = base_points(1, :);
point3 = base_points(2, :);
point4 = base_points(3, :);

% Check if test point is in the pyramid
if is_point_in_pyramid(test_point, tip_position, base_points)
    fprintf('Test point [%.3f, %.3f, %.3f] is inside the pyramid\n', test_point);
else
    fprintf('Test point [%.3f, %.3f, %.3f] is outside the pyramid boundaries\n', test_point);
end

% Create 3D plot
figure('Position', [100, 100, 1000, 800]);

% Plot pyramid edges
% Base triangle edges
hold on;
% for i = 1:3
%     p1 = base_points(i, :);
%     p2 = base_points(mod(i, 3) + 1, :);
%     plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], ...
%           'b-', 'LineWidth', 2, 'Color', [0, 0, 1, 0.5], 'HandleVisibility', 'off');
% end
% 
% % Edges from base vertices to tip
% for i = 1:3
%     base_point = base_points(i, :);
%     plot3([base_point(1), tip_position(1)], ...
%           [base_point(2), tip_position(2)], ...
%           [base_point(3), tip_position(3)], ...
%           'b-', 'LineWidth', 2, 'Color', [0, 0, 1, 0.5], 'HandleVisibility', 'off');
% end

% Plot intermediate triangular cross-sections at different heights
% heights = [0.005, 0.010];
% for h = heights
%     % Linear interpolation parameter
%     t = h / 0.0195;  % Updated to use actual tip height
%     intermediate_points = zeros(3, 3);
%     
%     for i = 1:3
%         base_point = base_points(i, :);
%         % Interpolate between base point and tip
%         interp_point = base_point + t * (tip_position - base_point);
%         intermediate_points(i, :) = interp_point;
%     end
%     
%     % Draw the triangle at this height
%     for i = 1:3
%         p1 = intermediate_points(i, :);
%         p2 = intermediate_points(mod(i, 3) + 1, :);
%         plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], ...
%               'b--', 'LineWidth', 1, 'Color', [0, 0, 1, 0.3], 'HandleVisibility', 'off');
%     end
% end

% Plot the test point
% scatter3(test_point(1), test_point(2), test_point(3), 150, 'cyan', 'd', 'filled', 'DisplayName', 'Test Point');

% Draw a circle with 20mm (0.02m) radius on the xy plane
circle_radius = 0.02;
theta = linspace(0, 2*pi, 100);
circle_x = circle_radius * cos(theta);
circle_y = circle_radius * sin(theta);
circle_z = zeros(size(theta));
plot3(circle_x, circle_y, circle_z, 'r-', 'LineWidth', 2, 'Color', [1, 0, 0, 0.3], 'HandleVisibility', 'off');

% Create secondary pyramid with base points 20% closer to tip
% Each base point moves 20% of the distance from base to tip
secondary_base_points = zeros(3, 3);
for i = 1:3
    base_point = base_points(i, :);
    % Move 20% closer to tip
    secondary_base_points(i, :) = base_point + 0.2 * (tip_position - base_point);
end

% Plot secondary pyramid edges (without legend)
% Base triangle edges
for i = 1:3
    p1 = secondary_base_points(i, :);
    p2 = secondary_base_points(mod(i, 3) + 1, :);
    plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], ...
          'g-', 'LineWidth', 2, 'Color', [0, 1, 0, 0.5], 'HandleVisibility', 'off');
end

% Edges from secondary base vertices to tip
for i = 1:3
    sec_base_point = secondary_base_points(i, :);
    plot3([sec_base_point(1), tip_position(1)], ...
          [sec_base_point(2), tip_position(2)], ...
          [sec_base_point(3), tip_position(3)], ...
          'g-', 'LineWidth', 2, 'Color', [0, 1, 0, 0.5], 'HandleVisibility', 'off');
end

% Now plot the legend items
% Plot tip
scatter3(point1(1), point1(2), point1(3), 100, 'red', 'o', 'filled', 'DisplayName', 'Initial End Effector Location');
% Plot one base point for legend
scatter3(point2(1), point2(2), point2(3), 100, 'blue', 's', 'filled', 'DisplayName', 'Tendon Location');
% Plot other base points without legend
scatter3(point3(1), point3(2), point3(3), 100, 'blue', 's', 'filled', 'HandleVisibility', 'off');
scatter3(point4(1), point4(2), point4(3), 100, 'blue', 's', 'filled', 'HandleVisibility', 'off');
% Plot one workspace edge for legend
p1 = secondary_base_points(1, :);
p2 = secondary_base_points(2, :);
plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], ...
      'g-', 'LineWidth', 2, 'Color', [0, 1, 0, 0.5], 'DisplayName', 'Workspace');

% Draw lines from each base point to the test point and calculate distances
base_point_list = [point2; point3; point4];
distances = zeros(3, 1);

for i = 1:3
    base_point = base_point_list(i, :);
    
    % Draw line
%     line_x = [base_point(1), test_point(1)];
%     line_y = [base_point(2), test_point(2)];
%     line_z = [base_point(3), test_point(3)];
%     plot3(line_x, line_y, line_z, 'k--', 'LineWidth', 2, 'Color', [0, 0, 0, 0.6], 'HandleVisibility', 'off');
    
    % Calculate distance
    distance = norm(test_point - base_point);
    distances(i) = distance;
    fprintf('Distance from Base Point %d to Test Point: %.6f m\n', i, distance);
end

% Calculate rest lengths (distances from base points to tip)
rest_lengths = zeros(3, 1);
fprintf('\nRest lengths (distance to tip at rest):\n');
for i = 1:3
    base_point = base_point_list(i, :);
    rest_length = norm(tip_position - base_point);
    rest_lengths(i) = rest_length;
    fprintf('Base Point %d rest length: %.6f m\n', i, rest_length);
end

% Calculate percentages relative to rest lengths
fprintf('\nPercentages relative to rest length:\n');
for i = 1:3
    percentage = (distances(i) / rest_lengths(i)) * 100;
    fprintf('Base Point %d: %.2f%%\n', i, percentage);
end

% Labels and formatting
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Triangular Pyramid with Test Point and Distance Lines');
legend('Location', 'best');

% Set equal aspect ratio
max_range = 0.04;  % Updated to show the full 40mm circle
xlim([-max_range, max_range]);
ylim([-max_range, max_range]);
zlim([0, 0.025]);  % Updated to accommodate taller pyramid

% Set viewing angle
view(3);
grid on;
axis equal;
pbaspect([1 1 1]);  % Force equal aspect ratio for all axes
hold off;

% Function to check if a point is inside the triangular pyramid
function result = is_point_in_pyramid(point, tip, base_points)
    % Check if a point is inside a triangular pyramid (tetrahedron).
    % Uses barycentric coordinates method.
    
    % The four vertices of the tetrahedron
    v0 = tip;
    v1 = base_points(1, :);
    v2 = base_points(2, :);
    v3 = base_points(3, :);
    
    % Create matrix for barycentric coordinate calculation
    mat = [v0(1) - v3(1), v1(1) - v3(1), v2(1) - v3(1);
           v0(2) - v3(2), v1(2) - v3(2), v2(2) - v3(2);
           v0(3) - v3(3), v1(3) - v3(3), v2(3) - v3(3)];
    
    % Vector from v3 to the point
    vec = point - v3;
    
    % Check if matrix is singular
    if abs(det(mat)) < 1e-10
        result = false;
        return;
    end
    
    % Solve for barycentric coordinates
    bary = mat \ vec';
    
    % Check if point is inside: all coordinates non-negative and sum <= 1
    if all(bary >= -1e-10) && sum(bary) <= 1 + 1e-10
        result = true;
    else
        result = false;
    end
end