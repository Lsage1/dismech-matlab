% NURBS 3-Branch Cone Structure Generator
% Generates evenly spaced points along a NURBS curve and rotates to create
% 3 branches at 120-degree intervals, starting from origin

%% 1. Define control points for one branch (starting from elevated position)
control_points = [
    0,      0,       .022;   % Start at elevated position
    
    .006, .00340,  .02166;
    .012, .00596,  .01978;
    .016, .00565,  .01551;
    .015, .00449,  .01235;
    
    .0088, .00339,  .00931;
    .00747,.00239,  .00656;
    .01,   .00171,   .0047;
    .0144, .00171,   .0047;
    .017,  .00082,  .00226;
    .017,   0,       0;      % End at base
];

%% 2. Create NURBS curve
n_ctrl = size(control_points, 1);
degree = min(3, n_ctrl - 1);

% Create CLAMPED knot vector
n_internal = n_ctrl - degree - 1;
if n_internal > 0
    internal_knots = linspace(0, 1, n_internal + 2);
    internal_knots = internal_knots(2:end-1);
else
    internal_knots = [];
end
knots = [zeros(1, degree+1), internal_knots, ones(1, degree+1)];

% Weights (all equal to 1 for standard B-spline)
weights = ones(n_ctrl, 1);

%% 3. Generate evenly spaced points along the curve
n_spaced = 31;  % Number of points per branch

% First evaluate at many points for arc length calculation
n_eval = 1000;
u_eval = linspace(0, 1, n_eval);
curve_points = evaluateNURBS(control_points, weights, knots, degree, u_eval);

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
    idx = find(arc_lengths >= target_lengths(i), 1, 'first');
    if isempty(idx)
        u_spaced(i) = 1;
    else
        u_spaced(i) = u_eval(idx);
    end
end

% Evaluate curve at evenly spaced parameters
branch_template = evaluateNURBS(control_points, weights, knots, degree, u_spaced);

% Force first and last points to match control points exactly
branch_template(1, :) = control_points(1, :);
branch_template(end, :) = control_points(end, :);

%% 4. Create 3 branches by rotating around Z-axis
n_branches = 3;
angles = [0, 120, 240];  % Degrees from x-axis

% Origin node (shared by all branches)
origin_node = [0, 0, 0.022];

% Initialize with origin node
nodes = origin_node;  % Node 1: origin
edges = [];

% Number of nodes per branch excluding the origin (which is shared)
n_nodes_per_branch = n_spaced - 1;

for branch = 1:n_branches
    % Rotation angle in radians
    angle_rad = deg2rad(angles(branch));
    
    % Create rotation matrix around Z-axis
    R_z = [cos(angle_rad), -sin(angle_rad), 0;
           sin(angle_rad),  cos(angle_rad), 0;
           0,               0,              1];
    
    % Rotate branch points (excluding the first point which is the origin)
    branch_nodes = zeros(n_nodes_per_branch, 3);
    for j = 1:n_nodes_per_branch
        point_rotated = R_z * branch_template(j+1, :)';
        branch_nodes(j, :) = point_rotated';
    end
    
    % Add nodes from this branch
    nodes = [nodes; branch_nodes];
    
    % Create edges for this branch
    start_idx = 1 + (branch-1)*n_nodes_per_branch;
    for i = 0:n_nodes_per_branch-1
        if i == 0
            % Connect origin to first node of branch
            edges = [edges; 1, start_idx + i + 1];
        else
            % Connect consecutive nodes within branch
            edges = [edges; start_idx + i, start_idx + i + 1];
        end
    end
end

%% 5. Write to file
script_path = fileparts(mfilename('fullpath'));
filename = fullfile(script_path, 'nurbs_cone_structure_1tendon.txt');

% Write nodes
fid = fopen(filename,'w');
if fid ~= -1
    fprintf(fid,'*Nodes\n');
    fclose(fid);
end
writematrix(nodes, filename, 'WriteMode', 'append');

% Write edges
fid = fopen(filename, 'at');
if fid ~= -1
    fprintf(fid,'*Edges\n');
    fclose(fid);
end
writematrix(edges, filename, 'WriteMode', 'append');

disp(['Created NURBS structure with ' num2str(size(nodes,1)) ' nodes and ' num2str(size(edges,1)) ' edges']);
disp(['Output file: ' filename]);
fprintf('Total curve length per branch: %.4f\n', total_length);
fprintf('Average spacing between points: %.4f\n', total_length / (n_spaced - 1));

%% 6. Visualize the structure
figure('Color', 'w');
hold on; grid on; axis equal;

% Plot all nodes
plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

% Plot edges
for i = 1:size(edges, 1)
    n1 = edges(i, 1);
    n2 = edges(i, 2);
    plot3([nodes(n1,1), nodes(n2,1)], ...
          [nodes(n1,2), nodes(n2,2)], ...
          [nodes(n1,3), nodes(n2,3)], 'r-', 'LineWidth', 2);
end

% Highlight origin node
plot3(origin_node(1), origin_node(2), origin_node(3), ...
    'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

xlabel('X'); ylabel('Y'); zlabel('Z');
title('NURBS-Based 3-Branch Cone Structure');
view(3);
rotate3d on;

%% Helper Functions

function points = evaluateNURBS(ctrl_pts, weights, knots, degree, u)
    n_eval = length(u);
    n_ctrl = size(ctrl_pts, 1);
    points = zeros(n_eval, 3);
    
    for i = 1:n_eval
        span = findSpan(n_ctrl - 1, degree, u(i), knots);
        N = basisFunctions(span, u(i), degree, knots);
        
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