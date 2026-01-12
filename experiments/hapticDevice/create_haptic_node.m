% NURBS 3-Branch Cone Structure Generator
% Generates evenly spaced points along a NURBS curve and rotates to create
% 3 branches at 120-degree intervals

%% 1. Define control points for one branch
control_points = [
     0,      0,       .015;
    .007,    .005,       .015;
    .019,    .0040,       .007;
    .007,    .0011,       .005;
    .004,    .0007,       .003;
    .004,    .0005,   .002;
    .009,    .0005,   .002;
    .009,    0,     0;
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
n_spaced = 31;  % Number of points per branch (matching your original n_nodes_per_branch)

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

% Central node (first point of the branch)
central_node = branch_template(1, :);

hardcoded_vertices = [
     0.011 * cosd(-30),  0.011 * sind(-30),    0;   % Vertex 1
     0.011 * cosd(90),   0.011 * sind(90),     0;   % Vertex 2
     0.011 * cosd(210),  0.011 * sind(210),    0;   % Vertex 3
];
nodes = [hardcoded_vertices; central_node];  % Vertices 1-3, central node is node 4
edges = [];

% Create edges from hardcoded vertices to central node
for i = 1:3
    edges = [edges; i, 4];
end

for branch = 1:n_branches
    % Rotation angle in radians
    angle_rad = deg2rad(angles(branch));
    
    % Create rotation matrix around Z-axis
    R_z = [cos(angle_rad), -sin(angle_rad), 0;
           sin(angle_rad),  cos(angle_rad), 0;
           0,               0,              1];
    
    % Rotate all points in the branch template
    % Subtract central node, rotate, then add back
    branch_nodes = zeros(size(branch_template));
    for j = 1:n_spaced
        point_centered = branch_template(j, :)' - central_node';
        point_rotated = R_z * point_centered;
        branch_nodes(j, :) = (point_rotated + central_node')';
    end
    
    % Add nodes (skip first node as it's the shared central node)
    nodes = [nodes; branch_nodes(2:end, :)];
    
    % Create edges for this branch
    start_idx = 4 + (branch-1)*(n_spaced-1);
    for i = 0:n_spaced-2
        if i == 0
            edges = [edges; 4, start_idx + i + 1];
        else
            edges = [edges; start_idx + i, start_idx + i + 1];
        end
    end
end

%% 5. Write to file in the same format as cone_structure
script_path = fileparts(mfilename('fullpath'));
filename = fullfile(script_path, 'nurbs_cone_structure.txt');

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
fprintf('Spacing between points: %.4f\n', total_length / (n_spaced - 1));

%% 6. Optional: Visualize the structure
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

% Highlight hardcoded vertices
plot3(hardcoded_vertices(:,1), hardcoded_vertices(:,2), hardcoded_vertices(:,3), ...
    'mo', 'MarkerSize', 12, 'MarkerFaceColor', 'm');

% Highlight central node
plot3(central_node(1), central_node(2), central_node(3), ...
    'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');

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