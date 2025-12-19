% create 3-branch cone structure with semicircular curves
n_nodes_per_branch = 10;
n_branches = 3;

% Central node
central_node = [0, 0, 0.015];

% Define the three branch endpoints (120 degrees apart when viewed from top)
% First branch: diagonal down to roughly -.015, -.015, 0
branch1_end = [-0.015, -0.015, 0];

% Calculate branch length
branch_length = norm(branch1_end - central_node);

% For 120-degree spacing in top view:
% Branch 1: -45 degrees from x-axis (SW direction)
% Branch 2: -45 + 120 = 75 degrees from x-axis
% Branch 3: -45 + 240 = 195 degrees from x-axis

angles = [-45, 75, 195]; % degrees from x-axis
branch_ends = zeros(3, 3);

for i = 1:n_branches
    angle_rad = deg2rad(angles(i));
    % Project onto xy-plane
    xy_distance = branch_length * cos(atan2(0.015, branch_length * sqrt(2)));
    branch_ends(i, 1) = central_node(1) + xy_distance * cosd(angles(i));
    branch_ends(i, 2) = central_node(2) + xy_distance * sind(angles(i));
    branch_ends(i, 3) = 0;
end

% Create nodes for all branches
nodes = central_node; % Start with central node
edges = [];

for branch = 1:n_branches
    % Create semicircular path from central_node to branch_ends(branch,:)
    
    % Vector from central to end
    vec = branch_ends(branch, :) - central_node;
    
    % Create parameter t from 0 to pi (semicircle)
    t = linspace(0, pi, n_nodes_per_branch);
    
    % Create semicircle in a 2D plane, then rotate to proper orientation
    % Semicircle goes from (0,R) down through (R,0) to (0,-R) when concave down
    % But we want it to go from central_node to branch_end
    
    % Distance in xy-plane and total 3D distance
    xy_dist = norm(vec(1:2));
    z_dist = vec(3);
    total_dist = norm(vec);
    
    % Radius of semicircle (half the chord length for a semicircle)
    radius = total_dist / 2;
    
    % Midpoint between central and end
    midpoint = central_node + vec/2;
    
    % Create semicircle centered at midpoint
    % Parametric semicircle: going from start to end, curving downward
    branch_nodes = zeros(n_nodes_per_branch, 3);
    
    for j = 1:n_nodes_per_branch
        % Angle parameter from 0 to pi
        theta = t(j);
        
        % Local coordinates: semicircle from (-R, 0) to (R, 0) with max at (0, -R)
        local_x = radius * cos(theta - pi/2);  % Goes from -R to +R
        local_z = -radius * sin(theta - pi/2); % Goes from 0 down to -R and back to 0
        
        % Direction vector in xy-plane
        xy_dir = vec(1:2) / xy_dist;
        
        % Transform to global coordinates
        % Move along the xy direction
        branch_nodes(j, 1) = midpoint(1) + local_x * xy_dir(1);
        branch_nodes(j, 2) = midpoint(2) + local_x * xy_dir(2);
        % z-coordinate: linear drop + semicircular bulge downward
        branch_nodes(j, 3) = midpoint(3) + local_z;
    end
    
    % Skip first node (it's the shared central node)
    nodes = [nodes; branch_nodes(2:end, :)];
    
    % Create edges for this branch
    start_idx = 1 + (branch-1)*(n_nodes_per_branch-1);
    for i = 0:n_nodes_per_branch-2
        if i == 0
            edges = [edges; 1, start_idx + i + 1];
        else
            edges = [edges; start_idx + i, start_idx + i + 1];
        end
    end
end

% Write to file
filename = 'cone_structure_n8.txt';
fid = fopen(filename,'w');
if fid ~= -1
    fprintf(fid,'*Nodes\n');
    fclose(fid);
end
writematrix(nodes, filename, 'WriteMode', 'append')

fid = fopen(filename, 'at');
if fid ~= -1
    fprintf(fid,'*Edges\n');
    fclose(fid);
end
writematrix(edges, filename, 'WriteMode', 'append')

disp(['Created structure with ' num2str(size(nodes,1)) ' nodes and ' num2str(size(edges,1)) ' edges']);