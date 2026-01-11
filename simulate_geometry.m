function [tilt_x, tilt_y, convergence_failed] = simulate_geometry(edges_to_actuate, actuation_factors, top_verts_ind, ...
    material, geom, env, sim_params, nodes, edges, face_nodes, fixed_node_indices, fixed_edge_indices)
% SIMULATE_GEOMETRY Run physics simulation and return plane tilt in two axes
%
% Inputs:
%   edges_to_actuate - Vector of edge indices to actuate (e.g., [1, 2, 3])
%   actuation_factors - Factors to multiply initial rest length (e.g., [1, 0.5, 1])
%   top_verts_ind - Indices of top vertices to track (e.g., [5, 35, 65])
%   material - Material properties struct from robotDescriptionHapticDevice
%   geom - Geometry struct from robotDescriptionHapticDevice
%   env - Environment struct from robotDescriptionHapticDevice
%   sim_params - Simulation parameters struct from robotDescriptionHapticDevice
%   nodes - Node array from robotDescriptionHapticDevice
%   edges - Edge array from robotDescriptionHapticDevice
%   face_nodes - Face nodes from robotDescriptionHapticDevice
%   fixed_node_indices - Fixed node indices from robotDescriptionHapticDevice
%   fixed_edge_indices - Fixed edge indices from robotDescriptionHapticDevice
%
% Outputs:
%   tilt_x - Tilt angle about X-axis (pitch) in degrees
%   tilt_y - Tilt angle about Y-axis (roll) in degrees
%   convergence_failed - Flag: 1 if convergence error occurred, 0 otherwise
%
% Example:
%   robotDescriptionHapticDevice;  % Run this first to set up variables
%   [tx, ty] = simulate_geometry([1,2,3], [1, 0.5, 1], [5, 35, 65], ...
%       material, geom, env, sim_params, nodes, edges, face_nodes, ...
%       fixed_node_indices, fixed_edge_indices);

% Initialize convergence flag
convergence_failed = 0;

% add to path
addpath springs/
addpath util_functions/
addpath contact_functions/
addpath rod_dynamics/
addpath shell_dynamics/
addpath external_forces/
addpath adaptive_stepping/
addpath logging/
addpath(genpath('experiments')); 
addpath experiments/hapticDevice

n_actuated_edges = length(edges_to_actuate);

% create geometry
[nodes, edges, rod_edges, shell_edges, rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, face_edges, face_shell_edges, ...
    elStretchRod, elStretchShell, elBendRod, elBendSign, elBendShell, sign_faces, face_unit_norms]...
    = createGeometry(nodes, edges, face_nodes);

% intialize twist angles for rod-edges to 0
twist_angles=zeros(size(rod_edges,1)+size(rod_shell_joint_total_edges,1),1);

% create environment and imc structs
[environment,imc] = createEnvironmentAndIMCStructs(env,geom,material,sim_params);

%% Create the soft robot structure
softRobot = MultiRod(geom, material, twist_angles,...
    nodes, edges, rod_edges, shell_edges, rod_shell_joint_edges, rod_shell_joint_total_edges, ...
    face_nodes, sign_faces, face_edges, face_shell_edges, sim_params, environment);

%% Creating stretching, bending, twisting and hinge springs

n_stretch = size(elStretchRod,1) + size(elStretchShell,1);
n_bend_twist = size(elBendRod,1);

% stretching spring
if(n_stretch==0)
    stretch_springs = [];
else
    for s=1:n_stretch
        if (s <= size(elStretchRod,1)) % rod
            stretch_springs (s) = stretchSpring (...
                softRobot.refLen(s), elStretchRod(s,:),softRobot);
        else % shell
            stretch_springs (s) = stretchSpring (...
                softRobot.refLen(s), ...
                elStretchShell(s-(size(elStretchRod,1)),:), ...
                softRobot, softRobot.ks(s));
        end
    end
end

% bending and twisting spring
if(n_bend_twist==0)
    bend_twist_springs = [];
else
    for b=1:n_bend_twist
        bend_twist_springs(b) = bendTwistSpring ( ...
            elBendRod(b,:), elBendSign(b,:), [0 0], 0, softRobot);
    end
end

%% Store initial rest length for actuation
initial_rest_length = zeros(n_actuated_edges, 1);
for i = 1:n_actuated_edges
    edge_idx = edges_to_actuate(i);
    initial_rest_length(i) = stretch_springs(edge_idx).refLen;
end

% Calculate final rest lengths using provided factors
final_rest_length = initial_rest_length .* actuation_factors(:);

% shell bending spring
n_hinge = size(elBendShell,1);
n_triangle = softRobot.n_faces;
if(n_triangle==0)
    hinge_springs = [];
    triangle_springs = [];
else
    if(~sim_params.use_midedge)
        triangle_springs = [];
        for h=1:n_hinge
            hinge_springs(h) = hingeSpring (...
                0, elBendShell(h,:), softRobot, softRobot.kb);
        end
        hinge_springs = setThetaBar(hinge_springs, softRobot);
    else
        hinge_springs = [];
        for t=1:n_triangle
            triangle_springs(t) = triangleSpring(softRobot.face_nodes_shell(t,:), softRobot.face_edges(t,:), softRobot.face_shell_edges(t,:), softRobot.sign_faces(t,:), softRobot);
        end
    end
end

%% Prepare system
% Reference frame (Space parallel transport at t=0)
softRobot = computeSpaceParallel(softRobot);

% Material frame from reference frame and twist angle
theta = softRobot.q0(3*softRobot.n_nodes+1:3*softRobot.n_nodes+softRobot.n_edges_dof);
[softRobot.m1, softRobot.m2] = computeMaterialDirectors(softRobot.a1,softRobot.a2,theta);

% Set rod natural curvature
bend_twist_springs = setkappa(softRobot, bend_twist_springs);

% Reference twist
softRobot.undef_refTwist = computeRefTwist_bend_twist_spring ...
    (bend_twist_springs, softRobot.a1, softRobot.tangent, ...
    zeros(n_bend_twist,1));
softRobot.refTwist = computeRefTwist_bend_twist_spring ...
    (bend_twist_springs, softRobot.a1, softRobot.tangent, ...
    softRobot.undef_refTwist);

%% Boundary Conditions
softRobot.fixed_nodes = fixed_node_indices;
for i=1:size(softRobot.Edges,1)
    if ( ismember(softRobot.Edges(i,1),fixed_node_indices) && ismember(softRobot.Edges(i,2),fixed_node_indices) )
        fixed_edge_indices = [fixed_edge_indices, i];
    end
end
if(sim_params.TwoDsim)
    fixed_edge_indices = [fixed_edge_indices, 1:softRobot.n_edges_dof];
end
softRobot.fixed_edges = fixed_edge_indices;
[softRobot.fixedDOF, softRobot.freeDOF] = FindFixedFreeDOF(softRobot.fixed_nodes, softRobot.fixed_edges, softRobot.n_DOF, softRobot.n_nodes);

%% Time stepping scheme
Nsteps = round(sim_params.totalTime/sim_params.dt);
ctime = 0;

for timeStep = 1:Nsteps
    if(sim_params.static_sim)
        environment.g = timeStep*environment.static_g/Nsteps;
        
        for i = 1:n_actuated_edges
            edge_idx = edges_to_actuate(i);
            current_rest_length = initial_rest_length(i) + ...
                timeStep * (final_rest_length(i) - initial_rest_length(i)) / Nsteps;
            
            stretch_springs(edge_idx).refLen = current_rest_length;
            softRobot.refLen(edge_idx) = current_rest_length;
        end
    end
    
    if(sim_params.use_midedge)
        tau_0 = updatePreComp_without_sign(softRobot.q, softRobot);
    else
        tau_0 = [];
    end
    
    [softRobot, stretch_springs, bend_twist_springs, hinge_springs, step_convergence_failed] = ...
        timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, triangle_springs, tau_0, environment, imc, sim_params);
    
    if step_convergence_failed
        convergence_failed = 1;
        break; % Exit loop
    end

    ctime = ctime + sim_params.dt;
    softRobot.q0 = softRobot.q;
end

%% Extract final positions of top nodes and compute tilt
top_node_positions = zeros(length(top_verts_ind), 3);
for i = 1:length(top_verts_ind)
    node_idx = top_verts_ind(i);
    top_node_positions(i, :) = [softRobot.q0(3*node_idx-2), ...
                                 softRobot.q0(3*node_idx-1), ...
                                 softRobot.q0(3*node_idx)];
end

%% Compute tilt of the plane formed by top nodes
[tilt_x, tilt_y] = compute_plane_tilt(top_node_positions);

%% Plot final geometry
plot_MultiRod(softRobot, ctime, sim_params, environment, imc, top_verts_ind);
if convergence_failed
    title(sprintf('Final Geometry (CONVERGENCE FAILED) - Tilt: X=%.2f°, Y=%.2f°', tilt_x, tilt_y));
else
    title(sprintf('Final Geometry - Tilt: X=%.2f°, Y=%.2f°', tilt_x, tilt_y));
end

end

function [tilt_x, tilt_y] = compute_plane_tilt(points)
% COMPUTE_PLANE_TILT Calculate tilt angles of plane from horizontal
%
% Input: points - Nx3 matrix of [x, y, z] coordinates
% Output: 
%   tilt_x - tilt about X-axis (pitch) in degrees
%   tilt_y - tilt about Y-axis (roll) in degrees

if size(points, 1) < 3
    error('Need at least 3 points to define a plane');
end

% Fit plane using first 3 points
p1 = points(1, :);
p2 = points(2, :);
p3 = points(3, :);

% Calculate normal vector to the plane
v1 = p2 - p1;
v2 = p3 - p1;
normal = cross(v1, v2);
normal = normal / norm(normal);  % normalize

% Ensure normal points upward (positive z-component)
if normal(3) < 0
    normal = -normal;
end

% Extract normal components
nx = normal(1);
ny = normal(2);
nz = normal(3);

% Tilt about X-axis (pitch): rotation in Y-Z plane
% When plane tilts forward/back, ny changes
tilt_x = atand(ny / nz);

% Tilt about Y-axis (roll): rotation in X-Z plane  
% When plane tilts left/right, nx changes
tilt_y = atand(-nx / nz);

end