function robotDescriptionHapticDevice_1tendon_parametric(applied_force)
% Parametric robot description that accepts applied force as input
% Input: applied_force - force to apply in Newtons

%% Simulation parameters
sim_params.static_sim = true;
sim_params.TwoDsim = false;
sim_params.use_midedge = false;
sim_params.use_lineSearch = false;
sim_params.showFrames = false;
sim_params.logStep = 1;
sim_params.log_data = true;
sim_params.bergou_DER = 0;
sim_params.FDM = 0;

% Time step
sim_params.dt = 1e-4;
% Maximum number of iterations in Newton Solver
sim_params.maximum_iter = 100;
% Total simulation time
if(sim_params.static_sim)
    sim_params.totalTime = sim_params.dt*5;
else
    sim_params.totalTime = 0.0015; % sec
end
% How often the plot should be saved?
sim_params.plotStep = 1;

%% Input text file
inputFileName = 'experiments/hapticDevice/nurbs_cone_structure_noTendons.txt';
% reading the input text file
[nodes, edges, face_nodes] = inputProcessorNew(inputFileName);

%% Input parameters
% geometry parameters
geom.shell_h = 0;
geom.rod_r0 = 0.0036/2; % for contact

% geom cross section of rod
ro = .0015;
geom.Axs = pi * ro ^ 2;
geom.Ixs1 = pi * ro^4 / 4;
geom.Ixs2 = pi * ro^4 / 4;
geom.Jxs = pi * ro^4 / 2;

% material parameters
material.density = 1200;
material.youngs_rod = 6e7;
material.youngs_shell = 0;
material.poisson_rod = 0.5;
material.poisson_shell = 0;

%% external force list
env.ext_force_list = ["gravity", "pointForce"];

% environment parameters
env.g = [0, 0, 0]';

% Apply the parametric force (upward in Z direction)
env.ptForce = [0; 0; applied_force];  % Force in Z direction
% Define which node to apply the force to
env.ptForce_node = 1;  % Apply to node 1 (origin node)

%% Tolerance on force function
sim_params.tol = 1e-4;
sim_params.ftol = 1e-4;
sim_params.dtol = 1e-2;

%% Boundary conditions
fixed_node_indices = [30, 31, 60, 61, 90, 91];
fixed_edge_indices = [];

%% logging
input_log_node = size(nodes,1);

%% Plot dimensions
sim_params.plot_x = [-0.02,0.02];
sim_params.plot_y = [-0.02,0.02];
sim_params.plot_z = [-0.01,0.03];

%% Assign to base workspace so main script can access them
assignin('base', 'sim_params', sim_params);
assignin('base', 'nodes', nodes);
assignin('base', 'edges', edges);
assignin('base', 'face_nodes', face_nodes);
assignin('base', 'geom', geom);
assignin('base', 'material', material);
assignin('base', 'env', env);
assignin('base', 'fixed_node_indices', fixed_node_indices);
assignin('base', 'fixed_edge_indices', fixed_edge_indices);
assignin('base', 'input_log_node', input_log_node);

end