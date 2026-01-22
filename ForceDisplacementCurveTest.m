clc
clear all
close all
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

%% Force sweep parameters
force_values = linspace(-2, -6, 4);  % Forces from 1N to 10N in 10 steps
n_forces = length(force_values);

% Arrays to store results
displacement_results = zeros(n_forces, 1);
force_results = zeros(n_forces, 1);

%% Run simulation for each force value
for force_idx = 1:n_forces
    fprintf('\n========================================\n');
    fprintf('Running simulation %d/%d with Force = %.2f N\n', force_idx, n_forces, force_values(force_idx));
    fprintf('========================================\n');
    
    % Clear previous simulation data
    clear nodes edges face_nodes softRobot stretch_springs bend_twist_springs hinge_springs
    
    % Set the force value for this iteration
    applied_force = force_values(force_idx);
    
    % Call robot description with the current force value
    robotDescriptionHapticDevice_1tendon(applied_force);
    
    % Define nodes to highlight in green
    top_verts_ind = [1];
    
    edges_to_actuate = [];  % No tendon actuation
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
    
    % Store initial position
    initial_pos_z = softRobot.q0(3*1);
    
    %% Time stepping scheme
    Nsteps = round(sim_params.totalTime/sim_params.dt);
    ctime = 0;
    
    for timeStep = 1:Nsteps
        if(sim_params.static_sim)
            environment.g = timeStep*environment.static_g/Nsteps; % ramp gravity
        end
        
        %% Precomputation at each timeStep
        if(sim_params.use_midedge)
            tau_0 = updatePreComp_without_sign(softRobot.q, softRobot);
        else
            tau_0 = [];
        end
        
        %%  Implicit stepping error iteration
        [softRobot, stretch_springs, bend_twist_springs, hinge_springs] = ...
            timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, triangle_springs, tau_0,environment,imc, sim_params);
        
        ctime = ctime + sim_params.dt;
        
        % Update q
        softRobot.q0 = softRobot.q;
    end
    
    %% Extract results for this force value
    final_pos_z = softRobot.q0(3*1);
    displacement = final_pos_z - initial_pos_z;
    
    force_results(force_idx) = applied_force;
    displacement_results(force_idx) = displacement;
    
    fprintf('Applied Force: %.2f N\n', applied_force);
    fprintf('Initial Z position: %.6f m\n', initial_pos_z);
    fprintf('Final Z position: %.6f m\n', final_pos_z);
    fprintf('Displacement: %.6f m\n', displacement);
    
    %% Visualize final configuration
    figure(1);
    clf;
    plot_MultiRod(softRobot, 1.0, sim_params, environment, imc, top_verts_ind);
    title(sprintf('Final Configuration - Force = %.2f N, Displacement = %.4f mm', ...
        applied_force, displacement*1000), 'FontSize', 12, 'FontWeight', 'bold');
    
    % Pause to show the plot (optional - comment out for faster execution)
    pause(0.5);
    drawnow;
end

%% Plot Force-Displacement Curve
figure(2);
clf;
plot(displacement_results * 1000, force_results, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('Displacement (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Applied Force (N)', 'FontSize', 12, 'FontWeight', 'bold');
title('Force-Displacement Curve', 'FontSize', 14, 'FontWeight', 'bold');

% Add data labels
for i = 1:n_forces
    text(displacement_results(i)*1000, force_results(i), ...
        sprintf('  %.1fN', force_results(i)), ...
        'FontSize', 8, 'VerticalAlignment', 'bottom');
end

%% Display results table
fprintf('\n========================================\n');
fprintf('FORCE-DISPLACEMENT RESULTS\n');
fprintf('========================================\n');
fprintf('Force (N)\tDisplacement (mm)\n');
fprintf('----------------------------------------\n');
for i = 1:n_forces
    fprintf('%.2f\t\t%.4f\n', force_results(i), displacement_results(i)*1000);
end

%% Save results
results.force = force_results;
results.displacement = displacement_results;
save('force_displacement_results.mat', 'results');
fprintf('\nResults saved to force_displacement_results.mat\n');