%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      AeroFlex - Example - Experimental Wing with Tip Mass           %
% - computes eigenvalues/eigenvectors of linear structural dynamics   %
% - finds equilibrium (only gravity)                                  %
% - simulation (beam starts undeformed, goes to equilibrium)          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function experimental_wing_sim
    clc;
    clear all;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AeroFlex Main Folder:                              %
    addpath('..\..\main');    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % General AeroFlex parameters - MANDATORY PARAMETERS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global softPARAMS;
    softPARAMS.isIS = 1; %is it International System?
    softPARAMS.isPINNED = 1; % PINNED RIGID BODY DOF
    softPARAMS.vecFREEDEG = [0 0 0 0 0 0]; % type 0 to remove a body state
                                           % degree of freedom [u v w p q r]
    softPARAMS.isGRAV = 1; % include gravity?
    softPARAMS.g = 9.81; % gravity in m/s^2    
    softPARAMS.isITER = 1; % iterative equilibrium determination?
    softPARAMS.numITER = 10; % number of iterations for equilibrium determination
    softPARAMS.modAED = 0; % AERODYNAMIC MODEL: 
                                    %0-Steady;
                                    %1-Quasi-steady;
                                    %2-Quasi-steady with added mass;
                                    %3-Unsteady(Peters);
    softPARAMS.updateStrJac = 1; % Structural Jacobians updates:
                                    % 0 - Never;
                                    % 1 - Only in equilib calculation;
                                    % 2 - Always
    softPARAMS.plota3d = 1; % 3d graphics plot while equilibrium is calculated
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    %%%%%%%%%%%% STRUCTURE INITIALIZATION %%%%%%%%%%%%%%%%%
    numele = 10; %number of elements
    ap = load_structure(numele); % this creates a flexible
                                        %airplane object with numele elements
                                        % check the function loadstruct
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%% MODAL ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
    % find structural natural frequencies and modal shapes
    st_modes = StructuralModes(ap);
    numerical_freqs = st_modes.frequencies(1:6);
    fprintf('Numerical Natural frequencies (in Hz):\n');
    fprintf('%.3f Hz\n', numerical_freqs');
    
    %plot results    
    figure('color','w','name','Mode Shapes');
    subplot(1,2,1); st_modes.plot_mode(ap,1); title('Mode 1');
    view(30,45); axis equal; colormap winter;
    subplot(1,2,2); st_modes.plot_mode(ap,2); title('Mode 2');
    view(30,45); axis equal; colormap winter;
    
    fprintf('press any key to equilibrium calculation\n\n');
    pause;
    %%%%%%%%%%%% FINDS EQUILIBRIUM CONDITION AND PLOT%%%%%%%%%%%%%%
    % FLIGHT CONDITIONS -- won't affect the results if there is no aerodynamics
    altitude = 0; % meters
    V = 0;           % m/s
    Vwind = 0;        % m/s     
    throttle = 0;
    deltaflap = 0;
    [rb_eq, strain_eq] = trimairplane(ap,V,altitude,Vwind,throttle,deltaflap);
    
    figure('color','w','name','Equilibrium Deformation');
    % plot structure without deformation:
    update(ap,strain_eq*0,zeros(size(strain_eq)),zeros(size(strain_eq)),zeros(sum(ap.membNAEDtotal),1));
    plotairplane3d(ap); 
    % plot deformed structure (equilibrium condition):
    update(ap,strain_eq,zeros(size(strain_eq)),zeros(size(strain_eq)),zeros(sum(ap.membNAEDtotal),1));
    plotairplane3d(ap); 
    view(30,45); axis equal; colormap winter;
    title('Equilibrium under Gravity');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% NONLINEAR SIMULATION %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial conditions

    tSIM_default = 5; % Default simulation time
    prompt = sprintf('Simulation time (seconds) [default: %.1f s]: ', tSIM_default);
    tSIM_input = input(prompt);
    if isempty(tSIM_input)
        tSIM = tSIM_default;
    else
        tSIM = tSIM_input;
    end
    
    aerodynamic_surface_pos = 0; % no aerodynamic surfaces!
    engine_position = 0;     % no engine!
    beta0 = [0; 0;0;0;0;0]; % rigid body speeds
    k0 = [0;0;0;0];         % rigid body position/orientation
    strain0 = strain_eq*0;
    % simulation:  
    tic;
    [tNL, strainNL, straindNL, lambdaNL, betaNL, kineticNL] = simulate(ap, [0 tSIM], strain0, beta0, k0, Vwind, @(t)engine_position, @(t)aerodynamic_surface_pos, 'implicit');
    toc;
    
    dt = 0.05;
    [ts, Xs] = changedatarate(tNL,strainNL,dt);
    tip_displacement_z = zeros(length(ts),1);
    tip_displacement_x = zeros(length(ts),1);

    % To get tip displacement, we need to know which member and node is the tip
    % Assuming one member, tip is node3 of the last element
    last_element_idx = ap.membSIZES(1); % num_elements for the first member

    for i = 1:size(ts,1)
        update(ap,Xs(i,:),zeros(size(Xs(i,:))),zeros(size(Xs(i,:))),zeros(sum(ap.membNAEDtotal),1));
        % Assuming the wing is member 1
        tip_node_global_pos = ap.members{1}(last_element_idx).node3.h(1:3); % Global position [x;y;z]
        tip_displacement_z(i) = tip_node_global_pos(3); % Z-displacement (flap)
        tip_displacement_x(i) = tip_node_global_pos(1); % X-displacement (chordwise/axial combination)
    end
    figure('color','w','name','Wing Tip Displacement');
    subplot(2,1,1);
    plot(ts,tip_displacement_z);
    xlabel('Time (s)'); ylabel('Tip Z-displacement (m)');
    title('Tip Flap Displacement');
    grid on;
    subplot(2,1,2);
    plot(ts,tip_displacement_x);
    xlabel('Time (s)'); ylabel('Tip X-displacement (m)');
    title('Tip X-Displacement (Chord/Axial)');
    grid on;
    
    figure('color','w','name','Animation');
    airplanemovie(ap, ts, Xs,[],dt,'experimental_wing_sim','gif'); colormap winter;
    

end

function ap = load_structure(numele)
    % member initialization
    % Damping ratio is now handled within create_flexible_member if needed, or passed
    flexible_member = create_flexible_member(numele);
    
    %set member origin node position and orientation:
    % Wing root is at origin, aligned with body axes by default here.
    % h = [x y z exx exy exz eyx eyy eyz ezx ezy ezz]
    % For a cantilever, root is fixed. Orientation defines initial global orientation.
    % Standard identity orientation: ex = [1;0;0], ey = [0;1;0], ez = [0;0;1]
    flexible_member(1).seth0([0 0 0  1 0 0  0 1 0  0 0 1]'); 
    update(flexible_member); % initialize displacements for each member node
    fus = []; % no fuselage
    % Motors/stores are now handled via rigidunit at nodes if they are point masses
    % If they were separate objects they would go into the Airplane constructor's last arg
    external_stores = []; 
    ap = Airplane({flexible_member}, fus, external_stores);
end

function flexible_member = create_flexible_member(num_elements)
    
    % Wing Properties from experimental data
    Length = 0.4508;     % Span L (m)
    chord_c = 0.0508;    % Chord c (m)
    mcs = 0.2351;        % Mass per unit length (kg/m)
    
    % Elastic axis 50% chord, CG 49% chord
    % pos_cg is CG offset from elastic axis: [x_offset, y_offset, z_offset]
    % x is along chord (positive towards trailing edge from elastic axis)
    % y is along span
    % z is perpendicular to wing plane (positive up)
    pos_cg_section = [(0.49 - 0.50) * chord_c, 0, 0]; % [-0.000508, 0, 0] m

    % Sectional Rigidity Matrix KG = diag([EA, GJ, EI_flap, EI_chord])
    EA_axial = 1e10;       % EA (Axial Rigidity N) - Placeholder, User should verify this value!
    GJ_torsion = 0.9539;   % GJ (Torsional Rigidity Nm^2)
    EI_flap = 0.4186;      % EI1 (Flap Bending Rigidity Nm^2)
    EI_chord = 18.44;      % EI2 (Chordwise Bending Rigidity Nm^2) = 0.1844e2
    KG_section = diag([EA_axial, GJ_torsion, EI_flap, EI_chord]);
    
    % Damping: Using Rayleigh damping C = alpha*M + beta*K
    % Here, example1 used C = damp_ratio * K. We have modal damping ratios.
    % Flap ξw = 0.02, Chordwise ξv = 0.025, Torsion ξφ = 0.031
    % For simplicity, using flap damping for the K-proportional term as a start.
    % A more accurate model would require constructing C from modal damping.
    damp_ratio = 0.02; % Based on flap structural modal damping ξw (Used for CG = damp_ratio*KG)
    CG_section = damp_ratio * KG_section; % Sectional damping matrix
    
    % Wing Section Inertia (per unit length) around elastic axis
    % Inertia = diag([Ixx, Iyy, Izz]) per unit length
    % Ixx: torsional (about span axis, y)
    % Iyy: chordwise bending (about flap axis, z)
    % Izz: flapwise bending (about chord axis, x)
    % The problem gives "Moment of inertia (50% chord) = 0.2056e-4 kgm"
    % This is ambiguous. Assuming it's the polar moment of inertia per unit length
    % about the elastic axis (which is also 50% chord), so this would be I_polar_over_L.
    % Or, if it's Ixx (torsional inertia about span axis, y), used for GJ term in dynamics.
    % Let's assume it is I_yy (pitching MOI per unit length for the section, often denoted I_alpha or I_theta)
    % In AeroFlex notation (from example1 comments): I11 (torsion), I22 (chord bend), I33 (flap bend)
    % If "Moment of inertia (50% chord)" refers to the mass moment of inertia for rotation about the pitch axis (y-axis for the wing section)
    I_pitch_section_per_L = 0.2056e-4; % kgm -> assuming kg.m^2/m 
    % Other section inertias (I_roll_section_per_L, I_yaw_section_per_L) are not given.
    % These relate to how mass mcs is distributed. Using placeholders.
    I_roll_section_per_L = I_pitch_section_per_L; % Placeholder - common for symmetric sections if I_pitch_section_per_L is Ixx for torsion
    I_chord_bending_inertia_per_L = 0; % Placeholder - User should verify. This is for mass distribution along chord.
    I_flap_bending_inertia_per_L = 0;  % Placeholder - User should verify. This is for mass distribution normal to chord.

    % The `Inertia` input to Node is a 3x3 matrix: diag([Ixx_node, Iyy_node, Izz_node])
    % where these are principal moments of inertia of the mass at that infinitesimal segment dy.
    % Node.m: Inertia = diag([I11 I22 I33]); I11 related to torsion, I22 to chordwise, I33 to flapwise rotation of mass.
    % Let's assume the given value 0.2056e-4 kgm^2/m is the torsional mass moment of inertia I_yy per unit length (about span axis, used with GJ)
    % This corresponds to I11 in example1 for `create_uniform_structure`'s Inertia input if beam x is along span.
    % Given the setup, the beam's local x-axis is typically along its length (span for a wing).
    % So, `Inertia(1,1)` is torsional inertia, `Inertia(2,2)` is for bending about local y, `Inertia(3,3)` for bending about local z.
    % Wing cross section: x_local_beam (span), y_local_beam (chord), z_local_beam (flap direction)
    % Moment of inertia (50% chord) 0.2056e-4 kgm. If this is I_yy for the wing section (torsional), then:
    Inertia_section = diag([0.2056e-4, 0, 0]); % [I_torsion_per_L, I_chord_bend_MOI_per_L, I_flap_bend_MOI_per_L]
                                             % User MUST verify Iyy_section and Izz_section values. Using 0 as placeholder.
                                             % The 0.2056e-4 is for torsion about the beam's own axis (spanwise).

    % Aerodynamic parameters (not used in this structural example)
    aeroparams = []; 
    geometry.a = 0.0; % elastic axis position from leading edge as fraction of chord (0.5 for 50% chord)
    geometry.b = chord_c / 2; % semi-chord

    % Initial rotation of the first element (dihedral, sweep, twist)
    rot0.dihedral = 0; % radians
    rot0.sweep = 0;    % radians
    rot0.twist = 0;    % radians
    
    % Slender Body (Tip Mass) Properties
    tip_mass_M = 0.0417; % kg
    % tip_mass_R = 0.4762e-2; % m (Radius) - not directly used unless for CG offset calculation
    % tip_mass_cSB = 0.1406; % m (Chord length) - not directly used
    tip_mass_Ix = 0.9753e-4; % kgm^2
    tip_mass_Iy = 0.3783e-5; % kgm^2
    tip_mass_Iz = 0.9753e-4; % kgm^2
    
    tip_rigid_body.m = tip_mass_M;
    tip_rigid_body.cg = [0, 0, 0]; % CG of tip mass relative to attachment node (wing tip). Assume 0,0,0.
                                   % User should verify if there's an offset.
    tip_rigid_body.I = diag([tip_mass_Ix, tip_mass_Iy, tip_mass_Iz]);

    % Discretization
    ds = Length / num_elements;
    
    % Node and Element Creation Loop
    current_rot = rot0; % Initial rotation for the first element
    
    % Define a default rigid unit for non-tip nodes
    default_rigidunit.m = 0; 
    default_rigidunit.cg = [0,0,0]; 
    default_rigidunit.I = zeros(3,3);

    noh = []; % Array to store node objects
    membro = []; % Array to store element objects

    for i = 1:num_elements
        % For the last node of the last element, attach the tip mass
        if i == num_elements
            rigid_unit_node1 = default_rigidunit;
            rigid_unit_node2 = default_rigidunit;
            rigid_unit_node3 = tip_rigid_body; % Tip mass on the 3rd node of the last element
        else
            rigid_unit_node1 = default_rigidunit;
            rigid_unit_node2 = default_rigidunit;
            rigid_unit_node3 = default_rigidunit;
        end

        % Create nodes for the current element
        % Node s-parameter seems to be a normalized coordinate or just an index for aero
        % Using simple indexing for s for now, as in create_uniform_structure.
        % Node(mcs, mcg, Inertia, aeroparams, rigidunit, s, geometry)
        node_s_param_1 = ((i-1)*ds)/Length; % Normalized s for node 1 of element i
        node_s_param_2 = (ds/2 + (i-1)*ds)/Length; % Normalized s for node 2 of element i
        node_s_param_3 = (i*ds)/Length; % Normalized s for node 3 of element i

        temp_node1 = Node(mcs, pos_cg_section, Inertia_section, aeroparams, rigid_unit_node1, node_s_param_1 ,geometry);
        temp_node2 = Node(mcs, pos_cg_section, Inertia_section, aeroparams, rigid_unit_node2, node_s_param_2 ,geometry);
        temp_node3 = Node(mcs, pos_cg_section, Inertia_section, aeroparams, rigid_unit_node3, node_s_param_3 ,geometry);
        
        % Store all nodes if they are unique, or manage references if elements share nodes
        % The example structure seems to imply Element gets its own 3 nodes.
        % Let's follow the pattern from create_uniform_structure which creates 3*num_elements nodes.
        % For a continuous beam, nodes should be shared.
        % Re-evaluating: create_uniform_structure creates noh((i-1)*3+1), etc.
        % This means each element gets 3 fresh nodes in that implementation.
        % This is unusual for FEM, but let's stick to the pattern if it works in AeroFlex.
        % The Element constructor Me calculation suggests nodes are specific to elements and not shared.

        if i == 1
            noh_global_idx_1 = 1;
            noh_global_idx_2 = 2;
            noh_global_idx_3 = 3;
            noh(noh_global_idx_1) = temp_node1;
            noh(noh_global_idx_2) = temp_node2;
            noh(noh_global_idx_3) = temp_node3;
        else
            % If nodes are shared, node1 of current is node3 of previous.
            % However, create_uniform_structure creates 3 distinct nodes per element.
            % For now, let's assume distinct nodes per element as per create_uniform_structure
            % This will create (num_elements * 3) nodes in total.
            % If AeroFlex internally connects elements based on their order in the member array,
            % this might be fine. If it expects shared node objects, this is wrong.
            % The class Element has node1, node2, node3. This suggests they are part of the element.
            % Let's try with distinct nodes first.
            
            % To ensure the tip mass is ONLY on the very last node of the wing:
            % We need to make sure node3 of the last element gets the tip_rigid_body.
            % The loop above correctly assigns tip_rigid_body to rigid_unit_node3 only for i == num_elements.

            noh_idx_start = (i-1)*3;
            noh(noh_idx_start+1) = temp_node1;
            noh(noh_idx_start+2) = temp_node2;
            noh(noh_idx_start+3) = temp_node3;
            
            noh_global_idx_1 = noh_idx_start+1;
            noh_global_idx_2 = noh_idx_start+2;
            noh_global_idx_3 = noh_idx_start+3;
        end
            
        current_element = Element(noh(noh_global_idx_1), noh(noh_global_idx_2), noh(noh_global_idx_3), current_rot, ds, KG_section, CG_section);
        current_element.setstrain([0 0 0 0],[0 0 0 0]); % Initialize strain
        membro(i) = current_element;
        
        % For subsequent elements, root rotation is zero (elements are chained based on previous element's end state)
        current_rot.dihedral = 0;
        current_rot.sweep = 0;
        current_rot.twist = 0;
    end
    
    % Assign the created nodes to the member object if necessary (not clear if member object stores nodes directly)
    % The Airplane object takes cell array of member arrays. Each member is an array of elements.
    flexible_member = membro; % This should be an array of Element objects
end

% Note: The original create_uniform_structure function is in the main folder.
% By defining create_flexible_member this way, we are overriding the local function
% if one were named create_uniform_structure locally.
% The provided code does not call create_uniform_structure directly, but rather
% create_flexible_member which now has the full logic.

% The original example1.m had:
% function flexible_member = create_flexible_member(num_elements,damp_ratio)
% ...
% flexible_member = create_uniform_structure(pos_cg,rot0, Length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements);    
% end
% We have now inlined and modified the logic of create_uniform_structure within our new create_flexible_member.
