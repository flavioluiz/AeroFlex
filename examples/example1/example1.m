%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      AeroFlex - Example 1  - Structural Dynamics                    %
% - computes eigenvalues/eigenvectors of linear structural dynamics   %
% - finds equilibrium (only gravity)                                  %
% - simulation (beam starts undeformed, goes to equilibrium)          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function example1
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
    softPARAMS.g = 9.8; % gravity in m/s^2    
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
    damping = 0.04; %damping coefficient (damping proportional to rigidity matrix)
    ap = load_structure(numele,damping); % this creates a flexible
                                        %airplane object with numele elements
                                        % check the function loadstruct
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%% MODAL ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
    % find structural natural frequencies and modal shapes
    st_modes = StructuralModes(ap);
    numerical_freqs = st_modes.frequencies(1:6);
    exact_freqs = sort(exact(4,1,1));
    fprintf('Natural frequencies (in Hz):\n');
    fprintf('| Numerical | Exact \n');
    fprintf('| %.3f | %.3f \n', [numerical_freqs, exact_freqs]')
    
    %plot results    
    figure('color','w','name','mode 1');
    subplot(1,2,1); st_modes.plot_mode(ap,1); % plot modal shape 1
    view(30,45); axis equal; colormap winter;
    subplot(1,2,2); st_modes.plot_mode(ap,2); % plot modal shape 1
    view(30,45); axis equal; colormap winter;
    
    fprintf('press any key to equilibrium calculation\n\n');
    pause;
    %%%%%%%%%%%% FINDS EQUILIBRIUM CONDITION AND PLOT%%%%%%%%%%%%%%
    % FLIGHT CONDITIONS -- won't affect the results if there is no aerodynamics
    altitude = -1; % meters
    V = -1;           % m/s
    Vwind = 0;        % m/s     
    throttle = 0;
    deltaflap = 0;
    [rb_eq, strain_eq] = trimairplane(ap,V,altitude,Vwind,throttle,deltaflap);
    
    figure('color','w');
    % plot structure without deformation:
    update(ap,strain_eq*0,zeros(size(strain_eq)),zeros(size(strain_eq)),zeros(sum(ap.membNAEDtotal),1));
    plotairplane3d(ap); 
    % plot deformed structure (equilibrium condition):
    update(ap,strain_eq,zeros(size(strain_eq)),zeros(size(strain_eq)),zeros(sum(ap.membNAEDtotal),1));
    plotairplane3d(ap); 
    view(30,45); axis equal; colormap winter;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% NONLINEAR SIMULATION %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial conditions

    tSIM = input('Simulation time: (seconds)');    
    aerodynamic_surface_pos = 0; % no aerodynamic surfaces!
    engine_position = 0;     % no engine!
    beta0 = [0; 0;0;0;0;0]; % rigid body speeds
    k0 = [0;0;0;0];         % rigid body position/orientation
    strain0 = strain_eq*0;
    % simulation:  
    tic;
    [tNL, strainNL, straindNL, lambdaNL, betaNL, kineticNL] = simulate(ap, [0 tSIM], strain0, beta0, k0, Vwind, @(t)engine_position, @(t)aerodynamic_surface_pos, 'implicit');
    toc;
    
    dt = 0.1;
    [ts, Xs] = changedatarate(tNL,strainNL,dt);
    tip_displacement = zeros(length(ts),1);
    for i = 1:size(ts,1)
        update(ap,Xs(i,:),zeros(size(Xs(i,:))),zeros(size(Xs(i,:))),zeros(sum(ap.membNAEDtotal),1));
        tip_displacement(i) = ap.members{1}(numele).node3.h(3);
    end
    figure('color','w','name','Wing tip displacement');
    plot(ts,tip_displacement);
    xlabel('Time (s)'); ylabel('Tip displacement (m)');
    grid on;
    
    figure('color','w');
    airplanemovie(ap, ts, Xs,[],dt,'test','gif'); colormap winter;
    

end

function ap = load_structure(numele, damp_ratio)
    % member initialization
    flexible_member = create_flexible_member(numele,damp_ratio);
    
    %set member origin node position and orientation:
    flexible_member(1).seth0([0 -0.0 0 1 0 0 0 1 0 0 0 1]'); 
    update(flexible_member); % initialize displacements for each member node
    fus = []; % no fuselage
    motor1 = []; % no engines
    ap = Airplane({flexible_member}, fus, [motor1]);
end

function flexible_member = create_flexible_member(num_elements,damp_ratio)
    

    % beam length
    Length = 16;
    
    % sectional rigidity matrix
    K11 = 1e10; %EA
    K22 = 1e4; %GJ
    K33 = 2e4; %flat bend: EI
    K44 = 4e6; %chord bend: EI
    KG = diag([K11 K22 K33 K44]);
    
    % sectional damping matrix
    CG = damp_ratio*diag([K11 K22 K33 K44]);
    
    % aerodynamic data
    aeroparams = [];
    
    pos_cg = [0 0.0 0]; % position of section gravity center
                        % relative to elastic axis
    geometry.a = 0.0;
    geometry.b = 0.5;    
    I22 = 0.0;
    I33 = 0.1;
    I11 = 0.1;
    mcs = 0.75; %mass per unit length (kg/m)
    Inertia = diag([I11 I22 I33]);
    
    % rotation of initial node with respect to body frame (allow
    % creating uniform beams with a twist, dihedral or sweep angle)
    rot0.dihedral = 0; %angles in RADIANS
    rot0.sweep = 0;
    rot0.twist = 0;
    % the following function creates a uniform structure automatically; if
    % you need a more complicate wing (with non-uniform parameters, check
    % how the following function creates the structure. you should modify
    % this function to define the correct parameters for each structural
    % node)
    flexible_member = create_uniform_structure(pos_cg,rot0, Length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements);    
    
end
