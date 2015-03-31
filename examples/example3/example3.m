%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      AeroFlex - Example 3  - Aileron Reversals                      %
% - include aerodynamic control surfaces                              %
% - include second semi-wing                                          %
% - perform simulations to verify airleron reversals phenomenon       %
% - one rigid body degree of freedom (rotations only)                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function example3
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
    softPARAMS.vecFREEDEG = [0 0 0 0 1 0]; % type 0 to remove a body state
                                           % degree of freedom [u v w p q r]
    softPARAMS.isGRAV = 1; % include gravity?
    softPARAMS.g = 0; % gravity in m/s^2    
    softPARAMS.isITER = 1; % iterative equilibrium determination?
    softPARAMS.numITER = 10; % number of iterations for equilibrium determination
    softPARAMS.modAED = 1; % AERODYNAMIC MODEL: 
                                    %0-Steady;
                                    %1-Quasi-steady;
                                    %2-Quasi-steady with added mass;
                                    %3-Unsteady(Peters);
    softPARAMS.updateStrJac = 0; % Structural Jacobians updates:
                                    % 0 - Never;
                                    % 1 - Only in equilib calculation;
                                    % 2 - Always
    softPARAMS.plota3d = 0; % 3d graphics plot while equilibrium is calculated
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    %%%%%%%%%%%% AIRPLANE INITIALIZATION %%%%%%%%%%%%%%%%%
    numele = 3; %number of elements
    damping = 0.04; %damping coefficient (damping proportional to rigidity matrix)
    ap = load_structure(numele,damping); % this creates a flexible
                                        %airplane object with numele elements
                                        % check the function loadstruct
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    %%%%%%%%%%%% FINDS EQUILIBRIUM CONDITION %%%%%%%%%%%%%%
    % FLIGHT CONDITIONS
    altitude = 19931.7; % meters
    V = 0;           % rigid body speed m/s
    throttle = 0;
    deltaflap = 0;
    Vwind = 10; % wind speed
    %%%%%% LINEARIZATION OF EQUATIONS OF MOTION %%%%%%%%%%
    [rb_eq, strain_eq] = trimairplane(ap,V,altitude,Vwind,throttle,deltaflap);        
      
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% NONLINEAR SIMULATION %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial conditions

    tSIM = input('Simulation time: (seconds)');    
    T=[0.00 0.1 0.11 4.00 4.01 5 10.0 100];
    elev=[0 0 1 1 1 1 1 1]*20;
    aerodynamic_surface_pos = @(t) [interp1(T,elev,t)]; % no aerodynamic surfaces!
    engine_position = 0;     % no engine!
    beta0 = [0; 0;0;0;0;0]; % rigid body speeds
    k0 = [0;0;0;20000];         % rigid body position/orientation
    strain0 = strain_eq*0;
    Vwind = 10;
    % simulation:  
    tic;
    [tNL, strainNL, straindNL, lambdaNL, betaNL, kineticNL] = simulate(ap, [0 tSIM], strain0, beta0, k0, Vwind, @(t)engine_position, @(t)aerodynamic_surface_pos, 'implicit');
    toc;
    
    dt = 0.1;
    [ts, Xs] = changedatarate(tNL,strainNL,dt);
    [ts, kinetics] = changedatarate(tNL,kineticNL,dt);
    tip_displacement = zeros(length(ts),1);
    for i = 1:length(ts)
        update(ap,Xs(i,:),zeros(size(Xs(i,:))),zeros(size(Xs(i,:))),zeros(sum(ap.membNAEDtotal),1));
        tip_displacement(i) = ap.members{1}(numele).node3.h(3);
    end
    figure('color','w','name','Wing tip displacement');
    plot(ts,tip_displacement);
    xlabel('Time (s)'); ylabel('Tip displacement (m)');
    grid on;
    
    figure('color','w');
    airplanemovie(ap, ts', Xs,kinetics,dt,'test','gif'); colormap winter;
    

end

function ap = load_structure(numele, damp_ratio)
    % member initialization
    [right_wing, left_wing] = create_flexible_member(numele,damp_ratio);
    
    %set member origin node position and orientation:
    right_wing(1).seth0([0 -0.0 0 1 0 0 0 1 0 0 0 1]');
    left_wing(1).seth0([0 -0.0 0 1 0 0 0 1 0 0 0 1]');
    update(right_wing); % initialize displacements for each member node
    update(left_wing); % initialize displacements for each member node
    fus = []; % no fuselage
    motor1 = []; % no engines
    ap = Airplane({right_wing, left_wing}, fus, [motor1]);
end

function [right_wing, left_wing] = create_flexible_member(num_elements,damp_ratio)   
    % beam length
    Length = 16;
    
    % sectional rigidity matrix
    K11 = 1e10; %EA
    K22 = 1e4; %GJ
    K33 = 2e4; %flat bend: EI
    K44 = 4e6; %chord bend: EI
    KG = 1000*diag([K11 K22 K33 K44]);
    
    % sectional damping matrix
    CG = damp_ratio*diag([K11 K22 K33 K44]);
    
    % aerodynamic data
    c = 1; % chord
    aeroparams.b = c/2; %semi-chord
    aeroparams.N = 0; %number of lag states (Peter's Unsteady model)
    aeroparams.a = 0.0; % position of aerodynamic center relative to elastic axis
                        % relative to elastic axis (in terms of semi-chord)
    aeroparams.alpha0 = 0*pi/180; % alpha_0 (in radians)
    aeroparams.clalpha = 2*pi; % cl_alpha  lift coeff/rad
    aeroparams.cm0 = 0;        % cm_0      moment coeff
    aeroparams.cd0 = 0.02;     % cd_0
    
    % aerodynamic data for flap/aileron, if exists
    aeroparams.ndelta = 1;   % Identification of the flap (1,2,3,...)
    aeroparams.cldelta = 0.01; % cl_delta
    aeroparams.cmdelta = -0.1; % cm_delta
            
    % cg position, mass and inertia data
    
    pos_cg = [0 0 0]; % position of section gravity center
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
    %right wing
    right_wing = create_uniform_structure(pos_cg, rot0, Length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements);    
    
    
    % rotation of initial node with respect to body frame (allow
    % creating uniform beams with a twist, dihedral or sweep angle)
    rot0.dihedral = pi; %angles in RADIANS
    rot0.sweep = 0;
    rot0.twist = 0;
    % left wing
    left_wing = create_uniform_structure(pos_cg, rot0, Length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements);    
    
end
