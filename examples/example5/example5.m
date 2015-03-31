%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      AeroFlex - Example 5  - Flexible flying wing stability         %
% to do: find the instability speed;                                  %
%        plot unstable mode (it should be an unstable                 %
%           rigid body/flexible one)                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function example5
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
    softPARAMS.isPINNED = 0; % PINNED RIGID BODY DOF
    softPARAMS.vecFREEDEG = [1 1 1 1 1 1]; % type 0 to remove a body state
                                           % degree of freedom [u v w p q r]
    softPARAMS.isGRAV = 1; % include gravity?
    softPARAMS.g = 9.8; % gravity in m/s^2    
    softPARAMS.isITER = 1; % iterative equilibrium determination?
    softPARAMS.numITER = 10; % number of iterations for equilibrium determination
    softPARAMS.modAED = 1; % AERODYNAMIC MODEL: 
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
    V = 10;           % rigid body speed m/s
    throttle = 0;
    deltaflap = 0;
    Vwind = 0; % wind speed
    [rb_eq, strain_eq] = trimairplane(ap,V,altitude,Vwind,throttle,deltaflap);                   
    
    theta = rb_eq(1);
    deltaflap = rb_eq(2);
    engine_position = rb_eq(3);

   %%%%%% LINEARIZATION OF EQUATIONS OF MOTION %%%%%%%%%%
    tic;
    manete = rb_eq(3);    deltaflap = rb_eq(2);
    thetaeq = rb_eq(1);    straineq = strain_eq;
    betaeq = [0 V*cos(thetaeq) -V*sin(thetaeq) 0 0 0]';
    keq = [thetaeq 0 0 altitude]';
    [Alin Aaeroelast Abody] = linearize(ap, straineq, betaeq, keq, manete, deltaflap, Vwind);
    toc;
    
    fprintf('Linearização numérica:\n');
    autoval=eig(Alin);
    fprintf('\nAutovalores da dinâmica completa com parte real maior que -10:');
    autoval(find(autoval>-10))
    
    fprintf('\nAutovalores da aeroelasticidade com parte real maior que -10 (asa engastada):');
    autoaeroelast = eig(Aaeroelast);
    autoaeroelast(find(autoaeroelast>-10))
    
    fprintf('\nAutovalores da matriz do corpo:');
    autobody = eig(Abody)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% NONLINEAR SIMULATION %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial conditions

    tSIM = input('Simulation time: (seconds)');    
    T=[0.00 0.49 0.50 1.99 2.00 3.49 3.5 100];
    elev=[0 0 1 1 -1 -1 0 0]/5;
    aerodynamic_surface_pos = @(t) (deltaflap+interp1(T,elev,t));
    beta0 = [0; V*cos(theta);-V*sin(theta);0;0;0]; % rigid body speeds
    k0 = [theta;0;0;altitude];         % rigid body position/orientation
    strain0 = strain_eq;
    Vwind = 0;
    % simulation:  
    [tNL, strainNL, straindNL, lambdaNL, betaNL, kineticNL] = simulate(ap, [0 tSIM], strain0, beta0, k0, Vwind, @(t)engine_position, aerodynamic_surface_pos, 'implicit');
    
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
        
    longfig = figure('color','w');
    subplot(2,2,1); plot(tNL,betaNL(:,2),'r'); xlabel('t'); ylabel('v (m/s)'); hold all;%velocidade eixo y
    subplot(2,2,2); plot(tNL,betaNL(:,3),'r'); xlabel('t'); ylabel('w (m/s)'); hold all;%velocidade eixo z
    subplot(2,2,3); plot(tNL,betaNL(:,4),'r'); xlabel('t'); ylabel('q (rad/s)');hold all;%q
    subplot(2,2,4); plot(tNL,kineticNL(:,4),'r'); xlabel('t'); ylabel('Altitude (m)');hold all; %H
    deltaalfa = 0; deltau = 0; deltav = 0; deltaw = 0;
    dinamicarigida(V,altitude,longfig, tSIM, deltav, deltaw, deltaalfa,@(t)0, @(t)interp1(T,elev,t));
    
    airplanemovie(ap, ts', Xs,kinetics,dt,'test','gif'); colormap winter;
end

function ap = load_structure(numele, damp_ratio)
    % member initialization
    [right_wing, left_wing] = create_flexible_member(numele,damp_ratio);
    
    %set member origin node position and orientation:
    right_wing(1).seth0([0 -.3 0 1 0 0 0 1 0 0 0 1]');
    left_wing(1).seth0([0 -.3 0 1 0 0 0 1 0 0 0 1]');
    update(right_wing); % initialize displacements for each member node
    update(left_wing); % initialize displacements for each member node
    fus = RigidFuselage(10, [0 0 0], zeros(3,3) + 0*[0.2^2*10 0 0;0 0 0; 0 0 0.2^2*10]);
    engparams.Fmax = 1; engparams.V0 = 1; engparams.rho0 = 1;
    engparams.nv= -1; engparams.nrho = 0; engparams.alphaf = 0;
    engparams.betaf = 0;
    motor1 = Engine(1, [1 1 1], engparams); %numManete, posicao do motor[MEMB,ELM,ND], params

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
    KG = diag([K11 K22 K33 K44]);
    
    % sectional damping matrix
    CG = damp_ratio*diag([K11 K22 K33 K44]);
    
    % aerodynamic data
    c = 1; % chord
    aeroparams.b = c/2; %semi-chord
    aeroparams.N = 0; %number of lag states (Peter's Unsteady model)
    aeroparams.a = 0.0; % position of aerodynamic center relative to elastic axis
                        % relative to elastic axis (in terms of semi-chord)
    aeroparams.alpha0 = -5*pi/180; % alpha_0 (in radians)
    aeroparams.clalpha = 2*pi; % cl_alpha  lift coeff/rad
    aeroparams.cm0 = 0;        % cm_0      moment coeff
    aeroparams.cd0 = 0.02;     % cd_0
    
    % aerodynamic data for flap/aileron, if exists
    aeroparams.ndelta = 1;   % Identification of the flap (1,2,3,...)
    aeroparams.cldelta = 0.01; % cl_delta
    aeroparams.cmdelta = -0.1; % cm_delta
            
    % cg position, mass and inertia data
    
    pos_cg = [0 0.3 0]; % position of section gravity center
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
    
    aeroparams.alpha0 = 5*pi/180; % alpha_0 (in radians)
    aeroparams.cldelta = -0.01; % cl_delta
    aeroparams.cmdelta = 0.1; % cm_delta
    left_wing = create_uniform_structure(pos_cg, rot0, Length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements);
end
