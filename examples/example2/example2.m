%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      AeroFlex - Example 2  - Aeroelasticity                         %
% - finds equilibrium (aerodynamics + gravity)                        %
% - linearizes to verify local stability                              %
% - find flutter speed:    linear vs nonlinear beam deflections       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function example2
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
    softPARAMS.modAED = 2; % AERODYNAMIC MODEL: 
                                    %0-Steady;
                                    %1-Quasi-steady;
                                    %2-Quasi-steady with added mass;
                                    %3-Unsteady(Peters);
    softPARAMS.updateStrJac = 0; % Structural Jacobians updates:
                                    % 0 - Never;
                                    % 1 - Only in equilib calculation;
                                    % 2 - Always
    softPARAMS.plota3d = 1; % 3d graphics plot while equilibrium is calculated
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    %%%%%%%%%%%% STRUCTURE INITIALIZATION %%%%%%%%%%%%%%%%%
    numele = 3; %number of elements
    damping = 0.01; %damping coefficient (damping proportional to rigidity matrix)
    ap = load_structure(numele,damping); % this creates a flexible
                                        %airplane object with numele elements
                                        % check the function loadstruct
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    Vwind = 10;
    
    %%%%%%%%%%%% FINDS EQUILIBRIUM CONDITION %%%%%%%%%%%%%%
    % FLIGHT CONDITIONS -- won't affect the results if there is no aerodynamics
    altitude = 20000; % meters
    V = 0;           % m/s
        throttle = 0;
        deltaflap = 0;
        [rb_eq, strain_eq] = trimairplane(ap,V,altitude,Vwind,throttle,deltaflap);    

    
    while Vwind>0
          Vwind = input('speed: ');        % m/s    
 

          %%%%%% LINEARIZATION OF EQUATIONS OF MOTION %%%%%%%%%%

    tic;
    betaeq = [0 0 0 0 0 0]';
    keq = [0 0 0 altitude]';
    [Alin, Aaeroelast, Abody] = linearize(ap, strain_eq*0, betaeq, keq, throttle, deltaflap, Vwind);
    toc;        
    fprintf('Linearização numérica:\n');  
    fprintf('\nAutovalores da aeroelasticidade com parte real maior que -10 (asa engastada):');
    autoaeroelast = eig(Aaeroelast);
    autoaeroelast(find(autoaeroelast>0))
    end
    
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
        tip_displacement(i) = ap.membros{1}(numele).node3.h(3);
    end
    figure('color','w','name','Wing tip displacement');
    plot(ts,tip_displacement);
    xlabel('Time (s)'); ylabel('Tip displacement (m)');
    grid on;
    
    figure('color','w');
    airplanemovie(ap, ts, Xs,dt,'test','gif'); colormap winter;
    

end

function ap = load_structure(numele, damp_ratio)
    % member initialization
    flexible_member = create_flexible_member(numele,damp_ratio);
    
    %set member origin node position and orientation:
    flexible_member(1).seth0([0 -0.0 0 1 0 0 0 1 0 0 0 1]'); 
    update(flexible_member); % initialize displacements for each member node
    fus = []; % no fuselage
    motor1 = []; % no engines
    ap = airplane({flexible_member}, fus, [motor1]);
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
    aeroparams.n = 2; aeroparams.m = 2;
    c = 1;
    aeroparams.b = c/2;
    aeroparams.N = 4;
    aeroparams.a = 0.5;

    aeroparams.alpha0 = -5*pi/180*0;
    aeroparams.clalpha = 2*pi;
    aeroparams.cm0 = 0;
    aeroparams.cd0 = 0.02;
    
    aeroparams.cldelta = 0.01*0;
    aeroparams.cmdelta = -0.1*0;
    aeroparams.ndelta = 0; %numero da superficie de controle ativada
    
    
    % cg position, mass and inertia data
    
    pos_cg = [0 0.25 0]; % position of section gravity center
                        % relative to elastic axis
    geometry.a = 0.0;
    geometry.b = 0.5;    
    I22 = 0.0;
    I33 = 0.1;
    I11 = 0.1;
    mcs = 0.75; %mass per unit length (kg/m)
    Inertia = diag([I11 I22 I33]);
    
    % the following function creates a uniform structure automatically; if
    % you need a more complicate wing (with non-uniform parameters, check
    % how the following function creates the structure. you should modify
    % this function to define the correct parameters for each structural
    % node)
    flexible_member = create_uniform_structure(pos_cg, Length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements);    
    
end
