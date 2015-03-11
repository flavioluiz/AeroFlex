%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      AeroFlex - Example 1  - Structural Dynamics                    %
% - computes eigenvalues/eigenvectors of linear structural dynamics   %
% - finds equilibrium (only gravity)                                  %
% - simulation (beam starts undeformed, goes to equilibrium)          %
%                                                                     %
% Several things not working properly yet:                            %
%     - graphical results  (problems with video rate)                 %
%     - should validate with exact beam results (at least eigenvalues)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function example1
    clc
    
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
    numele = 3; %number of elements
    damping = 0.04; %damping coefficient (eg.: 0.0001)
    rigmult = 1; %multiplier for the rigidity matrix (eg.: 1)
    ap = load_structure(numele,damping,rigmult); % this creates a flexible
                                        %airplane object with numele elements
                                        % check the function loadstruct
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%% FINDS TRIM CONDITION %%%%%%%%%%%%%%%%%%%%
    %%%% FLIGHT CONDITIONS -- won't affect the results if there is no aerodynamics
    altitude = 20000; % meters
    V = 15;           % m/s
    Vwind = 0;        % m/s 
    tic;
    tracao = 0;
    deltaflap = 0;
    [vecequilibrio, Xeq] = trimairplane(ap,V,altitude,Vwind,tracao,deltaflap);
    toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% LINEARIZATION OF EQUATIONS OF MOTION %%%%%%%%%%
    tic;
    manete = vecequilibrio(3);    deltaflap = vecequilibrio(2);
    thetaeq = vecequilibrio(1);    straineq = Xeq;
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
    deltau = 0;
    deltav = 0;
    deltaw = 0;
    deltaalfa = 0;
    tSIM = input('Tempo de simulação: ');
    
    theta = vecequilibrio(1);
    deltaflap = vecequilibrio(2);
    manete = vecequilibrio(3);
    
    %%%%%%%%%% Nonlinear implicit simulation %%%%%%%%%%%%%
    T=[0.00 0.49 0.50 0.74 0.75 0.99 1.00 100];
    elev=[0 0 1 1 -1 -1 0 0]*0;
    doublet = @(t) (deltaflap+interp1(T,elev,t));
    strain0 = Xeq*0;
    betaeq = [0;V*cos(theta);-V*sin(theta);0;0;0];
    keq = [theta;0;0;altitude];
    beta0 = [deltau; V*cos(theta+deltaalfa)+deltav; -V*sin(theta+deltaalfa)+deltaw;0;0;0];
    k0 = [theta+deltaalfa;0;0;altitude];
    tic;
    [tNL xNL xpNL lambdaNL betaNL kineticNL] = simulate(ap, [0 tSIM], strain0, beta0, k0, Vwind, @(t)manete, doublet, 'implicit');
    toc;
    
    for i = 1:size(tNL,1)
        dxNL(i,1:ap.NUMele*4) =  xNL(i,1:ap.NUMele*4) - Xeq(1:ap.NUMele*4);
    end
    longfig = figure;
    subplot(3,2,1); plot(tNL,dxNL,'r'); legend('strain');hold all;
    subplot(3,2,2); plot(tNL,xpNL,'r'); legend('strainp');hold all;
    subplot(3,2,3); plot(tNL,betaNL(:,2)-V*cos(theta),'r'); xlabel('t'); ylabel('beta(2)'); hold all;%velocidade eixo y
    subplot(3,2,4); plot(tNL,betaNL(:,3)+V*sin(theta),'r'); xlabel('t'); ylabel('beta(3)'); hold all;%velocidade eixo z
    subplot(3,2,5); plot(tNL,betaNL(:,4),'r'); xlabel('t'); ylabel('q');hold all;%q
    subplot(3,2,6); plot(tNL,kineticNL(:,4),'r'); xlabel('t'); ylabel('H');hold all; %H
    
    dt = 0.1;
    [ts Xs] = changedatarate(tNL,xNL,dt);
    ts = tNL; Xs = xNL;
    for i = 1:size(ts,1)
        update(ap,Xs(i,:),zeros(size(Xs(i,:))),zeros(size(Xs(i,:))),zeros(sum(ap.membNAEDtotal),1));
        desloctip(i) = ap.membros{1}(numele).node3.h(3);
    end
    figure('name','Wing tip displacement');
    plot(ts,desloctip); xlabel('Time (s)'); ylabel('Tip displacement (m)');
    
    
        dt = 1e-4;
        %[ts Xs] = changedatarate(tNL,xNL,dt);
        for i = 1:size(ts,1)
            Xs(i,:) =  Xs(i,:) + Xeq(1:ap.NUMele*4)*0;
        end
        airplanemovie(ap, ts, Xs,dt); axis equal;
    

end

function ap = load_structure(numele, amort, rigidez)
    % member initialization
    membro = create_flexible_member(numele,amort,rigidez);
    
    %set member origin node position and orientation:
    membro(1).seth0([0 -0.0 0 1 0 0 0 1 0 0 0 1]'); 
    update(membro); % initialize displacements for each member node
    fus = []; % no fuselage
    motor1 = []; % no engines
    ap = airplane({membro}, fus, [motor1]);
end

function membro = create_flexible_member(num_elements,amort, rigidez)
    

    % beam length
    Length = 16;
    
    % sectional rigidity matrix
    K11 = 1e10; %EA
    K22 = 1e4; %GJ
    K33 = 2e4; %flat bend: EI
    K44 = 4e6; %chord bend: EI
    KG = rigidez*diag([K11 K22 K33 K44]);
    
    % sectional damping matrix
    CG = diag([amort*K11 amort*K22 amort*K33 amort*K44]);    % tese: 0.00005*KG;
    
    % aerodynamic data
    aeroparams = [];
    
    pos_cg = [0 0.3 0]; % position of section gravity center
                        % relative to elastic axis
    geometry.a = 0.5;
    geometry.b = 0.5;    
    I22 = 0.0;
    I33 = 0.1;
    I11 = 0.1;
    mcs = 0.75;
    Inertia = diag([I11 I22 I33]);
    
    % the following function creates a uniform structure automatically; if
    % you need a more complicate wing (with non-uniform parameters, check
    % how the following function creates the structure. you should modify
    % this function to define the correct parameters for each structural
    % node)
    membro = create_uniform_structure(pos_cg, Length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements);
    
    
end
 