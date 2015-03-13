%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      AeroFlex Beta Ver1 - Example - FlyingWing Airplane        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%

function FlyingWingExample
    clc
    clear all
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AeroFlex Main Folder:                              %
    mainFolder = strcat(pwd,'\main');
    addpath(genpath(mainFolder));    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % General AeroFlex parameters - MANDATORY PARAMETERS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global softPARAMS;
    softPARAMS.isIS = 1; %is it International System?
    softPARAMS.isPINNED = 0; % engastar a asa para cálculo do equilíbrio? (estudos aeroelásticos)
    softPARAMS.vecFREEDEG = [1 1 1 1 1 1]; % type 0 to remove a body state degree of freedom [u v w p q r]
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
    softPARAMS.plota3d = 0; %plota graficos 3d no cálculo do equilíbrio
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    %%%%%%%%%%%% AIRPLANE INITIALIZATION %%%%%%%%%%%%%%%%%
    numele = 3; %number of elements
    amort = 0.04; %damping coefficient (eg.: 0.0001)
    rigidez = 1; %multiplier for the rigidity matrix (eg.: 1)
    ap = carregaasavoadora(numele,amort,rigidez); % this creates a flexible
                                        %airplane object with numele elements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % the following variables are script specific FLAGs (they are specific for
    % this file and won't affect anything else in AeroFlex)
    makemovie = 0;       % record a movie from simulation results
    
    
    
    %%%%%%%%%%%% FINDS TRIM CONDITION %%%%%%%%%%%%%%%%%%%%
    %%%% FLIGHT CONDITIONS
    altitude = 20000; % meters
    V = 15;           % m/s
    Vwind = 0;        % m/s (Take care! By now Vwind is aligned with y direction/body frame!)
    tic;
    tracao = 0;
    deltaflap = 0;
    [vecequilibrio Xeq] = trimairplane(ap,V,altitude,Vwind,tracao,deltaflap);
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
    elev=[0 0 1 1 -1 -1 0 0];
    doublet = @(t) (deltaflap+interp1(T,elev,t));
    strain0 = Xeq;
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
    
    dt = 0.05;
    [ts Xs] = changedatarate(tNL,xNL,dt);
    for i = 1:size(ts,1)
        update(ap,Xs(i,:),zeros(size(Xs(i,:))),zeros(size(Xs(i,:))),zeros(sum(ap.membNAEDtotal),1));
        desloctip(i) = ap.membros{1}(numele).node3.h(3);
    end
    save('results\flyingwingdoubletinput.mat','tNL','xNL','betaNL','kineticNL','ts','desloctip','betaeq','keq');
    figure('name','Wing tip displacement');
    plot(ts,desloctip); xlabel('Time (s)'); ylabel('Tip displacement (m)');
    
    if makemovie
        dt = 0.1;
        [ts Xs] = changedatarate(tNL,xNL,dt);
        for i = 1:size(ts,1)
            Xs(i,:) =  Xs(i,:) + Xeq(1:ap.NUMele*4);
        end
        airplanemovie(ap, ts, Xs,dt);
    end
    simulaRIGID = input('\nDigite 1 para apresentar resultados de simulação da aeronave rigida:');

    if simulaRIGID
        dinamicarigida(V,altitude,longfig, tSIM, deltav, deltaw, deltaalfa,@(t)0, @(t)interp1(T,elev,t));
    end
        
end

function ap = carregaasavoadora(numele, amort, rigidez)
    % inicializacao dos membros
    isRIGHT = 1;
    rot.sweep = 0;
    rot.dihedral = 0;
    rot.twist = 0;
    membroA = geradorhighlyflex(numele,0,amort,rigidez, rot, isRIGHT);
    rot.dihedral = pi;
    membroB = geradorhighlyflex(numele,0,amort,rigidez, rot, ~isRIGHT);
    
    fus = rigidfus(0, [0 0.0 0], zeros(3,3) + 0*[0.2^2*10 0 0;0 0 0; 0 0 0.2^2*10]);
    
    membroA(1).seth0([0 -.3 0 1 0 0 0 1 0 0 0 1]');
    membroB(1).seth0([0 -.3 0 1 0 0 0 1 0 0 0 1]');
    update(membroA); %seta deslocamentos
    update(membroB); %seta deslocamentos
    
    
    engparams.Fmax = 1; engparams.V0 = 1; engparams.rho0 = 1;
    engparams.nv= -1; engparams.nrho = 0; engparams.alphaf = 0;
    engparams.betaf = 0;
    motor1 = engine(1, [1 1 1], engparams); %numManete, posicao do motor[MEMB,ELM,ND], params
    ap = airplane({membroA membroB}, fus, [motor1]);
end

function membro = geradorhighlyflex(n,lixo,amort, rigidez, rot, isRIGHT)
    if isRIGHT
        rightMTP = 1;
    else
        rightMTP = -1;
    end
    K11 = 1e10; %EA??
    K22 = 1e4; %GJ
    K33 = 2e4; %flat bend: EI
    K44 = 4e6; %chord bend: EI
    KG = rigidez*diag([K11 K22 K33 K44]);
    
    %e = 0.25;
    %qd = (pi/2)^2*(K22)/16*(1/(16*1*e*2*pi))
    %Vd = sqrt(2*qd/0.0889)
    
    CG = diag([K11 amort*K22 amort*K33 amort*K44]);    % tese: 0.00005*KG;
    ds = 16/n;
    
    I22 = 0.0;
    I33 = 0.1;
    mcs = 0.75;
    I11 = 0.1;
    
    aeroparams.n = 2;
    c = 1;
    aeroparams.b = c/2;
    aeroparams.N = 0;
    aeroparams.a = 0;
    aeroparams.m = 2;
    aeroparams.alpha0 = -5*pi/180*rightMTP;
    aeroparams.clalpha = 2*pi;
    aeroparams.cldelta = 0.01*rightMTP;
    aeroparams.cm0 = 0;
    aeroparams.cmdelta = -0.1*rightMTP;
    aeroparams.ndelta = 1; %numero da superficie de controle ativada
    aeroparams.cd0 = 0.02;

    ycg = 0.3;
    for i = 1:n
        rigidunit.m = 0; rigidunit.cg = [0 0.3 0]; rigidunit.I = zeros(3);
        if i == 1
            rigidunit.m = 5;
        end
        noh((i-1)*3+1) = node(mcs, [0 ycg 0], diag([I11 I22 I33]),aeroparams, rigidunit,((i-1)*ds)/20);
        rigidunit.m = 0;
        noh((i-1)*3+2) = node(mcs, [0 ycg 0], diag([I11 I22 I33]),aeroparams,rigidunit,(ds/2+(i-1)*ds)/20);
        noh((i-1)*3+3) = node(mcs, [0 ycg 0], diag([I11 I22 I33]),aeroparams,rigidunit,i*ds/20);
        if i > 1
            rot.dihedral = 0;
            rot.sweep = 0;
            rot.twist = 0;
        end
        membro(i) = element(noh((i-1)*3+1),noh((i-1)*3+2),noh((i-1)*3+3),rot,ds,KG, CG);
        membro(i).setstrain([0 0 0 0],[0 0 0 0]);
    end
    membro(1).node1.aero.setmembermatrices(membro);
end

function xp = dinamicalinearnumerica(t,x,A)
    xp = A*x;
end