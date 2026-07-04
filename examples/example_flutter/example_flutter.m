%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      AeroFlex - Example Flutter - Cálculo da velocidade de flutter     %
%  Baseado em example2.m, mas focado apenas no cálculo de flutter        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function example_flutter
    clc;
    clear all;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AeroFlex Main Folder:                              %
    addpath('../../main');    % ajuste conforme necessário para o caminho correto
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % General AeroFlex parameters - MANDATORY PARAMETERS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global softPARAMS;
    softPARAMS.isIS = 1; %is it International System?
    softPARAMS.isPINNED = 1; % PINNED RIGID BODY DOF
    softPARAMS.vecFREEDEG = [0 0 0 0 0 0];
    softPARAMS.isGRAV = 1; % inclui gravidade
    softPARAMS.g = 9.8; % gravity in m/s^2    
    softPARAMS.isITER = 1; % iterative equilibrium determination?
    softPARAMS.numITER = 10; % number of iterations for equilibrium determination
    softPARAMS.modAED = 3; % AERODYNAMIC MODEL: Unsteady(Peters)
    softPARAMS.updateStrJac = 1; % Structural Jacobians updates
    softPARAMS.plota3d = 1; % 3d graphics plot while equilibrium is calculated
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%% STRUCTURE INITIALIZATION %%%%%%%%%%%%%%%%%
    numele = 3; % número de elementos
    damping = 0.0001; % coeficiente de amortecimento
    ap = load_structure_flutter(numele, damping);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%% CONDIÇÃO DE EQUILÍBRIO %%%%%%%%%%%%%%%%%%%
    altitude = 0; % metros (ajuste conforme necessário)
    V = 0;           % velocidade inicial m/s
    throttle = 0;
    deltaflap = 0;
    Vwind = 0; % sem vento
%    [rb_eq, strain_eq] = trimairplane(ap,V,altitude,Vwind,throttle,deltaflap);
    strain_eq = zeros(4*numele,1)';
    %rb_eq = zeros(4,1)
    %%%%%%%%%%%% CÁLCULO DA VELOCIDADE DE FLUTTER %%%%%%%%%
    [flut_speed, flut_eig_val, flut_eig_vec] = flutter_speed(1,25,0.01,ap, strain_eq*0, altitude);
    fprintf('\nVelocidade de flutter encontrada: %.3f m/s\n', flut_speed);
    disp('Autovalores associados ao flutter:');
    disp(flut_eig_val);
end

function ap = load_structure_flutter(numele, damp_ratio)
    %% Propriedades Estruturais    
    a = 0;      % Localização do eixo elástico [m] 
    Ssc = 0;    % Static mass moment per unit span about the elastic axis
    ch = 0;     % Amortecimento estrutural de elevação
    ctheta = 0; % Amortecimento estrutural de arfagem
    est = struct('a',a,'Ssc',Ssc,'ctheta',ctheta,'ch',ch);
    %==========================================================================
    %% INPUTS de asa
    % Geométricos da asa (Viga)
    L = 0.35;      % Comprimento [m]
    t = 0.000813;  % Espessura [m]
    c = 0.04;      % Corda [m]
    b = c/2;       % Meia-corda [m]
    % secao homogenea -> eixo elastico no centro da asa
    %-----------------------------------------
    % Material da asa
    E = 73.1e9;  % Módulo de Young [Pa]
    G = 28e9;    % Módulo de Cisalhamento [Pa]
    rho = 2697;  % Densidade [kg/m³]
    %-----------------------------------------
    % Massa na ponta da asa
    d = -0.005;   % Offset [m] 
    m = 0.03458;  % Massa [kg]

    % Montagem do membro flexível
    flexible_member = create_flexible_member_flutter(numele, damp_ratio, L, t, c, b, E, G, rho, m, d, est);
    flexible_member(1).seth0([0 -0.0 0 1 0 0 0 1 0 0 0 1]'); 
    update(flexible_member); % inicializa deslocamentos
    fus = [];
    motor1 = [];
    ap = Airplane({flexible_member}, fus, [motor1]);
end

function flexible_member = create_flexible_member_flutter(num_elements, damp_ratio, L, t, c, b, E, G, rho, m, d, est)
    % Propriedades de seção
    Area = t * c;
    Iyy = (c * t^3) / 12;
    Izz = (t * c^3) / 12;
    J = (1/3) * (c * t^3); % Aproximação para seção retangular fina
    mcs = rho * Area; % massa por unidade de comprimento
    
    % Rigidez
    K11 = E * Area; % EA
    K22 = G * J;    % GJ
    K33 = E * Iyy;  % EI (flexão)
    K44 = E * Izz;  % EI (flexão)
    KG = diag([K11 K22 K33 K44]);
    
    % Amortecimento
    CG = damp_ratio*diag([K11 K22 K33 K44]);
    
    % Dados aerodinâmicos
    aeroparams.b = b;
    aeroparams.N = 4;
    aeroparams.a = est.a;
    aeroparams.alpha0 = 0;
    aeroparams.clalpha = 2*pi;
    aeroparams.cm0 = 0;
    aeroparams.cd0 = 0.02;
    aeroparams.ndelta = 0;
    aeroparams.cldelta = 0;
    aeroparams.cmdelta = 0;
    
    % Geometria
    pos_cg = [0 0 0];
    geometry.a = est.a;
    geometry.b = b;
    I22 = 0.0;
    I33 = 0.1;
    I11 = 0.1;
    Inertia = diag([I11 I22 I33]);
    
    rot0.dihedral = 0;
    rot0.sweep = 0;
    rot0.twist = 0;
    
    flexible_member = create_uniform_wing(pos_cg, rot0, L, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements, m, d);
    
end 

% creates a structure with constant parameters (mass / unit length,
% sectional Inertia, aerodynamic parameters, etc.)

function membro = create_uniform_wing(pos_cg, initial_node_rot, length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements, m, d)
    ds = length/num_elements;
    
    rot = initial_node_rot;
    rigidunit.m = m; rigidunit.cg = [0,d,0]; rigidunit.I = zeros(3,3); % no rigid units
    
    for i = 1:num_elements
            noh((i-1)*3+1) = Node(mcs, pos_cg, Inertia, aeroparams, rigidunit,((i-1)*ds)/20 ,geometry);
            noh((i-1)*3+2) = Node(mcs, pos_cg, Inertia, aeroparams,rigidunit,(ds/2+(i-1)*ds)/20, geometry);
            noh((i-1)*3+3) = Node(mcs, pos_cg, Inertia, aeroparams,rigidunit,i*ds/20, geometry);
            membro(i) = Element(noh((i-1)*3+1),noh((i-1)*3+2),noh((i-1)*3+3),rot,ds,KG, CG);
            membro(i).setstrain([0 0 0 0],[0 0 0 0]);
            rot.dihedral = 0;
            rot.sweep = 0;
            rot.twist = 0;
    end
%    membro(1).node1.aero.setmembermatrices(membro);
end