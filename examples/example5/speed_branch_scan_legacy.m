function speed_branch_scan_legacy
    clc;
    addpath('/Users/flavioribeiro/github/AeroFlex/main');

    global softPARAMS;
    softPARAMS.isIS = 1;
    softPARAMS.isPINNED = 0;
    softPARAMS.vecFREEDEG = [1 1 1 1 1 1];
    softPARAMS.isGRAV = 1;
    softPARAMS.g = 9.8;
    softPARAMS.isITER = 1;
    softPARAMS.numITER = 10;
    softPARAMS.modAED = 1;
    softPARAMS.updateStrJac = 1;
    softPARAMS.plota3d = 0;

    ap = load_structure_legacy(3, 0.04);
    H = 19931.7;
    Vwind = 0;

    Vvec = 8:0.5:12;

    fprintf('Legacy speed branch scan\n');
    fprintf('V, theta, delta, tracao, bp2, bp3, bp4\n');
    for i = 1:length(Vvec)
        V = Vvec(i);
        [rb_eq, strain_eq] = trimairplane(ap, V, H, Vwind, 0, 0);
        theta = rb_eq(1);
        delta = rb_eq(2);
        tracao = rb_eq(3);
        bp = eval_bp(ap, strain_eq, theta, delta, tracao, V, H, Vwind);
        fprintf('%5.2f, % .6f, % .6f, % .6f, % .6e, % .6e, % .6e\n', V, theta, delta, tracao, bp(2), bp(3), bp(4));
    end
end

function bp = eval_bp(ap, strain, theta, delta, tracao, V, H, Vwind)
    strainp = zeros(size(strain));
    strainpp = zeros(size(strain));
    lambda = zeros(sum(ap.membNAEDtotal),1);
    beta = [0; V*cos(theta); -V*sin(theta); 0; 0; 0];
    betap = zeros(6,1);
    kinetic = [theta; 0; 0; H];
    [~, bp, ~, ~] = dinamicaflex(0, strain, strainp, strainpp, lambda, beta, betap, kinetic, ap, Vwind, tracao, [delta 0 0], 0);
end

function ap = load_structure_legacy(numele, damp_ratio)
    [right_wing, left_wing] = create_flexible_member(numele, damp_ratio);
    right_wing(1).seth0([0 -.3 0 1 0 0 0 1 0 0 0 1]');
    left_wing(1).seth0([0 -.3 0 1 0 0 0 1 0 0 0 1]');
    update(right_wing);
    update(left_wing);
    fus = RigidFuselage(10, [0 0 0], zeros(3,3));
    engparams.Fmax = 1; engparams.V0 = 1; engparams.rho0 = 1;
    engparams.nv= -1; engparams.nrho = 0; engparams.alphaf = 0;
    engparams.betaf = 0;
    motor1 = Engine(1, [1 1 1], engparams);
    ap = Airplane({right_wing, left_wing}, fus, [motor1]);
end

function [right_wing, left_wing] = create_flexible_member(num_elements,damp_ratio)
    Length = 16;
    K11 = 1e10; K22 = 1e4; K33 = 2e4; K44 = 4e6;
    KG = diag([K11 K22 K33 K44]);
    CG = damp_ratio*diag([K11 K22 K33 K44]);

    c = 1;
    aeroparams.b = c/2;
    aeroparams.N = 0;
    aeroparams.a = 0.0;
    aeroparams.alpha0 = -5*pi/180;
    aeroparams.clalpha = 2*pi;
    aeroparams.cm0 = 0;
    aeroparams.cd0 = 0.02;
    aeroparams.ndelta = 1;
    aeroparams.cldelta = 0.01;
    aeroparams.cmdelta = -0.1;

    pos_cg = [0 0.3 0];
    geometry.a = 0.0;
    geometry.b = 0.5;
    Inertia = diag([0.1 0.0 0.1]);
    mcs = 0.75;

    rot0.dihedral = 0;
    rot0.sweep = 0;
    rot0.twist = 0;
    right_wing = create_uniform_structure(pos_cg, rot0, Length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements);

    rot0.dihedral = pi;
    aeroparams.alpha0 = 5*pi/180;
    aeroparams.cldelta = -0.01;
    aeroparams.cmdelta = 0.1;
    left_wing = create_uniform_structure(pos_cg, rot0, Length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements);
end
