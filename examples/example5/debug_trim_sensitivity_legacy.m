function debug_trim_sensitivity_legacy
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

    numele = 3;
    damping = 0.04;
    ap = load_structure_legacy(numele, damping);

    altitude = 19931.7;
    V = 10;
    Vwind = 0;

    [rb_eq, strain_eq] = trimairplane(ap, V, altitude, Vwind, 0, 0);
    theta = rb_eq(1);
    delta = rb_eq(2);
    tracao = rb_eq(3);

    [~, bp0] = eval_bp_legacy(ap, strain_eq, theta, delta, tracao, V, altitude, Vwind);

    h = 1e-4;
    [~, bp_dp] = eval_bp_legacy(ap, strain_eq, theta, delta + h, tracao, V, altitude, Vwind);
    [~, bp_dm] = eval_bp_legacy(ap, strain_eq, theta, delta - h, tracao, V, altitude, Vwind);
    dbp_ddelta = (bp_dp - bp_dm)/(2*h);

    [~, bp_tp] = eval_bp_legacy(ap, strain_eq, theta + h, delta, tracao, V, altitude, Vwind);
    [~, bp_tm] = eval_bp_legacy(ap, strain_eq, theta - h, delta, tracao, V, altitude, Vwind);
    dbp_dtheta = (bp_tp - bp_tm)/(2*h);

    fprintf('Legacy AeroFlex example5 sensitivity:\n');
    fprintf('theta=%.6f rad delta=%.6f rad tracao=%.6f\n', theta, delta, tracao);
    fprintf('bp at trim = [% .6e % .6e % .6e % .6e % .6e % .6e]\n', bp0(1), bp0(2), bp0(3), bp0(4), bp0(5), bp0(6));
    fprintf('d(bp)/d(delta) = [% .6e % .6e % .6e % .6e % .6e % .6e]\n', dbp_ddelta(1), dbp_ddelta(2), dbp_ddelta(3), dbp_ddelta(4), dbp_ddelta(5), dbp_ddelta(6));
    fprintf('d(bp)/d(theta) = [% .6e % .6e % .6e % .6e % .6e % .6e]\n', dbp_dtheta(1), dbp_dtheta(2), dbp_dtheta(3), dbp_dtheta(4), dbp_dtheta(5), dbp_dtheta(6));
end

function [xp, bp] = eval_bp_legacy(ap, strain, theta, delta, tracao, V, H, Vwind)
    strainp = zeros(size(strain));
    strainpp = zeros(size(strain));
    lambda = zeros(sum(ap.membNAEDtotal),1);
    beta = [0; V*cos(theta); -V*sin(theta); 0; 0; 0];
    betap = zeros(6,1);
    kinetic = [theta; 0; 0; H];
    [xp, bp] = dinamicaflex(0, strain, strainp, strainpp, lambda, beta, betap, kinetic, ap, Vwind, tracao, [delta 0 0], 0);
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
