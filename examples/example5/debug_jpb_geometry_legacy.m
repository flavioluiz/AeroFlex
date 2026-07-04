function debug_jpb_geometry_legacy
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
    V = 10;
    Vwind = 0;

    [rb_eq, strain_eq] = trimairplane(ap, V, H, Vwind, 0, 0);
    theta = rb_eq(1);

    strainp = zeros(size(strain_eq));
    strainpp = zeros(size(strain_eq));
    lambda = zeros(sum(ap.membNAEDtotal),1);
    update(ap, strain_eq, strainp, strainpp, lambda);
    updateStrJac(ap);

    h = ap.members{1}(1).node2.h;
    p = h(1:3);
    blk = 4:6;
    S = ap.Jpb(blk,4:6);
    S_hb = ap.members{1}(1).memberJhb(13:15,4:6);
    p_from_jpb = [S(2,3); S(3,1); S(1,2)];

    fprintf('theta=%.9f\n', theta);
    fprintf('node2.h(1:3)      = [% .9e % .9e % .9e]\n', p(1), p(2), p(3));
    fprintf('Jpb block rows4:6 =\n');
    disp(S);
    fprintf('memberJhb rows13:15=\n');
    disp(S_hb);
    fprintf('p_from_Jpb        = [% .9e % .9e % .9e]\n', p_from_jpb(1), p_from_jpb(2), p_from_jpb(3));
    fprintf('matrixcross(p)    =\n');
    disp(matrixcross(p));
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
