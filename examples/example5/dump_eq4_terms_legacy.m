function dump_eq4_terms_legacy
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

    [rb_eq, strain] = trimairplane(ap, V, H, Vwind, 0, 0);
    theta = rb_eq(1);
    delta = rb_eq(2);
    tracao = rb_eq(3);
    strainp = zeros(size(strain));
    strainpp = zeros(size(strain));
    lambda = zeros(sum(ap.membNAEDtotal),1);
    beta = [0; V*cos(theta); -V*sin(theta); 0; 0; 0];
    betap = zeros(6,1);
    rho = atmosfera(H);

    update(ap, strain, strainp, strainpp, lambda);
    updateStrJac(ap);

    FAERO = [];
    MAERO = [];
    for i = 1:ap.NUMmembers
        [FAEROm, MAEROm, ~] = aeroforceandmoment(ap.members{i}(1).strainm, ap.members{i}(1).strainpm, ap.members{i}(1).strainppm, ap.members{i}(1).lambdam, beta, betap, ap.members{i}, Vwind, rho, [delta 0 0]);
        FAERO = [FAERO; FAEROm]; %#ok<AGROW>
        MAERO = [MAERO; MAEROm]; %#ok<AGROW>
    end

    Fdist = ap.B * FAERO;

    nRight = ap.membSIZES(1);
    nLeft = ap.membSIZES(2);
    ne = nRight + nLeft;
    nNodes = ne * 3;
    p_nodes = zeros(3, nNodes);
    for ie = 1:nRight
        h1 = ap.members{1}(ie).node1.h;
        h2 = ap.members{1}(ie).node2.h;
        h3 = ap.members{1}(ie).node3.h;
        base = (ie-1)*3;
        p_nodes(:, base+1) = h1(1:3);
        p_nodes(:, base+2) = h2(1:3);
        p_nodes(:, base+3) = h3(1:3);
    end
    offset = nRight * 3;
    for ie = 1:nLeft
        h1 = ap.members{2}(ie).node1.h;
        h2 = ap.members{2}(ie).node2.h;
        h3 = ap.members{2}(ie).node3.h;
        base = offset + (ie-1)*3;
        p_nodes(:, base+1) = h1(1:3);
        p_nodes(:, base+2) = h2(1:3);
        p_nodes(:, base+3) = h3(1:3);
    end

    Jpb = ap.Jpb;
    Bmat = ap.B;
    eq4_force_total = (Jpb' * Bmat * FAERO);

    save('eq4_dump_legacy.mat', 'theta', 'delta', 'tracao', 'strain', 'FAERO', 'MAERO', 'Fdist', 'p_nodes', 'beta', 'Jpb', 'Bmat', 'eq4_force_total');
    fprintf('Saved eq4_dump_legacy.mat\n');
    fprintf('theta=%.9f delta=%.9f tracao=%.9f\n', theta, delta, tracao);
    fprintf('eq4 force-path total = % .9e\n', eq4_force_total(4));
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
