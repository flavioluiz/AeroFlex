function debug_bp_components_legacy
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
    theta_legacy = rb_eq(1);
    delta_legacy = rb_eq(2);
    tracao_legacy = rb_eq(3);

    theta_af2 = 0.333710;
    delta_af2 = -0.097817;
    tracao_af2 = 23.662489;

    fprintf('--- Legacy trim point components ---\n');
    comp_leg = eval_components(ap, strain_eq, theta_legacy, delta_legacy, tracao_legacy, V, H, Vwind);
    report_components(comp_leg);

    fprintf('\n--- AF2 point in legacy model components ---\n');
    comp_af2 = eval_components(ap, strain_eq, theta_af2, delta_af2, tracao_af2, V, H, Vwind);
    report_components(comp_af2);
    report_force_node_breakdown(ap, comp_af2.FAERO, 'AF2 point in legacy (force path, eq4)');
end

function comp = eval_components(ap, strain, theta, delta, tracao, V, H, Vwind)
    global softPARAMS;

    strainp = zeros(size(strain));
    strainpp = zeros(size(strain));
    lambda = zeros(sum(ap.membNAEDtotal),1);
    beta = [0; V*cos(theta); -V*sin(theta); 0; 0; 0];
    betap = zeros(6,1);
    kinetic = [theta; 0; 0; H];

    rho = atmosfera(H);
    g = softPARAMS.g;

    update(ap, strain, strainp, strainpp, lambda);
    updateStrJac(ap);

    FAERO = [];
    MAERO = [];
    for i = 1:ap.NUMmembers
        if ~isempty(ap.members{i}(1).node1.aero)
            [FAEROm, MAEROm, ~] = aeroforceandmoment(ap.members{i}(1).strainm, ap.members{i}(1).strainpm, ap.members{i}(1).strainppm, ap.members{i}(1).lambdam, beta, betap, ap.members{i}, Vwind, rho, [delta 0 0]);
            FAERO = [FAERO; FAEROm]; %#ok<AGROW>
            MAERO = [MAERO; MAEROm]; %#ok<AGROW>
        else
            FAERO = [FAERO; zeros(ap.membSIZES(i)*9,1)]; %#ok<AGROW>
            MAERO = [MAERO; zeros(ap.membSIZES(i)*9,1)]; %#ok<AGROW>
        end
    end

    FPROP = ap.prop.getFPROP(ap, tracao, rho, beta(2));

    GX = g*sin(kinetic(2))*cos(theta);
    GY = -g*sin(theta);
    GZ = -g*cos(kinetic(2))*cos(theta);
    GRAVITY = [GX; GY; GZ];

    R_aero_force = ap.Jpb' * ap.B * FAERO;
    R_aero_moment = ap.Jthetab' * ap.B * MAERO;
    R_aero = R_aero_force + R_aero_moment;
    R_prop = ap.Jpb' * FPROP;
    R_grav = ap.Jhb' * ap.N * GRAVITY + ap.fus.N * GRAVITY;
    R_total = R_aero + R_prop + R_grav;

    [~, bp, ~, ~] = dinamicaflex(0, strain, strainp, strainpp, lambda, beta, betap, kinetic, ap, Vwind, tracao, [delta 0 0], 0);

    comp.theta = theta;
    comp.delta = delta;
    comp.tracao = tracao;
    comp.R_aero = R_aero;
    comp.R_aero_force = R_aero_force;
    comp.R_aero_moment = R_aero_moment;
    comp.R_prop = R_prop;
    comp.R_grav = R_grav;
    comp.R_total = R_total;
    comp.bp = bp;
    comp.FAERO = FAERO;
end

function report_force_node_breakdown(ap, FAERO, title_text)
    nNodes = ap.NUMele * 3;
    contrib = zeros(nNodes,1);
    for inode = 1:nNodes
        Fonly = zeros(size(FAERO));
        idx = (inode-1)*3 + (1:3);
        Fonly(idx) = FAERO(idx);
        R = ap.Jpb' * ap.B * Fonly;
        contrib(inode) = R(4);
    end

    [~, order] = sort(abs(contrib), 'descend');
    topn = min(8, nNodes);
    fprintf('\n%s\n', title_text);
    fprintf('Total eq4 force-path contribution: % .6e\n', sum(contrib));
    fprintf('Top %d node contributions (eq4):\n', topn);
    for k = 1:topn
        inode = order(k);
        ielem = ceil(inode/3);
        ilocal = inode - (ielem-1)*3;
        fprintf('  elem %d node %d -> % .6e\n', ielem, ilocal, contrib(inode));
    end
end

function report_components(comp)
    fprintf('theta=%.6f delta=%.6f tracao=%.6f\n', comp.theta, comp.delta, comp.tracao);
    fprintf('R_aero(2:4) = [% .6e % .6e % .6e]\n', comp.R_aero(2), comp.R_aero(3), comp.R_aero(4));
    fprintf('  from force  = [% .6e % .6e % .6e]\n', comp.R_aero_force(2), comp.R_aero_force(3), comp.R_aero_force(4));
    fprintf('  from moment = [% .6e % .6e % .6e]\n', comp.R_aero_moment(2), comp.R_aero_moment(3), comp.R_aero_moment(4));
    fprintf('R_prop(2:4) = [% .6e % .6e % .6e]\n', comp.R_prop(2), comp.R_prop(3), comp.R_prop(4));
    fprintf('R_grav(2:4) = [% .6e % .6e % .6e]\n', comp.R_grav(2), comp.R_grav(3), comp.R_grav(4));
    fprintf('R_total(2:4)= [% .6e % .6e % .6e]\n', comp.R_total(2), comp.R_total(3), comp.R_total(4));
    fprintf('bp(2:4)     = [% .6e % .6e % .6e]\n', comp.bp(2), comp.bp(3), comp.bp(4));
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
