function dump_trim_modes_legacy
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
    altitude = 19931.7;
    V = 10;
    Vwind = 0;

    [rb_eq, strain_eq] = trimairplane(ap, V, altitude, Vwind, 0, 0);
    thetaeq = rb_eq(1);
    deltaeq = rb_eq(2);
    throttleeq = rb_eq(3);
    betaeq = [0; V*cos(thetaeq); -V*sin(thetaeq); 0; 0; 0];
    keq = [thetaeq; 0; 0; altitude];

    [Alin, ~, Abody] = linearize(ap, strain_eq, betaeq, keq, throttleeq, deltaeq, Vwind);
    fullEig = eig(Alin);
    bodyEig = eig(Abody);

    st_modes = StructuralModes(ap);
    nf = min(6, numel(st_modes.frequencies));

    fprintf('RESULT_TRIM %.9f %.9f %.9f\n', thetaeq, deltaeq, throttleeq);
    fprintf('RESULT_STRUCT_FREQ');
    for i = 1:nf
        fprintf(' %.9f', st_modes.frequencies(i));
    end
    fprintf('\n');

    full_sorted = sort_complex(fullEig);
    body_sorted = sort_complex(bodyEig);
    nout = min(6, numel(full_sorted));
    fprintf('RESULT_FULL_EIG_REAL_IMAG');
    for i = 1:nout
        fprintf(' %.9f %.9f', real(full_sorted(i)), imag(full_sorted(i)));
    end
    fprintf('\n');
    noutb = min(6, numel(body_sorted));
    fprintf('RESULT_BODY_EIG_REAL_IMAG');
    for i = 1:noutb
        fprintf(' %.9f %.9f', real(body_sorted(i)), imag(body_sorted(i)));
    end
    fprintf('\n');
end

function vals = sort_complex(vals)
    key = [real(vals), imag(vals)];
    [~, idx] = sortrows(key, [1 2]);
    vals = vals(idx);
end

function ap = load_structure_legacy(numele, damp_ratio)
    [right_wing, left_wing] = create_flexible_member(numele, damp_ratio);
    right_wing(1).seth0([0 -.3 0 1 0 0 0 1 0 0 0 1]');
    left_wing(1).seth0([0 -.3 0 1 0 0 0 1 0 0 0 1]');
    update(right_wing);
    update(left_wing);
    fus = RigidFuselage(10, [0 0 0], zeros(3,3));
    engparams.Fmax = 1; engparams.V0 = 1; engparams.rho0 = 1;
    engparams.nv = -1; engparams.nrho = 0; engparams.alphaf = 0;
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
