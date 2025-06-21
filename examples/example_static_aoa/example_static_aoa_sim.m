%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      AeroFlex - Example - Static Equilibrium at Angle of Attack     %
% - finds equilibrium (aerodynamics + gravity) at a given AoA         %
% - extracts and plots spanwise vertical deflection and torsion       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function example_static_aoa_sim
    clc;
    clear all;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AeroFlex Main Folder:                              %
    addpath('../../main');    
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
    softPARAMS.numITER = 20; % number of iterations for equilibrium determination
    softPARAMS.modAED = 1; % AERODYNAMIC MODEL: 
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
    
    % Define Angle of Attack for this simulation
    AoA_deg_val = 1.0; % MODIFIED - Moved AoA definition here
    
    %%%%%%%%%%%% STRUCTURE INITIALIZATION %%%%%%%%%%%%%%%%%
    numele = 20; %number of elements
    damping = 0.0001; %damping coefficient (damping proportional to rigidity matrix)
    ap = load_structure(numele, damping, AoA_deg_val); % MODIFIED - Pass AoA_deg_val
                                        %airplane object with numele elements
                                        % check the function loadstruct
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    %%%%%%%%%%%% FINDS EQUILIBRIUM CONDITION %%%%%%%%%%%%%%
    % FLIGHT CONDITIONS
    altitude = 19931.7; % meters (from example2, for air density)
    V = 0;           % rigid body speed m/s (cantilevered wing)
    throttle = 0;
    deltaflap = 0; 
    Vwind = 30; % wind speed m/s
    
    fprintf('Calculating equilibrium at Vwind = %.1f m/s and AoA = %.1f deg...\\n', Vwind, AoA_deg_val); % MODIFIED - Use AoA_deg_val

    [rb_eq, strain_eq] = trimairplane(ap,V,altitude,Vwind,throttle,deltaflap);        
    
    fprintf('Equilibrium calculation complete.\\n');

    figure('color','w','name','Equilibrium Deformation');
    % plot structure without deformation:
    % update(ap,strain_eq*0,zeros(size(strain_eq)),zeros(size(strain_eq)),zeros(sum(ap.membNAEDtotal),1)); % Optional: undeformed
    % plotairplane3d(ap); 
    % plot deformed structure (equilibrium condition):
    update(ap,strain_eq,zeros(size(strain_eq)),zeros(size(strain_eq)),zeros(sum(ap.membNAEDtotal),1));
    plotairplane3d(ap); 
    view(30,45); axis equal; colormap winter;
    title(sprintf('Deformed Wing at Vwind = %.1f m/s, AoA = %.1f deg', Vwind, AoA_deg_val)); % MODIFIED - Use AoA_deg_val

    % Extract and Plot Spanwise Deflection and Torsion
    num_elements_on_member = ap.membSIZES(1); % Assuming first member is the wing
    span_loc_X = zeros(num_elements_on_member + 1, 1); 
    deflection_Z = zeros(num_elements_on_member + 1, 1); 
    elastic_torsion_deg = zeros(num_elements_on_member + 1, 1);

    phi_initial_rad = AoA_deg_val * pi/180; % Initial twist due to AoA setting

    % Get root position and torsion (node1 of first element)
    h_root = ap.members{1}(1).node1.h;
    span_loc_X(1) = h_root(1); 
    deflection_Z(1) = h_root(3); 
    ey_local_root = h_root(7:9); % Chordwise versor
    current_total_angle_root = atan2(ey_local_root(3), ey_local_root(2));
    elastic_torsion_rad_root = current_total_angle_root - phi_initial_rad;
    elastic_torsion_deg(1) = atan2(sin(elastic_torsion_rad_root), cos(elastic_torsion_rad_root)) * 180/pi;
    
    % Get positions and torsion of node3 of each element (element tips)
    for i = 1:num_elements_on_member
        h_node3 = ap.members{1}(i).node3.h;
        span_loc_X(i+1) = h_node3(1); 
        deflection_Z(i+1) = h_node3(3); 
        ey_local_node3 = h_node3(7:9); % Chordwise versor
        current_total_angle_node3 = atan2(ey_local_node3(3), ey_local_node3(2));
        elastic_torsion_rad_node3 = current_total_angle_node3 - phi_initial_rad;
        elastic_torsion_deg(i+1) = atan2(sin(elastic_torsion_rad_node3), cos(elastic_torsion_rad_node3)) * 180/pi;
    end

    figure('color','w','name','Spanwise Vertical Deflection');
    plot(span_loc_X, deflection_Z, '-o');
    xlabel('Spanwise Location (X global coordinate, m)');
    ylabel('Vertical Deflection (Z global coordinate, m)');
    title(sprintf('Spanwise Deflection at Vwind = %.1f m/s, AoA = %.1f deg', Vwind, AoA_deg_val)); % MODIFIED - Use AoA_deg_val
    grid on;
    %axis equal; % User commented out

    figure('color','w','name','Spanwise Elastic Torsion'); % NEW PLOT
    plot(span_loc_X, elastic_torsion_deg, '-o');
    xlabel('Spanwise Location (X global coordinate, m)');
    ylabel('Elastic Torsion Angle (degrees)');
    title(sprintf('Spanwise Elastic Torsion at Vwind = %.1f m/s, AoA = %.1f deg', Vwind, AoA_deg_val));
    grid on;

    fprintf('Maximum tip deflection (Z): %.4f m\\n', deflection_Z(end));
    fprintf('Maximum tip elastic torsion: %.4f deg\\n', elastic_torsion_deg(end));
    fprintf('Program finished.\\n');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      REMOVED FLUTTER SPEED CALCULATION             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      REMOVED NONLINEAR SIMULATION                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

end

function ap = load_structure(numele, damp_ratio, AoA_deg_val_input) % MODIFIED - Accept AoA
    % member initialization
    % Angle of Attack is now passed from main function
    flexible_member = create_flexible_member(numele, damp_ratio, AoA_deg_val_input); % MODIFIED - Pass AoA
    
    %set member origin node position and orientation:
    % Default: wing span along X, chord in XY or XZ, airflow along Y.
    % h0 = [x_root, y_root, z_root, exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz]
    % This sets the member's root frame. Local element rotations (twist, sweep, dihedral)
    % are then applied relative to this.
    flexible_member(1).seth0([0 0 0  1 0 0  0 1 0  0 0 1]'); % Root at origin, member x along global X.
    update(flexible_member); % initialize displacements for each member node
    fus = []; % no fuselage
    motor1 = []; % no engines
    ap = Airplane({flexible_member}, fus, [motor1]);
    
    % Store AoA for plotting/reference, if not already in a standard place
    % This is a bit of a hack, ideally the object structure supports this.
    % We are adding it to the geometry struct of the first element's first node.
    % This assumes all elements of this member will have this initial twist.
    %if ~isempty(ap.members{1}) && ~isempty(ap.members{1}(1).node1.geometry)
    %    ap.members{1}(1).geometry.AoA_deg = AoA_deg;
    %else % Fallback if geometry struct is not initialized as expected
    %    ap.members{1}(1).node1.geometry = struct('AoA_deg', AoA_deg);
    %end

end

function flexible_member = create_flexible_member(num_elements,damp_ratio, AoA_deg) % MODIFIED - Accept AoA_deg
    

    % beam length
    Length = 16;
    
    % sectional rigidity matrix
    K11 = 1e10; %EA
    K22 = 1e4; %GJ
    K33 = 2e4; %flat bend: EI
    K44 = 4e6; %chord bend: EI
    KG = diag([K11 K22 K33 K44]);
    
    % sectional damping matrix
    CG = damp_ratio*KG;
    
    % aerodynamic data
    c = 1; % chord
    aeroparams.b = c/2; %semi-chord
    aeroparams.N = 4; %number of lag states (Peter's Unsteady model)
    aeroparams.a = 0.0; % position of aerodynamic center relative to elastic axis
                        % relative to elastic axis (in terms of semi-chord)
    aeroparams.alpha0 = 0*pi/180; % alpha_L0 (zero-lift AoA in radians)
    aeroparams.clalpha = 2*pi; % cl_alpha  lift coeff/rad
    aeroparams.cm0 = 0;        % cm_ac     moment coeff about AC
    aeroparams.cd0 = 0.00;     % cd_0
    
    % aerodynamic data for flap/aileron, if exists
    aeroparams.ndelta = 0;   % Identification of the flap (1,2,3,...)
    aeroparams.cldelta = 0;  % cl_delta
    aeroparams.cmdelta = 0;  % cm_delta
            
    % cg position, mass and inertia data
    
    pos_cg = [0 0 0]; % position of section gravity center
                        % relative to elastic axis
    geometry.a = -0.5;
    geometry.b = c/2;    
    I22 = 0.0;
    I33 = 0.1;
    I11 = 0.1;
    mcs = 0.75; %mass per unit length (kg/m)
    Inertia = diag([I11 I22 I33]);
    
    % rotation of initial node with respect to body frame (allow
    % creating uniform beams with a twist, dihedral or sweep angle)
    rot0.dihedral = 0; %angles in RADIANS
    rot0.sweep = 0;
    rot0.twist = -AoA_deg * pi/180;
    
    % the following function creates a uniform structure automatically; if
    % you need a more complicate wing (with non-uniform parameters, check
    % how the following function creates the structure. you should modify
    % this function to define the correct parameters for each structural
    % node)
    flexible_member = create_uniform_structure(pos_cg, rot0, Length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements);    
    load data_cl_alpha.mat
    for i=1:num_elements
        flexible_member(i).node1.aero.clalpha = interp1(data_cl_alpha(:,1),data_cl_alpha(:,2),flexible_member(i).node1.s );
        flexible_member(i).node2.aero.clalpha = interp1(data_cl_alpha(:,1),data_cl_alpha(:,2),flexible_member(i).node2.s );
        flexible_member(i).node3.aero.clalpha = interp1(data_cl_alpha(:,1),data_cl_alpha(:,2),flexible_member(i).node3.s );
    end
end
