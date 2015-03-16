Example 2 - Linear and nonlinear Aeroelasticity
==========

In this example, we use the same structure presented in example 1.
The only difference is that we include an aerodynamic model, allowing
aeroelastic studies.

This initialization is almost unchanged. First, the main folder is included to path:

    addpath('..\..\main');    
	
Then, the general parameters are defined:

    global softPARAMS;
    softPARAMS.isIS = 1; %is it International System?
    softPARAMS.isPINNED = 1; % PINNED RIGID BODY DOF
    softPARAMS.vecFREEDEG = [0 0 0 0 0 0]; % type 0 to remove a body state
                                           % degree of freedom [u v w p q r]
    softPARAMS.isGRAV = 1; % include gravity?
    softPARAMS.g = 9.8; % gravity in m/s^2    
    softPARAMS.isITER = 1; % iterative equilibrium determination?
    softPARAMS.numITER = 10; % number of iterations for equilibrium determination
    softPARAMS.modAED = 3; % AERODYNAMIC MODEL: 
                                    %0-Steady;
                                    %1-Quasi-steady;
                                    %2-Quasi-steady with added mass;
                                    %3-Unsteady(Peters);
    softPARAMS.updateStrJac = 1; % Structural Jacobians updates:
                                    % 0 - Never;
                                    % 1 - Only in equilib calculation;
                                    % 2 - Always
    softPARAMS.plota3d = 1; % 3d graphics plot while equilibrium is calculated

The only difference above is that now we chose the Peters aerodynamic model. This is
the best aerodynamic model implemented in AeroFlex, taking into account unsteady effect of aerodynamics.
(it is also considerably slower, since it adds new dynamic states relative to each structural node).

The, the airplane object is initialized by calling the function 'load_structure':

    numele = 3; %number of elements
    damping = 0.0001; %damping coefficient (damping proportional to rigidity matrix)
    ap = load_structure(numele,damping); % this creates a flexible
                                        %airplane object with numele elements
                                        % check the function loadstruct


The function 'load_structure' is exactly the same as that from example 1:
										
	function ap = load_structure(numele, damp_ratio)
		% member initialization
		flexible_member = create_flexible_member(numele,damp_ratio);
		
		%set member origin node position and orientation:
		flexible_member(1).seth0([0 -0.0 0 1 0 0 0 1 0 0 0 1]'); 
		update(flexible_member); % initialize displacements for each member node
		fus = []; % no fuselage
		motor1 = []; % no engines
		ap = airplane({flexible_member}, fus, [motor1]);
	end

The difference is on the function 'create_flexible_member', which now includes
the aerodynamic data:

	function flexible_member = create_flexible_member(num_elements,damp_ratio)  
		% beam length
		Length = 16;
		
		% sectional rigidity matrix
		K11 = 1e10; %EA
		K22 = 1e4; %GJ
		K33 = 2e4; %flat bend: EI
		K44 = 4e6; %chord bend: EI
		KG = diag([K11 K22 K33 K44]);
		
		% sectional damping matrix
		CG = damp_ratio*diag([K11 K22 K33 K44]);
		
		% aerodynamic data
		c = 1; % chord
		aeroparams.b = c/2; %semi-chord
		aeroparams.N = 4; %number of lag states (Peter's Unsteady model)
		aeroparams.a = 0.0; % position of aerodynamic center relative to elastic axis
							% relative to elastic axis (in terms of semi-chord)
		aeroparams.alpha0 = 0; % alpha_0 (in radians)
		aeroparams.clalpha = 2*pi; % cl_alpha  lift coeff/rad
		aeroparams.cm0 = 0;        % cm_0      moment coeff
		aeroparams.cd0 = 0.02;     % cd_0
		
		% aerodynamic data for flap/aileron, if exists
		aeroparams.ndelta = 0;   % Identification of the flap (1,2,3,...)
		aeroparams.cldelta = 0;  % cl_delta
		aeroparams.cmdelta = 0;  % cm_delta
				
		% cg position, mass and inertia data
		
		pos_cg = [0 0 0]; % position of section gravity center
							% relative to elastic axis
		geometry.a = 0.0;
		geometry.b = 0.5;    
		I22 = 0.0;
		I33 = 0.1;
		I11 = 0.1;
		mcs = 0.75; %mass per unit length (kg/m)
		Inertia = diag([I11 I22 I33]);
		
		% the following function creates a uniform structure automatically; if
		% you need a more complicate wing (with non-uniform parameters, check
		% how the following function creates the structure. you should modify
		% this function to define the correct parameters for each structural
		% node)
		flexible_member = create_uniform_structure(pos_cg, Length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements);        
	end


Notice that we now define a structure 'aeroparams', which includes all the aerodynamic
relevant that for each node (cl_alpha, alpha_0, cm_0, cd_0, etc.). In this example, the
wing is uniform. If you are working with a more complicated wing, you need all these parameters
for each node.


Now, the wing is defined! We can use several functions to study this system (like 'trimairplane'
to find the equilibrium, 'linearize' to study the linear behavior, 'simulate', 'flutter_speed', etc.)

Linear aeroelasticity
----------

32.56 m/s | 22.55 rad/s

Nonlinear aeroelasticity
---------

24.03 m/s | 12.17 rad/s
