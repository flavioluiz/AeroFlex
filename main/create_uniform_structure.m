% creates a structure with constant parameters (mass / unit length,
% sectional Inertia, aerodynamic parameters, etc.)

function membro = create_uniform_structure(pos_cg, initial_node_rot, length, Inertia, mcs, KG, CG, aeroparams, geometry, num_elements)
    ds = length/num_elements;
    
    rot = initial_node_rot;
    rigidunit.m = 0; rigidunit.cg = [0,0,0]; rigidunit.I = zeros(3,3); % no rigid units
    
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