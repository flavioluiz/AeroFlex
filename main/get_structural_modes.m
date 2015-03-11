% this function finds the structural dynamic modes of the flexible airplane
% (without taking into account rigid body degrees of freedom!)
%  in the future this function could be become a class, with several
%  methods for dealing with the results
function [eigen_vectors, eigen_values] = get_structural_modes(airplane)
    ap = airplane;
    
    KFF = ap.KG;
    updateStrJac(ap);
    MFF = ap.Jhep'*ap.Me*ap.Jhep;
    MFB = ap.Jhep'*ap.Me*ap.Jhb;
    MBF = ap.Jhb'*ap.Me*ap.Jhep;
    MRB = ap.fus.MRB;
    MBB = ap.Jhb'*ap.Me*ap.Jhb + MRB;
            
            
    [autovet, autovalD] = eig([KFF],[MFF]);
    %[autovet autovalD] = eig([MFF MFB; MBF MBB],[KFF zeros(size(MFB)); zeros(size([MBF MBB]))]);
    autoval = diag(autovalD);
    [valores, indexsort] = sort(autoval);
    
    sqrt(valores(1:6))/(2*pi)
    
    eigen_vectors = autovet(:, indexsort);
    eigen_values = autoval(indexsort);
end
