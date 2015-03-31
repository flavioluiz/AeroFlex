% this class finds the structural dynamic modes of the flexible airplane
% (without taking into account rigid body degrees of freedom!)
classdef StructuralModes < handle
    properties
        frequencies;
        eigen_vectors;
        num_of_modes;
    end
    methods
        function st_modes = StructuralModes(airplane)
            ap = airplane;

            KFF = ap.KG;
            updateStrJac(ap);
            MFF = ap.Jhep'*ap.Me*ap.Jhep;
            %MFB = ap.Jhep'*ap.Me*ap.Jhb;
            %MBF = ap.Jhb'*ap.Me*ap.Jhep;
            %MRB = ap.fus.MRB;
            %MBB = ap.Jhb'*ap.Me*ap.Jhb + MRB;

            [autovet, autovalD] = eig(KFF,MFF);
            %[autovet autovalD] = eig([MFF MFB; MBF MBB],[KFF zeros(size(MFB)); zeros(size([MBF MBB]))]);
            autoval = diag(autovalD);
            [valores, indexsort] = sort(autoval);

            st_modes.eigen_vectors = autovet(:, indexsort);
            st_modes.frequencies = sqrt(autoval(indexsort))/2/pi;
            st_modes.num_of_modes = length(st_modes.frequencies);
        end
        function plot_mode(st_modes, ap, num_mode)
            strain = st_modes.eigen_vectors(:,num_mode);
            update(ap,strain/25,zeros(size(strain)),zeros(size(strain)),zeros(sum(ap.membNAEDtotal),1));
            plotairplane3d(ap);
        end
        
    end
end