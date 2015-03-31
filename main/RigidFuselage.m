%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2011- Flávio Luiz C. Ribeiro
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef RigidFuselage < handle
    properties
        m;
        pcm;
        I;
        MRB;
        ttpcm;
        tpcm;
        N;
    end
    methods
        function fus = RigidFuselage(m, pcm, inertia)
          fus.m = m;
          fus.pcm = pcm;
          fus.I = inertia;
          fus.ttpcm = matrixcross(pcm);
          fus.tpcm = fus.ttpcm';          
          fus.MRB = [eye(3)*m, m*fus.ttpcm;
                     m*fus.tpcm, inertia];                 
          %fus.N = [eye(3) * m;
          %         ttpcm*m];    -- the line below is equivalent to this
          fus.N = fus.MRB(:,1:3);
        end
        function CRB = getCRB(rigidfus, omega)
            CRB = [rigidfus.m*omega, rigidfus.m*omega*rigidfus.ttpcm;
                   rigidfus.m*rigidfus.tpcm*omega, rigidfus.I*omega];
        end
    end
end