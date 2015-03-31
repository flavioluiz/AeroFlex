%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2011- Flavio Luiz C. Ribeiro
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Node < handle
    properties
        mcs; %massa por unidade de comprimento
        mcg; %posicao do CG da sescao (rx,ry,rz)
        Inertia; %diadica de inercia da secao (3x3)
        Mnode; %matriz de inercia do no (por unidade de comprimento)
        Mrb; %matriz de inercia de unidade rigida ligada ao no
        h; %vetor com 9 elementos com posicao e orientacao do noh
        aero;
        s;
        geometry; % structure with two elements: b (semi-chord) and a
    end
    methods
        function nd = Node(mcs, mcg, Inertia,petersinp, rigid, s, geometry)  %inicializa o noh
            if nargin > 6
                nd.geometry = geometry;
            else 
                nd.geometry = [];
            end
            nd.mcs = mcs;
            nd.mcg = mcg;
            nd.Inertia = Inertia;
            nd.s = s;
            if isempty(petersinp)
                nd.aero = [];
            else
                N = petersinp.N;
                b = petersinp.b;
                a = petersinp.a;
                nd.geometry.a = a;
                nd.geometry.b = b;
                nd.aero = Peter(N,b,a);
                nd.aero.alpha0 = petersinp.alpha0;
                nd.aero.clalpha = petersinp.clalpha;
                nd.aero.cldelta = petersinp.cldelta; 
                nd.aero.cm0 = petersinp.cm0;    
                nd.aero.cmdelta = petersinp.cmdelta;
                nd.aero.cd0 = petersinp.cd0;
                nd.aero.ndelta = petersinp.ndelta;
            end
            calculateM(nd);
            calculateMrb(nd, rigid);
        end
        function calculateM(nd) %calcula matriz de massa Mi do noh
            mcs = nd.mcs*eye(3);
            mcg = nd.mcg;
            Inertia = nd.Inertia;
            rx = mcg(1)*eye(3);  ry = mcg(2)*eye(3);  rz = mcg(3)*eye(3);
            Ixx = Inertia(1,1)*eye(3);  Ixy = Inertia(1,2)*eye(3);  Ixz = Inertia(1,3)*eye(3);
            Iyx = Inertia(2,1)*eye(3);  Iyy = Inertia(2,2)*eye(3);  Iyz = Inertia(2,3)*eye(3);
            Izx = Inertia(3,1)*eye(3);  Izy = Inertia(3,2)*eye(3);  Izz = Inertia(3,3)*eye(3);
            nd.Mnode = [mcs, mcs*rx, mcs*ry, mcs*rz;
                mcs*rx, 1/2*(Iyy + Izz - Ixx), Ixy, Ixz;
                mcs*ry, Iyx, 1/2*(Ixx + Izz - Iyy), Iyz;
                mcs*rz, Izx, Izy, 1/2*(Ixx+Iyy-Izz)];
        end
        function calculateMrb(nd,rigid)
            m = rigid.m*eye(3);
            cg = rigid.cg;
            Inertia = rigid.I;
            rx = cg(1)*eye(3);  ry = cg(2)*eye(3);  rz = cg(3)*eye(3);
            Ixx = Inertia(1,1)*eye(3);  Ixy = Inertia(1,2)*eye(3);  Ixz = Inertia(1,3)*eye(3);
            Iyx = Inertia(2,1)*eye(3);  Iyy = Inertia(2,2)*eye(3);  Iyz = Inertia(2,3)*eye(3);
            Izx = Inertia(3,1)*eye(3);  Izy = Inertia(3,2)*eye(3);  Izz = Inertia(3,3)*eye(3);
%            nd.Mrb = [m, m*rx, m*ry, m*rz;
%                m*rx, Ixx + m*rx^2, Ixy + m*rx*ry, Ixz + m*rx*rz;
%                m*ry, Iyx + m*ry*rx, Iyy + m*ry^2, Iyz + m*ry*rz;
%                m*rz, Izx + m*rz*rx, Izy + m*rz*ry, Izz + m*rz^2];       
            nd.Mrb = [m, m*rx, m*ry, m*rz;
                m*rx, (Iyy+Izz-Ixx)/2 + m*rx^2, Ixy + m*rx*ry, Ixz + m*rx*rz;
                m*ry, Iyx + m*ry*rx, (Izz+Ixx-Iyy)/2 + m*ry^2, Iyz + m*ry*rz;
                m*rz, Izx + m*rz*rx, Izy + m*rz*ry, (Ixx+Iyy-Izz)/2 + m*rz^2];       

        end
    end
end
