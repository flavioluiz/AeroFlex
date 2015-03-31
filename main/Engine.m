%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2011- Flávio Luiz C. Ribeiro
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Engine < handle
    properties
        numPI; %numero da manete que controla o motor
        numMEMB; %numero do membro onde está posicionado
        numELEM; %numero do elemento onde está posicionado
        numNODE; %numero do nó onde está posicionado
        params; %params.Fmax, params.V0, params.rho0, params.nrho, params.nv
        NODEpos; %número do nó onde se encontra o elemento
    end
    methods
        function prop = Engine(numPI, pos, params)
            prop.numPI = numPI;
            prop.numMEMB = pos(1);
            prop.numELEM = pos(2);
            prop.numNODE = pos(3);
            prop.params = params;
        end
        function FPROP = getFPROP(prop,ap, manete, rho, U)
            FPROP = zeros(ap.NUMele*9,1);
            if isempty(prop.params)
                Fmax = 0; rho0 = 1; V0 = 1; nrho = 0; nv = 0;
                alphaf = 0; betaf = 0;
                
            else
                Fmax = prop.params.Fmax; rho0 = prop.params.rho0;
                V0 = prop.params.V0; nrho = prop.params.nrho;
                nv = prop.params.nv; alphaf = prop.params.alphaf;
                betaf = prop.params.betaf;
            end
            for i = 1:size(prop,2)
                 Cwb = [1 0 0;0 cos(alphaf) -sin(alphaf); 0 sin(alphaf) cos(alphaf)]*...
                     [cos(betaf) -sin(betaf) 0; sin(betaf) cos(betaf) 0; 0 0 1];
                 CBb = eye(3)*Cwb;
                FPROP((1+3*(prop(i).NODEpos-1)):(3+3*(prop(i).NODEpos-1))) = CBb*[0; manete(prop(i).numPI)*Fmax*(U/V0)^nv*(rho/rho0)^nrho;0];
            end
            %FPROP(3+9*(6-1)) = -150;            
        end
    end
end

%% adicionar: follower loads!
%                 if ~isempty(prop(i).params)
%                     if prop(i).params.follower
%                         switch prop(i).numNODE
%                             case 1
%                                 hprop = ap.membros{prop(i).numMEMB}(prop(i).numELEM).node1.h;
%                             case 2
%                                 hprop = ap.membros{prop(i).numMEMB}(prop(i).numELEM).node2.h;
%                             case 3
%                                 hprop = ap.membros{prop(i).numMEMB}(prop(i).numELEM).node3.h;
%                         end
%                         Cwb = [hprop(4:6) hprop(7:9) hprop (10:12)];
%                     else
%                         Cwb = eye(3);
%                     end
%                 else
%                     Cwb = eye(3);
%                 end
