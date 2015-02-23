%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2011- Flávio Luiz C. Ribeiro
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cinematica (obtencao de h)
% obtem, através de expressão analítica, exponencial matricial
% e^((s-s0)KEps)
function [expG dexpGdEps] = getexpG(strain, ds)
     KEps = getKEps(strain);
     KEps2 = getKEps2(strain);
     KEps3 = getKEps3(strain);
     ex = strain(1); kx = strain(2);
     ky = strain(3); kz = strain(4);
     
     lambda2 = kx^2 + ky^2 + kz^2;
     lambda = sqrt(lambda2);
     %hnum = expm(KEps*ds)
     if lambda2 > 1e-10
         a = (ds/lambda2 - sin(lambda*ds)/lambda^3);
         b = (1-cos(lambda*ds))/lambda2;
         c = ds;
     else
         a = ds^3/6;
         b = ds^2/2;
         c = ds;
     end
     %expGast = a*KEps^3 + b*KEps^2 + c*KEps + eye(4);
     
     expGast = a*KEps3 + b*KEps2 + c*KEps + eye(4);
     Th = zeros(12);  Th(1,1) = 1;  Th(2,5) = 1;  Th(3,9) = 1; Th(4,2) = 1;
     Th(5,6) = 1;  Th(6,10) = 1;  Th(7,3) = 1;  Th(8,7) = 1;  Th(9,11) = 1;
     Th(10,4) = 1;   Th(11,8) = 1;  Th(12,12) = 1;
     expG = Th*makediag(expGast)*Th';
     
     %calculo da derivada parcial d EXP[G] /d epsilon
     if lambda2>1e-10
         dbdEps = (lambda*ds*sin(lambda*ds)-2*(1-cos(lambda*ds)))/lambda^4*[0*eye(4);kx*eye(4);ky*eye(4); kz*eye(4)];
         dadEps = (3*sin(lambda*ds)/lambda^5 - (2*ds+ds*cos(lambda*ds))/lambda^4)*[0*eye(4);kx*eye(4);ky*eye(4);kz*eye(4)];
     else
         dbdEps = zeros(16,4);
         dadEps = zeros(16,4);
     end
     dKEpsdex = zeros(4); dKEpsdex(1,2) = 1;
     dKEpsdkx = zeros(4); dKEpsdkx(3,4) = 1; dKEpsdkx(4,3) = -1;
     dKEpsdky = zeros(4); dKEpsdky(4,2) = 1; dKEpsdky(2,4) = -1;
     dKEpsdkz = zeros(4); dKEpsdkz(2,3) = 1; dKEpsdkz(3,2) = -1; %aparentemente erraram aqui!!! utilizando dKEpsdkz(3,2) = 0, os resultados de frequencia natural ficam exatamente iguais aos das teses
     dKEpsdEps = [dKEpsdex;dKEpsdkx;dKEpsdky;dKEpsdkz];
     KEpsident = [KEps zeros(4) zeros(4) zeros(4); zeros(4) KEps zeros(4) zeros(4); zeros(4) zeros(4) KEps zeros(4);zeros(4) zeros(4) zeros(4) KEps];
     KEps2ident= [KEps2 zeros(4) zeros(4) zeros(4); zeros(4) KEps2 zeros(4) zeros(4); zeros(4) zeros(4) KEps2 zeros(4);zeros(4) zeros(4) zeros(4) KEps2];
     dexpGdEps = dadEps * KEps3 + a*(dKEpsdEps *KEps2 + KEpsident *dKEpsdEps*KEps + KEps2ident*dKEpsdEps) + ...
         dbdEps * KEps2 + b*(dKEpsdEps * KEps + KEpsident *dKEpsdEps) + ... 
         c*dKEpsdEps;
     
     dexpGdEps = [[Th]* makediag(dexpGdEps(1:4,1:4)) *[Th]', [Th]* makediag(dexpGdEps(5:8,1:4)) *[Th]',[Th]* makediag(dexpGdEps(9:12,1:4)) *[Th]',[Th]* makediag(dexpGdEps(13:16,1:4)) *[Th]'];
     
end


%obtem a matriz de deformacoes KEps (página 28 do schearer)
function [KEps] = getKEps(strain)
    ex = strain(1); kx = strain(2);
    ky = strain(3); kz = strain(4);
    exb = (1+ex);            
    KEps = [0, exb, 0, 0;
           0, 0, kz, -ky;
            0,-kz,0,kx;
           0, ky,-kx,0];
end

%obtem a matriz de deformacoes KEps^2
function [KEps2] = getKEps2(strain)
    ex = strain(1); kx = strain(2);
    ky = strain(3); kz = strain(4);
    exb = (1+ex);
    lambda2 = (kx^2+ky^2+kz^2);
    KEps2 = [0, 0, exb*kz, -exb*ky;
           0, -(lambda2-kx^2), ky*kx, kz*kx;
            0,ky*kx,-(lambda2-ky^2),kz*ky;
           0, kz*kx,kz*ky,-(lambda2-kz^2)];
end

%obtem a matriz de deformacoes KEps^3
function [KEps3] = getKEps3(strain)
    ex = strain(1); kx = strain(2);
    ky = strain(3); kz = strain(4);
    exb = (1+ex);
    lambda2 = (kx^2+ky^2+kz^2);
    KEps3 = [0, -exb*(kz^2+ky^2), exb*ky*kx, exb*kz*kx;
           0, 0, -kz*lambda2, ky*lambda2;
            0,kz*lambda2,0,-kx*lambda2;
           0, -ky*lambda2,kx*lambda2,0];
end

function X = makediag(matrix)
    X = [matrix zeros(4) zeros(4); zeros(4) matrix zeros(4); zeros(4) zeros(4) matrix];
end