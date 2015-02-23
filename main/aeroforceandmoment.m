%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2011- Flávio Luiz C. Ribeiro 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cálculo do momento e sustentação aerodinâmicos considerando modelo
%estacionário
function [F M lambdap] = aeroforceandmoment(strain, strainp, strainpp,lambda,beta,betap,membro,Vwindy,rho,deltaflap)
    Vwind = [0;Vwindy;0];
    membersize = size(membro,2);
    pp = membro(1).memberJpep*strainp+ membro(1).memberJpb*beta;
    ppp = membro(1).memberJpep*strainpp + membro(1).memberJpb*betap;
    thetap = membro(1).memberJthetaep*strainp+membro(1).memberJthetab*beta;
    thetapp = membro(1).memberJthetaep*strainpp+membro(1).memberJthetab*betap;
    tao = 15;
    %checar matrizes Fl0 e Ml0
    for i = 1:membersize
        N = membro(i).node1.aero.N;
        lambda0 = 0.5*lambda((1+(i-1)*N*3):(N+(i-1)*N*3))*membro(i).node1.aero.B;
        b = membro(i).node1.aero.b;
        d = (membro(i).node1.aero.a)*b;
     %   fcorr = (1-exp((membro(i).node1.s-1)*tao));
     %   fcorr = 1;
        [Forca Momento flambda] = getforceandmoment(membro(i).node1.h, pp((1+(i-1)*9):(3+(i-1)*9)),thetap((1+(i-1)*9):(3+(i-1)*9)), ppp((1+(i-1)*9):(3+(i-1)*9)), thetapp((1+(i-1)*9):(3+(i-1)*9)), lambda0, Vwind,b,rho,d,membro(i).node1.aero,lambda((1+(i-1)*N*3):(N+(i-1)*N*3)),deltaflap);
        F((1+(i-1)*9):(3+(i-1)*9),1) = Forca;
        M((1+(i-1)*9):(3+(i-1)*9),1) = Momento;
        lambdap((1+(i-1)*3*N):(N+(i-1)*3*N),1) = flambda;
        
        N = membro(i).node2.aero.N;
        lambda0 = 0.5*lambda((N+1+(i-1)*N*3):(2*N+(i-1)*N*3))*membro(i).node2.aero.B;
        b = membro(i).node2.aero.b;
        d = (membro(i).node2.aero.a)*b;
      %  fcorr = (1-exp((membro(i).node2.s-1)*tao));
     %   fcorr = 1;
        [Forca Momento flambda] = getforceandmoment(membro(i).node2.h, pp((4+(i-1)*9):(6+(i-1)*9)),thetap((4+(i-1)*9):(6+(i-1)*9)), ppp((4+(i-1)*9):(6+(i-1)*9)), thetapp((4+(i-1)*9):(6+(i-1)*9)), lambda0, Vwind,b,rho,d,membro(i).node2.aero,lambda((N+1+(i-1)*N*3):(2*N+(i-1)*N*3)),deltaflap);
        F((4+(i-1)*9):(6+(i-1)*9),1) = Forca;
        M((4+(i-1)*9):(6+(i-1)*9),1) = Momento;
        lambdap((N+1+(i-1)*3*N):(2*N+(i-1)*3*N),1) = flambda;
        
        N = membro(i).node3.aero.N;
        lambda0 = 0.5*lambda((2*N+1+(i-1)*N*3):(3*N+(i-1)*N*3))*membro(i).node3.aero.B;
        b = membro(i).node3.aero.b;
        d = (membro(i).node3.aero.a)*b;
      %  fcorr = (1-exp((membro(i).node3.s-1)*tao));
     %   fcorr = 1;                                 (              h,  dotp,                       dottheta,                        ddotp,                         ddottheta,                      lambda0,  Vwind, b, rho, d,aero,             lambda,                                      delta)
        [Forca Momento flambda] = getforceandmoment(membro(i).node3.h, pp((7+(i-1)*9):(9+(i-1)*9)),thetap((7+(i-1)*9):(9+(i-1)*9)), ppp((7+(i-1)*9):(9+(i-1)*9)), thetapp((7+(i-1)*9):(9+(i-1)*9)), lambda0, Vwind,b,rho,d,membro(i).node3.aero,lambda((2*N+1+(i-1)*N*3):(3*N+(i-1)*N*3)),deltaflap);
        F((7+(i-1)*9):(9+(i-1)*9),1) = Forca;
        M((7+(i-1)*9):(9+(i-1)*9),1) = Momento;
        lambdap((2*N+1+(i-1)*3*N):(3*N+(i-1)*3*N),1) = flambda;
    end

    
end

function [Forca Momento flambda] = getforceandmoment(h, dotp, dottheta, ddotp, ddottheta, lambda0, Vwind, b, rho, d,aero,lambda,deltaflap)
    global softPARAMS;
    e1 = [1 0 0]';
    e2 = [0 1 0]';
    e3 = [0 0 1]';
    %verificar matrizes de rotação aerodinâmicas
    %verificar ponto de aplicação de momento e se está coerente!!
    if aero.ndelta>0
        delta = deltaflap(aero.ndelta);
    else
        delta = 0;
    end
    % os parametros ci e gi dependem da geometria do FLAP 
    CBw = [h(4:6) h(7:9) h(10:12)]; % retirar [wx,wy,wz] de h(epsilon)
    Cwa0 = eye(3); % CONFIRMAR !!!!! ???? % matriz de transformação entre o sistema de coordenadas local e o sistema de sustentação nula (relacionado por alfa0).
    alpha0 = aero.alpha0;
    Cwa0 = [1 0 0;
             0 cos(alpha0) sin(alpha0);
             0 -sin(alpha0) cos(alpha0)];  
    
    CBa0 = CBw*Cwa0;
    doty = e2'*CBa0'*(dotp + Vwind);
    dotz = e3'*CBa0'*(dotp + Vwind);
    alphat = atan(-dotz/doty);
    Ca0a1 = [1 0 0;
             0 cos(alphat) sin(alphat);
             0 -sin(alphat) cos(alphat)];  %checar, tive que mudar essa matriz para que os resultados ficassem coerentes (Fy e Fz) com valores esperados.

    CBA = CBa0*Ca0a1;
    
    ddotz = e3'*CBa0'*ddotp;
    dotalpha = e1'*CBa0'*dottheta;
    ddotalpha = e1'*CBa0'*ddottheta;
    
    % forças e momentos para placa plana
    %L = pi*rho*b^2 *(-ddotz + doty*dotalpha + d*ddotalpha) + 2*pi *rho *b*doty^2*(-dotz/doty +(b/2-d)*dotalpha/doty-lambda0/doty);
    %M = 2*pi*rho*b^2 * (-0.5*doty*dotz-0.5*d*doty*dotalpha - 0.5 *doty*lambda0 - b^2 *ddotalpha /16);
    %D = -2*pi*rho*b * (dotz^2 + d^2*dotalpha^2 + lambda0^2 + 2*d*dotalpha*dotz + 2*lambda0*dotz +2 *d*dotalpha*lambda0);
    % forcas para perfil arqueado
    clalpha = aero.clalpha;
    cldelta = aero.cldelta;
    cd0 = aero.cd0;
    cm0 = aero.cm0;
    cmdelta = aero.cmdelta;
    
    switch softPARAMS.modAED
        case 0 %AED ESTACIONARIA - apenas ângulo de ataque
            L = clalpha *rho *b*doty^2*(-dotz/doty) + rho*b*doty^2*cldelta*delta;
            M = 2*rho*b^2*doty^2*(cm0 + cmdelta*delta);
            D = -rho*b*doty^2*cd0;
        case 1 %AED ESTACIONARIA  - eqs 5 e 6 do haddadpour
            %V = sqrt(doty^2+dotz^2);
            V = doty;
            L = clalpha *rho *b*V^2*(-dotz/V +(b/2-d)*dotalpha/V) + rho*b*doty^2*cldelta*delta;
            M = -0.5*pi*rho*b^3*V*dotalpha + 2*rho*b^2*doty^2*(cm0 + cmdelta*delta);
            D = -rho*b*doty^2*cd0;
        case 2 % AED ESTACIONARIA - eqs 7 e 8 do haddadpour mas com c(k) = 1
            L = pi*rho*b^2 *(-ddotz + doty*dotalpha - d*ddotalpha) + clalpha *rho *b*doty^2*(-dotz/doty +(b/2-d)*dotalpha/doty) + rho*b*doty^2*cldelta*delta;
            M = -pi*rho*b^3*(-0.5*ddotz + doty*dotalpha +(0.125*b - 0.5*d)*ddotalpha) + 2*rho*b^2*doty^2*(cm0 + cmdelta*delta);
            D = -rho*b*doty^2*cd0;
        case 3 % AED NÃO-ESTACIONÁRIA - Peters
            L = pi*rho*b^2 *(-ddotz + doty*dotalpha - d*ddotalpha) + clalpha *rho *b*doty^2*(-dotz/doty +(b/2-d)*dotalpha/doty-lambda0/doty) + rho*b*doty^2*cldelta*delta;
            M = -pi*rho*b^3*(-0.5*ddotz + doty*dotalpha +(0.125*b - 0.5*d)*ddotalpha) + 2*rho*b^2*doty^2*(cm0 + cmdelta*delta);
            D = -rho*b*doty^2*cd0;
    end
        
    % AED FULL DO WEIHUA    
    %L = pi*rho*b^2 *(-ddotz + doty*dotalpha - d*ddotalpha) + clalpha *rho *b*doty^2*(-dotz/doty +(b/2-d)*dotalpha/doty-lambda0/doty) + rho*b*doty^2*cldelta*delta;
    %M = 2*pi*rho*b^2*(-0.5*dotz*doty - 0.5*d*doty*dotalpha -0.5*doty*lambda0 - 0.125*b^2*ddotalpha) + 2*rho*b^2*doty^2*(cm0 + cmdelta*delta);
    %D = -rho*b*doty^2*cd0*0; 
       
    Forca = CBA*[0;D;L];
    Mx = M + (d+0.5*b)*(L*cos(alphat+alpha0)-D*sin(alphat+alpha0)); %estou so considerando momento devido às forças verticais! se houver um offset vertical na posiçao do CA em relação ao eixo elástico, o momento não estará totalmente modelado
    %Mx = M + (d+0.5*b)*(L); %estou so considerando momento devido às forças verticais! se houver um offset vertical na posiçao do CA em relação ao eixo elástico, o momento não estará totalmente modelado
    Momento = CBA*[Mx;0;0];
    
    Vwind = doty;
    flambda = aero.E1*lambda'*norm(Vwind)/b + aero.E2*(-ddotz/b) + aero.E3*ddotalpha + aero.E4*dotalpha*norm(Vwind)/b;    
    %flambda = aero.E1*lambda' + aero.E2*alphat + aero.E3*dotalpha;    
end