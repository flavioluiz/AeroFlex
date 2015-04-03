%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2011- Flávio Luiz C. Ribeiro
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xp, bp, lambdap, kineticp] = dinamicaflex(t,strain, strainp, strainpp,lambda,beta,betap,kinetic,ap,V, manete, deltaflap, FLAG)
    global softPARAMS;
    H = kinetic(4);
    
    if softPARAMS.isIS
        g = softPARAMS.g*(softPARAMS.isGRAV==1);
        rho = atmosfera(H);
    else
        g = softPARAMS.g*1/0.3048*(softPARAMS.isGRAV==1);
        rho = atmosfera(H*0.3048)*0.00194032033; %slug/ft^3
    end
    
    switch FLAG
        case 0
            isEQ = 1;
        case 1
            isEQ = 0;
            isIMPLICIT = 1;
        case 2
            isEQ = 0;
            isIMPLICIT = 0;
    end
    
    % define dados strain, strainp, strainpp, lambda dos membros
    % atualiza deflexoes estruturais:
    
    update(ap,strain,strainp,strainpp,lambda);

    
    %atualização dos jacobianos (nessa parte, pode ser adicionada uma FLAG,
    %mantendo os jacobianos constantes se desejado)
    if softPARAMS.updateStrJac ==2 || (softPARAMS.updateStrJac == 1 && isEQ)
        updateStrJac(ap);
    end   
     

    KFF = ap.KG;
    
    if ~isEQ   %se estiver calculando o equilibrio, nenhuma dessas matrizes é necessária!!
        MFF = ap.Jhep'*ap.Me*ap.Jhep;
        CFF = ap.CG;%+Jhep'*Me*Jhepp;
        MFB = ap.Jhep'*ap.Me*ap.Jhb;
        MBF = ap.Jhb'*ap.Me*ap.Jhep;
        MRB = ap.fus.MRB;
        MBB = ap.Jhb'*ap.Me*ap.Jhb + MRB;
% %        eig(-MFF \ KFF)
%         %MFF epp + CFF ep + KFF e = 0
%         % EYE (e)p = EYE ep
%         matrizest=[MFF zeros(size(CFF)); zeros(size(MFF)) eye(size(MFF))] \ [-CFF -KFF; eye(size(MFF)) zeros(size(MFF))];
%         eig(matrizest)
        
        
        CBF = ap.Jhb'*ap.Me*ap.Jhep*0; %cuidado, na realidade: = Jhb'*Me*Jhepp;
        omega = matrixcross(beta(4:6))';
        nodesnumber = ap.NUMele*3;
        Hhb = OMEGA(omega,nodesnumber) * ap.Jhb;

        hp = ap.Jhep*strainp'; % hp em relação à origem do sistema do corpo! (não inercial, por isso não aparece beta)
        Jhbp = getJhbp(hp, nodesnumber);
        CFB = ap.Jhep'*ap.Me*Hhb + 2*ap.Jhep'*ap.Me*Jhbp;
        CRB = getCRB(ap.fus, omega);
        CBB = ap.Jhb'*ap.Me*Hhb + 2*ap.Jhb'*ap.Me*Jhbp + CRB; % faltam termos aqui ainda!!!!! ver eq. 2.31 do shearer, eq 2.33 do Weihua e 2.11 weihua
    end
  
    FAERO = [];   MAERO = [];   FLAMBDA = [];
    for i = 1:ap.NUMmembers
        if ~isempty(ap.members{i}(1).node1.aero)
            [FAEROm, MAEROm, FLAMBDAm] = aeroforceandmoment(ap.members{i}(1).strainm, ap.members{i}(1).strainpm, ap.members{i}(1).strainppm,ap.members{i}(1).lambdam,beta,betap,ap.members{i},V,rho,deltaflap);
            FAERO = [FAERO; FAEROm];
            MAERO = [MAERO; MAEROm];
            FLAMBDA = [FLAMBDA; FLAMBDAm];
        else
            FAERO = [FAERO; zeros(ap.membSIZES(i)*9,1)];
            MAERO = [MAERO; zeros(ap.membSIZES(i)*9,1)];
            FLAMBDA = [FLAMBDA; []];
        end
    end
    U = beta(2);

    if ~isempty(ap.prop)
        FPROP = ap.prop.getFPROP(ap, manete, rho, U);
    else
        FPROP = zeros(ap.NUMele*9,1);
    end

    %cálculo da força gravitacional
  %  kinetic = [0 0 0 0];
    theta = kinetic(1);
    phi = kinetic(2);
    psi = kinetic(3);
    H = kinetic(4);
    

   % theta = atan(beta(3)/beta(2));
    GX = g*sin(phi)*cos(theta);
    GY = -g*sin(theta);
    GZ = -g*cos(phi)*cos(theta);
    GRAVITY = [GX; GY; GZ];

    RAEROF = ap.Jpep'*ap.B*FAERO + ap.Jpep'*FPROP + ap.Jthetaep'*ap.B*MAERO+ap.Jhep'*ap.N*GRAVITY;
    RAEROB = ap.Jpb'*ap.B*FAERO + ap.Jpb'*FPROP + ap.Jthetab'*ap.B*MAERO+ap.Jhb'*ap.N*GRAVITY + ap.fus.N*GRAVITY;
    
    if isEQ        
        xp = - KFF * strain' +  RAEROF;                
        bp = + RAEROB.*softPARAMS.vecFREEDEG';
    else
        if isIMPLICIT
            xp = MFF\(-MFB*betap -CFB*beta -CFF*strainp' - KFF*strain' + RAEROF);
            bp = MBB\(-MBF*strainpp' - CBF*strainp' - CBB*beta + RAEROB).*softPARAMS.vecFREEDEG';
        else  % EXPLICITA
            if softPARAMS.vecFREEDEG == [0 0 0 0 0 0]
                xp = MFF \ (-CFF * strainp' - KFF * strain' + RAEROF);
                bp = zeros(6,1);                
            else
                xpbp = [MFF MFB; MBF MBB] \ ([-CFF, -CFB; -CBF, -CBB]*[strainp'; beta] -[KFF*strain';zeros(6,1)]+ [RAEROF;RAEROB]);
                xp = xpbp(1:ap.NUMele*4);
                bp = xpbp((ap.NUMele*4+1):(ap.NUMele*4+6)).*softPARAMS.vecFREEDEG';
            end
        end
    end
    lambdap = FLAMBDA;
       
    P = beta(5);
    Q = beta(4);
    R = -beta(6);
    thetap = Q*cos(phi)- R*sin(phi);
    phip = P + tan(theta)*(Q*sin(phi)+R*cos(phi));
    psip = (Q*sin(phi) + R*cos(phi))/cos(theta);
    U = beta(2);
    V = beta(1);
    W = -beta(3);
    Hp = U*sin(theta) - V*sin(phi)*cos(theta) - W*cos(phi)*cos(theta);
    
    kineticp = [thetap; phip; psip; Hp];
end
