function dinamicarigida(V,H,longfig,tSIM, deltav, deltaw, deltaalfa,tracaot,deltapt)
    equil = fsolve(@equilibra, [0 0 0], [],V,H);
    alpha = equil(1);
    tracaoeq = equil(2);
    deltapeq = equil(3);
    equilibrio(1) = V*cos(alpha);
    equilibrio(2) = -V*sin(alpha);
    equilibrio(3) = 0;
    equilibrio(4) = alpha;
    equilibrio(5) = H;
    tracaoeq
    deltapeq
    alpha
    A = lineariza(equilibrio,[tracaoeq deltapeq]);
    eig(A)
   % alpha = 5*pi/180;
   % alpha = 5*pi/180;
    xi(1) = V*cos(alpha+deltaalfa) + deltav;
    xi(2) = -V*sin(alpha+deltaalfa) +deltaw;
    xi(3) = 0;
    xi(4) = alpha+deltaalfa;
    xi(5) = H;
    
    %xp = dinamica(0,xi,[tracaoeq, deltapeq]);
    [t,x] = ode45(@(t,x)dinamica(t,x,[tracaot(t)+tracaoeq,deltapt(t)+deltapeq]), [0 tSIM], xi);
    %figure;
    %plot(t,x); legend('u','w','q','theta','H');
    Alpha = -atan(x(:,2)./x(:,1));
    Veloc = sqrt(x(:,1).^2+x(:,2).^2);
    for i = 1:size(x,1)
        dx(i,:) = x(i,:);
    end
%    save('results\doubletrigida.mat','t','x','equilibrio','tracaoeq', 'deltapeq');
    figure(longfig);
    subplot(2,2,1); plot(t,dx(:,1)); %velocidade eixo y
    subplot(2,2,2); plot(t,dx(:,2)); %velocidade eixo z
    subplot(2,2,3); plot(t,dx(:,3)); %q
    subplot(2,2,4); plot(t,dx(:,5)); legend('Flexible', 'Rigid'); %H
    
    %plot(t,[dx (Alpha-alpha) Veloc-V]); legend('u','w','q','theta', 'H','alpha','V');
    
    
    
end
function zero = equilibra(var,V,H)
    x(1) = V*cos(var(1));
    x(2) = -V*sin(var(1));
    x(3) = 0;
    x(4) = var(1);
    x(5) = H;
    u(1) = var(2);
    u(2) = var(3);
    
    xp = dinamica(0,x,u);
    zero = [xp(1), xp(2), xp(3)];
end
function [A B] = lineariza(xeq,ueq)
    delta = 1e-7;   
    for i = 1:5
        A(:,i) = (dinamica(0,xeq+veclin(i,5)*delta,ueq) -dinamica(0,xeq-veclin(i,5)*delta,ueq))/(2*delta);
    end
    for i = 1:2
        B(:,i) = (dinamica(0,xeq,ueq+veclin(i,2)*delta) -dinamica(0,xeq,ueq-veclin(i,2)*delta))/(2*delta);
    end    
end
function vec = veclin(i,numstates)
    vec = zeros(1,numstates);
    vec(1,i) = 1;
end
function xp = dinamica(t,x,u)
global softPARAMS
    % observação! note que o cálculo de clq foi feito através de uma
    % modificação no w aparente. isso foi feito para manter a coerência com
    % o que acontece no modelo flexível!
    T = u(1);
    deltap = u(2);
    u = x(1);
    w = x(2);
    q = x(3);
    theta = x(4);
    H = x(5);
    
    %dados do aviao
    m = 24+10; %kg
    rho = atmosfera(H);
    CD0 = 0.02;
    ME = 0.05;
    l = (-0.25-ME);
    d = -ME;
    cldelta = 0.01;
    cmdelta = -0.1;
    S = 32; c = 1;
    Iyy = 1.04;
    g = 9.8;
    alpha0 = 9.8*0.75*32/(0.5*0.0889*15^2 * 32 * 2*pi);
    alpha0 = -5*pi/180;
    ua = -(w + l*q)*sin(alpha0) + u*cos(alpha0);
    wa = (w + l*q)*cos(alpha0) + u*sin(alpha0);
    alphat = atan(-wa/ua);
    CL = 2*pi*((-wa/ua)) + cldelta*deltap;
    Cm = cmdelta*deltap;
    V = sqrt(u^2 + w^2);
    V = ua;
    L = 0.5*rho*V^2*S*CL;
    D = -0.5*rho*V^2*S*(CD0);
    
    Cwa0 = [1 0 0;
             0 cos(alpha0) sin(alpha0);
             0 -sin(alpha0) cos(alpha0)];  
    Ca0a1 = [1 0 0;
             0 cos(alphat) sin(alphat);
             0 -sin(alphat) cos(alphat)];
         
    CBA = Cwa0*Ca0a1;
    

    b = c/2;
    doty = ua;
    
    if softPARAMS.modAED == 0        
        alfap = 0;
    else
        alfap = q;
    end
    
    L = L + 32*2*pi*rho*b*doty * (b/2) *alfap;
    Forcaaero = CBA*[0;D;L];
    mx = (L*cos(alphat+alpha0)-D*sin(alphat+alpha0))*d + 0.5*rho*V^2*S*c*Cm; %atencao, pode melhorar aqui! as forças e braços devem ser escritas no sistema de coordenadas LOCAL
    
    mx = mx -0.5*pi*rho*b^3*V*alfap*32;
    %mx = 0.5 rho V^2 S CL -ME
    % Cmx = CL -ME/c = 2pi alfa 
    % CL = CLalfa alfa + CLq q
    Fx = T*u^(-1)+Forcaaero(2) - m*g*sin(theta); %-D*cos(alpha) + L*sin(alpha) - m*g*sin(theta);
    Fz = Forcaaero(3) - m*g*cos(theta);%L*cos(alpha) + D*sin(alpha) - m*g*cos(theta);
    up = q*w + Fx/m;
    wp = -q*u + Fz/m;
    qp = mx/ Iyy;
    thetap = q;
    
    Hp = u*sin(theta) + w*cos(theta);
    
    xp = [up,wp,qp,thetap, Hp]';
end