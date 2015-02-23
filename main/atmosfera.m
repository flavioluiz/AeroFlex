%......................Modelo densidade...................................
function rho = atmosfera(H)
g0 = 9.81;  
rho0 = 1.225;
R = 287.053;L=-6.5e-3;
T0 = 288.15;

if H<= 11000
    rho = rho0*(temperatura(H)/T0)^-(1+g0/(R*L));
else
    rho = rho0*(temperatura(11000)/T0)^-(1+g0/(R*L))*exp(-g0*(H-11000)/(R*temperatura(11000)));
end

%......................Modelo Temperatura..................................
function T = temperatura(H)
T0 = 288.15;
L = -6.5e-3;
H0 = 0;
if H <= 11000
    T = T0 + L*(H-H0);
else
    T = temperatura(11000);
end