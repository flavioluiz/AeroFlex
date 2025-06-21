function exact_freqs = exact(num_bending_vert, num_bending_hor, num_torsion)

    % beam length
    Length = 0.55;
    c = 0.06;
    t = 0.00225;
    A = c*t;
    Ixx = c*t^3/12;
    Iyy = t*c^3/12;
    J = 0.33*c*t^3;
    E = 71.7e9;
    G = 26.9e9;
    % sectional rigidity matrix
    K11 = E*A; %EA
    K22 = G*J; %GJ
    K33 = E*Ixx; %flat bend: EI
    K44 = E*Iyy; %chord bend: EI
    KG = diag([K11 K22 K33 K44]);
    rho = 2810;
    I22 = Ixx*rho*0;
    I33 = Iyy*rho*0;
    I11 = rho*(Ixx+Iyy);
    mcs = rho*A; %mass per unit length (kg/m)

    exact_bending_vert = exact_bending(num_bending_vert, K33, mcs, Length);
    exact_bending_hor = exact_bending(num_bending_hor, K44, mcs, Length);
    exact_tor = exact_torsion(num_torsion,K22, I11, Length);
    
    exact_freqs = [exact_bending_vert(:);exact_bending_hor(:);exact_tor(:)];

end

function exactfreqs = exact_bending(NUMBEN, EI, mass_length, length)
%% Calculate exact natural frequencies of euler-bernoulli beam
% Function usage: exactfreqs = eb_exact(NUNBEN)
%   where NUNBEN is the number of bending modes;
%   exactfreqs is a vector of frequencies in Hertz
%
    
    %modal frequencies omegan
        betan0 = [0.1 3 6 8 10 13 15 17 19.5 22 24 26 29 31 33.3 35.6]*1.36;
                % betan0: starting values for the fsolve algorithm
        for i = 1:NUMBEN
            betan(i) = fsolve(@(x) cos(x)*cosh(x) + 1,betan0(i));            
            betan(i) = betan(i)/length;
        end
        omegan = (betan).^2 * sqrt(EI/(mass_length)); % this finds a VECTOR of omegan (natural frequencies)
    
     exactfreqs = omegan'/2/pi;
    
end

function exactfreqs = exact_torsion(NUMTOR, GJ, Ip, Length)
%% Calculate exact natural frequencies of torsion beam
% Function usage: exactfreqs = tb_exact(NUMTOR)
%   where NUMTOR is the number of torsion modes;
%   exactfreqs is a vector of frequencies in Hertz
%
    omegan = sqrt(GJ/Ip)/Length*[1:2:(NUMTOR*2-1)]/2  /2;
    exactfreqs = omegan(:);
    
end