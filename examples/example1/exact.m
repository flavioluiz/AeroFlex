function exact_freqs = exact(num_bending_vert, num_bending_hor, num_torsion)

    % beam length
    Length = 16;
    
    % sectional rigidity matrix    
    K22 = 1e4; %GJ
    K33 = 2e4; %flat bend: EI
    K44 = 4e6; %chord bend: EI
        
    I22 = 0.0;
    I33 = 0.1;
    I11 = 0.1;
    mcs = 0.75;

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