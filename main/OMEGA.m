function [saida]  = OMEGA(omega,nodesnumber)
    saida = zeros(4*nodesnumber);
    for i = 1:(nodesnumber*4)
        saida((1+(i-1)*3):(3+(i-1)*3),(1+(i-1)*3):(3+(i-1)*3)) = omega;
    end
end