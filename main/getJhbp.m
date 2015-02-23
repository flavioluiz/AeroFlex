%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2011- Flávio Luiz C. Ribeiro
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [saida]  = getJhbp(hp,nodesnumber)
    saida = zeros(12*nodesnumber,6);
    for i = 1:(nodesnumber)
        saida((1+(i-1)*12):(3+(i-1)*12),4:6) = matrixcross(hp((1 + (i-1)*12):(3 + (i-1)*12)));
        saida((4+(i-1)*12):(6+(i-1)*12),4:6) = matrixcross(hp((4 + (i-1)*12):(6 + (i-1)*12)));
        saida((7+(i-1)*12):(9+(i-1)*12),4:6) = matrixcross(hp((7 + (i-1)*12):(9 + (i-1)*12)));
        saida((10+(i-1)*12):(12+(i-1)*12),4:6) = matrixcross(hp((10 + (i-1)*12):(12 + (i-1)*12)));
    end
end