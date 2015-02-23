%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2011- Flávio Luiz C. Ribeiro
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotaest3d(membro)
    membersize = size(membro,2);
    j = 1;
    for i=1:membersize
        deflexao(:,j) = membro(i).node1.h; j = j+1;        
        deflexao(:,j) = membro(i).node2.h; j = j+1;
        deflexao(:,j) = membro(i).node3.h; j = j+1;
    end
    for i=1:size(deflexao,2)
        switch mod(i,3)
            case 0
                a = membro(floor((i-1)/3)+1).node3.aero.a;
                b = membro(floor((i-1)/3)+1).node3.aero.b;
            case 1
                a = membro(floor((i-1)/3)+1).node1.aero.a;
                b = membro(floor((i-1)/3)+1).node1.aero.b;
            case 2
                a = membro(floor((i-1)/3)+1).node2.aero.a;
                b = membro(floor((i-1)/3)+1).node2.aero.b;
        end
        X(i,1) = deflexao(1,i);% + deflexao(7,i);     
        X(i,2) = deflexao(1,i);
        X(i,3) = deflexao(1,i);
        
        Y(i,1) = deflexao(2,i) + deflexao(8,i)*(b*a+b);
        Y(i,2) = deflexao(2,i);
        Y(i,3) = deflexao(2,i)+deflexao(8,i)*(b*a-b);
        
        Z(i,1) = deflexao(3,i) + deflexao(9,i)*(b*a+b);        
        Z(i,2) = deflexao(3,i);              
        Z(i,3) = deflexao(3,i)+deflexao(9,i)*(b*a-b);              
    end
    
    surf(X,Y,Z);
    hold on; %plot3(X,Y,Z,'o');
    axis([-20 20 -20 20]);
end