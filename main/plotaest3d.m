%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2011- Flávio Luiz C. Ribeiro
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotest = plotaest3d(member, translate)
    if nargin < 2
        translate = [0,0,0];
    end
    if isempty(member(1).node3.geometry)
        fprintf('\n WARNING: Can`t plot structure without geometric data \n');
    else
        membersize = size(member,2);
        j = 1;
        for i=1:membersize
            deflexao(:,j) = member(i).node1.h; j = j+1;        
            deflexao(:,j) = member(i).node2.h; j = j+1;
            deflexao(:,j) = member(i).node3.h; j = j+1;
        end
        for i=1:size(deflexao,2)
            switch mod(i,3)
                case 0
                    a = member(floor((i-1)/3)+1).node3.geometry.a;
                    b = member(floor((i-1)/3)+1).node3.geometry.b;
                case 1
                    a = member(floor((i-1)/3)+1).node1.geometry.a;
                    b = member(floor((i-1)/3)+1).node1.geometry.b;
                case 2
                    a = member(floor((i-1)/3)+1).node2.geometry.a;
                    b = member(floor((i-1)/3)+1).node2.geometry.b;
            end
            X(i,1) = deflexao(1,i)+translate(1);% + deflexao(7,i);     
            X(i,2) = deflexao(1,i)+translate(1);
            X(i,3) = deflexao(1,i)+translate(1);

            Y(i,1) = deflexao(2,i) + deflexao(8,i)*(b*a+b)+translate(2);
            Y(i,2) = deflexao(2,i) + translate(2);
            Y(i,3) = deflexao(2,i)+deflexao(8,i)*(b*a-b) + translate(2);

            Z(i,1) = deflexao(3,i) + deflexao(9,i)*(b*a+b) +translate(3);        
            Z(i,2) = deflexao(3,i) + translate(3);              
            Z(i,3) = deflexao(3,i)+deflexao(9,i)*(b*a-b)+ translate(3); 
        end

        plotest = surf(X,Y,Z);
        hold on; %plot3(X,Y,Z,'o');
        axis([-20 20 -20 20]);
    end
end