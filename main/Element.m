%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2011- Flávio Luiz C. Ribeiro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Element < handle
% The Element class brings together all the information about the airplane
% flexible structures.
% Each element has the data about structural rigidity and damping;
% In addition, each element is linked to three Node objects. The node
% objects have the informations about distributed and concentrated masses,
% as well as aerodynamic data.
%
%
% Methods of this object computes the Structural jacobians as well as
% structural displacements, from the strain vector.
%

    properties
        memberJhep;
        memberJpep;
        memberJthetaep;
        memberB;
        memberJhb;
        memberJpb;
        memberJthetab;
        strainm;
        strainpm;
        strainppm;
        lambdam;
    end
    properties (SetAccess = private)
        Me; %matriz de inércia do elemento
        Ne; % matriz gravitacional do elemento (2.39 Brown)
        node1; %nohs associados ao elemento
        node2; %nohs associados ao elemento
        node3; %nohs associados ao elemento
        ds; %comprimento do elemento
        D; %matriz de rotacao entre o elemento atual e o anterior
        strainvec = [];
        strainpvec = [];
        h0;
        Jhepp
        KG;
        CG;
        expG;
        exp2G;
        dexpGdEps;
        dexp2GdEps;
    end
    methods
        function el = Element(node1,node2,node3,rot,ds,KG, CG) %inicializa o elemento de viga, associando a 3 nohs
            M1 = node1.Mnode;
            M2 = node2.Mnode;
            M3 = node3.Mnode;
            el.ds = ds;
            el.KG = ds*KG;
            el.CG = ds*CG;
            el.node1 = node1;
            el.node2 = node2;
            el.node3 = node3;
            el.Me = 0.5*ds*[(1/4*M1 + 1/12* M2), (1/12*M1+1/12*M2), zeros(12,12);
                     (1/12*M1+1/12*M2),(1/12*M1+1/2*M2+1/12*M3), (1/12*M2+1/12*M3);
                     zeros(12,12), (1/12*M2+1/12*M3), (1/12*M2+1/4*M3)]+...
                     [node1.Mrb zeros(12,12) zeros(12,12);
                      zeros(12,12) node2.Mrb zeros(12,12);
                      zeros(12,12) zeros(12,12) node3.Mrb];
            N1 = M1(:,1:3);
            N2 = M2(:,1:3);
            N3 = M3(:,1:3);
            el.Ne = 0.5*ds*[N1/3 + N2/6;
                N1/6 + 2/3*N2 + N3/6;
                N2/6 + N3/3] + ...
                [node1.Mrb(:,1:3);
                node2.Mrb(:,1:3);
                node3.Mrb(:,1:3)];
            
            c = [1 0 0; 0 cos(rot.twist) -sin(rot.twist); 0 sin(rot.twist) cos(rot.twist)]*...
                [cos(rot.sweep) -sin(rot.sweep) 0; sin(rot.sweep) cos(rot.sweep) 0; 0 0 1]*...
                [cos(rot.dihedral) 0 sin(rot.dihedral); 0 1 0; -sin(rot.dihedral) 0 cos(rot.dihedral)];
            el.D = [eye(3) zeros(3,9); ...
                zeros(3,3) eye(3)*c(1,1) eye(3)*c(1,2) eye(3)*c(1,3);
                zeros(3,3) eye(3)*c(2,1) eye(3)*c(2,2) eye(3)*c(2,3);
                zeros(3,3) eye(3)*c(3,1) eye(3)*c(3,2) eye(3)*c(3,3)];
        end
        function setstrain(el, strain, strainp)
            el.strainvec = strain;
            el.strainpvec = strainp;
        end
        function seth0(el, h0)
            el.h0 = el.D*h0;
        end
        function setexponentials(el)
            membersize = size(el,2);
            for i =1:membersize
                [el(i).expG, el(i).dexpGdEps] = getexpG(el(i).strainvec,el(i).ds/2);
                [el(i).exp2G, el(i).dexp2GdEps] = getexpG(el(i).strainvec,el(i).ds);
            end
        end
        function sethnodes(el)
            membersize = size(el,2);
            for i = 1:membersize
                if i ~= 1
                    el(i).h0 = el(i).D*el(i-1).node3.h;
                end
                el(i).node1.h = el(i).h0;
                el(i).node2.h = el(i).expG*el(i).h0;
                el(i).node3.h = el(i).exp2G*el(i).h0;
            end
        end
        function setJhep(el)
            membersize = size(el,2);
            for i = 1:membersize
                h0 = hdiag(el(i).h0);
                el(i).node1.Jhep = zeros(12,4);
                el(i).node2.Jhep = el(i).dexpGdEps*h0;
                el(i).node3.Jhep = el(i).dexp2GdEps*h0;
                el(i).Jhep = [el(i).node1.Jhep;el(i).node2.Jhep;el(i).node3.Jhep];
                %el(i).Jhep = getnumJhep(el(i).strainvec,el(i).h0,el(i).ds);
            end
        end
        function setJhepp(el)
            membersize = size(el,2);
            for i = 1:membersize
                Jhepep1 = getJhepep(el(i).strainvec, el(i).h0, 0);
                Jhepep2 = getJhepep(el(i).strainvec, el(i).h0, el(i).ds);
                Jhepep3 = getJhepep(el(i).strainvec, el(i).h0, 2*el(i).ds);
                el(i).Jhepp = product3dXvec([Jhepep1; Jhepep2; Jhepep3],el(i).strainpvec);
                %el.Jhepp = zeros(36,4);
            end
        end
        function [Jhep] = getmemberJhep(el)
            membersize = size(el,2);
            Jhep = zeros(membersize*[12*3,4]);
            for j = 1:membersize
                for i = 1:membersize
                    if i < j
                        %mantem-se zero
                    elseif i == j
                        h0 = hdiag(el(i).h0);                        
                        Jhep1 = zeros(12,4);
                        Jhep2 = el(i).dexpGdEps*h0;
                        Jhep3 = el(i).dexp2GdEps*h0;
                        Jhep(((i-1)*36+1):((i-1)*36+36),((j-1)*4+1):((j-1)*4+4)) = [Jhep1; Jhep2; Jhep3];
                    else
                        Jhep(((i-1)*36+1):((i-1)*36+12),((j-1)*4+1):((j-1)*4+4))=el(i).D*Jhep(((i-2)*36+25):((i-2)*36+36),((j-1)*4+1):((j-1)*4+4));
                        Jhep(((i-1)*36+13):((i-1)*36+24),((j-1)*4+1):((j-1)*4+4))=el(i).expG*Jhep(((i-1)*36+1):((i-1)*36+12),((j-1)*4+1):((j-1)*4+4));
                        Jhep(((i-1)*36+25):((i-1)*36+36),((j-1)*4+1):((j-1)*4+4))=el(i).expG*Jhep(((i-1)*36+13):((i-1)*36+24),((j-1)*4+1):((j-1)*4+4));
                    end
                end
            end
        end
        function [Jhb] = getmemberJhb(el)
            membersize = size(el,2);
            Jhb = zeros([membersize*12*3,6]);
            for i = 1:membersize
                %node1
                h = el(i).node1.h;
                p = h(1:3); wx = h(4:6); wy = h(7:9); wz = h(10:12);
                Jhb((1+3*12*(i-1)):(3+3*12*(i-1)),1:3) = eye(3);
                Jhb((1+3*12*(i-1)):(12+3*12*(i-1)),4:6) = [matrixcross(p);matrixcross(wx);matrixcross(wy);matrixcross(wz)];
                
                %node2
                h = el(i).node2.h;
                p = h(1:3); wx = h(4:6); wy = h(7:9); wz = h(10:12);
                Jhb((13+3*12*(i-1)):(15+3*12*(i-1)),1:3) = eye(3);
                Jhb((13+3*12*(i-1)):(24+3*12*(i-1)),4:6) = [matrixcross(p);matrixcross(wx);matrixcross(wy);matrixcross(wz)];                

                %node3
                h = el(i).node3.h;
                p = h(1:3); wx = h(4:6); wy = h(7:9); wz = h(10:12);
                Jhb((25+3*12*(i-1)):(27+3*12*(i-1)),1:3) = eye(3);
                Jhb((25+3*12*(i-1)):(36+3*12*(i-1)),4:6) = [matrixcross(p);matrixcross(wx);matrixcross(wy);matrixcross(wz)];                                
            end
        end
        function [Jpbeta] = getmemberJpb(el,Jhb)
            membersize = size(el,2);
            Jpbeta = zeros(9*membersize,6);
            for i = 1:(membersize*3)
                Jpbeta((1+3*(i-1)):(3+3*(i-1)),:) = Jhb((1+(i-1)*12):(3+(i-1)*12),:);
            end
        end
        function [Jthetabeta] = getmemberJthetab(el)
            membersize = size(el,2);
            Jthetabeta = zeros([membersize*3,6]);
            ident = eye(3);
            for i = 1:membersize
                Jthetabeta((1+9*(i-1)):(9+9*(i-1)),4:6) = [ident;ident;ident];
            end
        end
        function Jpep = getJpep(member,Jhep)
            membersize = size(Jhep,2)/4;
            Jpep = zeros(membersize*3,membersize*4);
            for i=1:membersize*3
                m = 1+(i-1)*12;
                n = 3+(i-1)*12;
                Jpep((1+(i-1)*3):(3+(i-1)*3),:) = Jhep(m:n,:);
            end
        end
        function [Me] = getmemberMe(el)
            membersize = size(el,2);
            Me = zeros(membersize*[36,36]);
            for i = 1:membersize
                Me(((i-1)*36+1):((i-1)*36+36),((i-1)*36+1):((i-1)*36+36)) = el(i).Me;
            end
        end
        function N = getmemberN(el)
            membersize = size(el,2);
            N = zeros([membersize*36,3]);
            for i = 1:membersize
                N(((i-1)*36+1):((i-1)*36+36),1:3) = el(i).Ne;
            end            
        end
        function [KG] = getmemberKG(el)
            membersize = size(el,2);
            KG = zeros(membersize*[4,4]);
            for i = 1:membersize
                KG(((i-1)*4+1):((i-1)*4+4),((i-1)*4+1):((i-1)*4+4)) = el(i).KG;
            end
        end
        function [CG] = getmemberCG(el)
            membersize = size(el,2);
            CG = zeros(membersize*[4,4]);
            for i = 1:membersize
                CG(((i-1)*4+1):((i-1)*4+4),((i-1)*4+1):((i-1)*4+4)) = el(i).CG;
            end
        end
        function Bf = getB(membro) %matriz B para forcas e momentos distribuidos
            membersize = size(membro,2);
            Bfe = [eye(3)*1/3, eye(3)*1/6, zeros(3);
                   eye(3)*1/6, eye(3)*2/3,eye(3)*1/6;
                   zeros(3), eye(3)*1/6,eye(3)*1/3];
            Bf = zeros(membersize*[9,9]);
            for i = 1:membersize
                Bf(((i-1)*9+1):((i-1)*9+9),((i-1)*9+1):((i-1)*9+9)) = 0.5*membro(i).ds*Bfe;
            end
        end
        function Jthetaep = getJthetaep(membro,Jhep) 
            % obtem o jacobiano J theta epsilon. Deve-se procurar otimizar
            % essa função. Ela faz apenas operações simples de
            % multiplicação e soma, buscando deslocamentos e o jacobiano
            % Jheps. Para uma viga de 10 elementos, a função leva até 0.1s
            % para obter o jacobiano.
            membersize = size(membro,2);
            Jthetaep = zeros(membersize*[9,4]);
            for i = 1:membersize
                for ii = 1:membersize
                    for j=1:4
                        if i<ii
                        % nao faz nada
                        else
                        % cuidado! estou considerando que o sistema do corpo NÃO
                        % GIRA!!! VER SHEARER 183, D.39
                        % noh 1
                        dtdepz = Jhep((4+(i-1)*36):(6+(i-1)*36),j+(ii-1)*4)'*membro(i).node1.h(7:9);
                        dtdepx = Jhep((7+(i-1)*36):(9+(i-1)*36),j+(ii-1)*4)'*membro(i).node1.h(10:12);
                        dtdepy = Jhep((10+(i-1)*36):(12+(i-1)*36),j+(ii-1)*4)'*membro(i).node1.h(4:6);
                        dtdex1 = [membro(i).node1.h(4:6),membro(i).node1.h(7:9),membro(i).node1.h(10:12)]*[dtdepx;dtdepy;dtdepz];

                        %noh2
                        dtdepz = Jhep((4+12+(i-1)*36):(6+12+(i-1)*36),j+(ii-1)*4)'*membro(i).node2.h(7:9);
                        dtdepx = Jhep((7+12+(i-1)*36):(9+12+(i-1)*36),j+(ii-1)*4)'*membro(i).node2.h(10:12);
                        dtdepy = Jhep((10+12+(i-1)*36):(12+12+(i-1)*36),j+(ii-1)*4)'*membro(i).node2.h(4:6);
                        dtdex2 = [membro(i).node2.h(4:6),membro(i).node2.h(7:9),membro(i).node2.h(10:12)]*[dtdepx;dtdepy;dtdepz];

                        %noh3
                        dtdepz = Jhep((4+24+(i-1)*36):(6+24+(i-1)*36),j+(ii-1)*4)'*membro(i).node3.h(7:9);
                        dtdepx = Jhep((7+24+(i-1)*36):(9+24+(i-1)*36),j+(ii-1)*4)'*membro(i).node3.h(10:12);
                        dtdepy = Jhep((10+24+(i-1)*36):(12+24+(i-1)*36),j+(ii-1)*4)'*membro(i).node3.h(4:6);
                        dtdex3 = [membro(i).node3.h(4:6),membro(i).node3.h(7:9),membro(i).node3.h(10:12)]*[dtdepx;dtdepy;dtdepz];
                        
                        Jthetaep(((i-1)*9+1):((i-1)*9+9),((ii-1)*4+j)) = [dtdex1;dtdex2;dtdex3];
                        end
                    end
                end
            end
        end
        function update(el)
            setexponentials(el);
            sethnodes(el);  
           % setJhep(el);
            %setJhepp(el)
        end
    end
end

function matrix = product3dXvec(A, v)
    matrix = zeros(size(A,1), size(A,2));
    for i=1:size(A,3)
        matrix(:,i) = reshape(A(:,i,:),36,4)*v';
    end
end

function Jhepep = getJhepep(strain, h0, xi)
    Jhepep = zeros(12,4,4);
    for i=1:4
        delta = 1e-5;
        JHepp = getJhep(strain + delta*veclin(i),h0,xi);
        JHepm = getJhep(strain - delta*veclin(i),h0,xi);
        Jhepep(:,:,i) = (JHepp-JHepm)/(2*delta);
    end
end

function [Jhep] = getJhep(strain,h0, ds)
    [~, dexpGdEps] = getexpG(strain,ds/2);
    %h = expG*h0;
    Jhep = dexpGdEps*[h0,zeros(12,3);zeros(12,1),h0,zeros(12,2);zeros(12,2),h0,zeros(12,1);zeros(12,3),h0];%[dexpGdEps(1:12,:)*h0, dexpGdEps(13:24,:)*h0,dexpGdEps(25:36,:)*h0,dexpGdEps(37:48,:)*h0];
end

function [Jhep] = getnumJhep(strain,h0,ds)
        for i=1:4
            delta = 1e-5;
            [expG, ~] = getnumexpG(strain + delta*veclin(i),ds/2);
            [expGm, ~] = getnumexpG(strain - delta*veclin(i),ds/2);
            [exp2G, ~] = getnumexpG(strain + delta*veclin(i),ds);
            [exp2Gm, ~] = getnumexpG(strain - delta*veclin(i),ds);
            h0p = h0;
            h0m = h0;
            h1p = expG*h0;
            h1m = expGm*h0;
            h2p = exp2G*h0;
            h2m = exp2Gm*h0;
            Jhep(:,i) = ([h0p;h1p;h2p]-[h0m;h1m;h2m])/(2*delta);
        end
end

function [expG, dexpGdeps] = getnumexpG(strain, ds)
          ex = strain(1); kx = strain(2)*eye(3);
          ky = strain(3)*eye(3); kz = strain(4)*eye(3);
          exb = (1+ex)*eye(3);            
          KEps = [zeros(3), exb, zeros(3), zeros(3);
                  zeros(3), zeros(3), kz, -ky;
                  zeros(3),-kz,zeros(3),kx;
                  zeros(3), ky,-kx,zeros(3)];
          expG = expm(KEps*ds);
          dexpGdeps = [];
end

function vec = veclin(i)
    vec = [0 0 0 0];
    vec(i) = 1;
end

function matrix = hdiag(h0)
    matrix = [h0,zeros(12,3);zeros(12,1),h0,zeros(12,2);zeros(12,2),h0,zeros(12,1);zeros(12,3),h0];
end