classdef Peter
   
    properties
        n;                  % No. of degrees of freedom    [1,2,3 or 4]
        m;                  % No. of 
        b;                  % Semi-cord length             [Units of lenght]
        N;                  % No. of aerodynamic inflow states 
        
        a;                  % Center of pressure position  [% of b]
        c;                  % Flap hinge position          [% of b]
        d;                  % Tab hinge position           [% of b]
        
        ndelta;             % No da superficie de controle que ativa o deltaflap
        

        
        alpha0;
        clalpha;
        cldelta;
        cm0;
        cmdelta;
        cd0;
        
        T;
        K;
        C;
        M;
        A;
        B1;
        B2;
        Bm;
        
        E1;
        E2;
        E3;
        E4;
        
        B;
        
        memberE1;
        memberE2;
        memberE3;
        memberE4;
        
        matrizE2Jpep;
        matrizE3Jthetaep;
        matrizE4Jthetaep;
    
    end 

    methods
       
        function obj = Peter(varargin)
            
            obj.n = 2;
            obj.m = 2;
            obj.N = varargin{1};
            obj.b = varargin{2};
            obj.a = varargin{3};
            

            if nargin == 6
                obj.c = varargin{6}; 
            end
            
            obj = calcmat(obj);
            obj.E1 = -inv(obj.A); %cuidado! *Vwind/b
            obj.E2 = inv(obj.A)*obj.B1(:,1);
            obj.E3 = inv(obj.A)*obj.B1(:,2);
            obj.E4 = inv(obj.A)*obj.B2(:,2);  %cuidado! *Vwind/b
            
          %   obj.N = 2;
          %   obj.E1 = [-0.6828 -0.0633;1 0];
          %   obj.E2 = [2*pi;0];
          %   obj.E3 = [2*pi*(1-2*obj.a)/4;0];    
          %   obj.E4 = inv(obj.A)*obj.B2(:,2)/obj.b*0;
            %produtos por Vwind são feitos dentro do cálculo aerodinâmico
            %(aerounsteady.m)
        
        end
        
        function setmembermatrices(aero,membro)
            membersize = size(membro,2);
            membro(1).node1.aero.memberE1 = zeros(aero.N*membersize*3,3*membersize*aero.N);
            membro(1).node1.aero.memberE2 = zeros(aero.N*membersize*3,3*membersize);
            membro(1).node1.aero.memberE3 = zeros(aero.N*membersize*3,3*membersize);
            membro(1).node1.aero.memberE4 = zeros(aero.N*membersize*3,3*membersize);
            membro(1).node1.aero.matrizE2Jpep = zeros(membersize*3,membersize*9);
            membro(1).node1.aero.matrizE3Jthetaep = zeros(membersize*3,membersize*9);
            membro(1).node1.aero.matrizE4Jthetaep = zeros(membersize*3,membersize*9);
            
            for i = 1:membersize
                membro(1).node1.aero.memberE1((1+(i-1)*aero.N*3):(aero.N+(i-1)*aero.N*3),(1+(i-1)*aero.N*3):(aero.N+(i-1)*aero.N*3))=membro(i).node1.aero.E1;
                membro(1).node1.aero.memberE1((aero.N+1+(i-1)*aero.N*3):(2*aero.N+(i-1)*aero.N*3),(aero.N+1+(i-1)*aero.N*3):(2*aero.N+(i-1)*aero.N*3))=membro(i).node2.aero.E1;
                membro(1).node1.aero.memberE1((2*aero.N+1+(i-1)*aero.N*3):(3*aero.N+(i-1)*aero.N*3),(2*aero.N+1+(i-1)*aero.N*3):(3*aero.N+(i-1)*aero.N*3))=membro(i).node3.aero.E1;
                
                membro(1).node1.aero.memberE2((1+(i-1)*aero.N*3):(aero.N+(i-1)*aero.N*3),(i-1)*3+1)=membro(i).node1.aero.E2;
                membro(1).node1.aero.memberE2((aero.N+1+(i-1)*aero.N*3):(2*aero.N+(i-1)*aero.N*3),(i-1)*3+2)=membro(i).node2.aero.E2;
                membro(1).node1.aero.memberE2((2*aero.N+1+(i-1)*aero.N*3):(3*aero.N+(i-1)*aero.N*3),(i-1)*3+3)=membro(i).node3.aero.E2;
                
                membro(1).node1.aero.memberE3((1+(i-1)*aero.N*3):(aero.N+(i-1)*aero.N*3),(i-1)*3+1)=membro(i).node1.aero.E3;
                membro(1).node1.aero.memberE3((aero.N+1+(i-1)*aero.N*3):(2*aero.N+(i-1)*aero.N*3),(i-1)*3+2)=membro(i).node2.aero.E3;
                membro(1).node1.aero.memberE3((2*aero.N+1+(i-1)*aero.N*3):(3*aero.N+(i-1)*aero.N*3),(i-1)*3+3)=membro(i).node3.aero.E3;
                
                membro(1).node1.aero.memberE4((1+(i-1)*aero.N*3):(aero.N+(i-1)*aero.N*3),(i-1)*3+1)=membro(i).node1.aero.E4;
                membro(1).node1.aero.memberE4((aero.N+1+(i-1)*aero.N*3):(2*aero.N+(i-1)*aero.N*3),(i-1)*3+2)=membro(i).node2.aero.E4;
                membro(1).node1.aero.memberE4((2*aero.N+1+(i-1)*aero.N*3):(3*aero.N+(i-1)*aero.N*3),(i-1)*3+3)=membro(i).node3.aero.E4;
                
                membro(1).node1.aero.matrizE2Jpep((1+(i-1)*3),(1+(i-1)*9):(3+(i-1)*9)) = [0 0 1];
                membro(1).node1.aero.matrizE2Jpep((2+(i-1)*3),(4+(i-1)*9):(6+(i-1)*9)) = [0 0 1];
                membro(1).node1.aero.matrizE2Jpep((3+(i-1)*3),(7+(i-1)*9):(9+(i-1)*9)) = [0 0 1];
                
                membro(1).node1.aero.matrizE3Jthetaep((1+(i-1)*3),(1+(i-1)*9):(3+(i-1)*9)) = [1 0 0];
                membro(1).node1.aero.matrizE3Jthetaep((2+(i-1)*3),(4+(i-1)*9):(6+(i-1)*9)) = [1 0 0];
                membro(1).node1.aero.matrizE3Jthetaep((3+(i-1)*3),(7+(i-1)*9):(9+(i-1)*9)) = [1 0 0];
                
                membro(1).node1.aero.matrizE4Jthetaep((1+(i-1)*3),(1+(i-1)*9):(3+(i-1)*9)) = [1 0 0];
                membro(1).node1.aero.matrizE4Jthetaep((2+(i-1)*3),(4+(i-1)*9):(6+(i-1)*9)) = [1 0 0];
                membro(1).node1.aero.matrizE4Jthetaep((3+(i-1)*3),(7+(i-1)*9):(9+(i-1)*9)) = [1 0 0];
                
            end
                                                
        end
    end
    
    % Private Methods
    methods (Access='protected')
            
        function obj = calcmat(obj)
           
            n = obj.n;
            m = obj.m;
            N = obj.N;
            b = obj.b;
            a = obj.a;

            %---------- Calculate T ------------
            obj.T = zeros(m,n);
            
            % First DOF (Plunge)
            obj.T(1,1) = b;
            
            % Second DOF (Pitch)
            if (n >= 2)
                obj.T(1,2) = -b*a;
                obj.T(2,2) = b;
            end
            
            % Third DOF (Flap) 
            if (n >=3)
                p = acos(obj.c);
                if m >= 1 
                    obj.T(1,3) = (b/pi)*(sin(p) - p*cos(p));
                end
                if m >= 2
                    obj.T(2,3) = (b/pi)*(p - sin(p)*cos(p));
                end
                for ith=2:m-1
                    Ti = ( (1/(ith+1))*sin((ith+1)*p) + ... 
                           (1/(ith-1))*sin((ith-1)*p) - ...
                           (2/ith)*cos(p)*sin(ith*p) );
                    obj.T(ith+1,3) = (b/pi)*Ti;
                end
            end
            
            %--------- Calculate K -------------
            obj.K = zeros(m,m);
            for i=1:m
                obj.K(1,i) = i - 1; 
                obj.K(i,i) = -(i - 1)/2;
            end
            
            %--------- Calculate C -------------
            obj.C = zeros(m,m);  
            if m >= 1
                obj.C(1,1) = 1;
            end
            if m >= 2
                obj.C(1,2) = 1;
                obj.C(2,1) = -1/2;
            end
            if m >= 3;
                obj.C(2:m,2:m) = diag(repmat(-1/2,[1 m-2]),-1)+diag(repmat(1/2,[1 m-2]),1);
            end
            
                      
            %--------- Calculate M -------------
            obj.M = zeros(m,m);  
            if m >= 1
                obj.M(1,1) = 1/2;
            end
            if m >= 2
                obj.M(2,2) = 1/16;
            end
            if m >= 3;
                obj.M(1,3) = -1/4;
                obj.M(3,1) = -1/4;
                
                for i=3:m
                   obj.M(i,i) = (i-1)/(4*(((i-1)^2)-1));
                   if i < m
                      obj.M(i-1,i+1) = -1/(8*(i-1));
                      obj.M(i+1,i-1) = -1/(8*(i-1));
                   end
                end
            end
            
            %--------- Calculate B1 and B2 ------
            c2 = zeros(N,1);
            for i=1:N
                c2(i) = 2/i;
            end
            s2 = zeros(m,1);
            s3 = s2;
            for i=1:m
                if i == 1
                    s2(i) = 1;
                elseif  i == 2
                    s2(i) = 1/2;
                end
                s3(i) = i - 1;
            end
            obj.B1 = c2*s2'*obj.T;
            obj.B2 = c2*s3'*obj.T;
            
            
            
            
            %--------- Calculate A -------------
            obj.A = zeros(N,N);
            D = zeros(N,N);
            b = zeros(N,1);
            c = zeros(N,1);
            d = zeros(N,1);
            
            for i=1:N
                for j=1:N
                    if (i == j + 1)
                        D(i,j) = 1/(2*i);
                    elseif (i == j - 1)
                        D(i,j) = -1/(2*i);
                    end
                end
                if i ~= N
                    b(i) = ((-1)^(i-1))*(factorial(N-1+i)/factorial(N-1-i))*(1/factorial(i)^2);
                end
                if i == N
                    b(i) = (-1)^(N+1);
                end

                c(i) = 2/i;
                
                if i == 1;
                    d(i) = 1/2;
                end
            end   
            obj.B = b;
            obj.Bm = obj.C*[1;zeros(m-1,1)]*(1/2)*b';
            obj.A = D + d*b' + c*d' + (1/2)*c*b';
            
        end
    end
end