classdef Polycrystal
    %POLYCRYSTAL Summary of this class goes here
    %   Detailed explanation goes here

    properties
    end

    methods (Static)
        %% Voigt Bound
        function [Cvoigt] = Voigt4(Cijklccs,texture)
            %VOIGT Summary of this function goes here
            %   Detailed explanation goes here
            Cvoigt = zeros(3,3,3,3);
            ngrains = length(texture);
            for igrain = 1:ngrains   %Loop over all orientations
                phi1 = texture(igrain,1);
                PHI  = texture(igrain,2);
                phi2 = texture(igrain,3);
                gmatrix = Polycrystal.Orimat(phi1,PHI,phi2);
                Cijklscs = Polycrystal.Cijklscsgen1(Cijklccs,gmatrix);
                Cvoigt = Cvoigt + Cijklscs;
            end
            Cvoigt = Cvoigt/ngrains;
        end
        function [Cvoigt] = Voigt2(Cijklccs,texture)
            %VOIGT Summary of this function goes here
            %   Detailed explanation goes here
            Cvoigt = zeros(6,6);
            ngrains = length(texture);
            for igrain = 1:ngrains   %Loop over all orientations
                phi1 = texture(igrain,1);
                PHI  = texture(igrain,2);
                phi2 = texture(igrain,3);
                gmatrix = Polycrystal.Orimat(phi1,PHI,phi2);
                Cijklscs = Polycrystal.Cijklscsgen1(Cijklccs,gmatrix);
                Cij = Polycrystal.ijkl2Voigt(Cijklscs);
                Cvoigt = Cvoigt + Cij;
            end
            Cvoigt = Cvoigt/ngrains;
        end
        function [Cvoigt] = Voigt2a(Cijklccs,texture)
            %VOIGT Summary of this function goes here
            %   Detailed explanation goes here
            Cvoigt = zeros(6,6);
            ngrains = length(texture);
            CCS = [1 0 0; 0 1 0; 0 0 1];
            Cijccs = Polycrystal.ijkl2Voigt(Cijklccs);
            for igrain = 1:ngrains   %Loop over all orientations
                phi1 = texture(igrain,1);
                PHI  = texture(igrain,2);
                phi2 = texture(igrain,3);
                gmatrix = Polycrystal.Orimat(phi1,PHI,phi2);
                Cijscs = Polycrystal.Cijklscsgen2(Cijccs,gmatrix,CCS);
                Cvoigt = Cvoigt + Cijscs;
            end
            Cvoigt = Cvoigt/ngrains;
        end
        %% Reuss bound
        function [Creuss] = Reuss2(Cijklccs,texture)
            %Reuss Summary of this function goes here
            %   Detailed explanation goes here
            Creuss = zeros(6,6);
            ngrains = length(texture);
            for igrain = 1:ngrains   %Loop over all orientations
                phi1 = texture(igrain,1);
                PHI  = texture(igrain,2);
                phi2 = texture(igrain,3);
                gmatrix = Polycrystal.Orimat(phi1,PHI,phi2);
                Cijklscs = Polycrystal.Cijklscsgen1(Cijklccs,gmatrix);
                CijV = Polycrystal.ijkl2Voigt(Cijklscs);
                CijM = Polycrystal.Voigt2Mandel(CijV);
                Creuss = Creuss + inv(CijM);
            end
            Creuss = Creuss/ngrains;
            Creuss = inv(Creuss);
            Creuss = Polycrystal.Mandel2Voigt(Creuss);
        end
        %% Self-Consistent Method
        function [Csc] = autoSelfCons(Cijklccs,texture,grainshape,tolerance)
            % initial guess is Hill bound by default
            arguments
                Cijklccs (3,3,3,3)
                texture
                grainshape (1,3) = [1 1 1]
                tolerance = 1
            end
            Cvoigt2 = Polycrystal.Voigt2(Cijklccs,texture);
            Creuss2 = Polycrystal.Reuss2(Cijklccs,texture);
            Chill2 = 0.5*(Cvoigt2+Creuss2);

            CguessHill = Polycrystal.Voigt2ijkl(Chill2);
            [Csc,~] = Polycrystal.SelfCons(CguessHill,Cijklccs,grainshape,texture,tolerance);
        end

        function [Csc,iter] = SelfCons(CguessHill,Cijklccs,grainshape,texture,tolerance)
            Csc = CguessHill;
            % Cscguess = zeros(3,3,3,3);
            ngrains = length(texture);
            iter = 1;
            error = 1e9;
            while error > tolerance
                Cscguess = Csc;
                Csc2M = zeros(6,6);
                EshTen4 = Polycrystal.EshelbyTensor(Cscguess,grainshape);
                for igrain = 1:ngrains   %Loop over all orientations
                    phi1 = texture(igrain,1);
                    PHI  = texture(igrain,2);
                    phi2 = texture(igrain,3);
                    gmatrix = Polycrystal.Orimat(phi1,PHI,phi2);
                    Cijklscs = Polycrystal.Cijklscsgen1(Cijklccs,gmatrix);
                    EshTen2M = Polycrystal.Voigt2Mandel(Polycrystal.ijkl2Voigt(EshTen4));
                    % ReacTen2M = (eye(6)-EshTen2M)/EshTen2M;
                    Cscs2M = Polycrystal.Voigt2Mandel(Polycrystal.ijkl2Voigt(Cijklscs));
                    Cscguess2M = Polycrystal.Voigt2Mandel(Polycrystal.ijkl2Voigt(Cscguess));
                    % Csc2M = Csc2M + (Cscs2M + Cscguess2M*ReacTen2M)\(Cscguess2M + Cscguess2M*ReacTen2M)*Cscguess2M;
                    Csc2M = Csc2M + Cscs2M/(eye(6)+EshTen2M/(Cscguess2M)*(Cscs2M-Cscguess2M));
                end
                Csc2M = Csc2M/ngrains;
                error = norm((abs(Csc2M-Cscguess2M)));
                fprintf('Iteration:%d --> error = %f\n',iter,error);
                iter = iter + 1;
                Csc = Polycrystal.Voigt2ijkl(Polycrystal.Mandel2Voigt(Csc2M));
                % sqrt(sum(abs(Csc-Cscguess).*abs(Csc-Cscguess),'all'))/1e9
            end
            Csc = Polycrystal.ijkl2Voigt(Csc);
        end
        function [Eijmn] = EshelbyTensor(Cguess,grainshape)
            Gamma_ikjl = zeros(3,3,3,3);
            % Numerical Integration for Gamma_ikjl
            theta = 0:2:180;
            phi = 0:2:360;
            for itheta = 2:length(theta)
                for iphi = 2:length(phi)
                    Gamma_ikjl = Gamma_ikjl + (pi/180)^2*(Polycrystal.Gamma_atpoint(theta(itheta),phi(iphi),Cguess,grainshape)+ ...
                        Polycrystal.Gamma_atpoint(theta(itheta-1),phi(iphi),Cguess,grainshape)+ ...
                        Polycrystal.Gamma_atpoint(theta(itheta),phi(iphi-1),Cguess,grainshape)+ ...
                        Polycrystal.Gamma_atpoint(theta(itheta-1),phi(iphi-1),Cguess,grainshape));
                end
            end
            Eijmn = zeros(3,3,3,3);
            for i = 1:3
                for j = 1:3
                    for m = 1:3
                        for n = 1:3
                            dummy = 0;
                            for k = 1:3
                                for l = 1:3
                                    dummy = dummy + 0.25*(Gamma_ikjl(i,k,j,l) + Gamma_ikjl(i,l,j,k) + Gamma_ikjl(j,k,i,l) + Gamma_ikjl(j,l,i,k))*Cguess(k,l,m,n);
                                end
                            end
                            Eijmn(i,j,m,n) = dummy;
                        end
                    end
                end
            end
        end
        function [Gamma_atpt] = Gamma_atpoint(itheta,iphi,Cguess,grainshape)
            a1 = grainshape(1);
            a2 = grainshape(2);
            a3 = grainshape(3);
            xi1 = sind(itheta)*cosd(iphi)/a1;
            xi2 = sind(itheta)*sind(iphi)/a2;
            xi3 = cosd(itheta)/a3;
            xi = [xi1 xi2 xi3];
            Kip = zeros(3,3);
            for i = 1:3
                for p = 1:3
                    Kipdummy = 0;
                    for j = 1:3
                        for l = 1:3
                            Kipdummy = Kipdummy + Cguess(i,j,p,l)*xi(j)*xi(l);
                        end
                    end
                    Kip(i,p) = Kipdummy;
                end
            end
            gamma_ikjl = zeros(3,3,3,3);
            Kik_ = inv(Kip);
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        for l = 1:3
                            gamma_ikjl(i,k,j,l) = Kik_(i,k)*xi(j)*xi(l);
                        end
                    end
                end
            end
            Gamma_atpt = (1/4/pi)*sind(itheta)*gamma_ikjl;
        end

        %% Suppooting functions
        % Cijkl generator for cubic ---------------------------------------
        function [Cijkl] = Cijklccsgen(C111244_cubic)
            C11 = C111244_cubic(1);
            C12 = C111244_cubic(2);
            C44 = C111244_cubic(3);
            Cijkl = zeros(3,3,3,3);
            Iijkl = zeros(3,3,3,3);
            Iijkl(1,1,1,1) = 1;
            Iijkl(2,2,2,2) = 1;
            Iijkl(3,3,3,3) = 1;
            Iij = [1,0,0;0,1,0;0,0,1];
            lambda = C12;
            mu = C44;
            mu_ = C11 - C12 - 2*C44;
            for i=1:3
                for j=1:3
                    for k=1:3
                        for l=1:3
                            Cijkl(i,j,k,l)= lambda*Iij(i,j)*Iij(k,l)+mu*Iij(i,k)*Iij(j,l)+mu*Iij(i,l)*Iij(j,k)+mu_*Iijkl(i,j,k,l);
                        end
                    end
                end
            end
        end

        % Cijkl: ccs to scs - method 1 using index notation ---------------
        function [Cijklscs] = Cijklscsgen1(Cijklccs,gmatrix)
            Cijklscs = zeros(3,3,3,3);
            gT = gmatrix';
            for i=1:3
                for j=1:3
                    for k=1:3
                        for l=1:3
                            dummy = 0;
                            for i1=1:3
                                for j1=1:3
                                    for k1=1:3
                                        for l1=1:3
                                            dummy = dummy + gT(i,i1)*gT(j,j1)*gT(k,k1)*gT(l,l1)*Cijklccs(i1,j1,k1,l1);
                                        end
                                    end
                                end
                            end
                            Cijklscs(i,j,k,l) = dummy;
                        end
                    end
                end
            end
        end

        % Cijkl: ccs to scs - method 2 using transformation matrix --------
        function [Cscs] = Cijklscsgen2 (Cccs,SCS,CCS)
            % HELP - orientationC (Stiffness Matrix)
            %   -- stiffness matrix of the crystal under the given loads in SCS
            %   -- effective Young's modulus in the direction of loading

            % Stiffness Matrix in SCS
            % Compliance matrix of the crystal in CCS
            Sccs = inv(Cccs);

            % a -- transformation matrix from SCS to CCS
            a = zeros(3,3);
            for i = 1:3
                for j = 1:3
                    a(i,j) = dot(CCS(:,i),SCS(:,j))/(norm(CCS(:,i))*norm(SCS(:,j)));
                end
            end
            % A -- transformation matrix from SCS to CCS in voigt notation
            % stress[CCS] = A*stress[SCS]
            A = [a(1,1)*a(1,1) a(1,2)*a(1,2) a(1,3)*a(1,3) 2*a(1,2)*a(1,3) 2*a(1,3)*a(1,1) 2*a(1,1)*a(1,2);
                a(2,1)*a(2,1) a(2,2)*a(2,2) a(2,3)*a(2,3) 2*a(2,2)*a(2,3) 2*a(2,3)*a(2,1) 2*a(2,1)*a(2,2);
                a(3,1)*a(3,1) a(3,2)*a(3,2) a(3,3)*a(3,3) 2*a(3,2)*a(3,3) 2*a(3,3)*a(3,1) 2*a(3,1)*a(3,2);
                a(2,1)*a(3,1) a(2,2)*a(3,2) a(2,3)*a(3,3) (a(2,2)*a(3,3)+a(2,3)*a(3,2)) (a(2,3)*a(3,1)+a(2,1)*a(3,3)) (a(2,1)*a(3,2)+a(2,2)*a(3,1));
                a(3,1)*a(1,1) a(3,2)*a(1,2) a(3,3)*a(1,3) (a(3,2)*a(1,3)+a(3,3)*a(1,2)) (a(3,3)*a(1,1)+a(3,1)*a(1,3)) (a(3,1)*a(1,2)+a(3,2)*a(1,1));
                a(1,1)*a(2,1) a(1,2)*a(2,2) a(1,3)*a(2,3) (a(1,2)*a(2,3)+a(1,3)*a(2,2)) (a(1,3)*a(2,1)+a(1,1)*a(2,3)) (a(1,1)*a(2,2)+a(1,2)*a(2,1))];

            % b -- transformation matrix from CCS to SCS
            b = zeros(3,3);
            for i = 1:3
                for j = 1:3
                    b(i,j) = dot(SCS(:,i),CCS(:,j))/(norm(SCS(:,i))*norm(CCS(:,j)));
                end
            end
            % B -- transformation matrix from CCS to SCS in voigt notation;
            % strain[SCS] = B*strain[CCS]
            B = [b(1,1)*b(1,1) b(1,2)*b(1,2) b(1,3)*b(1,3) b(1,2)*b(1,3) b(1,3)*b(1,1) b(1,1)*b(1,2);
                b(2,1)*b(2,1) b(2,2)*b(2,2) b(2,3)*b(2,3) b(2,2)*b(2,3) b(2,3)*b(2,1) b(2,1)*b(2,2);
                b(3,1)*b(3,1) b(3,2)*b(3,2) b(3,3)*b(3,3) b(3,2)*b(3,3) b(3,3)*b(3,1) b(3,1)*b(3,2);
                2*b(2,1)*b(3,1) 2*b(2,2)*b(3,2) 2*b(2,3)*b(3,3) (b(2,2)*b(3,3)+b(2,3)*b(3,2)) (b(2,3)*b(3,1)+b(2,1)*b(3,3)) (b(2,1)*b(3,2)+b(2,2)*b(3,1));
                2*b(3,1)*b(1,1) 2*b(3,2)*b(1,2) 2*b(3,3)*b(1,3) (b(3,2)*b(1,3)+b(3,3)*b(1,2)) (b(3,3)*b(1,1)+b(3,1)*b(1,3)) (b(3,1)*b(1,2)+b(3,2)*b(1,1));
                2*b(1,1)*b(2,1) 2*b(1,2)*b(2,2) 2*b(1,3)*b(2,3) (b(1,2)*b(2,3)+b(1,3)*b(2,2)) (b(1,3)*b(2,1)+b(1,1)*b(2,3)) (b(1,1)*b(2,2)+b(1,2)*b(2,1))];

            % Strain[SCS] = S_*Stress[SCS]
            % Strain[SCS] = B*Strain[CCS] = B*S*Stress[CCS] = B*S*A*Stress[SCS]
            % which gives S_ = B*S*A -- Compliance Matrix in SCS

            Sscs = B*Sccs*A;                                                            % Compliance matrix of the crystal in SCS
            Cscs = inv(Sscs);                                                           % Stiffness matrix of the crystal in SCS
        end

        % Cijkl to Cvoigt to Cmandel --------------------------------------
        function [Cv]=ijkl2Voigt(C)
            Cv = [C(1,1,1,1) C(1,1,2,2) C(1,1,3,3) C(1,1,2,3) C(1,1,3,1) C(1,1,1,2);
                C(2,2,1,1) C(2,2,2,2) C(2,2,3,3) C(2,2,2,3) C(2,2,3,1) C(2,2,1,2);
                C(3,3,1,1) C(3,3,2,2) C(3,3,3,3) C(3,3,2,3) C(3,3,3,1) C(3,3,1,2);
                C(2,3,1,1) C(2,3,2,2) C(2,3,3,3) C(2,3,2,3) C(2,3,3,1) C(2,3,1,2);
                C(3,1,1,1) C(3,1,2,2) C(3,1,3,3) C(3,1,2,3) C(3,1,3,1) C(3,1,1,2);
                C(1,2,1,1) C(1,2,2,2) C(1,2,3,3) C(1,2,2,3) C(1,2,3,1) C(1,2,1,2)];
        end
        function [Cijkl]=Voigt2ijkl(Cv)
            Cijkl = zeros(3,3,3,3);
            for i=1:3
                for j=1:3
                    for k=1:3
                        for l=1:3
                            if i==1 && j==1
                                m = 1;
                            elseif i==2 && j==2
                                m = 2;
                            elseif i==3 && j==3
                                m = 3;
                            elseif (i==2 && j==3) || (i==3 && j==2)
                                m = 4;
                            elseif (i==3 && j==1) || (i==1 && j==3)
                                m = 5;
                            elseif (i==1 && j==2) || (i==2 && j==1)
                                m = 6;
                            end
                            if k==1 && l==1
                                n = 1;
                            elseif k==2 && l==2
                                n = 2;
                            elseif k==3 && l==3
                                n = 3;
                            elseif (k==2 && l==3) || (k==3 && l==2)
                                n = 4;
                            elseif (k==3 && l==1) || (k==1 && l==3)
                                n = 5;
                            elseif (k==1 && l==2) || (k==2 && l==1)
                                n = 6;
                            end
                            Cijkl(i,j,k,l) = Cv(m,n);
                        end
                    end
                end
            end
        end
        function [Cm]=Voigt2Mandel(C)
            Cm = [C(1,1) C(1,2) C(1,3) sqrt(2)*C(1,4) sqrt(2)*C(1,5) sqrt(2)*C(1,6);
                C(2,1) C(2,2) C(2,3) sqrt(2)*C(2,4) sqrt(2)*C(2,5) sqrt(2)*C(2,6);
                C(3,1) C(3,2) C(3,3) sqrt(2)*C(3,4) sqrt(2)*C(3,5) sqrt(2)*C(3,6);
                sqrt(2)*C(4,1) sqrt(2)*C(4,2) sqrt(2)*C(4,3) 2*C(4,4) 2*C(4,5) 2*C(4,6);
                sqrt(2)*C(5,1) sqrt(2)*C(5,2) sqrt(2)*C(5,3) 2*C(5,4) 2*C(5,5) 2*C(5,6);
                sqrt(2)*C(6,1) sqrt(2)*C(6,2) sqrt(2)*C(6,3) 2*C(6,4) 2*C(6,5) 2*C(6,6)];
        end
        function [Cv]=Mandel2Voigt(C)
            Cv = [C(1,1) C(1,2) C(1,3) C(1,4)/sqrt(2) C(1,5)/sqrt(2) C(1,6)/sqrt(2);
                C(2,1) C(2,2) C(2,3) C(2,4)/sqrt(2) C(2,5)/sqrt(2) C(2,6)/sqrt(2);
                C(3,1) C(3,2) C(3,3) C(3,4)/sqrt(2) C(3,5)/sqrt(2) C(3,6)/sqrt(2);
                C(4,1)/sqrt(2) C(4,2)/sqrt(2) C(4,3)/sqrt(2) C(4,4)/2 C(4,5)/2 C(4,6)/2;
                C(5,1)/sqrt(2) C(5,2)/sqrt(2) C(5,3)/sqrt(2) C(5,4)/2 C(5,5)/2 C(5,6)/2;
                C(6,1)/sqrt(2) C(6,2)/sqrt(2) C(6,3)/sqrt(2) C(6,4)/2 C(6,5)/2 C(6,6)/2];
        end

        % Orientation Matrix ----------------------------------------------
        function [gmatrix] = Orimat(phi1,PHI,phi2)
            % bunge
            gmatrix = [cosd(phi1)*cosd(phi2)-sind(phi1)*sind(phi2)*cosd(PHI), sind(phi1)*cosd(phi2)+cosd(phi1)*sind(phi2)*cosd(PHI), sind(phi2)*sind(PHI);
                -cosd(phi1)*sind(phi2)-sind(phi1)*cosd(phi2)*cosd(PHI), -sind(phi1)*sind(phi2)+cosd(phi1)*cosd(phi2)*cosd(PHI), cosd(phi2)*sind(PHI);
                sind(phi1)*sind(PHI), -cosd(phi1)*sind(PHI), cosd(PHI)];
            % symmetric
            % gmatrix = [-sind(phi1)*sind(phi2)-cosd(phi1)*cosd(phi2)*cosd(PHI), cosd(phi1)*sind(phi2)-sind(phi1)*cosd(phi2)*cosd(PHI), cosd(phi2)*sind(PHI);
            %     sind(phi1)*cosd(phi2)-cosd(phi1)*sind(phi2)*cosd(PHI), -cosd(phi1)*cosd(phi2)-sind(phi1)*sind(phi2)*cosd(PHI), sind(phi2)*sind(PHI);
            %     cosd(phi1)*sind(PHI), sind(phi1)*sind(PHI), cosd(PHI)];
        end

        % Zener ratio -----------------------------------------------------
        function [ZenerR] = ZenerRat(C)
            if length(C)==3
                ZenerR = 2*C(2,3,2,3)/(C(1,1,1,1)-C(1,1,2,2));
            elseif length(C) == 6
                ZenerR = 2*C(4,4)/(C(1,1)-C(1,2));
            end
        end

        % show 6-components -----------------------------------------------
        function [out] = show9(C)
            if length(C)==3
                out = [C(1,1,1,1);C(2,2,2,2);C(3,3,3,3);C(1,1,2,2);C(1,1,3,3);C(2,2,3,3);C(2,3,2,3);C(3,1,3,1);C(1,2,1,2)]';
                out = join(string(out), '    ');
            elseif length(C) == 6
                out = [C(1,1);C(2,2);C(3,3);C(1,2);C(1,3);C(2,3);C(4,4);C(5,5);C(6,6)]';
                out = join(string(out), '    ');
            end
        end

        % Directional Young's modulus for rolling -------------------------
        function PlotDYMrol(Cscs)
            alpha = 0:90;
            DYM = zeros(length(alpha),1);
            if length(Cscs)==6
                Cscs = Polycrystal.Voigt2ijkl(Cscs);
            end
            for ialpha = alpha
                gmatrix = [cosd(ialpha) -sind(ialpha) 0;
                    sind(ialpha) cosd(ialpha) 0;
                    0 0 1];
                Cscs_ = Polycrystal.ijkl2Voigt(Polycrystal.Cijklscsgen1(Cscs,gmatrix));
                Sscs_ = inv(Cscs_);
                DYM(ialpha+1) = 1/(Sscs_(1,1));
            end
            plot(alpha,DYM)
        end

        % Elastic isotropic stiffness tensor ------------------------------
        function [C]=isoCij(E,nu)
            % E = (9*K*mu)/(3*K+mu);
            % nu = (3*K-2*mu)/2/(3*K+mu);
            C = (E/(1+nu)/(1-2*nu))*[1-nu nu nu 0 0 0;
                nu 1-nu nu 0 0 0;
                nu nu 1-nu 0 0 0;
                0 0 0 (1-2*nu)/2 0 0;
                0 0 0 0 (1-2*nu)/2 0;
                0 0 0 0 0 (1-2*nu)/2];
        end

        % Tensor Anisotropy Index -----------------------------------------
        function [TAI] = TenAniInd(C,Aa)
            % refernce -- https://doi.org/10.1007/s00707-018-2174-7
            arguments
                C (6,6)
                Aa = false
            end
            Cg1 = [C(1,1) C(2,2) C(3,3)];
            Cg2 = [C(4,4) C(5,5) C(6,6)];
            Cg3 = [C(1,2) C(1,3) C(2,3)];
            Cg4 = [C(3,4) C(4,5) C(5,6)];
            Cg5 = [C(1,4) C(2,5) C(3,6)];
            Cg6 = [C(2,4) C(3,5) C(4,6) C(1,5) C(2,6) C(1,6)];

            Aiz = 2*sum(Cg2)/(sum(Cg1)-sum(Cg3));
            Aiconv = std(Cg1)/mean(Cg1) + std(Cg2)/mean(Cg2) + std(Cg3)/mean(Cg3);
            n = length(find([Cg4 Cg5 Cg6]));
            if Aa
                Aam = (n/12)*mean([Cg4 Cg5 Cg6])/mean([Cg1 Cg2 Cg3]);
                Aaconv = (var(Cg4)+var(Cg5)+var(Cg6))/mean([Cg1 Cg2 Cg3]);
            else
                Aam = 0;
                Aaconv = 0;
            end
            TAI = Aiz + Aiconv + Aam + Aaconv;
        end

        % direction young's modulus 3D
        function PlotDYM3d2(Cscs)
            [xcoord,ycoord,zcoord] = sphere(200);
            DYM = zeros(size(zcoord));
            if length(Cscs)==6
                Cscs = Polycrystal.Voigt2ijkl(Cscs);
            end
            for i = 1:length(xcoord)
                for j = 1:length(xcoord)
                    [az,el,~] = cart2sph(xcoord(i,j),ycoord(i,j),zcoord(i,j));
                    gmatrixZ = [cos(az) -sin(az) 0;
                        sin(az) cos(az) 0;
                        0 0 1];
                    CscsND = Polycrystal.Cijklscsgen1(Cscs,gmatrixZ);

                    gmatrixY = [cos(el) 0 sin(el);
                        0 1 0;
                        -sin(el) 0 cos(el)];
                    CscsTD = Polycrystal.ijkl2Voigt(Polycrystal.Cijklscsgen1(CscsND,gmatrixY));

                    Sscs_ = inv(CscsTD);
                    DYM(i,j) = 1/(Sscs_(1,1));
                end
            end
            % figure('Name',"Directional Young's Modulus",'NumberTitle','off');
            set(gca,'FontSize',13)
            surf(xcoord,ycoord,zcoord,DYM/1e9,"LineStyle",'none')
            grid on;
            box on;
            axis on;
            xlabel x-axis/RD
            ylabel y-axis/TD
            zlabel z-axis/ND
            xticks([-1 -0.5 0 0.5 1])
            yticks([-1 -0.5 0 0.5 1])
            zticks([-1 -0.5 0 0.5 1])
            set(gca,'FontSize',12)
            set(gca, 'color', 'none');                                  % To remove background
            % set(gca, 'Projection','perspective')
            set(gcf,'units','pixels','position',[300 300 550 500]);     % To change the size of the figure
            view(35,30)
            axis equal
            colormap turbo
            c = colorbar;
            c.Label.String = "Directional Young's Modulus in GPa";
            % optional
            % t = get(c,'Limits');
            % T = linspace(t(1),t(2),7);
            % set(c,'Ticks',T)
            % TL = arrayfun(@(x) sprintf('%.1f',x),T,'un',0);
            % set(c,'TickLabels',TL)
        end

        % direction young's modulus 3D
        function PlotDYM3d(Cscs)
            [xcoord,ycoord,zcoord] = sphere(50);
            DYM = zeros(size(zcoord));
            azcoord = zeros(size(zcoord));
            elcoord = zeros(size(zcoord));
            if length(Cscs)==6
                Cscs = Polycrystal.Voigt2ijkl(Cscs);
            end
            for i = 1:length(xcoord)
                for j = 1:length(xcoord)
                    [az,el,~] = cart2sph(xcoord(i,j),ycoord(i,j),zcoord(i,j));
                    gmatrixZ = [cos(az) -sin(az) 0;
                        sin(az) cos(az) 0;
                        0 0 1];
                    CscsND = Polycrystal.Cijklscsgen1(Cscs,gmatrixZ);

                    gmatrixY = [cos(el) 0 sin(el);
                        0 1 0;
                        -sin(el) 0 cos(el)];
                    CscsTD = Polycrystal.ijkl2Voigt(Polycrystal.Cijklscsgen1(CscsND,gmatrixY));

                    Sscs_ = inv(CscsTD);
                    DYM(i,j) = 1/(Sscs_(1,1));
                    azcoord(i,j) = az;
                    elcoord(i,j) = el;
                end
            end
            [xcoord2,ycoord2,zcoord2] = sph2cart(azcoord,elcoord,DYM);
            % figure('Name',"Directional Young's Modulus",'NumberTitle','off');
            set(gca,'FontSize',13)
            % surf(xcoord2,ycoord2,zcoord2,DYM/1e9,"LineStyle",'none')
            surf(xcoord2,ycoord2,zcoord2,DYM/1e9,'LineWidth',0.25,'LineStyle',':')
            % hold on
            % xline(0,'-','y-axis/TD','LabelHorizontalAlignment','right','LineWidth',2)
            % yline(0,'-','x-axis/RD','LineWidth',2)
            % line([0,0],[0,0],[min(zcoord2,[],"all"),max(zcoord2,[],"all")])
            grid on;
            box on;
            axis on;
            xlabel x-axis/RD
            ylabel y-axis/TD
            zlabel z-axis/ND
            xticks(linspace(min(xcoord2,[],"all"),max(xcoord2,[],"all"),5))
            yticks(linspace(min(ycoord2,[],"all"),max(ycoord2,[],"all"),5))
            zticks(linspace(min(zcoord2,[],"all"),max(zcoord2,[],"all"),5))
            set(gca,'xticklabel',[])
            set(gca,'yticklabel',[])
            set(gca,'zticklabel',[])
            set(gca,'FontSize',12)
            set(gca, 'color', 'none');                                  % To remove background
            % set(gca, 'Projection','perspective')
            set(gcf,'units','pixels','position',[300 300 550 500]);     % To change the size of the figure
            view(35,30)
            axis equal
            colormap turbo
            c = colorbar;
            c.Label.String = "Directional Young's Modulus in GPa";
            % optional
            % t = get(c,'Limits');
            % T = linspace(t(1),t(2),7);
            % set(c,'Ticks',T)
            % TL = arrayfun(@(x) sprintf('%.1f',x),T,'un',0);
            % set(c,'TickLabels',TL)
        end

    end
end

