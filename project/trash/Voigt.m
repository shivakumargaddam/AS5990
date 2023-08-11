function [voigtbound] = Voigt(C111244_Cu,texture)
%VOIGT Summary of this function goes here
%   Detailed explanation goes here
voigtbound = zeros(3,3,3,3);
Cijklccs = Cijklccsgen(C111244_Cu);
ngrains = length(texture);
for igrain = 1:ngrains   %Loop over all orientations
    phi1 = texture(igrain,1);
    PHI  = texture(igrain,2);
    phi2 = texture(igrain,3);
    % g-matrix (Orientation Matrix)
    gmatrix = [cosd(phi1)*cosd(phi2)-sind(phi1)*sind(phi2)*cosd(PHI), sind(phi1)*cosd(phi2)+cosd(phi1)*sind(phi2)*cosd(PHI), sind(phi2)*sind(PHI);
        -cosd(phi1)*sind(phi2)-sind(phi1)*cosd(phi2)*cosd(PHI), -sind(phi1)*sind(phi2)+cosd(phi1)*cosd(phi2)*cosd(PHI), cosd(phi2)*sind(PHI);
        sind(phi1)*sind(PHI), -cosd(phi1)*sind(PHI), cosd(PHI)];
    Cijklscs = Cijklscsgen(Cijklccs,gmatrix);
    voigtbound = voigtbound + Cijklscs;
end
voigtbound = voigtbound/ngrains;
end

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

function [Cijklscs] = Cijklscsgen(Cijklccs,gmatrix)
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