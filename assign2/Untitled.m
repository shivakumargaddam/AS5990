clc
clear all
close all

% epoxy
K0 = 2.83e9;
mu0 = 1.31e9;
E0 = (9*K0*mu0)/(3*K0+mu0);
nu0 = (3*K0-2*mu0)/2/(3*K0+mu0);

% glass fibers
K1 = 51.2e9;
mu1 = 35.2e9;
E1 = (9*K1*mu1)/(3*K1+mu1);
nu1 = (3*K1-2*mu1)/2/(3*K1+mu1);

syms f

%% Problem 1 and 2
C0 = Voigt2Mandel(StiffnessTensor(K0,mu0))
C1 = Voigt2Mandel(StiffnessTensor(K1,mu1))
S0 = Voigt2Mandel(EshelbyTensorCylinder(K0,mu0))
A1 = inv(eye(6) + S0*inv(C0)*(C1-C0))
Cbar_mori = ((1-f)*C0+f*C1*A1)*inv((1-f)*eye(6)+f*A1);
Cbar_mori_voigt = Mandel2Voigt(Cbar_mori)
(Cbar_mori_voigt(1,1)-Cbar_mori_voigt(1,2))/2

%% Functions
function [S]=EshelbyTensorCylinder(K,mu)
E = (9*K*mu)/(3*K+mu);
nu = (3*K-2*mu)/2/(3*K+mu);
S = zeros(6,6);
S(1,1) = (5-4*nu)/8/(1-nu);
S(2,2) = S(1,1);
S(1,2) = (-1+4*nu)/8/(1-nu);
S(2,1) = S(1,2);
S(1,3) = nu/2/(1-nu);
S(2,3) = S(1,3);
S(4,4) = 0.25;
S(5,5) = 0.25;
S(6,6) = (3-4*nu)/8/(1-nu);
end
function [C]=StiffnessTensor(K,mu)
E = (9*K*mu)/(3*K+mu);
nu = (3*K-2*mu)/2/(3*K+mu);
C = (E/(1+nu)/(1-2*nu))*[1-nu nu nu 0 0 0;
    nu 1-nu nu 0 0 0;
    nu nu 1-nu 0 0 0;
    0 0 0 (1-2*nu)/2 0 0;
    0 0 0 0 (1-2*nu)/2 0;
    0 0 0 0 0 (1-2*nu)/2];
end
function [Cm]=Voigt2Mandel(C)
Cm = [C(1,1) C(1,2) C(1,3) sqrt(2)*C(1,4) sqrt(2)*C(1,5) sqrt(2)*C(1,6);
    C(2,1) C(2,2) C(2,3) sqrt(2)*C(2,4) sqrt(2)*C(2,5) sqrt(2)*C(2,6);
    C(3,1) C(3,2) C(3,3) sqrt(2)*C(3,4) sqrt(2)*C(3,5) sqrt(2)*C(3,6);
    sqrt(2)*C(4,1) sqrt(2)*C(4,2) sqrt(2)*C(4,3) 2*C(4,4) 2*C(4,5) 2*C(4,6);
    sqrt(2)*C(5,1) sqrt(2)*C(5,2) sqrt(2)*C(5,3) 2*C(5,4) 2*C(5,5) 2*C(5,6);
    sqrt(2)*C(6,1) sqrt(2)*C(6,2) sqrt(2)*C(6,3) 2*C(6,4) 2*C(6,5) 2*C(6,6)];
end
function [Cm]=Mandel2Voigt(C)
Cm = [C(1,1) C(1,2) C(1,3) C(1,4)/sqrt(2) C(1,5)/sqrt(2) C(1,6)/sqrt(2);
    C(2,1) C(2,2) C(2,3) C(2,4)/sqrt(2) C(2,5)/sqrt(2) C(2,6)/sqrt(2);
    C(3,1) C(3,2) C(3,3) C(3,4)/sqrt(2) C(3,5)/sqrt(2) C(3,6)/sqrt(2);
    C(4,1)/sqrt(2) C(4,2)/sqrt(2) C(4,3)/sqrt(2) C(4,4)/2 C(4,5)/2 C(4,6)/2;
    C(5,1)/sqrt(2) C(5,2)/sqrt(2) C(5,3)/sqrt(2) C(5,4)/2 C(5,5)/2 C(5,6)/2;
    C(6,1)/sqrt(2) C(6,2)/sqrt(2) C(6,3)/sqrt(2) C(6,4)/2 C(6,5)/2 C(6,6)/2];
end
function [S]=EshelbyTensorTisomatrix(C)
S = zeros(6,6);
S(1,1) = (5*C(1,1)+C(1,2))/8/C(1,1);
S(2,2) = S(1,1);
S(1,2) = (3*C(1,2)-C(1,1))/8/C(1,1);
S(2,1) = S(1,2);
S(1,3) = C(1,3)/2/C(1,1);
S(2,3) = S(1,3);
S(4,4) = 0.25;
S(5,5) = 0.25;
S(6,6) = (3*C(1,1)-C(1,2))/8/C(1,1);
end
function [Ecl,Ect,Gip,Gop,nulf] = TransverseElasticConstants(C)
Ecl = C(3,3)-((2*C(1,3)*C(1,3))/(C(1,1)+C(1,2)));
Ect = ((C(1,1)-C(1,2))*(C(1,1)*C(3,3)+C(1,2)*C(3,3)-2*C(1,3)*C(1,3)))/(C(1,1)*C(3,3)-C(1,3)*C(1,3));
Gip = C(6,6);
Gop = C(4,4);
nulf = C(1,3)/(C(1,1)+C(1,2));
end