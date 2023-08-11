clear all
close all
clc

% epoxy
K0 = 2.83e9;
mu0 = 1.31e9;
E0 = (9*K0*mu0)/(3*K0+mu0);
nu0 = (3*K0-2*mu0)/2/(3*K0+mu0);
lambda0 = K0 - 2*mu0/3;

% glass fibers
K1 = 51.2e9;
mu1 = 35.2e9;
E1 = (9*K1*mu1)/(3*K1+mu1);
nu1 = (3*K1-2*mu1)/2/(3*K1+mu1);
lambda1 = K1 - 2*mu1/3;

%% Problem 4,5,6
Gff_ipK = [];
Gmori_ipK = [];
Gself_ipK = [];
Greuss_ipK = [];
Gvoigt_ipK= [];


volfrac = 0:0.01:1;

for f = volfrac
    %% mori-tanaka
    C0 = Voigt2Mandel(StiffnessTensor(K0,mu0));
    C1 = Voigt2Mandel(StiffnessTensor(K1,mu1));
    S0 = Voigt2Mandel(EshelbyTensorCylinder(K0,mu0));
    A1 = inv(eye(6) + S0*inv(C0)*(C1-C0));
    Cbar_mori = ((1-f)*C0+f*C1*A1)*inv((1-f)*eye(6)+f*A1);
    Cbar_mori_voigt = Mandel2Voigt(Cbar_mori);
    mori_ipK = (Cbar_mori(1,1)+Cbar_mori(1,2))/2;

    %% self-consistent
    Cbar_self = Cbar_mori;
    Cbar_self_dummy = zeros(6,6);
    i=1;
    while norm(Cbar_self-Cbar_self_dummy) > 1e-4
        Cbar_self_dummy = Cbar_self;
        S_ = EshelbyTensorTisomatrix(Cbar_self_dummy);
        A1_ = inv(eye(6) + S_*inv(Cbar_self_dummy)*(C1-Cbar_self_dummy));
        Cbar_self = C0 + f*(C1-C0)*A1_;
        i = i+1;
    end
    C_self_voigt = Mandel2Voigt(Cbar_self);
    self_ipK = (Cbar_self(1,1)+Cbar_self(1,2))/2;
    
    %% reuss and voigt
    C_reussbound = inv((1-f)*inv(C0)+f*inv(C1));
    C_reussbound_voigt = Mandel2Voigt(C_reussbound);
    reuss_ipK = (C_reussbound(1,1)+C_reussbound(1,2))/2;
    
    C_voigtbound = (1-f)*C0+f*C1;
    C_voigtbound_voigt = Mandel2Voigt(C_voigtbound);
    voigt_ipK = (C_voigtbound(1,1)+C_voigtbound(1,2))/2;
    %% full field method
    ff_ipK = (mu0+lambda0) + (2*mu0+lambda0)*f/((mu0+lambda1+mu1)/(mu1+lambda1-mu0-lambda0) - f);
    
    %% ipK - inplane bulk modulus
    Gff_ipK = [Gff_ipK; ff_ipK];
    Gmori_ipK = [Gmori_ipK; mori_ipK];
    Gself_ipK = [Gself_ipK; self_ipK];
    Greuss_ipK = [Greuss_ipK; reuss_ipK];
    Gvoigt_ipK= [Gvoigt_ipK; voigt_ipK];
    
end

%% Plots
p1 = round(length(volfrac)*3/6);
p2 = round(length(volfrac)*4/6);
p3 = round(length(volfrac)*4.5/6);
p4 = round(length(volfrac)*5/6);
p5 = round(length(volfrac)*5/6);

figure
set(gca,'FontSize',13)
hold on
plot(volfrac',Gff_ipK,'g','MarkerSize',17,'LineWidth',5,'DisplayName','Full-Field')
plot(volfrac',Gmori_ipK,'MarkerSize',17,'LineWidth',1.5,'DisplayName','Mori-Tanaka')
plot(volfrac',Gself_ipK,'MarkerSize',17,'LineWidth',1.5,'DisplayName','Self-Consistent')
plot(volfrac',Greuss_ipK,'MarkerSize',17,'LineWidth',1.5,'DisplayName','Reuss')
plot(volfrac',Gvoigt_ipK,'MarkerSize',17,'LineWidth',1.5,'DisplayName','Voigt')
text(volfrac(p1),Gff_ipK(p1),'\leftarrow Full-Field')
text(volfrac(p2),Gmori_ipK(p2),'\leftarrow Mori-Tanaka')
text(volfrac(p3),Gself_ipK(p3),'\leftarrow Self-Consistent')
text(volfrac(p4),Greuss_ipK(p4),'\leftarrow Reuss')
text(volfrac(p5),Gvoigt_ipK(p5),'\leftarrow Voigt')
xlabel('volume fraction - f')
ylabel('Plane-Strain K [Pa]')
legend('Location','northwest')
box on
pbaspect([1 1 1]);
set(gcf,'units','pixels','position',[100 100 600 600]);
title("Plane-Strain Bulk Modulus")


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