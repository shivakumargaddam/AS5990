clear all
close all
clc

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

%% Problem 4,5,6
Ecl = [];
GmoriEcl = [];
GselfEcl = [];
GreussEcl = [];
GvoigtEcl = [];

Ect = [];
GmoriEct = [];
GselfEct = [];
GreussEct = [];
GvoigtEct = [];

Gip = [];
GmoriGip = [];
GselfGip = [];
GreussGip = [];
GvoigtGip = [];

Gop = [];
GmoriGop = [];
GselfGop = [];
GreussGop = [];
GvoigtGop = [];

nulf = [];
Gmorinulf = [];
Gselfnulf = [];
Greussnulf = [];
Gvoigtnulf = [];

volfrac = 0:0.01:1;

for f = volfrac
    %% mori-tanaka
    C0 = Voigt2Mandel(StiffnessTensor(K0,mu0));
    C1 = Voigt2Mandel(StiffnessTensor(K1,mu1));
    S0 = Voigt2Mandel(EshelbyTensorCylinder(K0,mu0));
    A1 = inv(eye(6) + S0*inv(C0)*(C1-C0));
    Cbar_mori = ((1-f)*C0+f*C1*A1)*inv((1-f)*eye(6)+f*A1);
    Cbar_mori_voigt = Mandel2Voigt(Cbar_mori);
    [moriEcl,moriEct,moriGip,moriGop,morinulf] = TransverseElasticConstants(Cbar_mori_voigt);
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
    [selfEcl,selfEct,selfGip,selfGop,selfnulf] = TransverseElasticConstants(C_self_voigt);
    %% reuss and voigt
    C_reussbound = inv((1-f)*inv(C0)+f*inv(C1));
    C_reussbound_voigt = Mandel2Voigt(C_reussbound);
    [reussEcl,reussEct,reussGip,reussGop,reussnulf] = TransverseElasticConstants(C_reussbound_voigt);
    C_voigtbound = (1-f)*C0+f*C1;
    C_voigtbound_voigt = Mandel2Voigt(C_voigtbound);
    [voigtEcl,voigtEct,voigtGip,voigtGop,voigtnulf] = TransverseElasticConstants(C_voigtbound_voigt);
    
    %% E - longitudinal
    Ecl = [Ecl; f*E1 + (1-f)*E0];
    GmoriEcl = [GmoriEcl;moriEcl];
    GselfEcl = [GselfEcl;selfEcl];
    GreussEcl = [GreussEcl;reussEcl];
    GvoigtEcl = [GvoigtEcl;voigtEcl];
    %% E - transverse
    Ect = [Ect; 1/(f/E1+(1-f)/E0)];
    GmoriEct = [GmoriEct;moriEct];
    GselfEct = [GselfEct;selfEct];
    GreussEct = [GreussEct;reussEct];
    GvoigtEct = [GvoigtEct;voigtEct];
    %% G - inplane
    Gip = [Gip; 1/(f/mu1+(1-f)/mu0)];
    GmoriGip = [GmoriGip; moriGip];
    GselfGip = [GselfGip; selfGip];
    GreussGip = [GreussGip; reussGip];
    GvoigtGip = [GvoigtGip; voigtGip];
    %% G - outplane
    Gop = [Gop; f*mu1 + (1-f)*mu0];
    GmoriGop = [GmoriGop; moriGop];
    GselfGop = [GselfGop; selfGop];
    GreussGop = [GreussGop; reussGop];
    GvoigtGop = [GvoigtGop; voigtGop];
    %% poisson's ration
    nulf = [nulf; f*nu1 + (1-f)*nu0];
    Gmorinulf = [Gmorinulf; morinulf];
    Gselfnulf = [Gselfnulf; selfnulf];
    Greussnulf = [Greussnulf; reussnulf];
    Gvoigtnulf = [Gvoigtnulf; voigtnulf];
end

%% Plots
p1 = round(length(volfrac)/6);
p2 = round(length(volfrac)*2/6);
p3 = round(length(volfrac)*3/6);
p4 = round(length(volfrac)*4/6);
p5 = round(length(volfrac)*5/6);

figure
set(gca,'FontSize',13)
hold on
plot(volfrac',Ecl,'MarkerSize',17,'LineWidth',1.5,'DisplayName','E_{cl}')
plot(volfrac',GmoriEcl,'MarkerSize',17,'LineWidth',1.5,'DisplayName','mori-tanaka')
plot(volfrac',GselfEcl,'MarkerSize',17,'LineWidth',1.5,'DisplayName','self-consistent')
plot(volfrac',GreussEcl,'MarkerSize',17,'LineWidth',1.5,'DisplayName','reuss')
plot(volfrac',GvoigtEcl,'MarkerSize',17,'LineWidth',1.5,'DisplayName','voigt')
text(volfrac(p1),Ecl(p1),'\leftarrow E_{cl}')
text(volfrac(p2),GmoriEcl(p2),'\leftarrow mori-tanaka')
text(volfrac(p3),GselfEcl(p3),'\leftarrow self-consistent')
text(volfrac(p4),GreussEcl(p4),'\leftarrow reuss')
text(volfrac(p5),GvoigtEcl(p5),'\leftarrow voigt')
xlabel('volume fraction - f')
ylabel('E-fibre [Pa]')
legend('Location','northwest')
box on
pbaspect([1 1 1]);
set(gcf,'units','pixels','position',[100 100 600 600]);
title("Young's modulus in fiber direction")

figure
set(gca,'FontSize',13)
hold on
plot(volfrac',Ect,'MarkerSize',17,'LineWidth',1.5,'DisplayName','E_{ct}')
plot(volfrac',GmoriEct,'MarkerSize',17,'LineWidth',1.5,'DisplayName','mori-tanaka')
plot(volfrac',GselfEct,'MarkerSize',17,'LineWidth',1.5,'DisplayName','self-consistent')
plot(volfrac',GreussEct,'MarkerSize',17,'LineWidth',1.5,'DisplayName','reuss')
plot(volfrac',GvoigtEct,'MarkerSize',17,'LineWidth',1.5,'DisplayName','voigt')
text(volfrac(p1),Ect(p1),'\leftarrow E_{ct}')
text(volfrac(p2),GmoriEct(p2),'\leftarrow mori-tanaka')
text(volfrac(p3),GselfEct(p3),'\leftarrow self-consistent')
text(volfrac(p4),GreussEct(p4),'\leftarrow reuss')
text(volfrac(p5),GvoigtEct(p5),'\leftarrow voigt')
xlabel('volume fraction - f')
ylabel('E-transverse [Pa]')
legend('Location','northwest')
box on
pbaspect([1 1 1]);
set(gcf,'units','pixels','position',[100 100 600 600]);
title("Young's modulus transverse to fiber direction")


figure
set(gca,'FontSize',13)
hold on
plot(volfrac',Gip,'MarkerSize',17,'LineWidth',1.5,'DisplayName','G_{ip}')
plot(volfrac',GmoriGip,'MarkerSize',17,'LineWidth',1.5,'DisplayName','mori-tanaka')
plot(volfrac',GselfGip,'MarkerSize',17,'LineWidth',1.5,'DisplayName','self-consistent')
plot(volfrac',GreussGip,'MarkerSize',17,'LineWidth',1.5,'DisplayName','reuss')
plot(volfrac',GvoigtGip,'MarkerSize',17,'LineWidth',1.5,'DisplayName','voigt')
text(volfrac(p1),Gip(p1),'\leftarrow G_{ip}')
text(volfrac(p2),GmoriGip(p2),'\leftarrow mori-tanaka')
text(volfrac(p3),GselfGip(p3),'\leftarrow self-consistent')
text(volfrac(p4),GreussGip(p4),'\leftarrow reuss')
text(volfrac(p5),GvoigtGip(p5),'\leftarrow voigt')
xlabel('volume fraction - f')
ylabel('G-inplane [Pa]')
legend('Location','northwest')
box on
pbaspect([1 1 1]);
set(gcf,'units','pixels','position',[100 100 600 600]);
title("Shear modulus in-plane")


figure
set(gca,'FontSize',13)
hold on
plot(volfrac',Gop,'MarkerSize',17,'LineWidth',1.5,'DisplayName','G_{op}')
plot(volfrac',GmoriGop,'MarkerSize',17,'LineWidth',1.5,'DisplayName','mori-tanaka')
plot(volfrac',GselfGop,'MarkerSize',17,'LineWidth',1.5,'DisplayName','self-consistent')
plot(volfrac',GreussGop,'MarkerSize',17,'LineWidth',1.5,'DisplayName','reuss')
plot(volfrac',GvoigtGop,'MarkerSize',17,'LineWidth',1.5,'DisplayName','voigt')
text(volfrac(p1),Gop(p1),'\leftarrow G_{op}')
text(volfrac(p2),GmoriGop(p2),'\leftarrow mori-tanaka')
text(volfrac(p3),GselfGop(p3),'\leftarrow self-consistent')
text(volfrac(p4),GreussGop(p4),'\leftarrow reuss')
text(volfrac(p5),GvoigtGop(p5),'\leftarrow voigt')
xlabel('volume fraction - f')
ylabel('G-outplane [Pa]')
legend('Location','northwest')
box on
pbaspect([1 1 1]);
set(gcf,'units','pixels','position',[100 100 600 600]);
title("Shear modulus out-plane")

figure
set(gca,'FontSize',13)
hold on
plot(volfrac',nulf,'MarkerSize',17,'LineWidth',1.5,'DisplayName','\nu_{lf}')
plot(volfrac',Gmorinulf,'MarkerSize',17,'LineWidth',1.5,'DisplayName','mori-tanaka')
plot(volfrac',Gselfnulf,'MarkerSize',17,'LineWidth',1.5,'DisplayName','self-consistent')
plot(volfrac',Greussnulf,'MarkerSize',17,'LineWidth',1.5,'DisplayName','reuss')
plot(volfrac',Gvoigtnulf,'MarkerSize',17,'LineWidth',1.5,'DisplayName','voigt')
text(volfrac(p1),nulf(p1),'\leftarrow \nu_{lf}')
text(volfrac(p2),Gmorinulf(p2),'\leftarrow mori-tanaka')
text(volfrac(p3),Gselfnulf(p3),'\leftarrow self-consistent')
text(volfrac(p4),Greussnulf(p4),'\leftarrow reuss')
text(volfrac(p5),Gvoigtnulf(p5),'\leftarrow voigt')
xlabel('volume fraction - f')
ylabel("poisson's ration")
legend('Location','northwest')
box on
pbaspect([1 1 1]);
set(gcf,'units','pixels','position',[100 100 600 600]);
title("Poisson's ratio")
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