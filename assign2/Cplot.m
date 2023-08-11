close all
clear
clc

% Material - Epoxy
E0 = 3.35;
nu0 = 0.33;
mu0 = E0/(2*(1+nu0));
K0 = E0/(3*(1-2*nu0));
% nu0 = (3*K0-2*mu0)/2/(3*K0+mu0);

% Material - Glass Fibre
E1 = 72.3;
nu1 = 0.21;
mu1 = E1/(2*(1+nu1));
K1 = E1/(3*(1-2*nu1));

f = 0:0.0001:1;
figure
hold on

% Hill bounds
reuss_mu = 1./((1-f)/mu0 + f/mu1);
voigt_mu = (1-f)*mu0 + f*mu1;

% Hashin-Shtrikman bounds
s10 = (8-10*nu0)/15/(1-nu0);
s20 = (1+nu0)/3/(1-nu0);
a10 = 1/(1+s10*(mu1/mu0-1));

mu_hsl = (((1-f)*mu0)+f*mu1*a10)./(1-f+f*a10);

s11 = (8-10*nu1)/15/(1-nu1);
s21 = (1+nu1)/3/(1-nu1);
a11 = 1/(1+s11*(mu0/mu1-1));
mu_hsu = (((1-f)*mu0*a11)+f*mu1)./((1-f)*a11+f);

% Mori-Tanaka Estimation

% Self-Consistent Estimation

plot(f,reuss_mu/mu0,'DisplayName','Reuss',LineWidth=2)
plot(f,voigt_mu/mu0,'DisplayName','Voigt',LineWidth=2)
plot(f,mu_hsl/mu0,'DisplayName','HS lower',LineWidth=2)
plot(f,mu_hsu/mu0,'DisplayName','HS upper',LineWidth=2)

pbaspect([1 1 1])
ylim([1 inf])
legend