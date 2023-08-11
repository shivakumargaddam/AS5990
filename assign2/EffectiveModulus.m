clear all;

% Epoxy
E0 = 3.5;
nu0 = 0.3;
mu0 = E0/2/(1+nu0);
k0 = E0/3/(1-2*nu0);

% Glass
E1 = 84;
nu1 = 0.2;
mu1 = E1/2/(1+nu1);
k1 = E1/3/(1-2*nu1);

f = 0:0.001:1;

% Voigt
muv = (1-f)*mu0+f*mu1;
kv = (1-f)*k0+f*k1;
nuv = (3*kv-2*muv)./(6*kv+2*muv);

% Reuss
mur = 1./((1-f)/mu0+f/mu1);
kr = 1./((1-f)/k0+f/k1);
nur = (3*kr-2*mur)./(6*kr+2*mur);

% Hashin-Shtrikman lower bound
s1 = (8-10*nu0)/15/(1-nu0);
s2 = (1+nu0)/3/(1-nu0);
a1 = 1/(1+s1*(mu1/mu0-1));
a2 = 1/(1+s2*(k1/k0-1));

muhsl = ((1-f)*mu0+f*mu1*a1)./(1-f+f*a1);
khsl = ((1-f)*k0+f*k1*a2)./(1-f+f*a2);
nuhsl = (3*khsl-2*muhsl)./(6*khsl+2*muhsl);

% Hashin-Shtrikman upper bound
s1 = (8-10*nu1)/15/(1-nu1);
s2 = (1+nu1)/3/(1-nu1);
a1 = 1/(1+s1*(mu0/mu1-1));
a2 = 1/(1+s2*(k0/k1-1));

muhsu = ((1-f)*mu0*a1+f*mu1)./((1-f)*a1+f);
khsu = ((1-f)*k0*a2+f*k1)./((1-f)*a2+f);
nuhsu = (3*khsu-2*muhsu)./(6*khsu+2*muhsu);

% Self-consistent
mubar = muhsl;
kbar = khsl;
nubar = (3*kbar-2*mubar)./(6*kbar+2*mubar);

for i=1:100
    mubarprev = mubar;
    kbarprev = kbar;

    s1bar = (8-10*nubar)/15./(1-nubar);
    s2bar = (1+nubar)/3./(1-nubar);
    a1bar = 1./(1+s1bar.*(mu1./mubar-1));
    a2bar = 1./(1+s2bar.*(k1./kbar-1));

    mubar = mu0 + f.*(mu1-mu0).*a1bar;
    kbar = k0 + f.*(k1-k0).*a2bar;
    nubar = (3*kbar-2*mubar)./(6*kbar+2*mubar);

    if (norm((mubar-mubarprev)/mu0) < 1.e-4)
        break
    end
end

figure;
hold on
plot(f,muv/mu0,'b.')
plot(f,mur/mu0,'r.')
plot(f,muhsu/mu0,'b-')
plot(f,muhsl/mu0,'r-')
plot(f,mubar/mu0,'k-')

figure;
hold on
plot(f,nuv,'b.')
plot(f,nur,'r.')
plot(f,nuhsu,'b-')
plot(f,nuhsl,'r-')
plot(f,nubar,'k-')