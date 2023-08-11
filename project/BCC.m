clear all
close all
clc

%% PROGRAM TITLE
disp('==================================================================');
disp('                        Project - AS5990                          ');
disp('      Effective elastic properties of textured polycrystals       ');
disp('                  Shiva Kumar Gaddam - MM22D014                   ');
disp('==================================================================');
pause

%% Input
C111244_Fe = [231.4e9,134.7e9,116.4e9];
grainshape = [5 1 0.2];
Cijklccs = Polycrystal.Cijklccsgen(C111244_Fe);

%% Calculate Effective Stiffness Tensors
% Rolling spherical
texture = importdata("textures/BCC_rolling25p.txt").data(:,1:3);
CscBccRol25p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/BCC_rolling50p.txt").data(:,1:3);
CscBccRol50p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/BCC_rolling75p.txt").data(:,1:3);
CscBccRol75p = Polycrystal.autoSelfCons(Cijklccs,texture);

% Rolling elongated and flattened
texture = importdata("textures/BCC_rolling25p.txt").data(:,1:3);
CscBccRol25pE = Polycrystal.autoSelfCons(Cijklccs,texture,grainshape);
texture = importdata("textures/BCC_rolling50p.txt").data(:,1:3);
CscBccRol50pE = Polycrystal.autoSelfCons(Cijklccs,texture,grainshape);
texture = importdata("textures/BCC_rolling75p.txt").data(:,1:3);
CscBccRol75pE = Polycrystal.autoSelfCons(Cijklccs,texture,grainshape);

% Uniaxial Tension
texture = importdata("textures/BCC_tension25p.txt").data(:,1:3);
CscBccTen25p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/BCC_tension50p.txt").data(:,1:3);
CscBccTen50p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/BCC_tension75p.txt").data(:,1:3);
CscBccTen75p = Polycrystal.autoSelfCons(Cijklccs,texture);

% Shear
texture = importdata("textures/BCC_shear25p.txt").data(:,1:3);
CscBccShr25p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/BCC_shear50p.txt").data(:,1:3);
CscBccShr50p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/BCC_shear75p.txt").data(:,1:3);
CscBccShr75p = Polycrystal.autoSelfCons(Cijklccs,texture);

save BCC_results.mat

%% Analysis - rolling with spherical grains
clear
clc
load BCC_results.mat

Ro75 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccRol75pE)
title('75% Rolling')
c = colorbar;
c.Label.String = "Directional Young's Modulus in GPa";
maxlimits = get(c,'Limits');
grid off
axis off
exportgraphics(Ro75,'CscBccRol75pE.png','Resolution',500);

Ro50 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccRol50pE)
title('50% Rolling')
clim(maxlimits)
colorbar('off')
grid off
axis off
exportgraphics(Ro50,'CscBccRol50pE.png','Resolution',500);

Ro25 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccRol25pE)
title('25% Rolling')
clim(maxlimits)
colorbar('off')
grid off
axis off
exportgraphics(Ro25,'CscBccRol25pE.png','Resolution',500);

%% Analysis - rolling with spherical grains
Ro75 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccRol75p)
title('75% Rolling')
clim(maxlimits)
c = colorbar;
c.Label.String = "Directional Young's Modulus in GPa";
grid off
axis off
exportgraphics(Ro75,'CscBccRol75p.png','Resolution',500);

Ro50 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccRol50p)
title('50% Rolling')
clim(maxlimits)
colorbar('off')
exportgraphics(Ro50,'CscBccRol50p.png','Resolution',500);

Ro25 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccRol25p)
title('25% Rolling')
clim(maxlimits)
colorbar('off')
exportgraphics(Ro25,'CscBccRol25p.png','Resolution',500);

%% Analysis - uniaxial tension with spherical grains
Te75 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccTen75p)
title('75% Uniaxial Tension')
c = colorbar;
c.Label.String = "Directional Young's Modulus in GPa";
maxlimits = get(c,'Limits');
grid off
axis off
exportgraphics(Te75,'CscBccTen75p.png','Resolution',500);

Te50 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccTen50p)
title('50% Uniaxial Tension')
clim(maxlimits)
colorbar('off')
exportgraphics(Te50,'CscBccTen50p.png','Resolution',500);

Te25 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccTen25p)
title('25% Uniaxial Tension')
clim(maxlimits)
colorbar('off')
exportgraphics(Te25,'CscBccTen25p.png','Resolution',500);

%% Analysis - shear with spherical grains
sh75 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccShr75p)
title('75% Simple shear')
c = colorbar;
c.Label.String = "Directional Young's Modulus in GPa";
maxlimits = get(c,'Limits');
grid off
axis off
exportgraphics(sh75,'CscBccShr75p.png','Resolution',500);

sh50 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccShr50p)
title('50% Shear')
clim(maxlimits)
colorbar('off')
exportgraphics(sh50,'CscBccShr50p.png','Resolution',500);

sh25 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscBccShr25p)
title('25% Shear')
clim(maxlimits)
colorbar('off')
exportgraphics(sh25,'CscBccShr25p.png','Resolution',500);




















