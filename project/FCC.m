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
C111244_Cu = [168e9,121.4e9,75.4e9];
grainshape = [5 1 0.2];
Cijklccs = Polycrystal.Cijklccsgen(C111244_Cu);

%% Calculate Effective Stiffness Tensors
% Rolling spherical
texture = importdata("textures/FCC_rolling25p.txt").data(:,1:3);
CscFccRol25p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/FCC_rolling50p.txt").data(:,1:3);
CscFccRol50p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/FCC_rolling75p.txt").data(:,1:3);
CscFccRol75p = Polycrystal.autoSelfCons(Cijklccs,texture);

% Rolling elongated and flattened
texture = importdata("textures/FCC_rolling25p.txt").data(:,1:3);
CscFccRol25pE = Polycrystal.autoSelfCons(Cijklccs,texture,grainshape);
texture = importdata("textures/FCC_rolling50p.txt").data(:,1:3);
CscFccRol50pE = Polycrystal.autoSelfCons(Cijklccs,texture,grainshape);
texture = importdata("textures/FCC_rolling75p.txt").data(:,1:3);
CscFccRol75pE = Polycrystal.autoSelfCons(Cijklccs,texture,grainshape);

% Uniaxial Tension
texture = importdata("textures/FCC_tension25p.txt").data(:,1:3);
CscFccTen25p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/FCC_tension50p.txt").data(:,1:3);
CscFccTen50p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/FCC_tension75p.txt").data(:,1:3);
CscFccTen75p = Polycrystal.autoSelfCons(Cijklccs,texture);

% Shear
texture = importdata("textures/FCC_shear25p.txt").data(:,1:3);
CscFccShr25p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/FCC_shear50p.txt").data(:,1:3);
CscFccShr50p = Polycrystal.autoSelfCons(Cijklccs,texture);
texture = importdata("textures/FCC_shear75p.txt").data(:,1:3);
CscFccShr75p = Polycrystal.autoSelfCons(Cijklccs,texture);

save FCC_results.mat

%% Analysis - rolling with spherical grains
clear
clc
load FCC_results.mat

Ro75 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccRol75pE)
title('75% Rolling')
c = colorbar;
c.Label.String = "Directional Young's Modulus in GPa";
maxlimits = get(c,'Limits');
grid off
axis off
exportgraphics(Ro75,'CscFccRol75pE.png','Resolution',500);

Ro50 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccRol50pE)
title('50% Rolling')
clim(maxlimits)
colorbar('off')
exportgraphics(Ro50,'CscFccRol50pE.png','Resolution',500);

Ro25 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccRol25pE)
title('25% Rolling')
clim(maxlimits)
colorbar('off')
exportgraphics(Ro25,'CscFccRol25pE.png','Resolution',500);

%% Analysis - rolling with spherical grains
Ro75 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccRol75p)
title('75% Rolling')
clim(maxlimits)
c = colorbar;
c.Label.String = "Directional Young's Modulus in GPa";
grid off
axis off
exportgraphics(Ro75,'CscFccRol75p.png','Resolution',500);

Ro50 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccRol50p)
title('50% Rolling')
clim(maxlimits)
colorbar('off')
exportgraphics(Ro50,'CscFccRol50p.png','Resolution',500);

Ro25 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccRol25p)
title('25% Rolling')
clim(maxlimits)
colorbar('off')
exportgraphics(Ro25,'CscFccRol25p.png','Resolution',500);

%% Analysis - uniaxial tension with spherical grains
Te75 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccTen75p)
title('75% Uniaxial Tension')
c = colorbar;
c.Label.String = "Directional Young's Modulus in GPa";
maxlimits = get(c,'Limits');
grid off
axis off
exportgraphics(Te75,'CscFccTen75p.png','Resolution',500);

Te50 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccTen50p)
title('50% Uniaxial Tension')
clim(maxlimits)
colorbar('off')
exportgraphics(Te50,'CscFccTen50p.png','Resolution',500);

Te25 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccTen25p)
title('25% Uniaxial Tension')
clim(maxlimits)
colorbar('off')
exportgraphics(Te25,'CscFccTen25p.png','Resolution',500);
%% Analysis - Shear with spherical grains
sh75 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccShr75p)
title('75% Simple shear')
c = colorbar;
c.Label.String = "Directional Young's Modulus in GPa";
maxlimits = get(c,'Limits');
grid off
axis off
exportgraphics(sh75,'CscFccShr75p.png','Resolution',500);

sh50 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccShr50p)
title('50% Shear')
clim(maxlimits)
colorbar('off')
exportgraphics(sh50,'CscFccShr50p.png','Resolution',500);

sh25 = figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(CscFccShr25p)
title('25% Shear')
clim(maxlimits)
colorbar('off')
exportgraphics(sh25,'CscFccShr25p.png','Resolution',500);





















