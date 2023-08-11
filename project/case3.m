clear all
close all
clc

%% PROGRAM TITLE
disp('==================================================================');
disp('                        Project - AS5990                          ');
disp('      Effective elastic properties of textured polycrystals       ');
disp('                  Shiva Kumar Gaddam - MM22D014                   ');
disp('==================================================================');

%% Input
texture = importdata("textures/ani_1000.txt").data(:,1:3);
C111244_Cu = [168e9,121.4e9,75.4e9];
grainshape = [1 1 1];
tolerance = 1;
Cijklccs = Polycrystal.Cijklccsgen(C111244_Cu);

Cvoigt = zeros(6,6);
Creuss = zeros(6,6);
Chill = zeros(6,6);
Csc = zeros(6,6);

%% main

% Upper bound
Cvoigt2 = Polycrystal.Voigt2(Cijklccs,texture);

% Lower bound
Creuss2 = Polycrystal.Reuss2(Cijklccs,texture);

% Hill bound
Chill2 = 0.5*(Cvoigt2+Creuss2);

% Self-consistent
CguessHill = Polycrystal.Voigt2ijkl(Chill2);
[Csc2,iter] = Polycrystal.SelfCons(CguessHill,Cijklccs,grainshape,texture,tolerance);

Cvoigt(:,:) = Cvoigt2;
Creuss(:,:) = Creuss2;
Chill(:,:) = Chill2;
Csc(:,:) = Csc2;

save case3_results.mat

%%  Analysis
load case3_results.mat

Polycrystal.PlotDYMrol(Csc/1e9)
xlim([0,90])
% ylim([100,200])
hold on
Polycrystal.PlotDYMrol(Creuss/1e9)
Polycrystal.PlotDYMrol(Cvoigt/1e9)









