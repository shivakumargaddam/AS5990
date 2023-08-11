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
texture = importdata("textures/random_30000.txt").data(:,1:3);
grainshape = [1 1 1];
% grainshape = [5 1 0.2];
tolerance = 1;
% C111244_Cu = [168e9,121.4e9,75.4e9];
% Cijklccs = Polycrystal.Cijklccsgen(C111244_Cu);
C111244_Fe = [231.4e9,134.7e9,116.4e9];
Cijklccs = Polycrystal.Cijklccsgen(C111244_Fe);

%% main
% Cvoigt4 = Polycrystal.Voigt4(Cijklccs,texture);
% fprintf('Voigt bound:  %s\n',Polycrystal.show9(Cvoigt4/1e9));


Cvoigt2 = Polycrystal.Voigt2(Cijklccs,texture);
Polycrystal.ZenerRat(Cvoigt2)
fprintf('Voigt bound:  %s\n\n',Polycrystal.show9(Cvoigt2/1e9));

Creuss2 = Polycrystal.Reuss2(Cijklccs,texture);
Polycrystal.ZenerRat(Creuss2)
fprintf('Reuss bound:  %s\n\n',Polycrystal.show9(Creuss2/1e9));

Chill2 = 0.5*(Cvoigt2+Creuss2);
Polycrystal.ZenerRat(Chill2)
fprintf('Hill bound:  %s\n\n',Polycrystal.show9(Chill2/1e9));

CguessHill = Polycrystal.Voigt2ijkl(Chill2);
[Csc,iter] = Polycrystal.SelfCons(CguessHill,Cijklccs,grainshape,texture,tolerance);
Polycrystal.ZenerRat(Csc)
fprintf('Self-consistent:  %s\n',Polycrystal.show9(Csc/1e9));

% Csc_ = Polycrystal.autoSelfCons(Cijklccs,texture);

Polycrystal.PlotDYMrol(Csc/1e9)
xlim([0,90])
ylim([100,200])
hold on
Polycrystal.PlotDYMrol(Creuss2/1e9)
Polycrystal.PlotDYMrol(Cvoigt2/1e9)



















