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
ntextures = ["textures/random_1000.txt" "textures/random_10000.txt" "textures/random_20000.txt" "textures/random_30000.txt" "textures/random_40000.txt" "textures/random_60000.txt" "textures/random_100000.txt"];
C111244_Cu = [168e9,121.4e9,75.4e9];
grainshape = [1 1 1];
tolerance = 1;
Cijklccs = Polycrystal.Cijklccsgen(C111244_Cu);

Cvoigt = zeros(6,6,length(ntextures));
Creuss = zeros(6,6,length(ntextures));
Chill = zeros(6,6,length(ntextures));
Csc = zeros(6,6,length(ntextures));

%% main
for itexture = 1:length(ntextures)
    texture = importdata(ntextures(itexture)).data(:,1:3);
    % Upper bound
    Cvoigt2 = Polycrystal.Voigt2(Cijklccs,texture);

    % Lower bound
    Creuss2 = Polycrystal.Reuss2(Cijklccs,texture);

    % Hill bound
    Chill2 = 0.5*(Cvoigt2+Creuss2);

    % Self-consistent
    CguessHill = Polycrystal.Voigt2ijkl(Chill2);
    [Csc2,iter] = Polycrystal.SelfCons(CguessHill,Cijklccs,grainshape,texture,tolerance);
    
    Cvoigt(:,:,itexture) = Cvoigt2;
    Creuss(:,:,itexture) = Creuss2;
    Chill(:,:,itexture) = Chill2;
    Csc(:,:,itexture) = Csc2;
end
save case1_results.mat

%%  Analysis
clear
clc
load case1_results.mat
for itexture = 1:length(ntextures)
    disp(' ');
    disp(ntextures(itexture))
    fprintf("%d: Tensor anisotropy index - %f\n",itexture,Polycrystal.TenAniInd(Csc(:,:,itexture)))
    evalues = eigs(Csc(:,:,itexture));
    fprintf("range: %f\n",max(evalues(2:end))-min(evalues(2:end)))
    disp(' ')
end





































