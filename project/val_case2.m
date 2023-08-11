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
ntextures = ["textures/random_1000.txt" "textures/random_10000.txt" "textures/random_20000.txt" "textures/random_30000.txt"  "textures/random_60000.txt"];
YoungsMod = 130e9;
PoissonR = 0.34;
grainshape = [1 1 1];
tolerance = 1;
Cijklccs = Polycrystal.Voigt2ijkl(Polycrystal.isoCij(130e9,0.34));

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
clear Csc2 Cvoigt2 Creuss2 CguessHill
save val_case2_results.mat

%%  Analysis
clear
clc
load val_case2_results.mat

TAI = zeros(length(ntextures),1);
grains = [1000 10000 20000 30000 60000];
lambda1 = zeros(5,length(ntextures));
Vlambda1 = zeros(5,length(ntextures));
Rlambda1 = zeros(5,length(ntextures));
for itexture = 1:length(ntextures)
    disp(' ');
    disp(ntextures(itexture))
    fprintf("The eigen values for the Upper/Voigt bound:\n")
    disp(eigs(Polycrystal.Voigt2Mandel(Cvoigt(:,:,itexture)))'/1e9)
    fprintf("The eigen values for the Lower/Reuss bound:\n")
    disp(eigs(Polycrystal.Voigt2Mandel(Creuss(:,:,itexture)))'/1e9)
    fprintf("The eigen values for the Hill bound:\n")
    disp(eigs(Polycrystal.Voigt2Mandel(Chill(:,:,itexture)))'/1e9)

    TAI(itexture) = Polycrystal.TenAniInd(Csc(:,:,itexture));
    fprintf("%d: Tensor anisotropy index - %f\n",itexture,TAI(itexture))
    evalues = eigs(Polycrystal.Voigt2Mandel(Csc(:,:,itexture)))/1e9;
    lambda1(:,itexture) = evalues(2:end);
    Vevalues = eigs(Polycrystal.Voigt2Mandel(Cvoigt(:,:,itexture)))/1e9;
    Vlambda1(:,itexture) = Vevalues(2:end);
    Revalues = eigs(Polycrystal.Voigt2Mandel(Creuss(:,:,itexture)))/1e9;
    Rlambda1(:,itexture) = Revalues(2:end);
    fprintf("range: %f\n",max(evalues(2:end))-min(evalues(2:end)))
    disp(' ')
end

%% Index Plot
indexplot = figure('Name','Tensor Anisotropy Index','NumberTitle','off');
plot(grains,TAI,'.-','LineWidth',1.5,'MarkerSize',20)
yline(1,'-','Isotropic','LabelHorizontalAlignment','left','LineWidth',1);
grid off
box on
xticks([1000 10000 20000 30000 40000 50000 60000])
xticklabels({'1000','10000','20000','30000','40000','50000','60000'})
xlim([0,65000])
ylim([0.9985,1.0255])
xlabel("no. of grains");
ylabel("TAI (isotropic - 1)");
pbaspect([1 1 1])
set(gca,'FontSize',12)
set(gca, 'color', 'none');                                                  % To remove background
set(gcf,'units','pixels','position',[1000 300 500 500]);                     % To change the size of the figure
title("Tensor Anisotropy Index (TAI) vs. no. of grains")

%% Boxplot
eigenplot = figure('Name','\lambda_1','NumberTitle','off');
boxchart(lambda1)
grid off
box on
xticklabels({'1000','10000','20000','30000','60000'})
xlabel("no. of grains");
ylabel("\lambda_1 in GPa");
pbaspect([1 1 1])
set(gca,'FontSize',12)
set(gca, 'color', 'none');                                                  % To remove background
set(gcf,'units','pixels','position',[1000 300 500 500]);                    % To change the size of the figure
title("Smallest eigen value (\lambda_1) vs. no. of grains")

%% Directional Young's Modulus
figure('Name',"Directional Young's Modulus",'NumberTitle','off');
Polycrystal.PlotDYM3d(Csc(:,:,1))
title('1000 grains')
c = colorbar;
c.Label.String = "Directional Young's Modulus in GPa";
maxlimits = get(c,'Limits');

figure
t = tiledlayout(2,2);
nexttile
Polycrystal.PlotDYM3d(Csc(:,:,2))
title('10000 grains')
clim(maxlimits)
colorbar('off')

nexttile
Polycrystal.PlotDYM3d(Csc(:,:,3))
title('20000 grains')
clim(maxlimits)
colorbar('off')

nexttile
Polycrystal.PlotDYM3d(Csc(:,:,4))
title('30000 grains')
clim(maxlimits)
colorbar('off')

nexttile
Polycrystal.PlotDYM3d(Csc(:,:,5))
title('60000 grains')
clim(maxlimits)
colorbar('off')

t.TileSpacing = 'compact';
t.Padding = 'compact';

%%
































