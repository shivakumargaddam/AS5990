clear
clc

cs = crystalSymmetry('cubic');
ss = specimenSymmetry('222');
npoints = 1000;
filepath = 'C:\Users\shiva\MATLAB Drive\AS5990\project\textures';

%% fibre texture
f = fibre.gamma(cs);
odf = fibreODF(f,'halfwidth',10*degree);
plot(odf,ss,'phi2',45*degree,'contourf','silent');
ori = discreteSample(odf,npoints);
% ori = calcOrientations(odf,npoints);
fname1 = fullfile(filepath,'goss_1000.txt');
export(ori,fname1,'Bunge')

%% random texture
h = Miller(1,1,1,cs);
ori = orientation.rand(npoints,cs);
odf = calcDensity(ori);
% figure('Name','random_40000');
plotPDF(odf,h,'antipodal')
mtexColorbar
fname2 = fullfile(filepath,'random_70000.txt');
export(ori,fname2,'Bunge')

%% copper texture
h = Miller(1,1,1,cs);
Cuori = orientation.byEuler(90*degree,35*degree,45*degree,'Bunge',cs);
odf = calcDensity(Cuori,'halfwidth',1*degree);
% ori = calcOrientations(odf,npoints);
ori = discreteSample(odf,npoints);
fname2 = fullfile(filepath,'copper_1000.txt');
export(ori,fname2,'Bunge')
plot(odf,ss,'phi2',45*degree,'contourf','silent');
plotPDF(odf,h,'antipodal')

%% modified rolling texture with halfwidth
h = Miller(1,1,1,cs);
roltex = importdata("rolling_1000.txt").data(:,1:3);
rolori = orientation.byEuler(roltex(:,1)*degree,roltex(:,2)*degree,roltex(:,3)*degree,'Bunge',cs);
odf = calcDensity(rolori,'halfwidth',1*degree);
% ori = calcOrientations(odf,npoints);
ori = discreteSample(odf,npoints);
fname2 = fullfile(filepath,'rollingmod_1000.txt');
export(ori,fname2,'Bunge')
plot(odf,ss,'phi2',45*degree,'contourf','silent');
plotPDF(odf,h,'antipodal')

%% modified rolling texture with intensity by multiplying
h = Miller(1,1,1,cs);
roltex = importdata("rolling_1000.txt").data(:,1:3);
rolori = orientation.byEuler(roltex(:,1)*degree,roltex(:,2)*degree,roltex(:,3)*degree,'Bunge',cs);
odf = calcDensity(rolori);
odf = odf*4;
% ori = calcOrientations(odf,npoints);
ori = discreteSample(odf,npoints);
fname2 = fullfile(filepath,'rollingmod_1000_x4.txt');
export(ori,fname2,'Bunge')
plot(odf,ss,'phi2',45*degree,'contourf','silent');
figure
plotPDF(odf,h,'antipodal')
%% Plot PDF
h = Miller(1,0,0,cs);
roltex = importdata("textures\FCC_rolling75p.txt").data(:,1:3);
rolori = orientation.byEuler(roltex(:,1)*degree,roltex(:,2)*degree,roltex(:,3)*degree,'Bunge',cs);
odf = calcDensity(rolori);
figure
plotPDF(odf,h,'antipodal')
mtexColorbar
%% Creating random texture for cubic
h = Miller(1,1,1,cs);
orientations = zeros(91^3,3);
iori = 1;
for iphi1 = 0:90
    for iPHI = 0:90
        for iphi2 = 0:90
            orientations(iori,:) = [iphi1,iPHI,iphi2];
            iori = iori+1;
        end
    end
end
random = orientation.byEuler(orientations(:,1)*degree,orientations(:,2)*degree,orientations(:,3)*degree,'Bunge',cs);
odf = calcDensity(random);
plotPDF(odf,h,'antipodal')
mtexColorbar
