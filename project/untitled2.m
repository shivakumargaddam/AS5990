clear
clc

% az = deg2rad(0:1:360);
% el = deg2rad(-90:1:90);
% len = length(az)*length(el);
% xcoord = zeros(len,1);
% ycoord = zeros(len,1);
% zcoord = zeros(len,1);
% i = 1;
% for iaz = az
%     for iel = el
%         [xcoord(i),ycoord(i),zcoord(i)] = sph2cart(iaz,iel,1);
%         i = i+1;
%     end
% end
% r = xcoord + ycoord + zcoord;
% surf(xcoord,ycoord,zcoord,r)

[xcoord,ycoord,zcoord] = sphere(100);
r = xcoord + ycoord + zcoord;
surf(xcoord,ycoord,zcoord,"LineStyle",'none')
xlabel x
ylabel y
zlabel z
axis equal
colormap jet
colorbar
