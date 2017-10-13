%% VisualizeSample.m
% Visualize Sample Data by ploting 6 surfaces & velocity fields
close all; clearvars -except u_org v_org w_org g;

sample = 'Data/sand500Sample'; % sample MAT file name
SampleConstants; % load sample properties & constants

edgeFr = 2;
edgeTo = gN-2;

solidColor = [0.3 0.3 0.3];

%% Define Surfaces
gS=NaN(gN,gN,gN);
gS(edgeFr,:,:)=g(edgeFr,:,:);
gS(edgeTo,:,:)=g(edgeTo,:,:);
gS(:,edgeFr,:)=g(:,edgeFr,:);
gS(:,edgeTo,:)=g(:,edgeTo,:);
gS(:,:,edgeFr)=g(:,:,edgeFr);
gS(:,:,edgeTo)=g(:,:,edgeTo);
[xP1, yP1, zP1] = ind2sub(size(gS), find(gS==1));

%% Define Velocity Planes
vel = sqrt(u_org.^2 + w_org.^2 + v_org.^2)/1000;
[planeX1,X1x,X1y,X1z] = coordinates(vel(edgeFr,:,:),edgeFr,gN);
[planeX2,X2x,X2y,X2z] = coordinates(vel(edgeTo,:,:),edgeTo,gN);
[planeY1,Y1y,Y1x,Y1z] = coordinates(vel(:,edgeFr,:),edgeFr,gN);
[planeY2,Y2y,Y2x,Y2z] = coordinates(vel(:,edgeTo,:),edgeTo,gN);
[planeZ1,Z1z,Z1x,Z1y] = coordinates(vel(:,:,edgeFr),edgeFr,gN);
[planeZ2,Z2z,Z2x,Z2y] = coordinates(vel(:,:,edgeTo),edgeTo,gN);

%% Plot Surfaces
figure('Name','Visualize Sample Data');
plot3(xP1, yP1, zP1, '.','Color',solidColor);
hold on;

%% Plot Velocity Fields
scatter3(X1x,X1y,X1z,[],reshape(planeX1(planeX1~=0),1,[]),'.');
scatter3(X2x,X2y,X2z,[],reshape(planeX2(planeX2~=0),1,[]),'.');
scatter3(Y1x,Y1y,Y1z,[],reshape(planeY1(planeY1~=0),1,[]),'.');
scatter3(Y2x,Y2y,Y2z,[],reshape(planeY2(planeY2~=0),1,[]),'.');
scatter3(Z1x,Z1y,Z1z,[],reshape(planeZ1(planeZ1~=0),1,[]),'.');
scatter3(Z2x,Z2y,Z2z,[],reshape(planeZ2(planeZ2~=0),1,[]),'.');

%% Figure Config
title(sprintf('%s \n(Velocity field layers from %d to %d)',sample,edgeFr,edgeTo));
axis square;
axis([0 gN 0 gN 0 gN]);

%% Nested Function
function [vel2,x,y,z] = coordinates(vel1,loc,N)
    vel2 = squeeze(vel1);
    yv = reshape(meshgrid(1:N)',1,[]);
    zv = reshape(meshgrid(1:N),1,[]);
    y = yv(vel2~=0);
    z = zv(vel2~=0);
    x = ones(size(y))*loc;
end
