%% Visualize a Particle Set
close all; clearvars -except u_org v_org w_org g;
FLOW_X = 1; FLOW_Y = 2; FLOW_Z = 3; % Flow direction
solidColor = [0.5 0.5 0.5];
spaceColor = [0.95 0.95 0.95];

sample = '../data/Berea'; % sample MAT file name
particleSet = [sample,'ParticleSetA']; % particle Set MAT file name

gN = 400;
gFlow = FLOW_Z;
dx = 5.3;

%% load
if ~(exist('g','var') && exist('u_org','var') && exist('v_org','var') && exist('w_org','var'))
    fprintf(sprintf('Load sample %s...',sample));   
    load(sample);
    fprintf(sprintf('Done!\n')); 
end

load(particleSet); % load Particles (xpt,ypt,zpt)

%% Grab data
ipt = ceil (xpt / dx);
jpt = ceil (ypt / dx);
kpt = ceil (zpt / dx);

iFr = min(ipt); iTo = max(ipt);
jFr = min(jpt); jTo = max(jpt);
kFr = min(kpt); kTo = max(kpt);
Nparticle = length(ipt);

%% Define Surfaces
gS=NaN(gN,gN,gN);
gS(iFr,jFr:jTo,kFr:kTo)=g(iFr,jFr:jTo,kFr:kTo);
gS(iTo,jFr:jTo,kFr:kTo)=g(iTo,jFr:jTo,kFr:kTo);
gS(iFr:iTo,jFr,kFr:kTo)=g(iFr:iTo,jFr,kFr:kTo);
gS(iFr:iTo,jTo,kFr:kTo)=g(iFr:iTo,jTo,kFr:kTo);
gS(iFr:iTo,jFr:jTo,kFr)=g(iFr:iTo,jFr:jTo,kFr);
gS(iFr:iTo,jFr:jTo,kTo)=g(iFr:iTo,jFr:jTo,kTo);
[xP1, yP1, zP1] = ind2sub(size(gS), find(gS==1));
[xP0, yP0, zP0] = ind2sub(size(gS), find(gS==0));

%% Identify Particles on Surfaces
ip = zeros(size(ipt)); jp = zeros(size(ipt)); kp = zeros(size(ipt));
counter = 0;
for i=1:Nparticle
    if ipt(i)==iFr || jpt(i)==jFr || kpt(i)==kFr || ipt(i)==iTo || jpt(i)==jTo || kpt(i)==kTo       
        counter = counter + 1;
        ip(counter) = ipt(i);
        jp(counter) = jpt(i);
        kp(counter) = kpt(i);          
    end
end
ip(counter+1:end)=[];jp(counter+1:end)=[];kp(counter+1:end)=[];

%% Plot 
figure('Name','Visualize Injected Particles');
plot3(xP1, yP1, zP1, '.','Color',solidColor);
hold on;
plot3(xP0, yP0, zP0, '.','Color',spaceColor);
plot3(ip, jp, kp, 'ro','MarkerSize',5);
title(sprintf('%s, %s \n(Nparticle=%d)',sample,particleSet,Nparticle))
axis square;
axis([0 gN 0 gN 0 gN]);