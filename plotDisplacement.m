%% plotDiplacement
close all; clearvars -except u_org v_org w_org g;

load('Result/test10');
ipt = ceil(displacement1/dx);

%%
figure('Name','Displacement');

% x axis: Displacement
iptMin = min(ipt); iptMax = max(ipt);
x_values = iptMin:1:iptMax;

% y axis (bar): Number of Particle
nbins = iptMax - iptMin + 1;
[N1,~] = histcounts(ipt,nbins);

% y axis (line): fit distribution line
pd1 = fitdist(ipt,'Normal');
y1 = pdf(pd1,x_values)*Nparticle;

% ploting
x_values = x_values*dx/1000;
bar(x_values,N1,'FaceColor',[0.5,0.5,1],'EdgeColor','none');
hold on;
plot(x_values,y1,'b','LineWidth',1);

% labels
title(sprintf('(Npe=%.2f) (average advection=%.2f) (random jump=%.2f)\n(number particle=%d) (sample=%s)',Npe,uavg*dt,zeta,Nparticle,sample));
xlabel('Displacement (mm)');
ylabel('Number of particles');
axis_xmin = min(x_values);
axis_xmax = max(x_values);
axis_ymax = max(N1);
axis ([axis_xmin axis_xmax 0 axis_ymax]);
