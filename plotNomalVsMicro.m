%% plotNormalVsMicro
close all; clearvars -except u_org v_org w_org g;

load('Result/test10');

%% plot
figure('Name','NormalVsMicro'); hold on;

yyaxis left;
h0 = plot(varianceTime,insolid,'m--');

yyaxis right;
h1 = plot(varianceTime,dispersion1/Dm,'b-');
h2 = plot(varianceTime,dispersion2/Dm,'r-');

%% legends
txt0 = 'Inside Solid Percentage';
txt1 = 'RDC (normal)';
txt2 = 'RDC (with microporosity)';
legend([h0 h1 h2],{txt0,txt1,txt2},'Location','South');
xlabel('Timestep');
yyaxis left; ylabel('Inside Solid Percentage','Color','k');
yyaxis right; ylabel('Reduced Dispersion Coefficient','Color','b');
title(sprintf('%s, Nparticle=%d\nNpe=%.2f, avgAdvec=%.4f, jump=%.4f',...
    sample, Nparticle, Npe, max([uavg,vavg,wavg])*dt, zeta));

