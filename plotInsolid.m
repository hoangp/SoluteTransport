%% plotInsolid
close all; clearvars -except u_org v_org w_org g;

load('Result/test10');

%%
figure('Name','Insolid');

yyaxis left;
h0 = plot(varianceTime,insolid,'k-'); hold on;

yyaxis right;
h1 = plot(varianceTime,dispersion1/Dm,'r');

txt0 = 'Inside Solid Percentage';
txt1 = 'RDC (with microporosity)';
legend([h0 h1],{txt0,txt1},'Location','NorthWest');

title(sprintf('(Npe=%.2f) (avg advection=%.2f) (random jump=%.2f)\n(number particle=%d) (sample=%s)',Npe,uavg*dt,zeta,Nparticle,sample));
xlabel('Time(s)');
yyaxis left; ylabel('Inside Solid Percentage','Color','k');
yyaxis right; ylabel('Reduced Dispersion Coefficient','Color','b');