%% plotVariance
close all; clearvars -except u_org v_org w_org g;

load('Result/test10');

%%
figure('Name','Variance');

yyaxis left;
h0 = plot(varianceTime,variance1,'k-'); hold on;

yyaxis right;
h1 = plot(varianceTime,dispersion1/Dm,'r');

title(sprintf('(Npe=%.2f) (avg advection=%.2f) (random jump=%.2f)\n(number particle=%d) (sample=%s)',Npe,uavg*dt,zeta,Nparticle,sample));
xlabel('Time(s)');
yyaxis left; ylabel('Variance of Displacement','Color','k');
yyaxis right; ylabel('Reduced Dispersion Coefficient','Color','b');
legend([h0 h1],{'Variance of Displacement',...
    sprintf('Reduced Dispersion Coefficient (slope=%d)',slope)},...
    'Location','northwest');

