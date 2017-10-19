%% plotNpeVsRDC
close all; clearvars -except u_org v_org w_org g;

resultName = 'Result/sandB_t40k_com';
run_list = 1:53;

maxSlope = 20;
fixTime = 0;

%% collect data
NPE = zeros(length(run_list),1);
RDC1 = zeros(length(run_list),1);
RDC2 = zeros(length(run_list),1);
c = 0;
for r = run_list
    c = c + 1;
    load([resultName,num2str(r)]);
    
    NPE(c) = Npe;
    if fixTime ~= 0
        tSlope = find(varianceTime >= fixTime,1,'first');
        RDC1(c) = dispersion1(tSlope)/Dm;
        RDC2(c) = dispersion2(tSlope)/Dm;
    elseif length(dispersion1)>=maxSlope
        RDC1(c) = dispersion1(maxSlope)/Dm;
        RDC2(c) = dispersion2(maxSlope)/Dm;
    else
        RDC1(c) = dispersion1(end)/Dm;
        RDC2(c) = dispersion2(end)/Dm;
    end
end

%% plot
figure('Name','NpeVsRDC');
xls = xlsread('SoluteTransport');
h0 = loglog(xls(:,1),xls(:,2),'k:','LineWidth',2); hold on;
h1 = loglog(NPE,RDC1,'bo-');
h2 = loglog(NPE,RDC2,'rx-');

%% legends
txt0 = 'Curve from the paper';
txt1 = 'Simulation results (normal)';
txt2 = 'Simulation results (with micro porosity)';
legend([h0 h1 h2],{txt0,txt1,txt2},'Location','NorthWest');
xlabel('Npe');
ylabel('Reduced Dispersion Coefficient');
title(sprintf('%s, Nparticle=%d',sample,Nparticle));
