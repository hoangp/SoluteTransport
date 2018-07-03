clear;

SAVE_TIMESTEP = 1000;
Nrun = 21;
Dm = 1000;

result_folder = '../result/';

result_name1 = 'BereaA0p'; 
%result_name2 = 'carbonateA1p'; 
%result_name3 = 'carbonateA0p'; 

result_list = 1:10;

Ntimestep = length(result_list)*SAVE_TIMESTEP;

fixTime = 10000;
fixSlope = 9000;

%% collect results
result_filename1 = [result_folder,result_name1];
load([result_filename1,num2str(result_list(1))]);

Nparticle = length(particleX(:,1));
NPE1 = run_Npe;

%% result 1
allVariance1 = zeros(Ntimestep,length(result_list));
for r = result_list
    load([result_filename1,num2str(r)]);
    for i = 1:Nrun
        allVariance1(((r-1)*SAVE_TIMESTEP + 1):(r*SAVE_TIMESTEP),i)=variance(:,i);
    end
end

RDC1 = zeros(1,Nrun);
for r = 1:Nrun
    RDC1(r) = 0.5*(allVariance1(fixTime,r)-allVariance1(fixTime-fixSlope,r))/(fixSlope*dt)/Dm;
end

%% result2
% result_filename2 = [result_folder,result_name2];
% load([result_filename2,num2str(result_list(1))]);
% 
% allVariance2 = zeros(Ntimestep,length(result_list));
% for r = result_list
%     load([result_filename2,num2str(r)]);
%     for i = 1:Nrun
%         allVariance2(((r-1)*SAVE_TIMESTEP + 1):(r*SAVE_TIMESTEP),i)=variance(:,i);
%     end
% end
% 
% RDC2 = zeros(1,Nrun);
% for r = 1:Nrun
%     RDC2(r) = 0.5*(allVariance2(fixTime,r)-allVariance2(fixTime-fixSlope,r))/(fixSlope*dt)/Dm;
% end
% 

%% result3
% result_filename3 = [result_folder,result_name3];
% load([result_filename3,num2str(result_list(1))]);
% 
% allVariance3 = zeros(Ntimestep,length(result_list));
% for r = result_list
%     load([result_filename3,num2str(r)]);
%     for i = 1:Nrun
%         allVariance3(((r-1)*SAVE_TIMESTEP + 1):(r*SAVE_TIMESTEP),i)=variance(:,i);
%     end
% end
% 
% RDC3 = zeros(1,Nrun);
% for r = 1:Nrun
%     RDC3(r) = 0.5*(allVariance3(fixTime,r)-allVariance3(fixTime-fixSlope,r))/(fixSlope*dt)/Dm;
% end

%% plot
close all;
figure('Name','NpeVsRDC'); grid on; 
%yyaxis left;
load([data_folder,'prev_paper']);
prev_paper(find(prev_paper(:,2)>max(NPE1),1,'first'):end,:)=[];
%xls(1:2,:)=NaN;
h0 = loglog(prev_paper(:,1),prev_paper(:,2),'k:','LineWidth',2); hold on;
h1 = loglog(NPE1,RDC1,'k-'); hold on;
%h2 = loglog(NPE1,RDC2,'b-'); hold on;
%h3 = loglog(NPE1,RDC3,'r-'); hold on;
%ht = loglog(NPE1,TIME,'ko');
ylabel('Reduced Dispersion Coefficient');


%% legends
txt0 = 'Mostaghimi et al. (2012)';
txt1 = result_name1;
%txt2 = result_name2;
%txt3 = result_name3;
txtt = 'Number of Timestep';

%legend([h0 h1 ht],{txt0,txt1,txtt},'Location','East');

legend([h0 h1],{txt0,txt1},'Location','North');
%legend([h0 h1 h2],{txt0,txt1,txt2},'Location','North');
%legend([h0 h1 h2 h3],{txt0,txt1,txt2,txt3},'Location','North');

xlabel('Npe');
title(sprintf('Nparticle=%d Ntimestep=%d slope=%d',Nparticle,Ntimestep,fixSlope));

xlim([0.001 1000]);
ylim([0.1 10000]);
%yyaxis right; ylim([0 1]);
grid on;
