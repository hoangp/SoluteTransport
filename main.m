%% MultiSimulation.m
% Simulate 
%
% Date: 13/10/2017
% Author: Hoang Pham
%
close all;
clearvars -except u_org v_org w_org g;

sample = 'Data/geometry2Sample'; % sample MAT file name
particleSet = 'Data/geoParticleSetA'; % particle Set MAT file name
resultFileName = 'Result/geoSetA_t1k_RefMicro'; % result MAT file name to be saved

MicroporosityFlag = true;

SampleConstants; % load sample properties & constants
SimulationParameters; % load preset of zeta, uavg and Npe

%% BIG FOR
for run_number = 1:length(run_zeta)
    fprintf(sprintf('Run %2d: ',run_number));
    
    u = u_org * run_coef(run_number); uavg = mean(u(u~=0));
    v = v_org * run_coef(run_number); vavg = mean(v(v~=0));
    w = w_org * run_coef(run_number); wavg = mean(w(w~=0));
    
    ReinjectionVariables; % inlet/outlet velocity for reinjection
    
    zeta = run_zeta(run_number);
    [Dm,Npe,slope] = CalculateNpe (dt,L,gFlow,zeta,uavg,vavg,wavg);
    
    load(particleSet); % load Particles (xpt,ypt,zpt)
    Nparticle = length (xpt); % get the number of particles
    Ntimestep = 100; % set number of timeStep
    
    % Simulation
    fprintf(sprintf('Nparticle=%d, Ntimestep=%d, Npe=%.2f, micro=%d.\n        ',...
        Nparticle,Ntimestep,Npe,MicroporosityFlag));
    tic; Simulation; toc;
    
    % Save results to file
    saveFile = sprintf('%s%d.mat',resultFileName,run_number);
    save(saveFile,'xpt','ypt','zpt','displacement','variance','sim','insolid');
    fprintf(sprintf('        Saved to file %s.\n',saveFile));
end
    
%% Nested Function

function [Dm,Npe,slope] = CalculateNpe (dt,L,gFlow,zeta,uavg,vavg,wavg)
    Dm = zeta^2 / (6*dt);
    if gFlow == 1 % FLOW_X
        Npe = uavg * L / Dm;
        slope = ceil(L/(uavg*dt));
    elseif gFlow == 2 % FLOW_Y
        Npe = vavg * L / Dm;
        slope = ceil(L/(vavg*dt));
    else % gFlow == 3 % FLOW_Z
        Npe = wavg * L / Dm;
        slope = ceil(L/(wavg*dt));
    end
end