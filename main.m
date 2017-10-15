%% SoluteTranport
%
% To study the effect of Peclet Number (Npe), 
% we vary the flow rate (u,v,w) & diffusive jump distance (zeta), and run simulation
% until we reach an asymptotic regime (computed dispersion coefficient no longer change with time)
%
% To study the significance of Microporosity,
% we run the simulation again (same parameters, but with microporosity)
%
close all; clearvars -except u_org v_org w_org g;

D_LIMIT = 0.01; % dispersion coefficient change limit

sample = 'Data/geometry2Sample'; % sample MAT file name
particleSet = 'Data/geoParticleSetB'; % particle et MAT file name
resultName = 'Result/geoBha'; % result MAT file name to be saved

SampleConstants; % load sample properties & constants
SimulationParameters; % load preset of Npe numbers (RUN parameters)

%% BIG FOR
for run_number = 1:10   
    % Vary the flow rate (u,v,w) & diffusive jump distance (zeta)
    u = u_org * run_coef(run_number); uavg = mean(u(u~=0));
    v = v_org * run_coef(run_number); vavg = mean(v(v~=0));
    w = w_org * run_coef(run_number); wavg = mean(w(w~=0));    
    zeta = run_zeta(run_number);
       
    ReinjectionVariables; % inlet/outlet velocity for reinjection
    
    load(particleSet); % load Particles (xpt,ypt,zpt)  
    Nparticle = length (xpt); % number of particles    
    [Npe,Dm,slope] = calculateNpe (uavg,vavg,wavg,zeta,gFlow,L,dt);
    
    % Simulation1 until reaching an asymptotic regime
    fprintf(sprintf('Run %2d: Npe=%.2f, Nparticle=%d, slope=%d, RDC=%.2f\n',...
        run_number, Npe, Nparticle, slope, run_RDC(run_number)));
    tic; Simulation1; toc;
    
    xptN = xpt; yptN = ypt; zptN = zpt;   
     
    % Simulation2 (microporosity) with Ntimestep
    load(particleSet);  % re-load particle
    tic; Simulation2; toc;
    
    % Save results to file
    saveFileName = sprintf('%s%d.mat',resultName,run_number);
    save(saveFileName,...
        'dt','dx','L','gN','gFlow','gL',... % geometry properties
        'sample','particleSet','Nparticle','Ntimestep',... % simulation parameters
        'zeta','uavg','vavg','wavg','Dm','Npe','slope',... % RUN parameters 
        'chance_into_solid','micro_zeta','inside_limit',... % microporosity parameters
        'xptN','yptN','zptN','displacement1','variance1','dispersion1','varianceTime',... % Simulation1 result
        'xpt','ypt','zpt','displacement2','variance2','dispersion2','insolid'... % Simulation2 result   
        );
    fprintf(sprintf('    Saved to file %s.\n',saveFileName));
end
    
%% Nested Function

function [Npe,Dm,slope] = calculateNpe (uavg,vavg,wavg,zeta,gFlow,L,dt)
    Dm = zeta^2 / (6*dt);
    if gFlow == 1 % FLOW_X
        Npe = uavg * L / Dm;
        slope = ceil (L / (uavg*dt));
    elseif gFlow == 2 % FLOW_Y
        Npe = vavg * L / Dm;
        slope = ceil (L / (vavg*dt));
    else % gFlow == 3 % FLOW_Z
        Npe = wavg * L / Dm;
        slope = ceil (L / (wavg*dt));
    end
end