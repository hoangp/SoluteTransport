%% SoluteTranport
%
% Simulation1: To study the effect of Peclet Number (Npe), 
%   we vary the flow rate (u,v,w) & diffusive jump distance (zeta), and run simulation
%   until we reach an asymptotic regime (computed dispersion coefficient no longer change with time)
%
% Simulation2: To study the significance of Microporosity,
%   we run the simulation again (same parameters, but with microporosity)
%
close all; clearvars -except u_org v_org w_org g;

D_LIMIT = 0.00; % dispersion coefficient change limit
                % D_LIMIT = 0 means: simulate until MAX_TIMESTEP
MAX_TIMESTEP = 100; % maximum timestep for a simulation
MIN_TIMESTEP = 10; % minimum timestep for a simulation
MIN_SLOPE = 10; % minimum timestep slope to calculate dispersion

sample = 'Data/geometry2Sample'; % sample MAT file name
particleSet = 'Data/geoParticleSetA'; % particle et MAT file name
resultName = 'Result/test'; % result MAT file name to be saved

SampleConstants; % load sample properties & constants
SimulationParameters; % load preset of Npe numbers (RUN parameters)

%% Micro-porosity parameters
chance_into_solid = 0.01;
micro_zeta = 0.1; 
inside_limit = 0.05; 

%% BIG FOR
for run_number = 10:length(run_zeta)  
    % Vary the flow rate (u,v,w) & diffusive jump distance (zeta)
    u = u_org * run_coef(run_number); uavg = mean(u(u~=0));
    v = v_org * run_coef(run_number); vavg = mean(v(v~=0));
    w = w_org * run_coef(run_number); wavg = mean(w(w~=0));    
    zeta = run_zeta(run_number);
    refRDC = run_RDC(run_number); % reference Reduced Dispersion Coefficient 
       
    ReinjectionVariables; % inlet/outlet velocity for reinjection
    SimulationEvents; % reset event structure  
    
    load(particleSet); % load Particles (xpt,ypt,zpt)
    xpt0 = xpt; ypt0 = ypt; zpt0 = zpt; % original injected particle (for ReinjectionDiffusion2)
    Nparticle = length (xpt); % number of particles    
    [Npe,Dm,slope] = calculateNpe (uavg,vavg,wavg,zeta,gFlow,L,dt,MAX_TIMESTEP,MIN_SLOPE);   
    
    % Simulation1 until reaching an asymptotic regime
    fprintf(sprintf('Run %d: Npe=%.2f, Nparticle=%d, slope=%d, RDC=%.2f\n',...
        run_number, Npe, Nparticle, slope, run_RDC(run_number)));
    tic; Simulation1; toc;    
    xpt1 = xpt; ypt1 = ypt; zpt1 = zpt; % save xpt,ypt,zpt of simulation1
     
    % Simulation2 (microporosity) with Ntimestep
    load(particleSet); % re-load particle
    tic; Simulation2; toc;
    xpt2 = xpt; ypt2 = ypt; zpt2 = zpt; % save xpt,ypt,zpt of simulation2
    
    % Save results to file
    saveFileName = sprintf('%s%d.mat',resultName,run_number);
    save(saveFileName,...
        'dt','dx','L','gN','gFlow','gL',... % geometry properties
        'sample','particleSet','Nparticle','Ntimestep',... % simulation parameters
        'zeta','uavg','vavg','wavg','Dm','Npe','slope','refRDC',... % RUN parameters 
        'chance_into_solid','micro_zeta','inside_limit',... % microporosity parameters
        'xpt1','ypt1','zpt1','displacement1','variance1','dispersion1','varianceTime','dispersionChange',... % Simulation1 result
        'xpt2','ypt2','zpt2','displacement2','variance2','dispersion2','insolid',... % Simulation2 result   
        'event'); % event structure
    fprintf(sprintf(' -> Saved to file %s.\n',saveFileName));
end
    
%% Nested Function

function [Npe,Dm,slope] = calculateNpe (uavg,vavg,wavg,zeta,gFlow,L,dt,MAX_TIMESTEP,MIN_SLOPE)
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
    
    if slope * 4 > MAX_TIMESTEP 
        slope = floor(MAX_TIMESTEP/4);
    elseif slope < MIN_SLOPE
        slope = MIN_SLOPE;
    end
end