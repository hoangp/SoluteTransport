clear;

%% constants
MAX_TIMESTEP = 2000; % limit timestep for a simulation
SAVE_TIMESTEP = 1000; % save result to hard drive every SAVE_TIMESTEP
FLOW_X = 1; FLOW_Y = 2; FLOW_Z = 3; % Flow direction
DM_LIQUIDS = 10^3; % free molecular diffusion in liquids (micro-m2.s-1)
VELOCITY_FACTORS = [ 0.0001 0.000316228 0.001 ...
    0.001467799 0.002154435 0.003162278 0.004641589 0.006812921 0.01 ...
    0.014677993 0.021544347 0.031622777 0.046415888 0.068129207 0.1 ...
    0.146779927 0.215443469 0.316227766 0.464158883 0.681292069 1 ]; 
Nsave = MAX_TIMESTEP/SAVE_TIMESTEP; % number of time to save data to hard drive
Nrun = length(VELOCITY_FACTORS); % number of run

%% simulation parameters
data_folder = '../data/';
result_folder = '../result/';
sample_name = 'carbonate'; % carbondate or sandpack
particle_set = 'A'; % A, B, C, D or E

continue_simulation = false; % continue a simulation?
continue_matname = 'carbonateA1p2';
continue_savenumber = 3;

microporosity_chance = 0.00;
micro_zeta = 0.1; 
inside_limit = 0.05; 

%% load sample geometry (g), velocity field (u,v,w)
if ~(exist('g','var') && exist('u_org','var') && exist('v_org','var') && exist('w_org','var'))
    fprintf(sprintf('Load sample %s...',sample_name));   
    load([data_folder,sample_name],'g','u_org','v_org','w_org');
    fprintf(sprintf('Done!\n')); 
end
g = g ~= 0;

%% sample properties (dt,zeta,dx,L,gN,gFlow,gLen)
dt = 0.0001; % timestep (seconds)
zeta = sqrt(6*DM_LIQUIDS*dt); % distance that particle diffusively jumps at a timestep (micro meter)

switch lower(sample_name)
    case 'carbonate'
        dx    = 3.8;     % dimension of the cell (dx=dy=dz) in micro meter
        L     = 432.2;   % characteristic length (micro meter)
        gN    = 400;     % number of voxels         
        gFlow = FLOW_Z;  % flow direction of the sample
    case 'sandpack'
        dx    = 5.345;   % dimension of the cell (dx=dy=dz) in micro meter
        L     = 75.72;   % characteristic length (micro meter)
        gN    = 300;     % number of voxels     
        gFlow = FLOW_X;  % flow direction of the sample     
    case 'sandstone'
        dx    = 2.32;    % dimension of the cell (dx=dy=dz) in micro meter
        L     = 54.93;   % characteristic length (micro meter)
        gN    = 500;     % number of voxels         
        gFlow = FLOW_Z;  % flow direction of the sample       
    case 'carbonate2'
        dx    = 3.8;     % dimension of the cell (dx=dy=dz) in micro meter
        L     = 236.33;  % characteristic length (micro meter) 390
        gN    = 636;     % number of voxels         
        gFlow = FLOW_Z;  % flow direction of the sample   
end
gLen = dx * gN - dx; % sample length (micro meter)

%% load initial particle location
fprintf(sprintf('Load initial particle location set %s...',particle_set));
load([data_folder,sample_name,'_particleset',particle_set],'xpt','ypt','zpt');
fprintf(sprintf('Done!\n'));

%% results to be saved (for every SAVE_TIMESTEP)
Nparticle = length (xpt); % number of particles
run_Npe = zeros(1,Nrun);
particleX = zeros(Nparticle,Nrun); % location of particles
particleY = zeros(Nparticle,Nrun);
particleZ = zeros(Nparticle,Nrun);
inside_solid = false (Nparticle,Nrun);  % status of each particle (in-solid or in-pore)
displacement = zeros(Nparticle,Nrun); % particle displacement in the main flow direction
variance = zeros(SAVE_TIMESTEP,Nrun); % variance of displacement

prev_particleX = zeros(Nparticle,Nrun); 
prev_particleY = zeros(Nparticle,Nrun); 
prev_particleZ = zeros(Nparticle,Nrun);
prev_inside_solid = false (Nparticle,Nrun); 
for rn = 1:Nrun
    prev_particleX(:,rn)=xpt;
    prev_particleY(:,rn)=ypt;
    prev_particleZ(:,rn)=zpt;
end

%% continue or new simulation?
if continue_simulation == true
    fprintf(sprintf('Load particle location from save file %s...',continue_matname));
    load([result_folder,continue_matname]);
    fprintf(sprintf('Done!\n'));
    savenumber_start = continue_savenumber;
else
    savenumber_start = 1;
end

%% SIMULATION
for sn = savenumber_start:Nsave
    tic
    fprintf(sprintf('Simulating %d particles, from %d to %d timestep, %.2f microporosity chance...',...
        Nparticle, (sn-1)*SAVE_TIMESTEP+1, sn*SAVE_TIMESTEP, microporosity_chance));
    parfor rn = 1:Nrun
        [run_Npe(rn),particleX(:,rn),particleY(:,rn),particleZ(:,rn),...
            displacement(:,rn),variance(:,rn),inside_solid(:,rn)] = ...
        particletracking(g,u_org*VELOCITY_FACTORS(rn),v_org*VELOCITY_FACTORS(rn),...
            w_org*VELOCITY_FACTORS(rn),dt,dx,L,gFlow,gLen,zeta,...
            micro_zeta,microporosity_chance,inside_limit,prev_inside_solid(:,rn),...
            prev_particleX(:,rn),prev_particleY(:,rn),prev_particleZ(:,rn),...
            xpt,ypt,zpt,SAVE_TIMESTEP);
    end
    prev_particleX = particleX;
    prev_particleY = particleY;
    prev_particleZ = particleZ;
    prev_inside_solid = inside_solid;
    
    saveFileName = sprintf('%s%s%s%d.mat',sample_name,particle_set,...
        [num2str(microporosity_chance*100),'p'],sn);
    save([result_folder,saveFileName],'run_Npe','particleX','particleY','particleZ',...
        'displacement','variance','inside_solid',...
        'prev_particleX','prev_particleY','prev_particleZ',...
        'data_folder','result_folder','sample_name','particle_set',... % simulation parameters
        'microporosity_chance','micro_zeta','inside_limit',... % microporosity parameters
        'dt','zeta','dx','L','gN','gFlow','gLen'); % sample properties       
    fprintf(sprintf('Done!\n -> Successfully saved to file %s! ',saveFileName));
    toc       
end
