%% constants
FLOW_X = 1; FLOW_Y = 2; FLOW_Z = 3; % Flow direction
DM_LIQUIDS = 10^3; % free molecular diffusion in liquids (microm2.s-1)
VELOCITY_FACTORS = [ 0.0001 0.000316228 0.001 ...
    0.001467799 0.002154435 0.003162278 0.004641589 0.006812921 0.01 ...
    0.014677993 0.021544347 0.031622777 0.046415888 0.068129207 0.1 ...
    0.146779927 0.215443469 0.316227766 0.464158883 0.681292069 1 ];

%% simulation parameters
MAX_TIMESTEP = 10000; % limit timestep for a simulation
SAVE_TIMESTEP = 1000; % save result to hard drive every SAVE_TIMESTEP

continue_simulation = true; % continue a simulation?
continue_matname = 'BereaA0p5';
continue_savenumber = 6; % must +1 from matName

microporosity_chance = 0.00; % chance
micro_zeta = 0.1;
inside_limit = 0.05;

data_folder = '../data/';
result_folder = '../result/';
sample_name = 'Berea'; % carbondate or sandpack
particle_set = 'A'; % A, B, C, D or E

gN    = 400;     % number of voxels
gFlow = FLOW_Z;  % flow direction of the sample
dx    = 5.3;     % dimension of the cell (dx=dy=dz) in micro meter
L     = 129.84;   % characteristic length (micro meter)
% geopack: L=150.21
% Berea: L=129.84
% bentheimer: L=154.47

gLen = dx * gN - dx; % sample length (micro meter)
dt = 0.0001; % timestep (seconds)
zeta = sqrt(6*DM_LIQUIDS*dt); % distance that particle diffusively jumps at a timestep (micro meter)

Nsave = MAX_TIMESTEP/SAVE_TIMESTEP; % number of time to save data to hard drive
Nrun = length(VELOCITY_FACTORS); % number of run

%% load data
if ~(exist('g','var') && exist('u_org','var') && exist('v_org','var') && exist('w_org','var'))
    fprintf(sprintf('Load sample %s...',sample_name));
    load([data_folder,sample_name]);
    fprintf(sprintf('Done!\n'));
end

fprintf(sprintf('Load initial particle location set %s...',particle_set));
load([data_folder,sample_name,'ParticleSet',particle_set],'xpt','ypt','zpt');
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
prev_displacement = zeros(Nparticle,Nrun);

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
    
    prev_particleX = particleX;
    prev_particleY = particleY;
    prev_particleZ = particleZ;
    prev_inside_solid = inside_solid;
    prev_displacement = displacement;
else
    savenumber_start = 1;
end

%% SIMULATION
for sn = savenumber_start:Nsave
    tic
    fprintf(sprintf('Simulating %d particles, from %d to %d timestep, %.2f microporosity chance...',...
        Nparticle, (sn-1)*SAVE_TIMESTEP+1, sn*SAVE_TIMESTEP, microporosity_chance));
    
    for rn = 1:Nrun
        u = u_org*VELOCITY_FACTORS(rn);
        v = v_org*VELOCITY_FACTORS(rn);
        w = w_org*VELOCITY_FACTORS(rn);
        
        Npe = [mean(u(u~=0)),mean(v(v~=0)),mean(w(w~=0))] * L / DM_LIQUIDS;
        run_Npe(rn) = Npe(gFlow);
        
        [particleX(:,rn),particleY(:,rn),particleZ(:,rn),...
            displacement(:,rn),variance(:,rn),inside_solid(:,rn)] = ...
            particletracking(g,u,v,w,dt,dx,gFlow,gLen,zeta,...
            micro_zeta,microporosity_chance,inside_limit,...
            prev_inside_solid(:,rn),prev_displacement(:,rn),...
            prev_particleX(:,rn),prev_particleY(:,rn),prev_particleZ(:,rn),...
            xpt,ypt,zpt,SAVE_TIMESTEP);
        
    end
    
    prev_particleX = particleX;
    prev_particleY = particleY;
    prev_particleZ = particleZ;
    prev_inside_solid = inside_solid;
    prev_displacement = displacement;
    
    saveFileName = sprintf('%s%s%s%d.mat',sample_name,particle_set,...
        [num2str(microporosity_chance*100),'p'],sn);
    
    save([result_folder,saveFileName],'run_Npe','particleX','particleY','particleZ',...
        'displacement','variance','inside_solid',...
        'data_folder','result_folder','sample_name','particle_set',... % simulation parameters
        'microporosity_chance','micro_zeta','inside_limit',... % microporosity parameters
        'dt','zeta','dx','L','gN','gFlow','gLen'); % sample properties
    fprintf(sprintf('Done!\n -> Successfully saved to file %s! ',saveFileName));
    toc
end
