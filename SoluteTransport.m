%% SoluteTranport: All-In-One script for better performance
%
% Simulation1: To study the effect of Peclet Number (Npe), 
%   we vary the flow rate (u,v,w) & diffusive jump distance (zeta), and run simulation
%   until we reach an asymptotic regime (computed dispersion coefficient no longer change with time)
%
% Simulation2: To study the significance of Microporosity,
%   we run the simulation again (same parameters, but with microporosity)

close all; clearvars -except u_org v_org w_org g;
xpt = zeros(1); ypt = zeros(1); zpt = zeros(1); % avoid Matlab style attention

D_LIMIT = 0.01; % dispersion coefficient change limit
                % D_LIMIT = 0 means: simulate until MAX_TIMESTEP
MAX_TIMESTEP = 40000; % maximum timestep for a simulation
MIN_TIMESTEP = 40000; % minimum timestep for a simulation
MIN_SLOPE = 1000; % minimum timestep slope to calculate dispersion

sample = 'Data/sand500Sample'; % sample MAT file name
particleSet = 'Data/sandParticleSetB'; % particle et MAT file name
resultName = 'Result/sandB_t40k_com'; % result MAT file name to be saved

%% SimulationParameters; % load preset of Npe numbers (RUN parameters)
% geometry2sample RUN Parameters Sets (preset Npe) 
% Combination
%run_zeta = [5.3         5.3         5.3         5.3         5.3         5.3         5.3         5.3         4.8837 4.0550 3.6123 3.0220 2.0730 1.4291 1.0360 0.6827 0.4089 0.2467 0.1277 ];
%run_RDC  = [0.32406     0.38463     0.45652     0.67704     1.1516      2.0504      3.448       6.537        14.46 20.723 28.861 48.257 125.97  319.6 641.55 1379.2 3419.9 9561   32832 ];
%run_coef = [0.004276384 0.011332418 0.023481625 0.046270475 0.084244766 0.165068424 0.278392601 0.534975643  1     1      1      1      1       1      1      1      1      1      1  ];
% sand500sample RUN Parameters Sets (preset Npe) 
% Combination
run_zeta = [ 2.32    2.32    2.32    2.32    2.32    2.32    2.32    2.32    2.32    2.32    2.32    2.32    2.32    2.32    2.32    2.32    2.32     2.2356 1.9665 1.7388 1.5142 1.3492 1.0923 0.9609 0.8202 0.7362 0.6476 0.6112 0.5784 0.5445 0.5127 0.4556 0.3689 0.3125 0.2694 0.2424 0.2154 0.1978 0.1784 0.1562 0.1384 0.1221 0.1164 0.1029 0.0842 0.0744 0.0616 0.0517 0.0449 0.0372 0.0193  0.027319 0.022306 ];
run_RDC  = [ 0.31136 0.30959 0.31136 0.31136 0.31136 0.30959 0.32406 0.34508 0.36328 0.38463 0.40492 0.41192 0.45652 0.51177 0.67704 0.8316  1.1516   1.641  2.0504 2.6362 3.448  4.3578 6.537  8.5012 11.639 14.46  18.698 20.723 24.596 28.861 34.451 48.257 82.552 125.97 180.53 238.84 319.6  397.07 490.51 641.55 801.62 1013.1 1078.8 1379.2 1965.2 2441.6 3419.9 4957.4 6671.7 9561   32832   17000    25000];
run_coef = [ 0.00359 0.00513 0.00690 0.01325 0.01801 0.02803 0.03640 0.05072 0.06720 0.09553 0.12981 0.14354 0.19800 0.27867 0.39024 0.49561 0.70970  1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1       1        1  ];
% Micro-porosity parameters
chance_into_solid = 0.01;
micro_zeta = 0.1; 
inside_limit = 0.05; 

%% SampleConstants; % load sample properties & constants
% Constants
FLOW_X = 1; FLOW_Y = 2; FLOW_Z = 3; % Flow direction
OUT_X1 = 1; OUT_X2 = 2; OUT_Y1 = 3; OUT_Y2 = 4; OUT_Z1 = 5; OUT_Z2 = 6; % out_system_location

% Sample Properties
% 
% Workspace Variables:
%   sample: sample file name
% 
% Output:
%   g,u,v,w
%   dx,gN,L,gL,gFlow
%

if ~(exist('g','var') && exist('u_org','var') && exist('v_org','var') && exist('w_org','var'))
    fprintf(sprintf('Load sample: %s...',sample)); tic
    load(sample);
    fprintf(sprintf('Done! ')); toc
end

dt = 0.0001; % time-step (seconds)

switch lower(sample)
    case 'data/geometry2sample'
        dx    = 5.345;   % dimension of the cell (dx=dy=dz) in micro meter
        L     = 131.13;  % characteristic length (micro meter)
        gN    = 300;     % number of voxels     
        gFlow = FLOW_X;  % sample main flow direction
        
    case 'data/sand500sample'
        dx    = 2.32;    % dimension of the cell (dx=dy=dz) in micro meter
        L     = 51.63;   % characteristic length (micro meter)
        gN    = 500;     % number of voxels         
        gFlow = FLOW_Z;  % sample main flow direction     
end

gL = dx * gN - dx; % sample length (micro meter)

%% BIG FOR
for run_number = 1:length(run_zeta)  
    %% Vary the flow rate (u,v,w) & diffusive jump distance (zeta)
    u = u_org * run_coef(run_number); uavg = mean(u(u~=0));
    v = v_org * run_coef(run_number); vavg = mean(v(v~=0));
    w = w_org * run_coef(run_number); wavg = mean(w(w~=0));    
    zeta = run_zeta(run_number);
    refRDC = run_RDC(run_number); % reference Reduced Dispersion Coefficient 
       
    %% ReinjectionVariables; % inlet/outlet velocity for reinjection
    % Prepare Velocity Fields for Reinjection
    %
    % Workspace variables
    %   u,v,w   : components of velocity field in x,y,z direction, respectively
    %   gN      : geometry number of voxel
    %
    % Output
    %   Uinlet  : velocity field of the inlet faces (X direction)
    %   Vinlet  : velocity field of the inlet faces (Y direction)
    %   Winlet  : velocity field of the inlet faces (Z direction)
    %   Uoutlet : velocity field of the outlet faces (X direction)
    %   Voutlet : velocity field of the outlet faces (Y direction)
    %   Woutlet : velocity field of the outlet faces (Z direction)
    %
    
    % Constants
    REINJECTION_INLET = 2;
    REINJECTION_OUTLET = gN - 1;
    REINJECTION_LIMIT = 0.1; % minimum velocity to re-inject
    
    % U velocity plane at inlet face
    Uinlet.face  = REINJECTION_INLET;
    Uinlet.plane = u(Uinlet.face,:,:);
    Uinlet.list  = reshape(Uinlet.plane,[],1);
    Uinlet.limit = mean(Uinlet.plane (Uinlet.plane~=0)) * REINJECTION_LIMIT;
    Uinlet.index = find(Uinlet.list > Uinlet.limit);
    Uinlet.list  = Uinlet.list(Uinlet.index);
    
    % U velocity plane at outlet face
    Uoutlet.face  = REINJECTION_OUTLET;
    Uoutlet.plane = u(Uoutlet.face,:,:);
    Uoutlet.list  = reshape(Uoutlet.plane,[],1);
    Uoutlet.limit = mean(Uoutlet.plane (Uoutlet.plane~=0)) * REINJECTION_LIMIT;
    Uoutlet.index = find(Uoutlet.list > Uoutlet.limit);
    Uoutlet.list  = Uoutlet.list(Uoutlet.index);
    
    % V velocity plane at inlet face
    Vinlet.face  = REINJECTION_INLET;
    Vinlet.plane = v(:,Vinlet.face,:);
    Vinlet.list  = reshape(Vinlet.plane,[],1);
    Vinlet.limit = mean(Vinlet.plane (Vinlet.plane~=0)) * REINJECTION_LIMIT;
    Vinlet.index = find(Vinlet.list > Vinlet.limit);
    Vinlet.list  = Vinlet.list(Vinlet.index);
    
    % V velocity plane at outlet face
    Voutlet.face  = REINJECTION_OUTLET;
    Voutlet.plane = v(:,Voutlet.face,:);
    Voutlet.list  = reshape(Voutlet.plane,[],1);
    Voutlet.limit = mean(Voutlet.plane (Voutlet.plane~=0)) * REINJECTION_LIMIT;
    Voutlet.index = find(Voutlet.list > Voutlet.limit);
    Voutlet.list  = Voutlet.list(Voutlet.index);
    
    % W velocity plane at inlet face
    Winlet.face  = REINJECTION_INLET;
    Winlet.plane = w(:,:,Winlet.face);
    Winlet.list  = reshape(Winlet.plane,[],1);
    Winlet.limit = mean(Winlet.plane (Winlet.plane~=0)) * REINJECTION_LIMIT;
    Winlet.index = find(Winlet.list > Winlet.limit);
    Winlet.list  = Winlet.list(Winlet.index);
    
    % W velocity plane at outlet face
    Woutlet.face  = REINJECTION_OUTLET;
    Woutlet.plane = w(:,:,Woutlet.face);
    Woutlet.list  = reshape(Woutlet.plane,[],1);
    Woutlet.limit = mean(Woutlet.plane (Woutlet.plane~=0)) * REINJECTION_LIMIT;
    Woutlet.index = find(Woutlet.list > Woutlet.limit);
    Woutlet.list  = Woutlet.list(Woutlet.index);
    
    %% SimulationEvents; % reset event structure
    % Simulation Events Structure
    event.ReinjectAdvectionX = 0;
    event.ReinjectAdvectionY = 0;
    event.ReinjectAdvectionZ = 0;
    event.ReinjectDiffusionX = 0;
    event.ReinjectDiffusionY = 0;
    event.ReinjectDiffusionZ = 0;
    event.AdvectionHitSolid = 0;
    event.DiffusionHitSolid = 0;
    event.GetInSolid = 0;
    event.GetOutSolid = 0;
    
    %% load Particles (xpt,ypt,zpt)
    load(particleSet);
    xpt0 = xpt; ypt0 = ypt; zpt0 = zpt; % original injected particle
    Nparticle = length (xpt); % number of particles    
    
    %% [Npe,Dm,slope] = calculateNpe (uavg,vavg,wavg,zeta,gFlow,L,dt,MAX_TIMESTEP);
    Dm = zeta^2 / (6*dt);
    
    if gFlow == 1 % FLOW_X
        Npe = uavg * L / Dm;
        slope = ceil (L / (uavg*dt));
        %slope = ceil (gN*dx / (uavg*dt));
    elseif gFlow == 2 % FLOW_Y
        Npe = vavg * L / Dm;
        slope = ceil (L / (vavg*dt));
        %slope = ceil (gN*dx / (vavg*dt));
    else % gFlow == 3 % FLOW_Z
        Npe = wavg * L / Dm;
        slope = ceil (L / (wavg*dt));
        %slope = ceil (gN*dx / (wavg*dt));
    end
    
    if slope * 4 > MAX_TIMESTEP 
        slope = floor(MAX_TIMESTEP/4);
    elseif slope < MIN_SLOPE
        slope = MIN_SLOPE;
    end
    
    %% Simulation1 until reaching an asymptotic regime
    fprintf(sprintf('Run %d: Npe=%.2f, Nparticle=%d, slope=%d, refRDC=%.2f\n',...
        run_number, Npe, Nparticle, slope, refRDC));
    tic; 
    % Simulation until we reach an asymptotic regime
    % (computed dispersion coefficient no longer change with time)
    %
    % Workspace Variables
    %   Nparticle      : number of particle
    %   slope          : number of timestep to calculate dispersion
    %   xpt,ypt,zpt    : particles locations
    %
    % Output
    %   xpt,ypt,zpt    : particles locations
    %   displacement1  : particles displacement in the main flow
    %   variance1      : variance of displacement (compute at each slope time)
    %   varianceTime   : timestep of the variance
    %   dispersion1    : longitudinal dispersion coefficient
    %
    % Constant used
    %   D_LIMIT        : dispersion coefficient change limit
    %   MAX_TIMESTEP   : maximum timestep for a simulation
    %   MIN_SLOPE      : minimum timestep slope to calculate dispersion
    
    maxnumSlope = ceil(MAX_TIMESTEP/MIN_SLOPE)+1; % Maximum number of slope
    
    % Simulation variables
    displacement1 = zeros(Nparticle,1);
    variance1 = zeros(maxnumSlope,1);
    varianceTime = zeros(maxnumSlope,1);
    dispersion1 = zeros(maxnumSlope,1);
    dispersionChange = zeros(maxnumSlope,1);
    
    % Simulation starts
    reverseStr = ''; % for diplay progress
    counter = 2; dispersionChange(counter-1) = 1; limitBeforeMinTime = 0;
    while (dispersionChange(counter-1) > D_LIMIT && varianceTime(counter-1) < MAX_TIMESTEP) || (varianceTime(counter-1) < MIN_TIMESTEP)
        for t = 1:slope
            for p = 1:Nparticle
                x = xpt(p); y = ypt(p); z = zpt(p); % short variables instead of array(index)
                
                %% AdvectionMovement
                % Avection Movement
                % Move particle along streamlines that follow velocity field
                %
                % Workspace variables
                %   g     : geometry
                %   dx    : dimension of the cell (dx=dy=dz)
                %   dt    : time-step
                %   gL    : sample length
                %   x,y,z : current location of particle
                %
                % Output
                %   x,y,z : Update location of particle
                %   advectionX, advectionY, advectionZ : Advection displacement
                %
                % Scipts Used:
                %   StreamlineTracing: Calculate ui,vi,wi
                %   ReinjectAdvection: Update x,y,z
                %
                % Functions Used:
                %   out_system : Check if particle out system or not
                %   get_ijk    : Voxel (i,j,k) of particle location (x,y,z)
                %
                
                % Particle location before advection
                xBA = x; yBA = y; zBA = z;
                
                % Streamline-tracing algorithm
                %% StreamlineTracing
                % Streamline-tracing algorithm
                % Semi-analytic formulation for the streamlines in voxels contain a solid boundary
                %
                % Workspace variables
                %   g        : geometry
                %   dx       : dimension of the cell (dx=dy=dz)
                %   x,y,z    : Location of particle
                %
                % Output
                %   ui,vi,wi : interpolated velocity between the normal velocities on opposite faces of gridblock
                %
                % Restriction:
                %   Only consider geometry from (2) --> (N-1)
                
                % Get particle boundary conditions
                
                % index location of the particle            
                i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);% [i,j,k] = get_ijk (x,y,z,dx);
                
                % location of the edges
                x1 = dx*(i-1); x2 = dx*i;
                y1 = dx*(j-1); y2 = dx*j;
                z1 = dx*(k-1); z2 = dx*k;
                
                % velocities on the faces normal to each direction
                u1 = u(i-1,j,k); u2 = u(i,j,k);
                v1 = v(i,j-1,k); v2 = v(i,j,k);
                w1 = w(i,j,k-1); w2 = w(i,j,k);
                
                % short variables names
                X1 = x-x1; X2 = x2-x;
                Y1 = y-y1; Y2 = y2-y;
                Z1 = z-z1; Z2 = z2-z;
                U = u2-u1; V = v2-v1; W = w2-w1;
                
                % geometry of the edge
                gx1 = g(i-1,j,k); gx2 = g(i+1,j,k);
                gy1 = g(i,j-1,k); gy2 = g(i,j+1,k);
                gz1 = g(i,j,k-1); gz2 = g(i,j,k+1);
                
                % short variables names
                gX = gx1+gx2; gY = gy1+gy2; gZ = gz1+gz2;
                gXY = gX+gY; gXZ = gX+gZ; gYZ = gY+gZ;
                
                % number of solid boundary arround the particle
                num_solid = gx1 + gx2 + gy1 + gy2 + gz1 + gz2;
                
                % Calculate interpolated velocities: ui,vi,wi
                
                if num_solid == 0
                    %% No solid boundaries
                    ui = U * X1 / dx + u1;
                    vi = V * Y1 / dx + v1;
                    wi = W * Z1 / dx + w1;
                    
                elseif num_solid == 1
                    %% One of the Neighbouring Voxels Is Solid (6 cases)
                    if gx2 == 1
                        % type = '1 solid face at x2'; SAMPLE FORMULATION
                        ui = u1 * X2^2 / dx^2;
                        vi = 2*v1 * X2 / dx + 2*V * X2*Y1 / dx^2;
                        wi = 2*w1 * X2 / dx + 2*W * X2*Z1 / dx^2;
                    elseif gx1 == 1
                        % type = '1 solid face at x1';
                        ui = u2 * X1^2 / dx^2;
                        vi = 2*v1 * X1 / dx + 2*V * X1*Y1 / dx^2;
                        wi = 2*w1 * X1 / dx + 2*W * X1*Z1 / dx^2;
                    elseif gy2 == 1
                        % type = '1 solid face at y2';
                        ui = 2*u1 * Y2 / dx + 2*U * Y2*X1 / dx^2;
                        vi = v1 * Y2^2 / dx^2;
                        wi = 2*w1 * Y2 / dx + 2*W * Y2*Z1 / dx^2;
                    elseif gy1 == 1
                        % type = '1 solid face at y1';
                        ui = 2*u1 * Y1 / dx + 2*U * Y1*X1 / dx^2;
                        vi = v2 * Y1^2 / dx^2;
                        wi = 2*w1 * Y1 / dx + 2*W * Y1*Z1 / dx^2;
                    elseif gz2 == 1
                        % type = '1 solid face at z2';
                        ui = 2*u1 * Z2 / dx + 2*U * Z2*X1 / dx^2;
                        vi = 2*v1 * Z2 / dx + 2*V * Z2*Y1 / dx^2;
                        wi = w1 * Z2^2 / dx^2;
                    elseif gz1 == 1
                        % type = '1 solid face at z1';
                        ui = 2*u1 * Z1 / dx + 2*U * Z1*X1 / dx^2;
                        vi = 2*v1 * Z1 / dx + 2*V * Z1*Y1 / dx^2;
                        wi = w2 * Z1^2 / dx^2;
                    end
                    
                elseif num_solid == 2
                    %% Two of the Neighbouring Voxels Are Solid (15 cases)
                    
                    if (gX == 2 || gY == 2 || gZ == 2)
                        %% There Are Two Opposing Solid Faces (3 cases)
                        if gX==2 % (gx2 == 1 && gx1 == 1)
                            % type = '2 opposing solid faces x direction'; % SAMPLE FORMULATION
                            ui = 0;
                            vi = 6*v1 * X2*X1 / dx^2 + 6*V * X2*X1 * Y1 / dx^3;
                            wi = 6*w1 * X2*X1 / dx^2 + 6*W * X2*X1 * Z1 / dx^3;
                        elseif gY==2 % (gy2 == 1 && gy1 == 1)
                            % type = '2 opposing solid faces y direction';
                            ui = 6*u1 * Y2*Y1 / dx^2 + 6*U * Y2*Y1 * X1 / dx^3;
                            vi = 0;
                            wi = 6*w1 * Y2*Y1 / dx^2 + 6*W * Y2*Y1 * Z1 / dx^3;
                        else % gz==2 % (gz2 == 1 && gz1 == 1)
                            % type = '2 opposing solid faces z direction';
                            ui = 6*u1 * Z2*Z1 / dx^2 + 6*U * Z2*Z1 * X1 / dx^3;
                            vi = 6*v1 * Z2*Z1 / dx^2 + 6*V * Z2*Z1 * Y1 / dx^3;
                            wi = 0;
                        end
                        
                    elseif gXY == 2
                        %% Two Adjacent Voxels (XY) Are Solid (4 cases), PAPER BUG
                        if (gx1 == 1 && gy1 == 1)
                            % type = '2 solid faces at x1y1'; % SAMPLE FORMULATION
                            ui = 2*u2 * X1^2 * Y1 / dx^3;
                            vi = 2*v2 * X1 * Y1^2 / dx^3;
                            wi = 4*w1 * X1*Y1 / dx^2 + 4*W * X1*Y1 * Z1 / dx^3;
                        elseif (gx2 == 1 && gy1 == 1)
                            % type = '2 solid faces at x2y1';
                            ui = 2*u1 * X2^2 * Y1 / dx^3;
                            vi = 2*v2 * X2 * Y1^2 / dx^3;
                            wi = 4*w1 * X2*Y1 / dx^2 + 4*W * X2*Y1 * Z1 / dx^3;
                        elseif (gx1 == 1 && gy2 == 1)
                            % type = '2 solid faces at x1y2';
                            ui = 2*u2 * X1^2 * Y2 / dx^3;
                            vi = 2*v1 * X1 * Y2^2 / dx^3;
                            wi = 4*w1 * X1*Y2 / dx^2 + 4*W * X1*Y2 * Z1 / dx^3;
                        else % (gx2 == 1 && gy2 == 1)
                            % type = '2 solid faces at x2y2';
                            ui = 2*u1 * X2^2 * Y2 / dx^3;
                            vi = 2*v1 * X2 * Y2^2 / dx^3;
                            wi = 4*w1 * X2*Y2 / dx^2 + 4*W * X2*Y2 * Z1 / dx^3;
                        end
                        
                    elseif gXZ == 2
                        %% Two Adjacent Voxels (XZ) Are Solid (4 cases)
                        if (gx1 == 1 && gz1 == 1)
                            % type = '2 solid faces at x1z1';
                            ui = 2*u2 * X1^2 * Z1 / dx^3;
                            vi = 4*v1 * X1*Z1 / dx^2 + 4*V * X1*Z1 * Y1 / dx^3;
                            wi = 2*w2 * X1 * Z1^2 / dx^3;
                        elseif (gx2 == 1 && gz1 == 1)
                            % type = '2 solid faces at x2z1';
                            ui = 2*u1 * X2^2 * Z1 / dx^3;
                            vi = 4*v1 * X2*Z1 / dx^2 + 4*V * X2*Z1 * Y1 / dx^3;
                            wi = 2*w2 * X2 * Z1^2 / dx^3;
                        elseif (gx1 == 1 && gz2 == 1)
                            % type = '2 solid faces at x1z2';
                            ui = 2*u2 * X1^2 * Z2 / dx^3;
                            vi = 4*v1 * X1*Z2 / dx^2 + 4*V * X1*Z2 * Y1 / dx^3;
                            wi = 2*w1 * X1 * Z2^2 / dx^3;
                        else % (gx2 == 1 && gz2 == 1)
                            % type = '2 solid faces at x2z2';
                            ui = 2*u1 * X2^2 * Z2 / dx^3;
                            vi = 4*v1 * X2*Z2 / dx^2 + 4*V * X2*Z2 * Y1 / dx^3;
                            wi = 2*w1 * X2 * Z2^2 / dx^3;
                        end
                        
                    else % gYZ == 2
                        %% Two Adjacent Voxels (YZ) Are Solid (4 cases)
                        if (gy1 == 1 && gz1 == 1)
                            % type = '2 solid faces at y1z1';
                            ui = 4*u1 * Y1*Z1 / dx^2 + 4*U * Y1*Z1 * X1 / dx^3 ;
                            vi = 2*v2 * Y1^2 * Z1 / dx^3;
                            wi = 2*w2 * Y1 *  Z1^2 / dx^3;
                        elseif (gy1 == 1 && gz2 == 1)
                            % type = '2 solid faces at y1z2';
                            ui = 4*u1 * Y1*Z2 / dx^2 + 4*U * Y1*Z2 * X1 / dx^3 ;
                            vi = 2*v2 * Y1^2 * Z2 / dx^3;
                            wi = 2*w1 * Y1 *  Z2^2 / dx^3;
                        elseif (gy2 == 1 && gz1 == 1)
                            % type = '2 solid faces at y2z1';
                            ui = 4*u1 * Y2*Z1 / dx^2 + 4*U * Y2*Z1 * X1 / dx^3 ;
                            vi = 2*v1 * Y2^2 * Z1 / dx^3;
                            wi = 2*w2 * Y2 *  Z1^2 / dx^3;
                        else % (gy2 == 1 && gz2 == 1)
                            % type = '2 solid faces at y2z2';
                            ui = 4*u1 * Y2*Z2 / dx^2 + 4*U * Y2*Z2 * X1 / dx^3 ;
                            vi = 2*v1 * Y2^2 * Z2 / dx^3;
                            wi = 2*w1 * Y2 *  Z2^2 / dx^3;
                        end
                        
                    end
                    
                elseif num_solid == 3
                    %% Three of the Neighbouring Voxels Are Solid (20 cases)
                    
                    if gX == 2
                        %% Two of Voxels Are in the Same X Direction (4 cases)
                        if gz1 == 1 % (gx2 == 1 && gx1 == 1 && gz1 == 1)
                            % type = '3 solid faces at x1x2-z1'; SAMPLE FORMULATION
                            ui =  0;
                            vi = 12*v1 * X1*X2*Z1 / dx^3 + 12*V * X1*X2*Z1 * Y1 / dx^4;
                            wi =  6*w2 * X1*X2*Z1^2 / dx^4;
                        elseif gz2 == 1 % (gx2 == 1 && gx1 == 1 && gz2 == 1)
                            % type = '3 solid faces at x1x2-z2'; %% MY BUG FOUND HERE
                            ui =  0;
                            vi = 12*v1 * X1*X2*Z2 / dx^3 + 12*V * X1*X2*Z2 * Y1 / dx^4;
                            wi =  6*w1 * X1*X2*Z2^2 / dx^4;
                        elseif gy1 == 1 % (gx2 == 1 && gx1 == 1 && gy1 == 1)
                            % type = '3 solid faces at x1x2-y1';
                            ui =  0;
                            vi =  6*v2 * X1*X2*Y1^2 / dx^4;
                            wi = 12*w1 * X1*X2*Y1 / dx^3 + 12*W * X1*X2*Y1 * Z1 / dx^4;
                        else % (gx2 == 1 && gx1 == 1 && gy2 == 1)
                            % type = '3 solid faces at x1x2-y2'; %% MY BUG FOUND HERE
                            ui =  0;
                            vi =  6*v1 * X1*X2*Y2^2 / dx^4;
                            wi = 12*w1 * X1*X2*Y2 / dx^3 + 12*W * X1*X2*Y2 * Z1 / dx^4;
                        end
                        
                    elseif gY == 2
                        %% Two of Voxels Are in the Same Y Direction (4 cases)
                        if gx1 == 1 % (gx1 == 1 && gy2 == 1 && gy1 == 1)
                            % type = '3 solid faces at y1y2-x1'; %% MY BUG FOUND HERE
                            ui =  6*u2 * Y1*Y2*X1^2 / dx^4;
                            vi =  0;
                            wi = 12*w1 * Y1*Y2*X1 / dx^3 + 12*W * Y1*Y2*X1 * Z1 / dx^4;
                        elseif gx2 == 1 % (gx2 == 1 && gy2 == 1 && gy1 == 1)
                            % type = '3 solid faces at y1y2-x2';
                            ui =  6*u1 * Y1*Y2*X2^2 / dx^4;
                            vi =  0;
                            wi = 12*w1 * Y1*Y2*X2 / dx^3 + 12*W * Y1*Y2*X2 * Z1 / dx^4;
                        elseif gz1 == 1 % (gy1 == 1 && gy2 == 1 && gz1 == 1)
                            % type = '3 solid faces at y1y2-z1'; %BUG FOUND HERE
                            ui = 12*u1 * Y1*Y2*Z1 / dx^3 + 12*U * Y1*Y2*Z1 * X1 / dx^4;
                            vi =  0;
                            wi =  6*w2 * Y1*Y2*Z1^2 / dx^4;
                        else % (gy2 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '3 solid faces at y1y2-z2';
                            ui = 12*u1 * Y1*Y2*Z2 / dx^3 + 12*U * Y1*Y2*Z2 * X1 / dx^4;
                            vi =  0;
                            wi =  6*w1 * Y1*Y2*Z2^2 / dx^4;
                        end
                        
                    elseif gZ == 2
                        %% Two of Voxels Are in the Same Z Direction (4 cases)
                        if gx1 == 1 % (gx1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '3 solid faces at z1z2-x1';
                            ui =  6*u2 * Z1*Z2*X1^2 / dx^4;
                            vi = 12*v1 * Z1*Z2*X1 / dx^3 + 12*V * Z1*Z2*X1 * Y1 / dx^4;
                            wi =  0;
                        elseif gx2 == 1 % (gx2 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '3 solid faces at z1z2-x2';
                            ui =  6*u1 * Z1*Z2*X2^2 / dx^4;
                            vi = 12*v1 * Z1*Z2*X2 / dx^3 + 12*V * Z1*Z2*X2 * Y1 / dx^4;
                            wi =  0;
                        elseif gy1 == 1 % (gy1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '3 solid faces at z1z2-y1';
                            ui = 12*u1 * Z1*Z2*Y1 / dx^3 + 12*U * Z1*Z2*Y1 * X1 / dx^4;
                            vi =  6*v2 * Z1*Z2*Y1^2 / dx^4;
                            wi =  0;
                        else % (gy2 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '3 solid faces at z1z2-y2';
                            ui = 12*u1 * Z1*Z2*Y2 / dx^3 + 12*U * Z1*Z2*Y2 * X1 / dx^4;
                            vi =  6*v1 * Z1*Z2*Y2^2 / dx^4;
                            wi =  0;
                        end
                        
                    else
                        %% Three block in diferent directions (8 cases)
                        if (gx2 == 1 && gy2 == 1 && gz1 == 1)
                            % type = '3 solid faces at x2y2z1'; SAMPLE FORMULATION
                            ui = 4*u1 * X2^2 * Y2 * Z1 / dx^4;
                            vi = 4*v1 * X2 * Y2^2 * Z1 / dx^4;
                            wi = 4*w2 * X2 * Y2 * Z1^2 / dx^4;
                        elseif (gx2 == 1 && gy2 == 1 && gz2 == 1)
                            % type = '3 solid faces at x2y2z2';
                            ui = 4*u1 * X2^2 * Y2 * Z2 / dx^4;
                            vi = 4*v1 * X2 * Y2^2 * Z2 / dx^4;
                            wi = 4*w1 * X2 * Y2 * Z2^2 / dx^4;
                        elseif (gx2 == 1 && gy1 == 1 && gz1 == 1)
                            % type = '3 solid faces at x2y1z1';
                            ui = 4*u1 * X2^2 * Y1 * Z1 / dx^4;
                            vi = 4*v2 * X2 * Y1^2 * Z1 / dx^4;
                            wi = 4*w2 * X2 * Y1 * Z1^2 / dx^4;
                        elseif (gx2 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '3 solid faces at x2y1z2';
                            ui = 4*u1 * X2^2 * Y1 * Z2 / dx^4;
                            vi = 4*v2 * X2 * Y1^2 * Z2 / dx^4;
                            wi = 4*w1 * X2 * Y1 * Z2^2 / dx^4;
                        elseif (gx1 == 1 && gy1 == 1 && gz1 == 1)
                            % type = '3 solid faces at x1y1z1';
                            ui = 4*u2 * X1^2 * Y1 * Z1 / dx^4;
                            vi = 4*v2 * X1 * Y1^2 * Z1 / dx^4;
                            wi = 4*w2 * X1 * Y1 * Z1^2 / dx^4;
                        elseif (gx1 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '3 solid faces at x1y1z2';
                            ui = 4*u2 * X1^2 * Y1 * Z2 / dx^4;
                            vi = 4*v2 * X1 * Y1^2 * Z2 / dx^4;
                            wi = 4*w1 * X1 * Y1 * Z2^2 / dx^4;
                        elseif (gx1 == 1 && gy2 == 1 && gz1 == 1)
                            % type = '3 solid faces at x1y2z1';
                            ui = 4*u2 * X1^2 * Y2 * Z1 / dx^4;
                            vi = 4*v1 * X1 * Y2^2 * Z1 / dx^4;
                            wi = 4*w2 * X1 * Y2 * Z1^2 / dx^4;
                        elseif (gx1 == 1 && gy2 == 1 && gz2 == 1)
                            % type = '3 solid faces at x1y2z2';
                            ui = 4*u2 * X1^2 * Y2 * Z2 / dx^4;
                            vi = 4*v1 * X1 * Y2^2 * Z2 / dx^4;
                            wi = 4*w1 * X1 * Y2 * Z2^2 / dx^4;
                        end
                        
                    end
                    
                elseif num_solid == 4
                    %% Four of the Neighbouring Voxels Are Solid (15 cases)
                    if (gX + gY == 4) || (gX + gZ == 4) || (gY + gZ == 4)
                        %% They block 2 directions (3 cases)
                        if (gX + gY == 4) % (gx2 == 1 && gx1 == 1 && gy2 == 1 && gy1 == 1)
                            % type = '4 solid faces at x1x2-y1y2'; % SAMPLE FORMULATION
                            ui = 0;
                            vi = 0;
                            wi = 36*w1 * X1*X2*Y1*Y2 / dx^4;
                        elseif (gX + gZ == 4) % (gx2 == 1 && gx1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at x1x2-z1z2';
                            ui = 0;
                            vi = 36*v1 * X1*X2*Z1*Z2 / dx^4;
                            wi = 0;
                        else % (gY + gZ == 4) % (gy2 == 1 && gy1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at y1y2-z1z2';
                            ui = 36*u1 * Y1*Y2*Z1*Z2 / dx^4;
                            vi = 0;
                            wi = 0;
                        end
                        
                    elseif gX == 2
                        %% Two in X Direction and Two in Different Directions (4 cases)
                        if (gy2 == 1 && gz1 == 1) % (gx2 == 1 && gx1 == 1 && gy2 == 1 && gz1 == 1)
                            % type = '4 solid faces at x1x2-y2z1'; SAMPLE FORMULATION
                            ui = 0;
                            vi = 12*v1 * X1*X2 * Y2^2 * Z1 / dx^5;
                            wi = 12*w2 * X1*X2 * Y2 * Z1^2 / dx^5;
                        elseif (gy2 == 1 && gz2 == 1) % (gx2 == 1 && gx1 == 1 && gy2 == 1 && gz2 == 1)
                            % type = '4 solid faces at x1x2-y2z2';
                            ui = 0;
                            vi = 12*v1 * X1*X2 * Y2^2 * Z2 / dx^5;
                            wi = 12*w1 * X1*X2 * Y2 * Z2^2 / dx^5;
                        elseif (gy1 == 1 && gz1 == 1) % (gx2 == 1 && gx1 == 1 && gy1 == 1 && gz1 == 1)
                            % type = '4 solid faces at x1x2-y1z1';
                            ui = 0;
                            vi = 12*v2 * X1*X2 * Y1^2 * Z1 / dx^5;
                            wi = 12*w2 * X1*X2 * Y1 * Z1^2 / dx^5;
                        else % (gx2 == 1 && gx1 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '4 solid faces at x1x2-y1z2';
                            ui = 0;
                            vi = 12*v2 * X1*X2 * Y1^2 * Z2 / dx^5;
                            wi = 12*w1 * X1*X2 * Y1 * Z2^2 / dx^5;
                        end
                        
                    elseif gY == 2
                        %% Two in Y Direction and Two in Different Directions (4 cases)
                        if (gx1 == 1 && gz1 == 1) % (gx1 == 1 && gy2 == 1 && gy1 == 1 && gz1 == 1)
                            % type = '4 solid faces at y1y2-x1z1';
                            ui = 12*u2 * Y1*Y2 * X1^2 * Z1 / dx^5;
                            vi = 0;
                            wi = 12*w2 * Y1*Y2 * X1 * Z1^2 / dx^5;
                        elseif (gx1 == 1 && gz2 == 1) % (gx1 == 1 && gy2 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '4 solid faces at y1y2-x1z2';
                            ui = 12*u2 * Y1*Y2 * X1^2 * Z2 / dx^5;
                            vi = 0;
                            wi = 12*w1 * Y1*Y2 * X1 * Z2^2 / dx^5;
                        elseif (gx2 == 1 && gz1 == 1) % (gx2 == 1 && gy2 == 1 && gy1 == 1 && gz1 == 1)
                            % type = '4 solid faces at y1y2-x2z1';
                            ui = 12*u1 * Y1*Y2 * X2^2 * Z1 / dx^5;
                            vi = 0;
                            wi = 12*w2 * Y1*Y2 * X2 * Z1^2 / dx^5;
                        else % (gx2 == 1 && gy2 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '4 solid faces at y1y2-x2z2';
                            ui = 12*u1 * Y1*Y2 * X2^2 * Z2 / dx^5;
                            vi = 0;
                            wi = 12*w1 * Y1*Y2 * X2 * Z2^2 / dx^5;
                        end
                        
                    else % gZ == 2
                        %% Two in Z Direction and Two in Different Directions (4 cases)
                        if (gx1 == 1 && gy1 == 1) % (gx1 == 1 && gy1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at z1z2-x1y1';
                            ui = 12*u2 * Z1*Z2 * X1^2 * Y1 / dx^5;
                            vi = 12*v2 * Z1*Z2 * X1 * Y1^2 / dx^5;
                            wi = 0;
                        elseif (gx1 == 1 && gy2 == 1) % (gx1 == 1 && gy2 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at z1z2-x1y2';
                            ui = 12*u2 * Z1*Z2 * X1^2 * Y2 / dx^5;
                            vi = 12*v1 * Z1*Z2 * X1 * Y2^2 / dx^5;
                            wi = 0;
                        elseif (gx2 == 1 && gy1 == 1) % (gx2 == 1 && gy1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at z1z2-x2y1';
                            ui = 12*u1 * Z1*Z2 * X2^2 * Y1 / dx^5;
                            vi = 12*v2 * Z1*Z2 * X2 * Y1^2 / dx^5;
                            wi = 0;
                        else % (gx2 == 1 && gy2 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at z1z2-x2y2';
                            ui = 12*u1 * Z1*Z2 * X2^2 * Y2 / dx^5;
                            vi = 12*v1 * Z1*Z2 * X2 * Y2^2 / dx^5;
                            wi = 0;
                        end
                        
                    end
                    
                else
                    %% type = '4, 5 or 6 solid faces';
                    ui = 0;
                    vi = 0;
                    wi = 0;
                end
                
                %% Advection displacement
                advectionX = ui * dt; x = x + advectionX;
                advectionY = vi * dt; y = y + advectionY;
                advectionZ = wi * dt; z = z + advectionZ;
                
                if x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                    % if particle moves out system -> reinject
                    %% ReinjectAdvection
                    % Reinject (during advection)
                    % If a particle exits the system during its advection,
                    % it will be reinject at the inlet face at a flow-weighted random location
                    %
                    % Workspace variables
                    %   x,y,z  : current location of particle
                    %   dx     : dimension of the cell (dx=dy=dz)
                    %   gL     : sample length
                    %   gFlow  : sample main flow direction
                    %   Uinlet : velocity field of the inlet faces (X direction)
                    %   Vinlet : velocity field of the inlet faces (Y direction)
                    %   Winlet : velocity field of the inlet faces (Z direction)
                    %
                    % Output
                    %   x,y,z  : Updated location of particle
                    %
                    % Functions used:
                    %   out_system : Check if particle out system or not
                    %
                    % Constants Used:
                    %   FLOW_X, FLOW_Y, FLOW_Z
                    %

                    if gFlow == FLOW_X
                        %% [x,y,z] = reinject_advection_X (Uinlet,x,y,z,dx,gL);
                        velField = Uinlet;
                        while x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                            random_index = randsample(velField.index, 1, true, velField.list);
                            [~,j,k] = ind2sub (size(velField.plane),random_index);
                            x = (velField.face-0.5) * dx;
                            y = (j-0.5) * dx;
                            z = (k-0.5) * dx;
                        end
                        event.ReinjectAdvectionX = event.ReinjectAdvectionX + 1;
                    elseif gFlow == FLOW_Y
                        %% [x,y,z] = reinject_advection_Y (Vinlet,x,y,z,dx,gL);
                        velField = Vinlet;
                        while x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                            random_index = randsample(velField.index, 1, true, velField.list);
                            [i,~,k] = ind2sub (size(velField.plane),random_index);
                            x = (i-0.5) * dx;
                            y = (velField.face-0.5) * dx;
                            z = (k-0.5) * dx;
                        end
                        event.ReinjectAdvectionY = event.ReinjectAdvectionY + 1;
                    else % gFlow == FLOW_Z
                        %% [x,y,z] = reinject_advection_Z (Winlet,x,y,z,dx,gL);
                        velField = Winlet;
                        while x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                            random_index = randsample(velField.index, 1, true, velField.list);
                            [i,j,~] = ind2sub(size(velField.plane),random_index);
                            x = (i-0.5) * dx;
                            y = (j-0.5) * dx;
                            z = (velField.face-0.5) * dx;
                        end
                        event.ReinjectAdvectionZ = event.ReinjectAdvectionZ + 1;
                    end
                    
                else
                    % advection movement within the system
                    i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);%[i,j,k] = get_ijk (x,y,z,dx);
                    if g(i,j,k) == 1
                        % if particle hits a solid surface -> NO MOVE
                        x = xBA; y = yBA; z = zBA; % restore location
                        advectionX = 0; advectionY = 0; advectionZ = 0; % update advection displacement
                        event.AdvectionHitSolid = event.AdvectionHitSolid + 1;
                    end
                end

                %% DiffusionMovement
                % Diffusion Movement
                % Move particle by using a random-walk method: the particle instantaneously
                % jumps over a fixed distance (zeta) in a random direction
                %
                % Workspace variables
                %   g     : geometry
                %   dx    : dimension of the cell (dx=dy=dz)
                %   dt    : time-step
                %   gL    : sample length
                %   zeta  : particle jump distance
                %   x,y,z : current location of particle
                %
                % Output
                %   x,y,z : Update location of particle
                %   diffusionX, diffusionY, diffusionZ : Diffusion displacement
                %
                % Scipts Used:
                %   ReinjectDiffusion : Update x,y,z
                %
                % Functions Used:
                %   out_system : Check if particle out system or not
                %   get_ijk    : Voxel (i,j,k) of particle location (x,y,z)
                %
                
                % Particle location before diffusion
                xBD = x; yBD = y; zBD = z;
                
                % Randomize theta & phi
                theta = abs (2*rand()*pi);
                phi = abs (rand()*pi);
                
                % Diffustion displacement
                diffusionX = zeta * sin(phi) * cos(theta); x = x + diffusionX;
                diffusionY = zeta * sin(phi) * sin(theta); y = y + diffusionY;
                diffusionZ = zeta * cos(phi);              z = z + diffusionZ;
                                            
                if x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                    % if particle jumps out system -> reinject
                    %% ReinjectDiffusion
                    random_index = ceil(Nparticle*rand());
                    x = xpt0(random_index);
                    y = ypt0(random_index);
                    z = zpt0(random_index);
                    if x >= dx || x <= gL
                        event.ReinjectDiffusionX = event.ReinjectDiffusionX + 1;
                    elseif y >= dx || y <= gL
                        event.ReinjectDiffusionY = event.ReinjectDiffusionY + 1;
                    elseif z >= dx || z <= gL
                        event.ReinjectDiffusionZ = event.ReinjectDiffusionZ + 1;
                    end
                    
                else
                    %% diffusion movement within the system
                    i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);%[i,j,k] = get_ijk (x,y,z,dx);
                    if g(i,j,k) == 1
                        % if hit a surface -> bounce movement
                        %[x,y,z] = bounce_movement(i,j,k,xpt1,ypt1,zpt1,...
                        %    x,y,z,diffusionX,diffusionY,diffusionZ,dx);
                        
                        i1=ceil(xBD/dx); j1=ceil(yBD/dx); k1=ceil(zBD/dx);%[i1,j1,k1] = get_ijk (xpt1,ypt1,zpt1,dx);

                        if (i ~= i1) %bounce X direction
                            if (diffusionX >= 0)
                                bounceX = - 2*abs(x - (dx*floor(x/dx)));
                            else
                                bounceX =   2*abs(x - (dx*ceil(x/dx)));
                            end
                            x = x + bounceX;

                        elseif (j ~= j1) %bounce Y direction
                            if (diffusionY >= 0)
                                bounceY = - 2*abs(y - (dx*floor(y/dx)));
                            else
                                bounceY =   2*abs(y - (dx*ceil(y/dx)));
                            end
                            y = y + bounceY;

                        elseif (k ~= k1) %bounce Z direction
                            if (diffusionZ >= 0)
                                bounceZ = - 2*abs(z - (dx*floor(z/dx)));
                            else
                                bounceZ =   2*abs(z - (dx*ceil(z/dx)));
                            end
                            z = z + bounceZ;

                        end
                        
                        diffusionX = x - xBD; diffusionY = y - zBD; diffusionZ = z - zBD; % update diffusion
                                                
                        if x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                            % after bouncing: if particle bounces out system -> reinject again
                            %% ReinjectDiffusion
                            random_index = ceil(Nparticle*rand());
                            x = xpt0(random_index);
                            y = ypt0(random_index);
                            z = zpt0(random_index);
                            if x >= dx || x <= gL
                                event.ReinjectDiffusionX = event.ReinjectDiffusionX + 1;
                            elseif y >= dx || y <= gL
                                event.ReinjectDiffusionY = event.ReinjectDiffusionY + 1;
                            elseif z >= dx || z <= gL
                                event.ReinjectDiffusionZ = event.ReinjectDiffusionZ + 1;
                            end
                            
                        else
                            %% bounce movement within system
                            i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);%[i,j,k] = get_ijk (x,y,z,dx);
                            if g(i,j,k) == 1
                                % if bounce movement into a solid surface -> NO MOVE
                                x = xBD; y = yBD; z = zBD; % restore location
                                diffusionX = 0; diffusionY = 0; diffusionZ = 0;% update diffusion
                            end
                        end
                        event.DiffusionHitSolid = event.DiffusionHitSolid + 1;
                    end
                end
                                                                         
                %% Update particle displacement
                if gFlow == FLOW_X
                    displacement1(p) = displacement1(p) + advectionX + diffusionX;
                elseif gFlow == FLOW_Y
                    displacement1(p) = displacement1(p) + advectionY + diffusionZ;
                else % gFlow == FLOW_Z
                    displacement1(p) = displacement1(p) + advectionZ + diffusionZ;
                end
                
                xpt(p) = x; ypt(p) = y; zpt(p) = z; % copy short variables back to array(index)
            end % for Nparticle
        end
        
        %% update variance & calculate dispersion coefficient
        variance1(counter) = var(displacement1); % Update variance
        varianceTime(counter) = varianceTime(counter-1) + slope;
        dispersion1(counter) = 0.5 * (variance1(counter) - variance1(counter-1)) / (slope*dt);
        dispersionChange(counter) = abs(dispersion1(counter)-dispersion1(counter-1))/dispersion1(counter);
        if dispersionChange(counter-1) <= D_LIMIT && limitBeforeMinTime == 0
            limitBeforeMinTime = counter;
        end
        
        %% display progress
        progressMsg = sprintf(' -> Sim1: time=%dx%d, RDC=%.2f, change=%.3f ',...
            counter-1, slope, dispersion1(counter)/Dm, dispersionChange(counter));
        fprintf([reverseStr, progressMsg]);
        reverseStr = repmat(sprintf('\b'), 1, length(progressMsg));  
        
        counter = counter + 1;
    end % while
    
    variance1(counter:end) = [];
    varianceTime(counter:end) = [];
    dispersion1(counter:end) = [];
    dispersionChange(counter:end) = [];
    
    if dispersionChange(end) > D_LIMIT
        fprintf(sprintf('* MAX_TIMESTEP * '));
    end
    
    if limitBeforeMinTime ~= 0
        fprintf(sprintf('* D_LIMIT at %dx%d * ',limitBeforeMinTime-1,slope));
    end
    
    toc;
    xpt1 = xpt; ypt1 = ypt; zpt1 = zpt; % save xpt,ypt,zpt of simulation1
     
    %% Simulation2 (microporosity) with Ntimestep
    load(particleSet);  % re-load particle
    tic; 
    % Simulation with Microporosity (Ntimestep)
    %
    % Workspace Variables
    %   Ntimestep      : number of timestep
    %   Nparticle      : number of particle
    %   xpt,ypt,zpt    : particles locations
    %   slope          : number of timestep to calculate dispersion
    %
    % Output
    %   xpt,ypt,zpt    : particles locations
    %   displacement2  : particles displacement in the main flow
    %   variance2      : variance of displacement (compute at each slope time)
    %   dispersion2    : longitudinal dispersion coefficient
    
    Ntimestep = varianceTime(end);
    
    % Simulation variables
    displacement2 = zeros(Nparticle,1); % particles displacement in the main flow
    variance2 = zeros(length(varianceTime),1); % variance of displacement
    dispersion2 = zeros(length(varianceTime),1); % variance of displacement
    insolid = zeros (length(varianceTime),1); % mean(inside_solid)
    
    % Microporosity variables
    inside_solid = false (Nparticle,1); % current status of each particle
    
    % Simulation starts
    reverseStr = ''; % for diplay progress
    counter = 1;
    for t = 1:Ntimestep
        for p = 1:Nparticle
            x = xpt(p); y = ypt(p); z = zpt(p); % short variables instead of array(index)
            
            if inside_solid(p) == false
                %% AdvectionMovement
                % Avection Movement
                % Move particle along streamlines that follow velocity field
                %
                % Workspace variables
                %   g     : geometry
                %   dx    : dimension of the cell (dx=dy=dz)
                %   dt    : time-step
                %   gL    : sample length
                %   x,y,z : current location of particle
                %
                % Output
                %   x,y,z : Update location of particle
                %   advectionX, advectionY, advectionZ : Advection displacement
                %
                % Scipts Used:
                %   StreamlineTracing: Calculate ui,vi,wi
                %   ReinjectAdvection: Update x,y,z
                %
                % Functions Used:
                %   out_system : Check if particle out system or not
                %   get_ijk    : Voxel (i,j,k) of particle location (x,y,z)
                %
                
                % Particle location before advection
                xBA = x; yBA = y; zBA = z;
                
                % Streamline-tracing algorithm
                %% StreamlineTracing
                % Streamline-tracing algorithm
                % Semi-analytic formulation for the streamlines in voxels contain a solid boundary
                %
                % Workspace variables
                %   g        : geometry
                %   dx       : dimension of the cell (dx=dy=dz)
                %   x,y,z    : Location of particle
                %
                % Output
                %   ui,vi,wi : interpolated velocity between the normal velocities on opposite faces of gridblock
                %
                % Restriction:
                %   Only consider geometry from (2) --> (N-1)
                
                % Get particle boundary conditions
                
                % index location of the particle
                i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);%[i,j,k] = get_ijk (x,y,z,dx);
                
                % location of the edges
                x1 = dx*(i-1); x2 = dx*i;
                y1 = dx*(j-1); y2 = dx*j;
                z1 = dx*(k-1); z2 = dx*k;
                
                % velocities on the faces normal to each direction
                u1 = u(i-1,j,k); u2 = u(i,j,k);
                v1 = v(i,j-1,k); v2 = v(i,j,k);
                w1 = w(i,j,k-1); w2 = w(i,j,k);
                
                % short variables names
                X1 = x-x1; X2 = x2-x;
                Y1 = y-y1; Y2 = y2-y;
                Z1 = z-z1; Z2 = z2-z;
                U = u2-u1; V = v2-v1; W = w2-w1;
                
                % geometry of the edge
                gx1 = g(i-1,j,k); gx2 = g(i+1,j,k);
                gy1 = g(i,j-1,k); gy2 = g(i,j+1,k);
                gz1 = g(i,j,k-1); gz2 = g(i,j,k+1);
                
                % short variables names
                gX = gx1+gx2; gY = gy1+gy2; gZ = gz1+gz2;
                gXY = gX+gY; gXZ = gX+gZ; gYZ = gY+gZ;
                
                % number of solid boundary arround the particle
                num_solid = gx1 + gx2 + gy1 + gy2 + gz1 + gz2;
                
                % Calculate interpolated velocities: ui,vi,wi
                
                if num_solid == 0
                    %% No solid boundaries
                    ui = U * X1 / dx + u1;
                    vi = V * Y1 / dx + v1;
                    wi = W * Z1 / dx + w1;
                    
                elseif num_solid == 1
                    %% One of the Neighbouring Voxels Is Solid (6 cases)
                    if gx2 == 1
                        % type = '1 solid face at x2'; SAMPLE FORMULATION
                        ui = u1 * X2^2 / dx^2;
                        vi = 2*v1 * X2 / dx + 2*V * X2*Y1 / dx^2;
                        wi = 2*w1 * X2 / dx + 2*W * X2*Z1 / dx^2;
                    elseif gx1 == 1
                        % type = '1 solid face at x1';
                        ui = u2 * X1^2 / dx^2;
                        vi = 2*v1 * X1 / dx + 2*V * X1*Y1 / dx^2;
                        wi = 2*w1 * X1 / dx + 2*W * X1*Z1 / dx^2;
                    elseif gy2 == 1
                        % type = '1 solid face at y2';
                        ui = 2*u1 * Y2 / dx + 2*U * Y2*X1 / dx^2;
                        vi = v1 * Y2^2 / dx^2;
                        wi = 2*w1 * Y2 / dx + 2*W * Y2*Z1 / dx^2;
                    elseif gy1 == 1
                        % type = '1 solid face at y1';
                        ui = 2*u1 * Y1 / dx + 2*U * Y1*X1 / dx^2;
                        vi = v2 * Y1^2 / dx^2;
                        wi = 2*w1 * Y1 / dx + 2*W * Y1*Z1 / dx^2;
                    elseif gz2 == 1
                        % type = '1 solid face at z2';
                        ui = 2*u1 * Z2 / dx + 2*U * Z2*X1 / dx^2;
                        vi = 2*v1 * Z2 / dx + 2*V * Z2*Y1 / dx^2;
                        wi = w1 * Z2^2 / dx^2;
                    elseif gz1 == 1
                        % type = '1 solid face at z1';
                        ui = 2*u1 * Z1 / dx + 2*U * Z1*X1 / dx^2;
                        vi = 2*v1 * Z1 / dx + 2*V * Z1*Y1 / dx^2;
                        wi = w2 * Z1^2 / dx^2;
                    end
                    
                elseif num_solid == 2
                    %% Two of the Neighbouring Voxels Are Solid (15 cases)
                    
                    if (gX == 2 || gY == 2 || gZ == 2)
                        %% There Are Two Opposing Solid Faces (3 cases)
                        if gX==2 % (gx2 == 1 && gx1 == 1)
                            % type = '2 opposing solid faces x direction'; % SAMPLE FORMULATION
                            ui = 0;
                            vi = 6*v1 * X2*X1 / dx^2 + 6*V * X2*X1 * Y1 / dx^3;
                            wi = 6*w1 * X2*X1 / dx^2 + 6*W * X2*X1 * Z1 / dx^3;
                        elseif gY==2 % (gy2 == 1 && gy1 == 1)
                            % type = '2 opposing solid faces y direction';
                            ui = 6*u1 * Y2*Y1 / dx^2 + 6*U * Y2*Y1 * X1 / dx^3;
                            vi = 0;
                            wi = 6*w1 * Y2*Y1 / dx^2 + 6*W * Y2*Y1 * Z1 / dx^3;
                        else % gz==2 % (gz2 == 1 && gz1 == 1)
                            % type = '2 opposing solid faces z direction';
                            ui = 6*u1 * Z2*Z1 / dx^2 + 6*U * Z2*Z1 * X1 / dx^3;
                            vi = 6*v1 * Z2*Z1 / dx^2 + 6*V * Z2*Z1 * Y1 / dx^3;
                            wi = 0;
                        end
                        
                    elseif gXY == 2
                        %% Two Adjacent Voxels (XY) Are Solid (4 cases), PAPER BUG
                        if (gx1 == 1 && gy1 == 1)
                            % type = '2 solid faces at x1y1'; % SAMPLE FORMULATION
                            ui = 2*u2 * X1^2 * Y1 / dx^3;
                            vi = 2*v2 * X1 * Y1^2 / dx^3;
                            wi = 4*w1 * X1*Y1 / dx^2 + 4*W * X1*Y1 * Z1 / dx^3;
                        elseif (gx2 == 1 && gy1 == 1)
                            % type = '2 solid faces at x2y1';
                            ui = 2*u1 * X2^2 * Y1 / dx^3;
                            vi = 2*v2 * X2 * Y1^2 / dx^3;
                            wi = 4*w1 * X2*Y1 / dx^2 + 4*W * X2*Y1 * Z1 / dx^3;
                        elseif (gx1 == 1 && gy2 == 1)
                            % type = '2 solid faces at x1y2';
                            ui = 2*u2 * X1^2 * Y2 / dx^3;
                            vi = 2*v1 * X1 * Y2^2 / dx^3;
                            wi = 4*w1 * X1*Y2 / dx^2 + 4*W * X1*Y2 * Z1 / dx^3;
                        else % (gx2 == 1 && gy2 == 1)
                            % type = '2 solid faces at x2y2';
                            ui = 2*u1 * X2^2 * Y2 / dx^3;
                            vi = 2*v1 * X2 * Y2^2 / dx^3;
                            wi = 4*w1 * X2*Y2 / dx^2 + 4*W * X2*Y2 * Z1 / dx^3;
                        end
                        
                    elseif gXZ == 2
                        %% Two Adjacent Voxels (XZ) Are Solid (4 cases)
                        if (gx1 == 1 && gz1 == 1)
                            % type = '2 solid faces at x1z1';
                            ui = 2*u2 * X1^2 * Z1 / dx^3;
                            vi = 4*v1 * X1*Z1 / dx^2 + 4*V * X1*Z1 * Y1 / dx^3;
                            wi = 2*w2 * X1 * Z1^2 / dx^3;
                        elseif (gx2 == 1 && gz1 == 1)
                            % type = '2 solid faces at x2z1';
                            ui = 2*u1 * X2^2 * Z1 / dx^3;
                            vi = 4*v1 * X2*Z1 / dx^2 + 4*V * X2*Z1 * Y1 / dx^3;
                            wi = 2*w2 * X2 * Z1^2 / dx^3;
                        elseif (gx1 == 1 && gz2 == 1)
                            % type = '2 solid faces at x1z2';
                            ui = 2*u2 * X1^2 * Z2 / dx^3;
                            vi = 4*v1 * X1*Z2 / dx^2 + 4*V * X1*Z2 * Y1 / dx^3;
                            wi = 2*w1 * X1 * Z2^2 / dx^3;
                        else % (gx2 == 1 && gz2 == 1)
                            % type = '2 solid faces at x2z2';
                            ui = 2*u1 * X2^2 * Z2 / dx^3;
                            vi = 4*v1 * X2*Z2 / dx^2 + 4*V * X2*Z2 * Y1 / dx^3;
                            wi = 2*w1 * X2 * Z2^2 / dx^3;
                        end
                        
                    else % gYZ == 2
                        %% Two Adjacent Voxels (YZ) Are Solid (4 cases)
                        if (gy1 == 1 && gz1 == 1)
                            % type = '2 solid faces at y1z1';
                            ui = 4*u1 * Y1*Z1 / dx^2 + 4*U * Y1*Z1 * X1 / dx^3 ;
                            vi = 2*v2 * Y1^2 * Z1 / dx^3;
                            wi = 2*w2 * Y1 *  Z1^2 / dx^3;
                        elseif (gy1 == 1 && gz2 == 1)
                            % type = '2 solid faces at y1z2';
                            ui = 4*u1 * Y1*Z2 / dx^2 + 4*U * Y1*Z2 * X1 / dx^3 ;
                            vi = 2*v2 * Y1^2 * Z2 / dx^3;
                            wi = 2*w1 * Y1 *  Z2^2 / dx^3;
                        elseif (gy2 == 1 && gz1 == 1)
                            % type = '2 solid faces at y2z1';
                            ui = 4*u1 * Y2*Z1 / dx^2 + 4*U * Y2*Z1 * X1 / dx^3 ;
                            vi = 2*v1 * Y2^2 * Z1 / dx^3;
                            wi = 2*w2 * Y2 *  Z1^2 / dx^3;
                        else % (gy2 == 1 && gz2 == 1)
                            % type = '2 solid faces at y2z2';
                            ui = 4*u1 * Y2*Z2 / dx^2 + 4*U * Y2*Z2 * X1 / dx^3 ;
                            vi = 2*v1 * Y2^2 * Z2 / dx^3;
                            wi = 2*w1 * Y2 *  Z2^2 / dx^3;
                        end
                        
                    end
                    
                elseif num_solid == 3
                    %% Three of the Neighbouring Voxels Are Solid (20 cases)
                    
                    if gX == 2
                        %% Two of Voxels Are in the Same X Direction (4 cases)
                        if gz1 == 1 % (gx2 == 1 && gx1 == 1 && gz1 == 1)
                            % type = '3 solid faces at x1x2-z1'; SAMPLE FORMULATION
                            ui =  0;
                            vi = 12*v1 * X1*X2*Z1 / dx^3 + 12*V * X1*X2*Z1 * Y1 / dx^4;
                            wi =  6*w2 * X1*X2*Z1^2 / dx^4;
                        elseif gz2 == 1 % (gx2 == 1 && gx1 == 1 && gz2 == 1)
                            % type = '3 solid faces at x1x2-z2'; %% MY BUG FOUND HERE
                            ui =  0;
                            vi = 12*v1 * X1*X2*Z2 / dx^3 + 12*V * X1*X2*Z2 * Y1 / dx^4;
                            wi =  6*w1 * X1*X2*Z2^2 / dx^4;
                        elseif gy1 == 1 % (gx2 == 1 && gx1 == 1 && gy1 == 1)
                            % type = '3 solid faces at x1x2-y1';
                            ui =  0;
                            vi =  6*v2 * X1*X2*Y1^2 / dx^4;
                            wi = 12*w1 * X1*X2*Y1 / dx^3 + 12*W * X1*X2*Y1 * Z1 / dx^4;
                        else % (gx2 == 1 && gx1 == 1 && gy2 == 1)
                            % type = '3 solid faces at x1x2-y2'; %% MY BUG FOUND HERE
                            ui =  0;
                            vi =  6*v1 * X1*X2*Y2^2 / dx^4;
                            wi = 12*w1 * X1*X2*Y2 / dx^3 + 12*W * X1*X2*Y2 * Z1 / dx^4;
                        end
                        
                    elseif gY == 2
                        %% Two of Voxels Are in the Same Y Direction (4 cases)
                        if gx1 == 1 % (gx1 == 1 && gy2 == 1 && gy1 == 1)
                            % type = '3 solid faces at y1y2-x1'; %% MY BUG FOUND HERE
                            ui =  6*u2 * Y1*Y2*X1^2 / dx^4;
                            vi =  0;
                            wi = 12*w1 * Y1*Y2*X1 / dx^3 + 12*W * Y1*Y2*X1 * Z1 / dx^4;
                        elseif gx2 == 1 % (gx2 == 1 && gy2 == 1 && gy1 == 1)
                            % type = '3 solid faces at y1y2-x2';
                            ui =  6*u1 * Y1*Y2*X2^2 / dx^4;
                            vi =  0;
                            wi = 12*w1 * Y1*Y2*X2 / dx^3 + 12*W * Y1*Y2*X2 * Z1 / dx^4;
                        elseif gz1 == 1 % (gy1 == 1 && gy2 == 1 && gz1 == 1)
                            % type = '3 solid faces at y1y2-z1'; %BUG FOUND HERE
                            ui = 12*u1 * Y1*Y2*Z1 / dx^3 + 12*U * Y1*Y2*Z1 * X1 / dx^4;
                            vi =  0;
                            wi =  6*w2 * Y1*Y2*Z1^2 / dx^4;
                        else % (gy2 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '3 solid faces at y1y2-z2';
                            ui = 12*u1 * Y1*Y2*Z2 / dx^3 + 12*U * Y1*Y2*Z2 * X1 / dx^4;
                            vi =  0;
                            wi =  6*w1 * Y1*Y2*Z2^2 / dx^4;
                        end
                        
                    elseif gZ == 2
                        %% Two of Voxels Are in the Same Z Direction (4 cases)
                        if gx1 == 1 % (gx1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '3 solid faces at z1z2-x1';
                            ui =  6*u2 * Z1*Z2*X1^2 / dx^4;
                            vi = 12*v1 * Z1*Z2*X1 / dx^3 + 12*V * Z1*Z2*X1 * Y1 / dx^4;
                            wi =  0;
                        elseif gx2 == 1 % (gx2 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '3 solid faces at z1z2-x2';
                            ui =  6*u1 * Z1*Z2*X2^2 / dx^4;
                            vi = 12*v1 * Z1*Z2*X2 / dx^3 + 12*V * Z1*Z2*X2 * Y1 / dx^4;
                            wi =  0;
                        elseif gy1 == 1 % (gy1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '3 solid faces at z1z2-y1';
                            ui = 12*u1 * Z1*Z2*Y1 / dx^3 + 12*U * Z1*Z2*Y1 * X1 / dx^4;
                            vi =  6*v2 * Z1*Z2*Y1^2 / dx^4;
                            wi =  0;
                        else % (gy2 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '3 solid faces at z1z2-y2';
                            ui = 12*u1 * Z1*Z2*Y2 / dx^3 + 12*U * Z1*Z2*Y2 * X1 / dx^4;
                            vi =  6*v1 * Z1*Z2*Y2^2 / dx^4;
                            wi =  0;
                        end
                        
                    else
                        %% Three block in diferent directions (8 cases)
                        if (gx2 == 1 && gy2 == 1 && gz1 == 1)
                            % type = '3 solid faces at x2y2z1'; SAMPLE FORMULATION
                            ui = 4*u1 * X2^2 * Y2 * Z1 / dx^4;
                            vi = 4*v1 * X2 * Y2^2 * Z1 / dx^4;
                            wi = 4*w2 * X2 * Y2 * Z1^2 / dx^4;
                        elseif (gx2 == 1 && gy2 == 1 && gz2 == 1)
                            % type = '3 solid faces at x2y2z2';
                            ui = 4*u1 * X2^2 * Y2 * Z2 / dx^4;
                            vi = 4*v1 * X2 * Y2^2 * Z2 / dx^4;
                            wi = 4*w1 * X2 * Y2 * Z2^2 / dx^4;
                        elseif (gx2 == 1 && gy1 == 1 && gz1 == 1)
                            % type = '3 solid faces at x2y1z1';
                            ui = 4*u1 * X2^2 * Y1 * Z1 / dx^4;
                            vi = 4*v2 * X2 * Y1^2 * Z1 / dx^4;
                            wi = 4*w2 * X2 * Y1 * Z1^2 / dx^4;
                        elseif (gx2 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '3 solid faces at x2y1z2';
                            ui = 4*u1 * X2^2 * Y1 * Z2 / dx^4;
                            vi = 4*v2 * X2 * Y1^2 * Z2 / dx^4;
                            wi = 4*w1 * X2 * Y1 * Z2^2 / dx^4;
                        elseif (gx1 == 1 && gy1 == 1 && gz1 == 1)
                            % type = '3 solid faces at x1y1z1';
                            ui = 4*u2 * X1^2 * Y1 * Z1 / dx^4;
                            vi = 4*v2 * X1 * Y1^2 * Z1 / dx^4;
                            wi = 4*w2 * X1 * Y1 * Z1^2 / dx^4;
                        elseif (gx1 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '3 solid faces at x1y1z2';
                            ui = 4*u2 * X1^2 * Y1 * Z2 / dx^4;
                            vi = 4*v2 * X1 * Y1^2 * Z2 / dx^4;
                            wi = 4*w1 * X1 * Y1 * Z2^2 / dx^4;
                        elseif (gx1 == 1 && gy2 == 1 && gz1 == 1)
                            % type = '3 solid faces at x1y2z1';
                            ui = 4*u2 * X1^2 * Y2 * Z1 / dx^4;
                            vi = 4*v1 * X1 * Y2^2 * Z1 / dx^4;
                            wi = 4*w2 * X1 * Y2 * Z1^2 / dx^4;
                        elseif (gx1 == 1 && gy2 == 1 && gz2 == 1)
                            % type = '3 solid faces at x1y2z2';
                            ui = 4*u2 * X1^2 * Y2 * Z2 / dx^4;
                            vi = 4*v1 * X1 * Y2^2 * Z2 / dx^4;
                            wi = 4*w1 * X1 * Y2 * Z2^2 / dx^4;
                        end
                        
                    end
                    
                elseif num_solid == 4
                    %% Four of the Neighbouring Voxels Are Solid (15 cases)
                    if (gX + gY == 4) || (gX + gZ == 4) || (gY + gZ == 4)
                        %% They block 2 directions (3 cases)
                        if (gX + gY == 4) % (gx2 == 1 && gx1 == 1 && gy2 == 1 && gy1 == 1)
                            % type = '4 solid faces at x1x2-y1y2'; % SAMPLE FORMULATION
                            ui = 0;
                            vi = 0;
                            wi = 36*w1 * X1*X2*Y1*Y2 / dx^4;
                        elseif (gX + gZ == 4) % (gx2 == 1 && gx1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at x1x2-z1z2';
                            ui = 0;
                            vi = 36*v1 * X1*X2*Z1*Z2 / dx^4;
                            wi = 0;
                        else % (gY + gZ == 4) % (gy2 == 1 && gy1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at y1y2-z1z2';
                            ui = 36*u1 * Y1*Y2*Z1*Z2 / dx^4;
                            vi = 0;
                            wi = 0;
                        end
                        
                    elseif gX == 2
                        %% Two in X Direction and Two in Different Directions (4 cases)
                        if (gy2 == 1 && gz1 == 1) % (gx2 == 1 && gx1 == 1 && gy2 == 1 && gz1 == 1)
                            % type = '4 solid faces at x1x2-y2z1'; SAMPLE FORMULATION
                            ui = 0;
                            vi = 12*v1 * X1*X2 * Y2^2 * Z1 / dx^5;
                            wi = 12*w2 * X1*X2 * Y2 * Z1^2 / dx^5;
                        elseif (gy2 == 1 && gz2 == 1) % (gx2 == 1 && gx1 == 1 && gy2 == 1 && gz2 == 1)
                            % type = '4 solid faces at x1x2-y2z2';
                            ui = 0;
                            vi = 12*v1 * X1*X2 * Y2^2 * Z2 / dx^5;
                            wi = 12*w1 * X1*X2 * Y2 * Z2^2 / dx^5;
                        elseif (gy1 == 1 && gz1 == 1) % (gx2 == 1 && gx1 == 1 && gy1 == 1 && gz1 == 1)
                            % type = '4 solid faces at x1x2-y1z1';
                            ui = 0;
                            vi = 12*v2 * X1*X2 * Y1^2 * Z1 / dx^5;
                            wi = 12*w2 * X1*X2 * Y1 * Z1^2 / dx^5;
                        else % (gx2 == 1 && gx1 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '4 solid faces at x1x2-y1z2';
                            ui = 0;
                            vi = 12*v2 * X1*X2 * Y1^2 * Z2 / dx^5;
                            wi = 12*w1 * X1*X2 * Y1 * Z2^2 / dx^5;
                        end
                        
                    elseif gY == 2
                        %% Two in Y Direction and Two in Different Directions (4 cases)
                        if (gx1 == 1 && gz1 == 1) % (gx1 == 1 && gy2 == 1 && gy1 == 1 && gz1 == 1)
                            % type = '4 solid faces at y1y2-x1z1';
                            ui = 12*u2 * Y1*Y2 * X1^2 * Z1 / dx^5;
                            vi = 0;
                            wi = 12*w2 * Y1*Y2 * X1 * Z1^2 / dx^5;
                        elseif (gx1 == 1 && gz2 == 1) % (gx1 == 1 && gy2 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '4 solid faces at y1y2-x1z2';
                            ui = 12*u2 * Y1*Y2 * X1^2 * Z2 / dx^5;
                            vi = 0;
                            wi = 12*w1 * Y1*Y2 * X1 * Z2^2 / dx^5;
                        elseif (gx2 == 1 && gz1 == 1) % (gx2 == 1 && gy2 == 1 && gy1 == 1 && gz1 == 1)
                            % type = '4 solid faces at y1y2-x2z1';
                            ui = 12*u1 * Y1*Y2 * X2^2 * Z1 / dx^5;
                            vi = 0;
                            wi = 12*w2 * Y1*Y2 * X2 * Z1^2 / dx^5;
                        else % (gx2 == 1 && gy2 == 1 && gy1 == 1 && gz2 == 1)
                            % type = '4 solid faces at y1y2-x2z2';
                            ui = 12*u1 * Y1*Y2 * X2^2 * Z2 / dx^5;
                            vi = 0;
                            wi = 12*w1 * Y1*Y2 * X2 * Z2^2 / dx^5;
                        end
                        
                    else % gZ == 2
                        %% Two in Z Direction and Two in Different Directions (4 cases)
                        if (gx1 == 1 && gy1 == 1) % (gx1 == 1 && gy1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at z1z2-x1y1';
                            ui = 12*u2 * Z1*Z2 * X1^2 * Y1 / dx^5;
                            vi = 12*v2 * Z1*Z2 * X1 * Y1^2 / dx^5;
                            wi = 0;
                        elseif (gx1 == 1 && gy2 == 1) % (gx1 == 1 && gy2 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at z1z2-x1y2';
                            ui = 12*u2 * Z1*Z2 * X1^2 * Y2 / dx^5;
                            vi = 12*v1 * Z1*Z2 * X1 * Y2^2 / dx^5;
                            wi = 0;
                        elseif (gx2 == 1 && gy1 == 1) % (gx2 == 1 && gy1 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at z1z2-x2y1';
                            ui = 12*u1 * Z1*Z2 * X2^2 * Y1 / dx^5;
                            vi = 12*v2 * Z1*Z2 * X2 * Y1^2 / dx^5;
                            wi = 0;
                        else % (gx2 == 1 && gy2 == 1 && gz2 == 1 && gz1 == 1)
                            % type = '4 solid faces at z1z2-x2y2';
                            ui = 12*u1 * Z1*Z2 * X2^2 * Y2 / dx^5;
                            vi = 12*v1 * Z1*Z2 * X2 * Y2^2 / dx^5;
                            wi = 0;
                        end
                        
                    end
                    
                else
                    %% type = '4, 5 or 6 solid faces';
                    ui = 0;
                    vi = 0;
                    wi = 0;
                end
                
                %% Advection displacement
                advectionX = ui * dt; x = x + advectionX;
                advectionY = vi * dt; y = y + advectionY;
                advectionZ = wi * dt; z = z + advectionZ;
                
                if x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                    % if particle moves out system -> reinject
                    %% ReinjectAdvection
                    % Reinject (during advection)
                    % If a particle exits the system during its advection,
                    % it will be reinject at the inlet face at a flow-weighted random location
                    %
                    % Workspace variables
                    %   x,y,z  : current location of particle
                    %   dx     : dimension of the cell (dx=dy=dz)
                    %   gL     : sample length
                    %   gFlow  : sample main flow direction
                    %   Uinlet : velocity field of the inlet faces (X direction)
                    %   Vinlet : velocity field of the inlet faces (Y direction)
                    %   Winlet : velocity field of the inlet faces (Z direction)
                    %
                    % Output
                    %   x,y,z  : Updated location of particle
                    %
                    % Functions used:
                    %   out_system : Check if particle out system or not
                    %
                    % Constants Used:
                    %   FLOW_X, FLOW_Y, FLOW_Z
                    %

                    if gFlow == FLOW_X
                        %% [x,y,z] = reinject_advection_X (Uinlet,x,y,z,dx,gL);
                        velField = Uinlet;
                        while x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                            random_index = randsample(velField.index, 1, true, velField.list);
                            [~,j,k] = ind2sub (size(velField.plane),random_index);
                            x = (velField.face-0.5) * dx;
                            y = (j-0.5) * dx;
                            z = (k-0.5) * dx;
                        end
                        event.ReinjectAdvectionX = event.ReinjectAdvectionX + 1;
                    elseif gFlow == FLOW_Y
                        %% [x,y,z] = reinject_advection_Y (Vinlet,x,y,z,dx,gL);
                        velField = Vinlet;
                        while x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                            random_index = randsample(velField.index, 1, true, velField.list);
                            [i,~,k] = ind2sub (size(velField.plane),random_index);
                            x = (i-0.5) * dx;
                            y = (velField.face-0.5) * dx;
                            z = (k-0.5) * dx;
                        end
                        event.ReinjectAdvectionY = event.ReinjectAdvectionY + 1;
                    else % gFlow == FLOW_Z
                        %% [x,y,z] = reinject_advection_Z (Winlet,x,y,z,dx,gL);
                        velField = Winlet;
                        while x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                            random_index = randsample(velField.index, 1, true, velField.list);
                            [i,j,~] = ind2sub(size(velField.plane),random_index);
                            x = (i-0.5) * dx;
                            y = (j-0.5) * dx;
                            z = (velField.face-0.5) * dx;
                        end
                        event.ReinjectAdvectionZ = event.ReinjectAdvectionZ + 1;
                    end
                    
                else
                    % advection movement within the system
                    i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);%[i,j,k] = get_ijk (x,y,z,dx);
                    if g(i,j,k) == 1
                        % if particle hits a solid surface -> NO MOVE
                        x = xBA; y = yBA; z = zBA; % restore location
                        advectionX = 0; advectionY = 0; advectionZ = 0; % update advection displacement
                        event.AdvectionHitSolid = event.AdvectionHitSolid + 1;
                    end
                end
            end
            
            %% DiffusionMovement_micro
            % Diffusion Movement (MicroPorosity)
            % Move particle by using a random-walk method: the particle instantaneously
            % jumps over a fixed distance (zeta) in a random direction
            %
            % Workspace variables
            %   g     : geometry
            %   dx    : dimension of the cell (dx=dy=dz)
            %   dt    : time-step
            %   gL    : sample length
            %   zeta  : particle jump distance
            %   x,y,z : current location of particle
            %
            %   inside_solid(p) : inside_solid status of particle (p)
            %   chance_into_solid, inside_limit, micro_zeta : microporosity parameters
            %
            % Output
            %   x,y,z : Update location of particle
            %   diffusionX, diffusionY, diffusionZ : Diffusion displacement
            %
            %   inside_solid(p) : Update inside_solid status of particle (p)
            %
            % Scipts Used:
            %   ReinjectDiffusion :  Update x,y,z
            %
            % Functions Used:
            %   out_system : Check if particle out system or not
            %   get_ijk    : Voxel (i,j,k) of particle location (x,y,z)
            %
            
            % Particle location before diffusion
            xBD = x; yBD = y; zBD = z;
            
            % Randomize theta & phi
            theta = abs(2*rand()*pi);
            phi = abs(rand()*pi);
            
            % Diffustion displacement
            if inside_solid(p) == true
                zetaTmp = micro_zeta;
            else
                zetaTmp = zeta;
            end
            
            diffusionX = zetaTmp * sin(phi) * cos(theta); x = x + diffusionX;
            diffusionY = zetaTmp * sin(phi) * sin(theta); y = y + diffusionY;
            diffusionZ = zetaTmp * cos(phi);              z = z + diffusionZ;
            
            if x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                % if particle jumps out system -> reinject
                %% ReinjectDiffusion
                random_index = ceil(Nparticle*rand());
                x = xpt0(random_index);
                y = ypt0(random_index);
                z = zpt0(random_index);
                if x >= dx || x <= gL
                    event.ReinjectDiffusionX = event.ReinjectDiffusionX + 1;
                elseif y >= dx || y <= gL
                    event.ReinjectDiffusionY = event.ReinjectDiffusionY + 1;
                elseif z >= dx || z <= gL
                    event.ReinjectDiffusionZ = event.ReinjectDiffusionZ + 1;
                end
                
            elseif inside_solid(p)
                % if inside_solid (and in system) -> update inside_solid status accordingly
                i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);%[i,j,k] = get_ijk (x,y,z,dx);
                if g(i,j,k)==0
                    inside_solid(p) = false;
                    event.GetOutSolid = event.GetOutSolid + 1;
                end
                
            else
                %% diffusion movement within the system (as normal)
                i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);%[i,j,k] = get_ijk (x,y,z,dx);
                if g(i,j,k) == 1
                    % if hit a surface -> ROLL a chance to STAY INSIDE
                    if (rand() <= chance_into_solid)
                        inside_solid(p) = true; % stay inside
                        %[x,y,z] = into_solid(i,j,k,xpt1,ypt1,zpt1,...
                        %    x,y,z,diffusionX,diffusionY,diffusionZ,dx,inside_limit);
                        
                        i1=ceil(xBD/dx); j1=ceil(yBD/dx); k1=ceil(zBD/dx);%[i1,j1,k1] = get_ijk (xpt1,ypt1,zpt1,dx);
    
                        if (i ~= i1) %bounce X direction
                            if (diffusionX >= 0)
                                insideX = abs(x - (dx*floor(x/dx)));
                                x = x - insideX + inside_limit;
                            else
                                insideX = abs(x - (dx*ceil(x/dx)));
                                x = x + insideX - inside_limit;
                            end

                        elseif (j ~= j1) %bounce Y direction
                            if (diffusionY >= 0)
                                insideY = abs(y - (dx*floor(y/dx)));
                                y = y - insideY + inside_limit;
                            else
                                insideY = abs(y - (dx*ceil(y/dx)));
                                y = y + insideY - inside_limit;
                            end

                        elseif (k ~= k1) %bounce Z direction
                            if (diffusionZ >= 0)
                                insideZ = abs(z - (dx*floor(z/dx)));
                                z = z - insideZ + inside_limit;
                            else
                                insideZ = abs(z - (dx*ceil(z/dx)));
                                z = z + insideZ - inside_limit;
                            end
                        end
                        diffusionX = x - xBD; diffusionY = y - zBD; diffusionZ = z - zBD; % update diffusion
                        event.GetInSolid = event.GetInSolid + 1;
                    else
                        %% if hit a surface -> bounce movement (as normal)
                        %[x,y,z] = bounce_movement(i,j,k,xpt1,ypt1,zpt1,...
                        %    x,y,z,diffusionX,diffusionY,diffusionZ,dx);
                        
                        i1=ceil(xBD/dx); j1=ceil(yBD/dx); k1=ceil(zBD/dx);%[i1,j1,k1] = get_ijk (xpt1,ypt1,zpt1,dx);

                        if (i ~= i1) %bounce X direction
                            if (diffusionX >= 0)
                                bounceX = - 2*abs(x - (dx*floor(x/dx)));
                            else
                                bounceX =   2*abs(x - (dx*ceil(x/dx)));
                            end
                            x = x + bounceX;

                        elseif (j ~= j1) %bounce Y direction
                            if (diffusionY >= 0)
                                bounceY = - 2*abs(y - (dx*floor(y/dx)));
                            else
                                bounceY =   2*abs(y - (dx*ceil(y/dx)));
                            end
                            y = y + bounceY;

                        elseif (k ~= k1) %bounce Z direction
                            if (diffusionZ >= 0)
                                bounceZ = - 2*abs(z - (dx*floor(z/dx)));
                            else
                                bounceZ =   2*abs(z - (dx*ceil(z/dx)));
                            end
                            z = z + bounceZ;

                        end
                        
                        diffusionX = x - xBD; diffusionY = y - zBD; diffusionZ = z - zBD; % update diffusion
                                               
                        if x<=dx || x>=gL || y<=dx || y>=gL || z<=dx || z>=gL
                            % after bouncing: if particle bounces out system -> reinject again
                            %% ReinjectDiffusion
                            random_index = ceil(Nparticle*rand());
                            x = xpt0(random_index);
                            y = ypt0(random_index);
                            z = zpt0(random_index);
                            if x >= dx || x <= gL
                                event.ReinjectDiffusionX = event.ReinjectDiffusionX + 1;
                            elseif y >= dx || y <= gL
                                event.ReinjectDiffusionY = event.ReinjectDiffusionY + 1;
                            elseif z >= dx || z <= gL
                                event.ReinjectDiffusionZ = event.ReinjectDiffusionZ + 1;
                            end
                            
                        else
                            %% bounce movement within system
                            i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);%[i,j,k] = get_ijk (x,y,z,dx);
                            if g(i,j,k) == 1
                                % if bounce movement into a solid surface -> NO MOVE
                                x = xBD; y = yBD; z = zBD; % restore location
                                diffusionX = 0; diffusionY = 0; diffusionZ = 0;% update diffusion
                            end
                        end
                    end
                end
            end
            
            %% Update particle displacement
            if gFlow == FLOW_X
                displacement2(p) = displacement2(p) + advectionX + diffusionX;
            elseif gFlow == FLOW_Y
                displacement2(p) = displacement2(p) + advectionY + diffusionZ;
            else % gFlow == FLOW_Z
                displacement2(p) = displacement2(p) + advectionZ + diffusionZ;
            end
            
            xpt(p) = x; ypt(p) = y; zpt(p) = z; % copy short variables back to array(index)
        end % for Nparticle
        
        if mod(t,slope) == 0
            % update variance & calculate dispersion coefficient
            counter = counter + 1;
            variance2(counter) = var(displacement2); % Update variance
            dispersion2(counter) = 0.5 * (variance2(counter) - variance2(counter-1)) / (slope*dt);
            insolid(counter) = mean(inside_solid); % mean(inside_solid) at each timestep
            
            % display progress
            progressMsg = sprintf(' -> Sim2: time=%dx%d, insolid=%.3f, RDC=%.2f, Diff=%.2f ',...
                counter-1, slope, insolid(counter), dispersion2(counter)/Dm,...
                (dispersion2(counter)-dispersion1(counter))/Dm);
            fprintf([reverseStr, progressMsg]);
            reverseStr = repmat(sprintf('\b'), 1, length(progressMsg));
        end
    end % for Ntimestep
    
    toc;
    xpt2 = xpt; ypt2 = ypt; zpt2 = zpt; % save xpt,ypt,zpt of simulation2
    
    %% Save results to file
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