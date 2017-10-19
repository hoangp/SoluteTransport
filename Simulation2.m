%% Simulation with Microporosity (Ntimestep)
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

%% Simulation variables
displacement2 = zeros(Nparticle,1); % particles displacement in the main flow
variance2 = zeros(length(varianceTime),1); % variance of displacement
dispersion2 = zeros(length(varianceTime),1); % variance of displacement
insolid = zeros (length(varianceTime),1); % mean(inside_solid)

% Microporosity variables
inside_solid = false (Nparticle,1); % current status of each particle

%% Simulation starts
reverseStr = ''; % for diplay progress
counter = 1;
for t = 1:Ntimestep
    for p = 1:Nparticle
        x = xpt(p); y = ypt(p); z = zpt(p); % short variables instead of array(index)
           
        if inside_solid(p) == false
            AdvectionMovement
        end
        
        DiffusionMovement_micro
             
        % Update particle displacement
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
