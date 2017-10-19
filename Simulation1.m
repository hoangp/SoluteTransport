%% Simulation until we reach an asymptotic regime
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
%   MIN_TIMESTEP   : minimum timestep for a simulation
%   MIN_SLOPE      : minimum timestep slope to calculate dispersion

maxnumSlope = ceil(MAX_TIMESTEP/MIN_SLOPE); % Maximum number of slope, assume min(slope) = 100

%% Simulation variables
displacement1 = zeros(Nparticle,1); 
variance1 = zeros(maxnumSlope,1); 
varianceTime = zeros(maxnumSlope,1);
dispersion1 = zeros(maxnumSlope,1);
dispersionChange = zeros(maxnumSlope,1);

%% Simulation starts
reverseStr = ''; % for diplay progress
counter = 2; dispersionChange(counter-1) = 1; limitBeforeMinTime = 0;
while (dispersionChange(counter-1) > D_LIMIT && varianceTime(counter-1) < MAX_TIMESTEP) || (varianceTime(counter-1) < MIN_TIMESTEP)
    for t = 1:slope
        for p = 1:Nparticle
            x = xpt(p); y = ypt(p); z = zpt(p); % short variables instead of array(index)

            AdvectionMovement          
            DiffusionMovement

            % Update particle displacement
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
    
    % update variance & calculate dispersion coefficient
    variance1(counter) = var(displacement1); % Update variance
    varianceTime(counter) = varianceTime(counter-1) + slope;
    dispersion1(counter) = 0.5 * (variance1(counter) - variance1(counter-1)) / (slope*dt);
    dispersionChange(counter) = abs(dispersion1(counter)-dispersion1(counter-1))/dispersion1(counter);
    if dispersionChange(counter-1) <= D_LIMIT && limitBeforeMinTime == 0
        limitBeforeMinTime = counter;
    end        
    
    % display progress
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