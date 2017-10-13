%% Simulation
%
% Workspace Variables
%   Ntimestep, Nparticle
%   xpt,ypt,zpt
%
% Output
%   xpt,ypt,zpt  : Final particles locations
%   displacement : Final particles displacement
%   variance     : Final variance
%   sim : Structure stores other simulation properties:
%       : dt,dx,Dm,Npe,zeta,gFlow,uavg,vavg,wavg,sample,Nparticle,Ntimestep

%% Microporosity variables
insolid = zeros (Ntimestep,1); % percent of insolid at a time
inside_solid = false (Nparticle,1);
[ip,jp,kp] = get_ijk(xpt,ypt,zpt,dx);
for p=1:Nparticle
    if g(ip(p),jp(p),kp(p)) == 1
        inside_solid(p) = true;
    end
end

%% Simulation properties
sim.dt = dt; sim.dx=dx; sim.Dm=Dm; sim.Npe=Npe; sim.zeta=zeta;
sim.gFlow=gFlow; sim.uavg=uavg; sim.uavg=vavg; sim.uavg=wavg;
sim.sample=sample; sim.Nparticle=Nparticle; sim.Ntimestep=Ntimestep;
sim.MicroporosityFlag=MicroporosityFlag;
sim.chance_into_solid = chance_into_solid; sim.micro_zeta=micro_zeta;
sim.inside_limit=inside_limit; sim.insolid = insolid; 

displacement = zeros(Nparticle,1); % particles displacement in the main flow
variance = zeros(Ntimestep,1); % variance of displacement

%% Simulation starts
reverseStr = ''; % for diplay progress
for t = 1:Ntimestep
    for p = 1:Nparticle
        x = xpt(p); y = ypt(p); z = zpt(p);
           
        if inside_solid(p) == false
            AdvectionMovement;
        end
        
        if MicroporosityFlag == true
            DiffusionMovement_micro;
        else
            DiffusionMovement;
        end
             
        % Update particle displacement
        if gFlow == FLOW_X
            displacement(p) = displacement(p) + advectionX + diffusionX;
        elseif gFlow == FLOW_Y
            displacement(p) = displacement(p) + advectionX + diffusionX;
        elseif gFlow == FLOW_Z
            displacement(p) = displacement(p) + advectionZ + diffusionZ;
        end
        
        xpt(p) = x; ypt(p) = y; zpt(p) = z;
        
    end
    
    variance(t) = var(displacement); % Update variance
    if MicroporosityFlag == true
        insolid(t) = mean(inside_solid); % percent of insolid at a time
    end
    
    % Display progress
    if mod(t,ceil(Ntimestep/100)) == 0
        ProgressMsg = sprintf('Percentage done: %.0f. ', 100 * t / Ntimestep);
        fprintf([reverseStr, ProgressMsg]);
        reverseStr = repmat(sprintf('\b'), 1, length(ProgressMsg));
    end
    
end
