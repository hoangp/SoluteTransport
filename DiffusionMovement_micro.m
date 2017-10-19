%% Diffusion Movement (MicroPorosity)
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

%% Diffustion displacement
if inside_solid(p) == true
    zetaTmp = micro_zeta; 
else
    zetaTmp = zeta;
end

diffusionX = zetaTmp * sin(phi) * cos(theta); x = x + diffusionX;
diffusionY = zetaTmp * sin(phi) * sin(theta); y = y + diffusionY;
diffusionZ = zetaTmp * cos(phi);              z = z + diffusionZ;

out_location = out_system_location (x,y,z,dx,gL); % for ReinjectionDiffusion

if out_location ~= 0
    % if particle jumps out system -> reinject
    ReinjectDiffusion2
    
elseif inside_solid(p)              
    % if inside_solid (and in system) -> update inside_solid status accordingly
    [i,j,k] = get_ijk (x,y,z,dx);
    if g(i,j,k)==0
        inside_solid(p) = false;
        event.GetOutSolid = event.GetOutSolid + 1;
    end
    
else
    %% diffusion movement within the system (as normal)
    [i,j,k] = get_ijk (x,y,z,dx);   
    if g(i,j,k) == 1
        % if hit a surface -> ROLL a chance to STAY INSIDE
        if roll_a_chance(chance_into_solid) == true
            inside_solid(p) = true; % stay inside
            
            [x,y,z] = into_solid_movement(i,j,k,xBD,yBD,zBD,...
                x,y,z,diffusionX,diffusionY,diffusionZ,dx,inside_limit);
            
            diffusionX = x - xBD; diffusionY = y - zBD; diffusionZ = z - zBD; % update diffusion
            
            event.GetInSolid = event.GetInSolid + 1;
        else     
            %% if hit a surface -> bounce movement (as normal)
            [x,y,z] = bounce_movement(i,j,k,xBD,yBD,zBD,...
                x,y,z,diffusionX,diffusionY,diffusionZ,dx);      
            
            diffusionX = x - xBD; diffusionY = y - zBD; diffusionZ = z - zBD; % update diffusion            
            
            out_location = out_system_location (x,y,z,dx,gL); % check out system again        
            if out_location ~= 0
                % after bouncing: if particle bounces out system -> reinject again
                ReinjectDiffusion2
                
            else
                %% bounce movement within system
                [i,j,k] = get_ijk (x,y,z,dx);                
                if g(i,j,k) == 1
                    % if bounce movement into a solid surface -> NO MOVE
                    x = xBD; y = yBD; z = zBD; % restore location
                    diffusionX = 0; diffusionY = 0; diffusionZ = 0;% update diffusion
                end               
            end           
        end        
    end    
end

%% Nested Fuctions

function roll = roll_a_chance(chance_into_solid)
    roll = (rand() <= chance_into_solid);
end

function [x,y,z] = into_solid_movement(i,j,k,xBD,yBD,zBD,...
    x,y,z,diffusionX,diffusionY,diffusionZ,dx,inside_limit)

    [iBD,jBD,kBD] = get_ijk (xBD,yBD,zBD,dx);
    
    if (i ~= iBD) %bounce X direction
        if (diffusionX >= 0)
            insideX = abs(x - (dx*floor(x/dx)));
            x = x - insideX + inside_limit;
        else
            insideX = abs(x - (dx*ceil(x/dx)));
            x = x + insideX - inside_limit;
        end
        
    elseif (j ~= jBD) %bounce Y direction
        if (diffusionY >= 0)
            insideY = abs(y - (dx*floor(y/dx)));
            y = y - insideY + inside_limit;
        else
            insideY = abs(y - (dx*ceil(y/dx)));
            y = y + insideY - inside_limit;
        end
        
    elseif (k ~= kBD) %bounce Z direction
        if (diffusionZ >= 0)
            insideZ = abs(z - (dx*floor(z/dx)));
            z = z - insideZ + inside_limit;
        else
            insideZ = abs(z - (dx*ceil(z/dx)));
            z = z + insideZ - inside_limit;
        end
    end
end