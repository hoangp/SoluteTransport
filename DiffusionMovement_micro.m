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
xpt1 = x; ypt1 = y; zpt1 = z;

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
    ReinjectDiffusion
    
elseif inside_solid(p)              
    % if inside_solid (and in system) -> update inside_solid status accordingly
    [i,j,k] = get_ijk (x,y,z,dx);
    if g(i,j,k)==0
        inside_solid(p) = false;
    end
    
else
    %% diffusion movement within the system (as normal)
    [i,j,k] = get_ijk (x,y,z,dx);   
    if g(i,j,k) == 1
        % if hit a surface -> ROLL a chance to STAY INSIDE
        if roll_a_chance(chance_into_solid) == true
            inside_solid(p) = true; % stay inside
            [x,y,z] = into_solid(i,j,k,xpt1,ypt1,zpt1,...
                x,y,z,diffusionX,diffusionY,diffusionZ,dx,inside_limit);
        else     
            %% if hit a surface -> bounce movement (as normal)
            [x,y,z] = bounce_movement(i,j,k,xpt1,ypt1,zpt1,...
                x,y,z,diffusionX,diffusionY,diffusionZ,dx);
            
            diffusionX = x - xpt1; diffusionY = y - zpt1; diffusionZ = z - zpt1; % update diffusion
            
            out_location = out_system_location (x,y,z,dx,gL); % check out system again
            
            if out_location ~= 0
                % after bouncing: if particle bounces out system -> reinject again
                ReinjectDiffusion
                
            else
                %% bounce movement within system
                [i,j,k] = get_ijk (x,y,z,dx);                
                if g(i,j,k) == 1
                    % if bounce movement into a solid surface -> NO MOVE
                    x = xpt1; y = ypt1; z = zpt1; % restore location
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

function [x,y,z] = into_solid(i,j,k,xpt1,ypt1,zpt1,...
    x,y,z,movementX,movementY,movementZ,dx,inside_limit)

    [i1,j1,k1] = get_ijk (xpt1,ypt1,zpt1,dx);
    
    if (i ~= i1) %bounce X direction
        if (movementX >= 0)
            insideX = abs(x - (dx*floor(x/dx)));
            x = x - insideX + inside_limit;
        else
            insideX = abs(x - (dx*ceil(x/dx)));
            x = x + insideX - inside_limit;
        end
        
    elseif (j ~= j1) %bounce Y direction
        if (movementY >= 0)
            insideY = abs(y - (dx*floor(y/dx)));
            y = y - insideY + inside_limit;
        else
            insideY = abs(y - (dx*ceil(y/dx)));
            y = y + insideY - inside_limit;
        end
        
    elseif (k ~= k1) %bounce Z direction
        if (movementZ >= 0)
            insideZ = abs(z - (dx*floor(z/dx)));
            z = z - insideZ + inside_limit;
        else
            insideZ = abs(z - (dx*ceil(z/dx)));
            z = z + insideZ - inside_limit;
        end
    end
end

function out_location = out_system_location (x,y,z,left,right)
    if (x <= left)
        out_location = 1; % OUT_X1;
    elseif (x >= right)
        out_location = 2; % OUT_X2;
    elseif (y <= left)
        out_location = 3; % OUT_Y1;
    elseif (y >= right)
        out_location = 4; % OUT_Y2;
    elseif (z <= left)
        out_location = 5; % OUT_Z1;
    elseif (z >= right)
        out_location = 6; % OUT_Z2;
    else
        out_location = 0;
    end
end

function [x,y,z] = bounce_movement(i,j,k,...
    xpt1,ypt1,zpt1,x,y,z,movementX,movementY,movementZ,dx)

    [i1,j1,k1] = get_ijk (xpt1,ypt1,zpt1,dx);
    
    if (i ~= i1) %bounce X direction
        if (movementX >= 0)
            bounceX = - 2*abs(x - (dx*floor(x/dx)));
        else
            bounceX =   2*abs(x - (dx*ceil(x/dx)));
        end
        x = x + bounceX;
        
    elseif (j ~= j1) %bounce Y direction
        if (movementY >= 0)
            bounceY = - 2*abs(y - (dx*floor(y/dx)));
        else
            bounceY =   2*abs(y - (dx*ceil(y/dx)));
        end
        y = y + bounceY;
        
    elseif (k ~= k1) %bounce Z direction
        if (movementZ >= 0)
            bounceZ = - 2*abs(z - (dx*floor(z/dx)));
        else
            bounceZ =   2*abs(z - (dx*ceil(z/dx)));
        end
        z = z + bounceZ;
        
    end                
end