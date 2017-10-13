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
%   inside_solid(p)
%   inside_limit
%   micro_zeta
%
% Output
%   x,y,z : Update location of particle
%   diffusionX, diffusionY, diffusionZ
%
% Scipts Used:
%   ReinjectDiffusion
%       Input : RV,x,y,z,dx,gL,out
%       Output: update x,y,z
%
% Functions Used:
%   out_system, get_ijk

% Particle location before diffusion
xpt1 = x; ypt1 = y; zpt1 = z;

% Randomize theta & phi
theta = abs(2*rand()*pi);
phi = abs(rand()*pi);

%% Diffustion movement
if inside_solid(p) == true
    p_zeta = micro_zeta; %flag_visualize(i) = true;
else
    p_zeta = zeta;
end

diffusionX = p_zeta * sin(phi) * cos(theta); x = x + diffusionX;
diffusionY = p_zeta * sin(phi) * sin(theta); y = y + diffusionY;
diffusionZ = p_zeta * cos(phi);              z = z + diffusionZ;

out = out_system_location (x,y,z,dx,gL); % need for ReinjectionDiffusion

if out ~= 0
    % if particle jumps out system -> reinject
    ReinjectDiffusion
    
elseif inside_solid(p)              
    % if inside_solid -> update inside_solid status accordingly
    if g(ip_move,jp_move,kp_move)==0
        inside_solid(p) = false;
    end
    
else
    %% diffusion movement within the system (as normal)
    [ip_move,jp_move,kp_move] = get_ijk (x,y,z,dx);
    
    if g(ip_move,jp_move,kp_move) == 1
        % if hit a surface -> ROLL a chance to STAY INSIDE
        if roll_a_chance(chance_into_solid) == true
            inside_solid(p) = true; % stay inside
            [x,y,z] = into_solid(ip_move,jp_move,kp_move,xpt1,ypt1,zpt1,...
                x,y,z,diffusionX,diffusionY,diffusionZ,dx,inside_limit);
        else     
            %% if hit a surface -> bounce movement (as normal)
            [x,y,z] = bounce_movement(ip_move,jp_move,kp_move,xpt1,ypt1,zpt1,...
                x,y,z,diffusionX,diffusionY,diffusionZ,dx);
            
            diffusionX = x - xpt1; diffusionY = y - zpt1; diffusionZ = z - zpt1; % update diffusion
            
            out = out_system_location (x,y,z,dx,gL); % check out system again
            
            if out ~= 0
                % after bouncing: if particle bounces out system -> reinject again
                ReinjectDiffusion
                
            else
                %% bounce movement within system
                [ip_move,jp_move,kp_move] = get_ijk (x,y,z,dx);
                
                if g(ip_move,jp_move,kp_move) == 1
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

function [x,y,z] = into_solid(ip_move,jp_move,kp_move,xpt1,ypt1,zpt1,...
    x,y,z,movementX,movementY,movementZ,dx,inside_limit)

    [ip_previous,jp_previous,kp_previous] = get_ijk (xpt1,ypt1,zpt1,dx);
    
    if (ip_move ~= ip_previous) %bounce X direction
        if (movementX >= 0)
            insideX = abs(x - (dx*floor(x/dx)));
            x = x - insideX + inside_limit;
        else
            insideX = abs(x - (dx*ceil(x/dx)));
            x = x + insideX - inside_limit;
        end
        
    elseif (jp_move ~= jp_previous) %bounce Y direction
        if (movementY >= 0)
            insideY = abs(y - (dx*floor(y/dx)));
            y = y - insideY + inside_limit;
        else
            insideY = abs(y - (dx*ceil(y/dx)));
            y = y + insideY - inside_limit;
        end
        
    elseif (kp_move ~= kp_previous) %bounce Z direction
        if (movementZ >= 0)
            insideZ = abs(z - (dx*floor(z/dx)));
            z = z - insideZ + inside_limit;
        else
            insideZ = abs(z - (dx*ceil(z/dx)));
            z = z + insideZ - inside_limit;
        end
    end
end

function out = out_system_location (x,y,z,left,right)
    if (x <= left)
        out = 1; % OUT_X1;
    elseif (x >= right)
        out = 2; % OUT_X2;
    elseif (y <= left)
        out = 3; % OUT_Y1;
    elseif (y >= right)
        out = 4; % OUT_Y2;
    elseif (z <= left)
        out = 5; % OUT_Z1;
    elseif (z >= right)
        out = 6; % OUT_Z2;
    else
        out = 0;
    end
end

function [x1,y1,z1] = bounce_movement(ip_move,jp_move,kp_move,...
    xpt1,ypt1,zpt1,x,y,z,movementX,movementY,movementZ,dx)

    [ip_previous,jp_previous,kp_previous] = get_ijk (xpt1,ypt1,zpt1,dx);

    bounceX = 0;
    bounceY = 0;
    bounceZ = 0;
    
    if (ip_move ~= ip_previous) %bounce X direction
        if (movementX >= 0)
            bounceX = - 2*abs(x - (dx*floor(x/dx)));
        else
            bounceX =   2*abs(x - (dx*ceil(x/dx)));
        end
        
    elseif (jp_move ~= jp_previous) %bounce Y direction
        if (movementY >= 0)
            bounceY = - 2*abs(y - (dx*floor(y/dx)));
        else
            bounceY =   2*abs(y - (dx*ceil(y/dx)));
        end
        
    elseif (kp_move ~= kp_previous) %bounce Z direction
        if (movementZ >= 0)
            bounceZ = - 2*abs(z - (dx*floor(z/dx)));
        else
            bounceZ =   2*abs(z - (dx*ceil(z/dx)));
        end
        
    end
    
    x1 = x + bounceX;
    y1 = y + bounceY;
    z1 = z + bounceZ;
end
