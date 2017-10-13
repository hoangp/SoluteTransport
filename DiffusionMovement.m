%% Diffusion Movement
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
diffusionX = zeta * sin(phi) * cos(theta); x = x + diffusionX;
diffusionY = zeta * sin(phi) * sin(theta); y = y + diffusionY;
diffusionZ = zeta * cos(phi);              z = z + diffusionZ;

out = out_system_location (x,y,z,dx,gL); % need for ReinjectionDiffusion

if out ~= 0
    % if particle jumps out system -> reinject
    ReinjectDiffusion
    
else
    %% diffusion movement within the system
    [ip_move,jp_move,kp_move] = get_ijk (x,y,z,dx);
    
    if g(ip_move,jp_move,kp_move) == 1
        % if hit a surface -> bounce movement
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

%% Nested Fuctions

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
