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

%% Diffustion displacement
diffusionX = zeta * sin(phi) * cos(theta); x = x + diffusionX;
diffusionY = zeta * sin(phi) * sin(theta); y = y + diffusionY;
diffusionZ = zeta * cos(phi);              z = z + diffusionZ;

out_location = out_system_location (x,y,z,dx,gL); % for ReinjectionDiffusion

if out_location ~= 0
    % if particle jumps out system -> reinject
    ReinjectDiffusion2
    
else
    %% diffusion movement within the system
    [i,j,k] = get_ijk (x,y,z,dx);    
    if g(i,j,k) == 1
        % if hit a surface -> bounce movement     
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
        event.DiffusionHitSolid = event.DiffusionHitSolid + 1;
    end    
end