%% Avection Movement
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
xpt0 = x; ypt0 = y; zpt0 = z;

% Streamline-tracing algorithm 
StreamlineTracing

% Advection displacement
advectionX = ui * dt; x = x + advectionX;
advectionY = vi * dt; y = y + advectionY;
advectionZ = wi * dt; z = z + advectionZ;

if out_system (x,y,z,dx,gL)
    % if particle moves out system -> reinject
    ReinjectAdvection
    
else
    % advection movement within the system
    [i,j,k] = get_ijk (x,y,z,dx);    
    if g(i,j,k) == 1
        % if particle hits a solid surface -> NO MOVE
        x = xpt0; y = ypt0; z = zpt0; % restore location
        advectionX = 0; advectionY = 0; advectionZ = 0; % update advection displacement
    end
end
