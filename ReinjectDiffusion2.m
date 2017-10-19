%% Reinject (during diffusion) - version 2
% If a particle exits the system during its diffusion,
% it will be reinject to a random original injected location
%
% Workspace variables
%   x,y,z          : current location of particle
%   dx             : dimension of the cell (dx=dy=dz)
%   gL             : sample length
%   xpt0,ypt0,zpt0 : original injected location
%
% Output
%   x,y,z           : Updated location of particle
%   event structure : ReinjectDiffusionX, ReinjectDiffusionY, ReinjectDiffusionZ
%

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