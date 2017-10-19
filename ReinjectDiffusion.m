%% Reinject (during diffusion)
% If a particle exits the system during its diffusion,
% it will be reinject by an area-weighted random location on inlet or outlet face
%
% Workspace variables
%   x,y,z        : current location of particle
%   dx           : dimension of the cell (dx=dy=dz)
%   gL           : sample length
%   out_location : out location of the particle
%   Uinlet       : velocity field of the inlet faces (X direction)
%   Vinlet       : velocity field of the inlet faces (Y direction)
%   Winlet       : velocity field of the inlet faces (Z direction)
%   Uoutlet      : velocity field of the outlet faces (X direction)
%   Voutlet      : velocity field of the outlet faces (Y direction)
%   Woutlet      : velocity field of the outlet faces (Z direction)
%
% Output
%   x,y,z        : Updated location of particle
%
% Functions used:
%   out_system   : Check if particle out system or not
%
% Constants Used:
%   OUT_X1, OUT_X2, OUT_Y1, OUT_Y2, OUT_Z1, OUT_Z2
%

if out_location == OUT_X2
    [x,y,z] = reinject_diffusion_X (Uinlet,x,y,z,dx,gL);
    event.ReinjectDiffusionX = event.ReinjectDiffusionX + 1;
elseif out_location == OUT_X1
    [x,y,z] = reinject_diffusion_X (Uoutlet,x,y,z,dx,gL);
    event.ReinjectDiffusionX = event.ReinjectDiffusionX + 1;
elseif out_location == OUT_Y2
    [x,y,z] = reinject_diffusion_Y (Vinlet,x,y,z,dx,gL);
    event.ReinjectDiffusionY = event.ReinjectDiffusionY + 1;
elseif out_location == OUT_Y1
    [x,y,z] = reinject_diffusion_Y (Voutlet,x,y,z,dx,gL);
    event.ReinjectDiffusionY = event.ReinjectDiffusionY + 1;
elseif out_location == OUT_Z2
    [x,y,z] = reinject_diffusion_Z (Winlet,x,y,z,dx,gL);
    event.ReinjectDiffusionZ = event.ReinjectDiffusionZ + 1;
elseif out_location == OUT_Z1
    [x,y,z] = reinject_diffusion_Z (Woutlet,x,y,z,dx,gL);
    event.ReinjectDiffusionZ = event.ReinjectDiffusionZ + 1;
end

%% Nested Functions

function [x,y,z] = reinject_diffusion_X (velField,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = velField.index(ceil(length(velField.index)*rand()));
        [~,j,k] = ind2sub (size(velField.plane),random_index);
        x = (velField.face-0.5) * dx;
        y = (j-0.5) * dx;
        z = (k-0.5) * dx;
    end
end

function [x,y,z] = reinject_diffusion_Y (velField,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = velField.index(ceil(length(velField.index)*rand()));
        [i,~,k] = ind2sub (size(velField.plane),random_index);
        x = (i-0.5) * dx;
        y = (velField.face-0.5) * dx;
        z = (k-0.5) * dx;
    end
end

function [x,y,z] = reinject_diffusion_Z (velField,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = velField.index(ceil(length(velField.index)*rand()));
        [i,j,~] = ind2sub(size(velField.plane),random_index);
        x = (i-0.5) * dx;
        y = (j-0.5) * dx;
        z = (velField.face-0.5) * dx;
    end
end