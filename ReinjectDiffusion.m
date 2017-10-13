%% Reinject (during diffusion)
% If a particle exits the system during its diffusion,
% it will be reinject by an area-weighted random location on inlet or outlet face
%
% Workspace variables
%   RV    : Reinjection Variables structure
%   x,y,z : current location of particle
%   dx    : dimension of the cell (dx=dy=dz)
%   gL    : sample length
%   out   : out location of the particle
%
% Output
%   x,y,z : Updated location of particle
%
% Functions used:
%   out_system
%
% Constants Used:
%   OUT_X1, OUT_X2, OUT_Y1, OUT_Y2, OUT_Z1, OUT_Z2
%

if out == OUT_X2
    [x,y,z] = reinject_diffusion_X (RV.Uinlet,x,y,z,dx,gL);
elseif out == OUT_X1
    [x,y,z] = reinject_diffusion_X (RV.Uoutlet,x,y,z,dx,gL);
elseif out == OUT_Y2
    [x,y,z] = reinject_diffusion_Y (RV.Vinlet,x,y,z,dx,gL);
elseif out == OUT_Y1
    [x,y,z] = reinject_diffusion_Y (RV.Voutlet,x,y,z,dx,gL);
elseif out == OUT_Z2
    [x,y,z] = reinject_diffusion_Z (RV.Winlet,x,y,z,dx,gL);
elseif out == OUT_Z1
    [x,y,z] = reinject_diffusion_Z (RV.Woutlet,x,y,z,dx,gL);
end

%% Nested Functions

function [x,y,z] = reinject_diffusion_X (velocityFace,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = velocityFace.index(ceil(length(velocityFace.index)*rand()));
        [~,j,k]=ind2sub(size(velocityFace.plane),random_index);
        x = (velocityFace.face-0.5) * dx;
        y = (j-0.5) * dx;
        z = (k-0.5) * dx;
    end
end

function [x,y,z] = reinject_diffusion_Y (velocityFace,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = velocityFace.index(ceil(length(velocityFace.index)*rand()));
        [i,~,k]=ind2sub(size(velocityFace.plane),random_index);
        x = (i-0.5) * dx;
        y = (velocityFace.face-0.5) * dx;
        z = (k-0.5) * dx;
    end
end

function [x,y,z] = reinject_diffusion_Z (velocityFace,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = velocityFace.index(ceil(length(velocityFace.index)*rand()));
        [i,j,~]=ind2sub(size(velocityFace.plane),random_index);
        x = (i-0.5) * dx;
        y = (j-0.5) * dx;
        z = (velocityFace.face-0.5) * dx;
    end
end