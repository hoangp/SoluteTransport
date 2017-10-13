%% Reinject (during advection)
% If a particle exits the system during its advection,
% it will be reinject at the inlet face at a flow-weighted random location
%
% Workspace variables
%   RV    : Reinjection Variables structure
%   x,y,z : current location of particle
%   dx    : dimension of the cell (dx=dy=dz)
%   gL    : sample length
%   gFlow : sample main flow direction
%
% Output
%   x,y,z : Updated location of particle
%
% Functions used:
%   out_system
%
% Constants Used:
%   FLOW_X, FLOW_Y, FLOW_Z
%

if gFlow == FLOW_X
    [x,y,z] = reinject_advection_X (RV.Uinlet,x,y,z,dx,gL);
elseif gFlow == FLOW_Y
    [x,y,z] = reinject_advection_Y (RV.Vinlet,x,y,z,dx,gL);
else % gFlow == FLOW_Z
    [x,y,z] = reinject_advection_Z (RV.Winlet,x,y,z,dx,gL);
end

%% Nested Functions

function [x,y,z] = reinject_advection_X (velocityFace,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = randsample(velocityFace.index, 1, true, velocityFace.list);
        [~,j,k]=ind2sub(size(velocityFace.plane),random_index);
        x = (velocityFace.face-0.5) * dx;
        y = (j-0.5) * dx;
        z = (k-0.5) * dx;
    end
end

function [x,y,z] = reinject_advection_Y (velocityFace,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = randsample(velocityFace.index, 1, true, velocityFace.list);
        [i,~,k]=ind2sub(size(velocityFace.plane),random_index);
        x = (i-0.5) * dx;
        y = (velocityFace.face-0.5) * dx;
        z = (k-0.5) * dx;
    end
end

function [x,y,z] = reinject_advection_Z (velocityFace,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = randsample(velocityFace.index, 1, true, velocityFace.list);
        [i,j,~]=ind2sub(size(velocityFace.plane),random_index);
        x = (i-0.5) * dx;
        y = (j-0.5) * dx;
        z = (velocityFace.face-0.5) * dx;
    end
end