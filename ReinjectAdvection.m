%% Reinject (during advection)
% If a particle exits the system during its advection,
% it will be reinject at the inlet face at a flow-weighted random location
%
% Workspace variables
%   x,y,z  : current location of particle
%   dx     : dimension of the cell (dx=dy=dz)
%   gL     : sample length
%   gFlow  : sample main flow direction
%   Uinlet : velocity field of the inlet faces (X direction)
%   Vinlet : velocity field of the inlet faces (Y direction)
%   Winlet : velocity field of the inlet faces (Z direction)
%
% Output
%   x,y,z  : Updated location of particle
%
% Functions used:
%   out_system : Check if particle out system or not
%
% Constants Used:
%   FLOW_X, FLOW_Y, FLOW_Z
%

if gFlow == FLOW_X
    [x,y,z] = reinject_advection_X (Uinlet,x,y,z,dx,gL);
    event.ReinjectAdvectionX = event.ReinjectAdvectionX + 1;
elseif gFlow == FLOW_Y
    [x,y,z] = reinject_advection_Y (Vinlet,x,y,z,dx,gL);
    event.ReinjectAdvectionY = event.ReinjectAdvectionY + 1;
else % gFlow == FLOW_Z
    [x,y,z] = reinject_advection_Z (Winlet,x,y,z,dx,gL);
    event.ReinjectAdvectionZ = event.ReinjectAdvectionZ + 1;
end

%% Nested Functions

function [x,y,z] = reinject_advection_X (velField,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = randsample(velField.index, 1, true, velField.list);
        [~,j,k] = ind2sub (size(velField.plane),random_index);
        x = (velField.face-0.5) * dx;
        y = (j-0.5) * dx;
        z = (k-0.5) * dx;
    end
end

function [x,y,z] = reinject_advection_Y (velField,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = randsample(velField.index, 1, true, velField.list);
        [i,~,k] = ind2sub (size(velField.plane),random_index);
        x = (i-0.5) * dx;
        y = (velField.face-0.5) * dx;
        z = (k-0.5) * dx;
    end
end

function [x,y,z] = reinject_advection_Z (velField,x,y,z,dx,gL)
    while out_system (x,y,z,dx,gL)
        random_index = randsample(velField.index, 1, true, velField.list);
        [i,j,~] = ind2sub(size(velField.plane),random_index);
        x = (i-0.5) * dx;
        y = (j-0.5) * dx;
        z = (velField.face-0.5) * dx;
    end
end