%% Prepare Velocity Fields for Reinjection
%
% Workspace variables
%   u,v,w   : components of velocity field in x,y,z direction, respectively
%   gN      : geometry number of voxel
%
% Output
%   Uinlet  : velocity field of the inlet faces (X direction)
%   Vinlet  : velocity field of the inlet faces (Y direction)
%   Winlet  : velocity field of the inlet faces (Z direction)
%   Uoutlet : velocity field of the outlet faces (X direction)
%   Voutlet : velocity field of the outlet faces (Y direction)
%   Woutlet : velocity field of the outlet faces (Z direction)
%

%% Constants
REINJECTION_INLET = 2;
REINJECTION_OUTLET = gN - 1;
REINJECTION_LIMIT = 0.1; % minimum velocity to re-inject

%% U velocity plane at inlet face
Uinlet.face  = REINJECTION_INLET;
Uinlet.plane = u(Uinlet.face,:,:);
Uinlet.list  = reshape(Uinlet.plane,[],1);
Uinlet.limit = mean(Uinlet.plane (Uinlet.plane~=0)) * REINJECTION_LIMIT;
Uinlet.index = find(Uinlet.list > Uinlet.limit);
Uinlet.list  = Uinlet.list(Uinlet.index);

%% U velocity plane at outlet face
Uoutlet.face  = REINJECTION_OUTLET;
Uoutlet.plane = u(Uoutlet.face,:,:);
Uoutlet.list  = reshape(Uoutlet.plane,[],1);
Uoutlet.limit = mean(Uoutlet.plane (Uoutlet.plane~=0)) * REINJECTION_LIMIT;
Uoutlet.index = find(Uoutlet.list > Uoutlet.limit);
Uoutlet.list  = Uoutlet.list(Uoutlet.index);

%% V velocity plane at inlet face
Vinlet.face  = REINJECTION_INLET;
Vinlet.plane = v(:,Vinlet.face,:);
Vinlet.list  = reshape(Vinlet.plane,[],1);
Vinlet.limit = mean(Vinlet.plane (Vinlet.plane~=0)) * REINJECTION_LIMIT;
Vinlet.index = find(Vinlet.list > Vinlet.limit);
Vinlet.list  = Vinlet.list(Vinlet.index);

%% V velocity plane at outlet face
Voutlet.face  = REINJECTION_OUTLET;
Voutlet.plane = v(:,Voutlet.face,:);
Voutlet.list  = reshape(Voutlet.plane,[],1);
Voutlet.limit = mean(Voutlet.plane (Voutlet.plane~=0)) * REINJECTION_LIMIT;
Voutlet.index = find(Voutlet.list > Voutlet.limit);
Voutlet.list  = Voutlet.list(Voutlet.index);

%% W velocity plane at inlet face
Winlet.face  = REINJECTION_INLET;
Winlet.plane = w(:,:,Winlet.face);
Winlet.list  = reshape(Winlet.plane,[],1);
Winlet.limit = mean(Winlet.plane (Winlet.plane~=0)) * REINJECTION_LIMIT;
Winlet.index = find(Winlet.list > Winlet.limit);
Winlet.list  = Winlet.list(Winlet.index);

%% W velocity plane at outlet face
Woutlet.face  = REINJECTION_OUTLET;
Woutlet.plane = w(:,:,Woutlet.face);
Woutlet.list  = reshape(Woutlet.plane,[],1);
Woutlet.limit = mean(Woutlet.plane (Woutlet.plane~=0)) * REINJECTION_LIMIT;
Woutlet.index = find(Woutlet.list > Woutlet.limit);
Woutlet.list  = Woutlet.list(Woutlet.index);

