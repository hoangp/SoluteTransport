%% Prepare Velocity Fields for Reinjection
%
% Workspace variables
%   u,v,w : components of velocity field in x,y,z direction, respectively
%   gN    : geometry number of voxel
%
% Output
%   RV    : Reinjection Variables Structure 
%

%% Constants
REINJECTION_INLET = 2;
REINJECTION_OUTLET = gN - 1;
REINJECTION_LIMIT = 0.1; % minimum velocity to re-inject

%% RV Structure
RV = [];

%% U velocity plane at inlet face
RV.Uinlet.face  = REINJECTION_INLET;
RV.Uinlet.plane = u(RV.Uinlet.face,:,:);
RV.Uinlet.list  = reshape(RV.Uinlet.plane,[],1);
RV.Uinlet.limit = mean(RV.Uinlet.plane (RV.Uinlet.plane~=0)) * REINJECTION_LIMIT;
RV.Uinlet.index = find(RV.Uinlet.list > RV.Uinlet.limit);
RV.Uinlet.list  = RV.Uinlet.list(RV.Uinlet.index);

%% U velocity plane at outlet face
RV.Uoutlet.face  = REINJECTION_OUTLET;
RV.Uoutlet.plane = u(RV.Uoutlet.face,:,:);
RV.Uoutlet.list  = reshape(RV.Uoutlet.plane,[],1);
RV.Uoutlet.limit = mean(RV.Uoutlet.plane (RV.Uoutlet.plane~=0)) * REINJECTION_LIMIT;
RV.Uoutlet.index = find(RV.Uoutlet.list > RV.Uoutlet.limit);
RV.Uoutlet.list  = RV.Uoutlet.list(RV.Uoutlet.index);

%% V velocity plane at inlet face
RV.Vinlet.face  = REINJECTION_INLET;
RV.Vinlet.plane = v(:,RV.Vinlet.face,:);
RV.Vinlet.list  = reshape(RV.Vinlet.plane,[],1);
RV.Vinlet.limit = mean(RV.Vinlet.plane (RV.Vinlet.plane~=0)) * REINJECTION_LIMIT;
RV.Vinlet.index = find(RV.Vinlet.list > RV.Vinlet.limit);
RV.Vinlet.list  = RV.Vinlet.list(RV.Vinlet.index);

%% V velocity plane at outlet face
RV.Voutlet.face  = REINJECTION_OUTLET;
RV.Voutlet.plane = v(:,RV.Voutlet.face,:);
RV.Voutlet.list  = reshape(RV.Voutlet.plane,[],1);
RV.Voutlet.limit = mean(RV.Voutlet.plane (RV.Voutlet.plane~=0)) * REINJECTION_LIMIT;
RV.Voutlet.index = find(RV.Voutlet.list > RV.Voutlet.limit);
RV.Voutlet.list  = RV.Voutlet.list(RV.Voutlet.index);

%% W velocity plane at inlet face
RV.Winlet.face  = REINJECTION_INLET;
RV.Winlet.plane = w(:,:,RV.Winlet.face);
RV.Winlet.list  = reshape(RV.Winlet.plane,[],1);
RV.Winlet.limit = mean(RV.Winlet.plane (RV.Winlet.plane~=0)) * REINJECTION_LIMIT;
RV.Winlet.index = find(RV.Winlet.list > RV.Winlet.limit);
RV.Winlet.list  = RV.Winlet.list(RV.Winlet.index);

%% W velocity plane at outlet face
RV.Woutlet.face  = REINJECTION_OUTLET;
RV.Woutlet.plane = w(:,:,RV.Woutlet.face);
RV.Woutlet.list  = reshape(RV.Woutlet.plane,[],1);
RV.Woutlet.limit = mean(RV.Woutlet.plane (RV.Woutlet.plane~=0)) * REINJECTION_LIMIT;
RV.Woutlet.index = find(RV.Woutlet.list > RV.Woutlet.limit);
RV.Woutlet.list  = RV.Woutlet.list(RV.Woutlet.index);

