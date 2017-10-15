%% Constants
FLOW_X = 1; FLOW_Y = 2; FLOW_Z = 3; % Flow direction
OUT_X1 = 1; OUT_X2 = 2; OUT_Y1 = 3; OUT_Y2 = 4; OUT_Z1 = 5; OUT_Z2 = 6; % out_system_location

%% Sample Properties
% 
% Workspace Variables:
%   sample: sample file name
% 
% Output:
%   g,u,v,w
%   dx,gN,L,gL,gFlow
%

if ~(exist('g','var') && exist('u_org','var') && exist('v_org','var') && exist('w_org','var'))
    fprintf(sprintf('Load sample: %s...',sample));
    load(sample);
    fprintf(sprintf('Done!\n'));
end

dt = 0.0001; % time-step (seconds)

switch lower(sample)
    case 'data/geometry2sample'
        dx    = 5.345;   % dimension of the cell (dx=dy=dz) in micro meter
        L     = 131.13;  % characteristic length (micro meter)
        gN    = 300;     % number of voxels     
        gFlow = FLOW_X;  % sample main flow direction
        
    case 'data/sand500sample'
        dx    = 2.32;    % dimension of the cell (dx=dy=dz) in micro meter
        L     = 51.63;   % characteristic length (micro meter)
        gN    = 500;     % number of voxels         
        gFlow = FLOW_Z;  % sample main flow direction     
end

gL = dx * gN - dx; % sample length (micro meter)


