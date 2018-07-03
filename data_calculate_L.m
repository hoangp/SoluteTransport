%% Calculate L
close all; clearvars -except u_org v_org w_org g;
FLOW_X = 1; FLOW_Y = 2; FLOW_Z = 3; % Flow direction

sample = '../data/Berea'; % sample MAT file name
gN = 400;
dx = 5.3;

% geopack: L=150.21
% Berea: L=129.84
% bentheimer: L=154.47

%% load sample
if ~(exist('g','var') && exist('u_org','var') && exist('v_org','var') && exist('w_org','var'))
    fprintf(sprintf('Load sample %s...',sample));   
    load(sample);
    fprintf(sprintf('Done!\n')); 
end

%% calculation
surface = 0;
edge = 0;

for i=2:gN-1 
    for j=2:gN-1
        for k=2:gN-1           
            if g(i,j,k)==1                               
                % geometry arround the solid
                gx1 = g(i-1,j,k); gx2 = g(i+1,j,k);
                gy1 = g(i,j-1,k); gy2 = g(i,j+1,k);
                gz1 = g(i,j,k-1); gz2 = g(i,j,k+1);
                
                % number of solid boundary arround the solid
                num_solid = gx1 + gx2 + gy1 + gy2 + gz1 + gz2;
                
                % number of pore arround the solid
                num_pore = 6 - num_solid;
                
                if num_pore > 0
                    surface = surface + num_pore;
                    edge = edge + 1;
                end
            end
        end
    end
end

S = surface * (dx^2);
V = ((gN-2)^3) * (dx^3);
L = pi*V/S;

fprintf('edges=%d surfaces=%d L=%.2f\n',edge,surface,L);
             
                
        
            