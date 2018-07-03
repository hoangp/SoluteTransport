%% Normalize velocity fields and save data to file
geoName = 'bentheimer'; eval(['g = ',geoName,';']); g = g ~= 0;
uName = 'xVel2'; eval(['u_org = ',uName,';']);
vName = 'yVel2'; eval(['v_org = ',vName,';']);
wName = 'zVel2'; eval(['w_org = ',wName,';']);
dt = 0.0001;
dx = 5.3;

%% maximize u,v,w (flow in z direcrion)
max_distance = max(max(max(w_org))) * dt; % assume flow unit is mico meter/second, flow in z direction

velmax_coef = dx / max_distance; % max distance is size of voxel
u_org = u_org * velmax_coef;
v_org = v_org * velmax_coef;
w_org = w_org * velmax_coef;

%% save data
data_folder = '../data/';
saveFileName = sprintf('%s.mat',geoName);
save([data_folder,saveFileName],'g','u_org','v_org','w_org')