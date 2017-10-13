%% Inject Partiles
% Inject particle with conditions:
%   ECTION_LIMIT : minimum velocity to inject
%   BOUNDARY_SPACE  : advoid boundary
%
close all; clearvars -except u_org v_org w_org g;
sample = 'Data/geometry2Sample'; % sample MAT file name
particleSet = 'Data/geo2ParticleSetD'; % particle Set MAT file name to be saved

INJECTION_LIMIT = 0.1; % minimum velocity to inject
BOUNDARY_SPACE = 1; % advoid boundary

%% sand500sample
% iFrom = 10; iTo = 490; jFrom = 10; jTo = 490; kFrom = 10; kTo = 490;
% Nparticle=1940; space=20; % A
% Nparticle=4400; space=15; % B 
% Nparticle=14580; space=10; % C 
% Nparticle=114489; space=5; % D

%% geometry2sample
iFrom = 10; iTo = 290; jFrom = 10; jTo = 290; kFrom = 10; kTo = 290;
%Nparticle=1155; space=20; % A
%Nparticle=2740; space=15; % B 
%Nparticle=8770; space=10; % C 
Nparticle=66970; space=5; % D

%%
SampleConstants; % load sample properties & constants
fprintf(sprintf('Injecting %d particles...',Nparticle));

uavg = mean(u_org(u_org~=0));
vavg = mean(v_org(v_org~=0));
wavg = mean(w_org(w_org~=0));

if gFlow == FLOW_X
    minVelocity = uavg * INJECTION_LIMIT;
    averageVelocity = uavg;
elseif gFlow == FLOW_Y
    minVelocity = vavg * INJECTION_LIMIT;
    averageVelocity = vavg;
else
    minVelocity = wavg * INJECTION_LIMIT;
    averageVelocity = wavg;
end

ipt = zeros (Nparticle,1);
jpt = zeros (Nparticle,1);
kpt = zeros (Nparticle,1);

velocityArray = zeros(Nparticle,1);

i = iFrom; j = jFrom; k = kFrom;

injectionCounter = 0;
smallVelCounter = 0;
boundaryCounter = 0;

while injectionCounter < Nparticle
    if (g(i,j,k) == 0)
        if any(g(i-BOUNDARY_SPACE:i+BOUNDARY_SPACE,j-BOUNDARY_SPACE:j+BOUNDARY_SPACE,k-BOUNDARY_SPACE:k+BOUNDARY_SPACE)==1)
            boundaryCounter = boundaryCounter + 1;
        else
            
            if gFlow == FLOW_X
                vel = u_org(i,j,k);
            elseif gFlow == FLOW_Y
                vel = v_org(i,j,k);
            else
                vel = w_org(i,j,k);
            end
            
            if vel > minVelocity
                injectionCounter = injectionCounter + 1;
                ipt(injectionCounter) = i;
                jpt(injectionCounter) = j;
                kpt(injectionCounter) = k;
                velocityArray(injectionCounter)=vel;
            else
                smallVelCounter = smallVelCounter + 1;
            end
        end
    end
    
    if k >= kTo
        k = kFrom;
        j = j + space;
        if j > jTo
            j = jTo;
        end
        
    elseif j >= jTo
        j = jFrom;
        k = kFrom;
        i = i + space;
        if i > iTo
            i = iTo;
            fprintf(sprintf('Full injection at Counter=%d \n',injectionCounter));
            return;
        end
    else
        k = k + space;
        if k > kTo
            k = kTo;
        end
    end
end

xpt = (ipt-0.5)*dx; ypt = (jpt-0.5)*dx; zpt = (kpt-0.5)*dx;

ratioVelocity = mean(velocityArray) / averageVelocity;

fprintf(sprintf('Done! (SmallVel=%d Boundary=%d Ratio=%.2f)\n',...
    smallVelCounter,boundaryCounter,ratioVelocity));

question = input(sprintf('Save to file: %s ? Please choose: Yes/No? ',particleSet),'s');
if lower(question) == 'y'
    save(sprintf('%s.mat',particleSet),'xpt','ypt','zpt','Nparticle',...
        'space','iFrom','iTo','jFrom','jTo','kFrom','kTo',...
        'smallVelCounter','boundaryCounter','minVelocity','averageVelocity','ratioVelocity');
    
end