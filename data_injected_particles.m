%% Inject Partiles
% Inject particle with conditions:
%   VELOCITY_LIMIT : minimum velocity to inject
%   BOUNDARY_SPACE  : advoid boundary
%
VELOCITY_LIMIT = 0.1; % minimum velocity to inject
BOUNDARY_SPACE = 1; % advoid boundary
FLOW_X = 1; FLOW_Y = 2; FLOW_Z = 3; % Flow direction

%% testblock2 (gN = 50)
%sample = '../data/testblock2'; % sample MAT file name
%iFrom = 5; iTo = 45; jFrom = iFrom; jTo = iTo; kFrom = iFrom; kTo = iTo;
%Nparticle=2900; space=3; ch='A'; % A
%Nparticle=8800; space=2; ch='B'; % B 
%Nparticle=67000; space=1; ch='C'; % C 
%gFlow = FLOW_Z;  % sample main flow direction  
%dx = 10;

%% geopack (gN = 150)
%sample = '../data/geopack'; % sample MAT file name
%iFrom = 5; iTo = 145; jFrom = iFrom; jTo = iTo; kFrom = iFrom; kTo = iTo;
%Nparticle=1050; space=10; ch='A'; % A
%Nparticle=7990; space=5; ch='B'; % B 
%Nparticle=36000; space=3; ch='C'; % C 
%gFlow = FLOW_Z;  % sample main flow direction  
%dx = 10;

%% Berea (gN = 400)
%sample = '../data/Berea'; % sample MAT file name
%iFrom = 5; iTo = 395; jFrom = iFrom; jTo = iTo; kFrom = iFrom; kTo = iTo;
%Nparticle=1150; space=20; ch='A'; % A
%Nparticle=8800; space=10; ch='B'; % B 
%Nparticle=40000; space=6; ch='C'; % C 
%gFlow = FLOW_Z;  % sample main flow direction  
%dx = 5.3;

%% bentheimer (gN = 500)
sample = '../data/bentheimer'; % sample MAT file name
iFrom = 5; iTo = 495; jFrom = iFrom; jTo = iTo; kFrom = iFrom; kTo = iTo;
%Nparticle=910; space=30; ch='A'; % A
%Nparticle=9800; space=13; ch='B'; % B 
Nparticle=42000; space=8; ch='C'; % C 
gFlow = FLOW_Z;  % sample main flow direction  
dx = 5.3;

%%
if ~(exist('g','var') && exist('u_org','var') && exist('v_org','var') && exist('w_org','var'))
    fprintf(sprintf('Load sample %s...',sample));   
    load(sample);
    fprintf(sprintf('Done!\n')); 
end

fprintf(sprintf('Injecting %d particles...',Nparticle));

uavg = mean(u_org(u_org~=0));
vavg = mean(v_org(v_org~=0));
wavg = mean(w_org(w_org~=0));

if gFlow == FLOW_X
    minVelocity = uavg * VELOCITY_LIMIT;
    averageVelocity = uavg;
elseif gFlow == FLOW_Y
    minVelocity = vavg * VELOCITY_LIMIT;
    averageVelocity = vavg;
else
    minVelocity = wavg * VELOCITY_LIMIT;
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


particleSet = [sample,'ParticleSet',ch]; % particle Set MAT file name to be saved


question = input(sprintf('Save to file: %s ? Please choose: Yes/No? ',particleSet),'s');
if lower(question) == 'y'
    
    save(sprintf('%s.mat',particleSet),'xpt','ypt','zpt','Nparticle',...
        'space','iFrom','iTo','jFrom','jTo','kFrom','kTo',...
        'smallVelCounter','boundaryCounter','minVelocity','averageVelocity','ratioVelocity');
    
end