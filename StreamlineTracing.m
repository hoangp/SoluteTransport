%% Streamline-tracing algorithm
% Semi-analytic formulation for the streamlines in voxels contain a solid boundary
%
% Workspace variables
%   g        : geometry
%   dx       : dimension of the cell (dx=dy=dz)
%   x,y,z    : Location of particle
%
% Output
%   ui,vi,wi : interpolated velocity between the normal velocities on opposite faces of gridblock
%
% Restriction:
%   Only consider geometry from (2) --> (N-1)

%% Get particle boundary conditions

% index location of the particle
ip = ceil(x/dx);
jp = ceil(y/dx);
kp = ceil(z/dx);

% location of the edges
x1 = dx*(ip-1); x2 = dx*ip;
y1 = dx*(jp-1); y2 = dx*jp;
z1 = dx*(kp-1); z2 = dx*kp;

% velocities on the faces normal to each direction
u1 = u(ip-1,jp,kp); u2 = u(ip,jp,kp);
v1 = v(ip,jp-1,kp); v2 = v(ip,jp,kp);
w1 = w(ip,jp,kp-1); w2 = w(ip,jp,kp);

% short variables names
X1 = x-x1; X2 = x2-x;
Y1 = y-y1; Y2 = y2-y;
Z1 = z-z1; Z2 = z2-z;
U = u2-u1; V = v2-v1; W = w2-w1;

% geometry of the edge
gx1 = g(ip-1,jp,kp); gx2 = g(ip+1,jp,kp);
gy1 = g(ip,jp-1,kp); gy2 = g(ip,jp+1,kp);
gz1 = g(ip,jp,kp-1); gz2 = g(ip,jp,kp+1);

% short variables names
gX = gx1+gx2; gY = gy1+gy2; gZ = gz1+gz2;
gXY = gX+gY; gXZ = gX+gZ; gYZ = gY+gZ;

% number of solid boundary arround the particle
num_solid = gx1 + gx2 + gy1 + gy2 + gz1 + gz2;

%% Calculate interpolated velocities: ui,vi,wi

if num_solid == 0
    %% No solid boundaries
    ui = U * X1 / dx + u1;
    vi = V * Y1 / dx + v1;
    wi = W * Z1 / dx + w1;
    
elseif num_solid == 1
    %% One of the Neighbouring Voxels Is Solid (6 cases)
    if gx2 == 1
        % type = '1 solid face at x2'; SAMPLE FORMULATION
        ui = u1 * X2^2 / dx^2;
        vi = 2*v1 * X2 / dx + 2*V * X2*Y1 / dx^2;
        wi = 2*w1 * X2 / dx + 2*W * X2*Z1 / dx^2;
    elseif gx1 == 1
        % type = '1 solid face at x1';
        ui = u2 * X1^2 / dx^2;
        vi = 2*v1 * X1 / dx + 2*V * X1*Y1 / dx^2;
        wi = 2*w1 * X1 / dx + 2*W * X1*Z1 / dx^2;
    elseif gy2 == 1
        % type = '1 solid face at y2';
        ui = 2*u1 * Y2 / dx + 2*U * Y2*X1 / dx^2;
        vi = v1 * Y2^2 / dx^2;
        wi = 2*w1 * Y2 / dx + 2*W * Y2*Z1 / dx^2;
    elseif gy1 == 1
        % type = '1 solid face at y1';
        ui = 2*u1 * Y1 / dx + 2*U * Y1*X1 / dx^2;
        vi = v2 * Y1^2 / dx^2;
        wi = 2*w1 * Y1 / dx + 2*W * Y1*Z1 / dx^2;
    elseif gz2 == 1
        % type = '1 solid face at z2';
        ui = 2*u1 * Z2 / dx + 2*U * Z2*X1 / dx^2;
        vi = 2*v1 * Z2 / dx + 2*V * Z2*Y1 / dx^2;
        wi = w1 * Z2^2 / dx^2;
    elseif gz1 == 1
        % type = '1 solid face at z1';
        ui = 2*u1 * Z1 / dx + 2*U * Z1*X1 / dx^2;
        vi = 2*v1 * Z1 / dx + 2*V * Z1*Y1 / dx^2;
        wi = w2 * Z1^2 / dx^2;
    end
    
elseif num_solid == 2
    %% Two of the Neighbouring Voxels Are Solid (15 cases)
    
    if (gX == 2 || gY == 2 || gZ == 2)
        %% There Are Two Opposing Solid Faces (3 cases)
        if gX==2 % (gx2 == 1 && gx1 == 1)
            % type = '2 opposing solid faces x direction'; % SAMPLE FORMULATION
            ui = 0;
            vi = 6*v1 * X2*X1 / dx^2 + 6*V * X2*X1 * Y1 / dx^3;
            wi = 6*w1 * X2*X1 / dx^2 + 6*W * X2*X1 * Z1 / dx^3;
        elseif gY==2 % (gy2 == 1 && gy1 == 1)
            % type = '2 opposing solid faces y direction';
            ui = 6*u1 * Y2*Y1 / dx^2 + 6*U * Y2*Y1 * X1 / dx^3;
            vi = 0;
            wi = 6*w1 * Y2*Y1 / dx^2 + 6*W * Y2*Y1 * Z1 / dx^3;
        else % gz==2 % (gz2 == 1 && gz1 == 1)
            % type = '2 opposing solid faces z direction';
            ui = 6*u1 * Z2*Z1 / dx^2 + 6*U * Z2*Z1 * X1 / dx^3;
            vi = 6*v1 * Z2*Z1 / dx^2 + 6*V * Z2*Z1 * Y1 / dx^3;
            wi = 0;
        end
        
    elseif gXY == 2
        %% Two Adjacent Voxels (XY) Are Solid (4 cases), PAPER BUG
        if (gx1 == 1 && gy1 == 1)
            % type = '2 solid faces at x1y1'; % SAMPLE FORMULATION
            ui = 2*u2 * X1^2 * Y1 / dx^3;
            vi = 2*v2 * X1 * Y1^2 / dx^3;
            wi = 4*w1 * X1*Y1 / dx^2 + 4*W * X1*Y1 * Z1 / dx^3;
        elseif (gx2 == 1 && gy1 == 1)
            % type = '2 solid faces at x2y1';
            ui = 2*u1 * X2^2 * Y1 / dx^3;
            vi = 2*v2 * X2 * Y1^2 / dx^3;
            wi = 4*w1 * X2*Y1 / dx^2 + 4*W * X2*Y1 * Z1 / dx^3;
        elseif (gx1 == 1 && gy2 == 1)
            % type = '2 solid faces at x1y2';
            ui = 2*u2 * X1^2 * Y2 / dx^3;
            vi = 2*v1 * X1 * Y2^2 / dx^3;
            wi = 4*w1 * X1*Y2 / dx^2 + 4*W * X1*Y2 * Z1 / dx^3;
        else % (gx2 == 1 && gy2 == 1)
            % type = '2 solid faces at x2y2';
            ui = 2*u1 * X2^2 * Y2 / dx^3;
            vi = 2*v1 * X2 * Y2^2 / dx^3;
            wi = 4*w1 * X2*Y2 / dx^2 + 4*W * X2*Y2 * Z1 / dx^3;
        end
        
    elseif gXZ == 2
        %% Two Adjacent Voxels (XZ) Are Solid (4 cases)
        if (gx1 == 1 && gz1 == 1)
            % type = '2 solid faces at x1z1';
            ui = 2*u2 * X1^2 * Z1 / dx^3;
            vi = 4*v1 * X1*Z1 / dx^2 + 4*V * X1*Z1 * Y1 / dx^3;
            wi = 2*w2 * X1 * Z1^2 / dx^3;
        elseif (gx2 == 1 && gz1 == 1)
            % type = '2 solid faces at x2z1';
            ui = 2*u1 * X2^2 * Z1 / dx^3;
            vi = 4*v1 * X2*Z1 / dx^2 + 4*V * X2*Z1 * Y1 / dx^3;
            wi = 2*w2 * X2 * Z1^2 / dx^3;
        elseif (gx1 == 1 && gz2 == 1)
            % type = '2 solid faces at x1z2';
            ui = 2*u2 * X1^2 * Z2 / dx^3;
            vi = 4*v1 * X1*Z2 / dx^2 + 4*V * X1*Z2 * Y1 / dx^3;
            wi = 2*w1 * X1 * Z2^2 / dx^3;
        else % (gx2 == 1 && gz2 == 1)
            % type = '2 solid faces at x2z2';
            ui = 2*u1 * X2^2 * Z2 / dx^3;
            vi = 4*v1 * X2*Z2 / dx^2 + 4*V * X2*Z2 * Y1 / dx^3;
            wi = 2*w1 * X2 * Z2^2 / dx^3;
        end
        
    else % gYZ == 2
        %% Two Adjacent Voxels (YZ) Are Solid (4 cases)
        if (gy1 == 1 && gz1 == 1)
            % type = '2 solid faces at y1z1';
            ui = 4*u1 * Y1*Z1 / dx^2 + 4*U * Y1*Z1 * X1 / dx^3 ;
            vi = 2*v2 * Y1^2 * Z1 / dx^3;
            wi = 2*w2 * Y1 *  Z1^2 / dx^3;
        elseif (gy1 == 1 && gz2 == 1)
            % type = '2 solid faces at y1z2';
            ui = 4*u1 * Y1*Z2 / dx^2 + 4*U * Y1*Z2 * X1 / dx^3 ;
            vi = 2*v2 * Y1^2 * Z2 / dx^3;
            wi = 2*w1 * Y1 *  Z2^2 / dx^3;
        elseif (gy2 == 1 && gz1 == 1)
            % type = '2 solid faces at y2z1';
            ui = 4*u1 * Y2*Z1 / dx^2 + 4*U * Y2*Z1 * X1 / dx^3 ;
            vi = 2*v1 * Y2^2 * Z1 / dx^3;
            wi = 2*w2 * Y2 *  Z1^2 / dx^3;
        else % (gy2 == 1 && gz2 == 1)
            % type = '2 solid faces at y2z2';
            ui = 4*u1 * Y2*Z2 / dx^2 + 4*U * Y2*Z2 * X1 / dx^3 ;
            vi = 2*v1 * Y2^2 * Z2 / dx^3;
            wi = 2*w1 * Y2 *  Z2^2 / dx^3;
        end
        
    end
    
elseif num_solid == 3
    %% Three of the Neighbouring Voxels Are Solid (20 cases)
    
    if gX == 2
        %% Two of Voxels Are in the Same X Direction (4 cases)
        if gz1 == 1 % (gx2 == 1 && gx1 == 1 && gz1 == 1)
            % type = '3 solid faces at x1x2-z1'; SAMPLE FORMULATION
            ui =  0;
            vi = 12*v1 * X1*X2*Z1 / dx^3 + 12*V * X1*X2*Z1 * Y1 / dx^4;
            wi =  6*w2 * X1*X2*Z1^2 / dx^4;
        elseif gz2 == 1 % (gx2 == 1 && gx1 == 1 && gz2 == 1)
            % type = '3 solid faces at x1x2-z2'; %% MY BUG FOUND HERE
            ui =  0;
            vi = 12*v1 * X1*X2*Z2 / dx^3 + 12*V * X1*X2*Z2 * Y1 / dx^4;
            wi =  6*w1 * X1*X2*Z2^2 / dx^4;
        elseif gy1 == 1 % (gx2 == 1 && gx1 == 1 && gy1 == 1)
            % type = '3 solid faces at x1x2-y1';
            ui =  0;
            vi =  6*v2 * X1*X2*Y1^2 / dx^4;
            wi = 12*w1 * X1*X2*Y1 / dx^3 + 12*W * X1*X2*Y1 * Z1 / dx^4;
        else % (gx2 == 1 && gx1 == 1 && gy2 == 1)
            % type = '3 solid faces at x1x2-y2'; %% MY BUG FOUND HERE
            ui =  0;
            vi =  6*v1 * X1*X2*Y2^2 / dx^4;
            wi = 12*w1 * X1*X2*Y2 / dx^3 + 12*W * X1*X2*Y2 * Z1 / dx^4;
        end
        
    elseif gY == 2
        %% Two of Voxels Are in the Same Y Direction (4 cases)
        if gx1 == 1 % (gx1 == 1 && gy2 == 1 && gy1 == 1)
            % type = '3 solid faces at y1y2-x1'; %% MY BUG FOUND HERE
            ui =  6*u2 * Y1*Y2*X1^2 / dx^4;
            vi =  0;
            wi = 12*w1 * Y1*Y2*X1 / dx^3 + 12*W * Y1*Y2*X1 * Z1 / dx^4;
        elseif gx2 == 1 % (gx2 == 1 && gy2 == 1 && gy1 == 1)
            % type = '3 solid faces at y1y2-x2';
            ui =  6*u1 * Y1*Y2*X2^2 / dx^4;
            vi =  0;
            wi = 12*w1 * Y1*Y2*X2 / dx^3 + 12*W * Y1*Y2*X2 * Z1 / dx^4;
        elseif gz1 == 1 % (gy1 == 1 && gy2 == 1 && gz1 == 1)
            % type = '3 solid faces at y1y2-z1'; %BUG FOUND HERE
            ui = 12*u1 * Y1*Y2*Z1 / dx^3 + 12*U * Y1*Y2*Z1 * X1 / dx^4;
            vi =  0;
            wi =  6*w2 * Y1*Y2*Z1^2 / dx^4;
        else % (gy2 == 1 && gy1 == 1 && gz2 == 1)
            % type = '3 solid faces at y1y2-z2';
            ui = 12*u1 * Y1*Y2*Z2 / dx^3 + 12*U * Y1*Y2*Z2 * X1 / dx^4;
            vi =  0;
            wi =  6*w1 * Y1*Y2*Z2^2 / dx^4;
        end
        
    elseif gZ == 2
        %% Two of Voxels Are in the Same Z Direction (4 cases)
        if gx1 == 1 % (gx1 == 1 && gz2 == 1 && gz1 == 1)
            % type = '3 solid faces at z1z2-x1';
            ui =  6*u2 * Z1*Z2*X1^2 / dx^4;
            vi = 12*v1 * Z1*Z2*X1 / dx^3 + 12*V * Z1*Z2*X1 * Y1 / dx^4;
            wi =  0;
        elseif gx2 == 1 % (gx2 == 1 && gz2 == 1 && gz1 == 1)
            % type = '3 solid faces at z1z2-x2';
            ui =  6*u1 * Z1*Z2*X2^2 / dx^4;
            vi = 12*v1 * Z1*Z2*X2 / dx^3 + 12*V * Z1*Z2*X2 * Y1 / dx^4;
            wi =  0;
        elseif gy1 == 1 % (gy1 == 1 && gz2 == 1 && gz1 == 1)
            % type = '3 solid faces at z1z2-y1';
            ui = 12*u1 * Z1*Z2*Y1 / dx^3 + 12*U * Z1*Z2*Y1 * X1 / dx^4;
            vi =  6*v2 * Z1*Z2*Y1^2 / dx^4;
            wi =  0;
        else % (gy2 == 1 && gz2 == 1 && gz1 == 1)
            % type = '3 solid faces at z1z2-y2';
            ui = 12*u1 * Z1*Z2*Y2 / dx^3 + 12*U * Z1*Z2*Y2 * X1 / dx^4;
            vi =  6*v1 * Z1*Z2*Y2^2 / dx^4;
            wi =  0;
        end
        
    else
        %% Three block in diferent directions (8 cases)
        if (gx2 == 1 && gy2 == 1 && gz1 == 1)
            % type = '3 solid faces at x2y2z1'; SAMPLE FORMULATION
            ui = 4*u1 * X2^2 * Y2 * Z1 / dx^4;
            vi = 4*v1 * X2 * Y2^2 * Z1 / dx^4;
            wi = 4*w2 * X2 * Y2 * Z1^2 / dx^4;
        elseif (gx2 == 1 && gy2 == 1 && gz2 == 1)
            % type = '3 solid faces at x2y2z2';
            ui = 4*u1 * X2^2 * Y2 * Z2 / dx^4;
            vi = 4*v1 * X2 * Y2^2 * Z2 / dx^4;
            wi = 4*w1 * X2 * Y2 * Z2^2 / dx^4;
        elseif (gx2 == 1 && gy1 == 1 && gz1 == 1)
            % type = '3 solid faces at x2y1z1';
            ui = 4*u1 * X2^2 * Y1 * Z1 / dx^4;
            vi = 4*v2 * X2 * Y1^2 * Z1 / dx^4;
            wi = 4*w2 * X2 * Y1 * Z1^2 / dx^4;
        elseif (gx2 == 1 && gy1 == 1 && gz2 == 1)
            % type = '3 solid faces at x2y1z2';
            ui = 4*u1 * X2^2 * Y1 * Z2 / dx^4;
            vi = 4*v2 * X2 * Y1^2 * Z2 / dx^4;
            wi = 4*w1 * X2 * Y1 * Z2^2 / dx^4;
        elseif (gx1 == 1 && gy1 == 1 && gz1 == 1)
            % type = '3 solid faces at x1y1z1';
            ui = 4*u2 * X1^2 * Y1 * Z1 / dx^4;
            vi = 4*v2 * X1 * Y1^2 * Z1 / dx^4;
            wi = 4*w2 * X1 * Y1 * Z1^2 / dx^4;
        elseif (gx1 == 1 && gy1 == 1 && gz2 == 1)
            % type = '3 solid faces at x1y1z2';
            ui = 4*u2 * X1^2 * Y1 * Z2 / dx^4;
            vi = 4*v2 * X1 * Y1^2 * Z2 / dx^4;
            wi = 4*w1 * X1 * Y1 * Z2^2 / dx^4;
        elseif (gx1 == 1 && gy2 == 1 && gz1 == 1)
            % type = '3 solid faces at x1y2z1';
            ui = 4*u2 * X1^2 * Y2 * Z1 / dx^4;
            vi = 4*v1 * X1 * Y2^2 * Z1 / dx^4;
            wi = 4*w2 * X1 * Y2 * Z1^2 / dx^4;
        elseif (gx1 == 1 && gy2 == 1 && gz2 == 1)
            % type = '3 solid faces at x1y2z2';
            ui = 4*u2 * X1^2 * Y2 * Z2 / dx^4;
            vi = 4*v1 * X1 * Y2^2 * Z2 / dx^4;
            wi = 4*w1 * X1 * Y2 * Z2^2 / dx^4;
        end
        
    end
    
elseif num_solid == 4
    %% Four of the Neighbouring Voxels Are Solid (15 cases)
    if (gX + gY == 4) || (gX + gZ == 4) || (gY + gZ == 4)
        %% They block 2 directions (3 cases)  
        if (gX + gY == 4) % (gx2 == 1 && gx1 == 1 && gy2 == 1 && gy1 == 1)
            % type = '4 solid faces at x1x2-y1y2'; % SAMPLE FORMULATION
            ui = 0;
            vi = 0;
            wi = 36*w1 * X1*X2*Y1*Y2 / dx^4;
        elseif (gX + gZ == 4) % (gx2 == 1 && gx1 == 1 && gz2 == 1 && gz1 == 1)
            % type = '4 solid faces at x1x2-z1z2';
            ui = 0;
            vi = 36*v1 * X1*X2*Z1*Z2 / dx^4;
            wi = 0;
        else % (gY + gZ == 4) % (gy2 == 1 && gy1 == 1 && gz2 == 1 && gz1 == 1)
            % type = '4 solid faces at y1y2-z1z2';
            ui = 36*u1 * Y1*Y2*Z1*Z2 / dx^4;
            vi = 0;
            wi = 0;
        end
        
    elseif gX == 2
        %% Two in X Direction and Two in Different Directions (4 cases)       
        if (gy2 == 1 && gz1 == 1) % (gx2 == 1 && gx1 == 1 && gy2 == 1 && gz1 == 1)
            % type = '4 solid faces at x1x2-y2z1'; SAMPLE FORMULATION
            ui = 0;
            vi = 12*v1 * X1*X2 * Y2^2 * Z1 / dx^5;
            wi = 12*w2 * X1*X2 * Y2 * Z1^2 / dx^5;
        elseif (gy2 == 1 && gz2 == 1) % (gx2 == 1 && gx1 == 1 && gy2 == 1 && gz2 == 1)
            % type = '4 solid faces at x1x2-y2z2';
            ui = 0;
            vi = 12*v1 * X1*X2 * Y2^2 * Z2 / dx^5;
            wi = 12*w1 * X1*X2 * Y2 * Z2^2 / dx^5;
        elseif (gy1 == 1 && gz1 == 1) % (gx2 == 1 && gx1 == 1 && gy1 == 1 && gz1 == 1)
            % type = '4 solid faces at x1x2-y1z1';
            ui = 0;
            vi = 12*v2 * X1*X2 * Y1^2 * Z1 / dx^5;
            wi = 12*w2 * X1*X2 * Y1 * Z1^2 / dx^5;
        else % (gx2 == 1 && gx1 == 1 && gy1 == 1 && gz2 == 1)
            % type = '4 solid faces at x1x2-y1z2';
            ui = 0;
            vi = 12*v2 * X1*X2 * Y1^2 * Z2 / dx^5;
            wi = 12*w1 * X1*X2 * Y1 * Z2^2 / dx^5;
        end
        
    elseif gY == 2
        %% Two in Y Direction and Two in Different Directions (4 cases)   
        if (gx1 == 1 && gz1 == 1) % (gx1 == 1 && gy2 == 1 && gy1 == 1 && gz1 == 1)
            % type = '4 solid faces at y1y2-x1z1';
            ui = 12*u2 * Y1*Y2 * X1^2 * Z1 / dx^5;
            vi = 0;
            wi = 12*w2 * Y1*Y2 * X1 * Z1^2 / dx^5;
        elseif (gx1 == 1 && gz2 == 1) % (gx1 == 1 && gy2 == 1 && gy1 == 1 && gz2 == 1)
            % type = '4 solid faces at y1y2-x1z2';
            ui = 12*u2 * Y1*Y2 * X1^2 * Z2 / dx^5;
            vi = 0;
            wi = 12*w1 * Y1*Y2 * X1 * Z2^2 / dx^5;
        elseif (gx2 == 1 && gz1 == 1) % (gx2 == 1 && gy2 == 1 && gy1 == 1 && gz1 == 1)
            % type = '4 solid faces at y1y2-x2z1';
            ui = 12*u1 * Y1*Y2 * X2^2 * Z1 / dx^5;
            vi = 0;
            wi = 12*w2 * Y1*Y2 * X2 * Z1^2 / dx^5;        
        else % (gx2 == 1 && gy2 == 1 && gy1 == 1 && gz2 == 1)
            % type = '4 solid faces at y1y2-x2z2';
            ui = 12*u1 * Y1*Y2 * X2^2 * Z2 / dx^5;
            vi = 0;
            wi = 12*w1 * Y1*Y2 * X2 * Z2^2 / dx^5; 
        end
        
    else % gZ == 2
        %% Two in Z Direction and Two in Different Directions (4 cases) 
        if (gx1 == 1 && gy1 == 1) % (gx1 == 1 && gy1 == 1 && gz2 == 1 && gz1 == 1)
            % type = '4 solid faces at z1z2-x1y1';
            ui = 12*u2 * Z1*Z2 * X1^2 * Y1 / dx^5;
            vi = 12*v2 * Z1*Z2 * X1 * Y1^2 / dx^5;
            wi = 0;
        elseif (gx1 == 1 && gy2 == 1) % (gx1 == 1 && gy2 == 1 && gz2 == 1 && gz1 == 1)
            % type = '4 solid faces at z1z2-x1y2';
            ui = 12*u2 * Z1*Z2 * X1^2 * Y2 / dx^5;
            vi = 12*v1 * Z1*Z2 * X1 * Y2^2 / dx^5;
            wi = 0;
        elseif (gx2 == 1 && gy1 == 1) % (gx2 == 1 && gy1 == 1 && gz2 == 1 && gz1 == 1)
            % type = '4 solid faces at z1z2-x2y1';
            ui = 12*u1 * Z1*Z2 * X2^2 * Y1 / dx^5;
            vi = 12*v2 * Z1*Z2 * X2 * Y1^2 / dx^5;
            wi = 0;        
        else % (gx2 == 1 && gy2 == 1 && gz2 == 1 && gz1 == 1)
            % type = '4 solid faces at z1z2-x2y2';
            ui = 12*u1 * Z1*Z2 * X2^2 * Y2 / dx^5;
            vi = 12*v1 * Z1*Z2 * X2 * Y2^2 / dx^5;
            wi = 0;  
        end
        
    end
    
else
    %% type = '4, 5 or 6 solid faces';
    ui = 0;
    vi = 0;
    wi = 0;
end