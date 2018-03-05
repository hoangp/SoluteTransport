function [run_Npe,xpt,ypt,zpt,flow_displacement,flow_variance,inside_solid] = ...
    particletracking(g,u,v,w,dt,dx,L,gFlow,gLen,geo_zeta,...
        micro_zeta,chance_into_solid,inside_limit,inside_solid,...
        xpt,ypt,zpt,xpt0,ypt0,zpt0,Ntimestep)

% constants
DM_LIQUIDS = 10^3; % free molecular diffusion in liquids (microm2.s-1)

% results
Npe = [mean(u(u~=0)),mean(v(v~=0)),mean(w(w~=0))] * L / DM_LIQUIDS;
Nparticle = length (xpt); % number of particles
displacement = zeros(Nparticle,3); % particle displacement (in 3 direction)
variance = zeros(Ntimestep,3); % variance of displacement (in 3 direction)

% simulation starts
for t = 1:Ntimestep
    for p = 1:Nparticle
        x = xpt(p); y = ypt(p); z = zpt(p); % shorten variables name
        
        %% advection movement
        if inside_solid(p) == false
            
            % Particle location before advection (BA)
            xBA = x; yBA = y; zBA = z;
            
            %% streamline tracing algorithm
   
            % index location of the particle
            i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);
            
            % location of the edges
            x1 = dx*(i-1); x2 = dx*i;
            y1 = dx*(j-1); y2 = dx*j;
            z1 = dx*(k-1); z2 = dx*k;
            
            % velocities on the faces normal to each direction
            u1 = u(i-1,j,k); u2 = u(i,j,k);
            v1 = v(i,j-1,k); v2 = v(i,j,k);
            w1 = w(i,j,k-1); w2 = w(i,j,k);
            
            % short variables names
            X1 = x-x1; X2 = x2-x;
            Y1 = y-y1; Y2 = y2-y;
            Z1 = z-z1; Z2 = z2-z;
            U = u2-u1; V = v2-v1; W = w2-w1;
            
            % geometry of the edge
            gx1 = g(i-1,j,k); gx2 = g(i+1,j,k);
            gy1 = g(i,j-1,k); gy2 = g(i,j+1,k);
            gz1 = g(i,j,k-1); gz2 = g(i,j,k+1);
            
            % short variables names
            gX = gx1+gx2; gY = gy1+gy2; gZ = gz1+gz2;
            gXY = gX+gY; gXZ = gX+gZ; % gYZ = gY+gZ;
            
            % number of solid boundary arround the particle
            num_solid = gx1 + gx2 + gy1 + gy2 + gz1 + gz2;
            
            % Calculate interpolated velocities: ui,vi,wi
            
            if num_solid == 0
                %% No solid boundaries
                ui = U * X1 / dx + u1;
                vi = V * Y1 / dx + v1;
                wi = W * Z1 / dx + w1;
                
            elseif num_solid == 1
                %% One of the Neighbouring Voxels Is Solid (6 cases)
                if gx2 == 1 % type = '1 solid face at x2'; SAMPLE FORMULATION
                    ui = u1 * X2^2 / dx^2;
                    vi = 2*v1 * X2 / dx + 2*V * X2*Y1 / dx^2;
                    wi = 2*w1 * X2 / dx + 2*W * X2*Z1 / dx^2;
                elseif gx1 == 1 % type = '1 solid face at x1';
                    ui = u2 * X1^2 / dx^2;
                    vi = 2*v1 * X1 / dx + 2*V * X1*Y1 / dx^2;
                    wi = 2*w1 * X1 / dx + 2*W * X1*Z1 / dx^2;
                elseif gy2 == 1 % type = '1 solid face at y2';
                    ui = 2*u1 * Y2 / dx + 2*U * Y2*X1 / dx^2;
                    vi = v1 * Y2^2 / dx^2;
                    wi = 2*w1 * Y2 / dx + 2*W * Y2*Z1 / dx^2;
                elseif gy1 == 1  % type = '1 solid face at y1';
                    ui = 2*u1 * Y1 / dx + 2*U * Y1*X1 / dx^2;
                    vi = v2 * Y1^2 / dx^2;
                    wi = 2*w1 * Y1 / dx + 2*W * Y1*Z1 / dx^2;
                elseif gz2 == 1 % type = '1 solid face at z2';
                    ui = 2*u1 * Z2 / dx + 2*U * Z2*X1 / dx^2;
                    vi = 2*v1 * Z2 / dx + 2*V * Z2*Y1 / dx^2;
                    wi = w1 * Z2^2 / dx^2;
                elseif gz1 == 1 % type = '1 solid face at z1';
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
            
            %% advection displacement
            advection = [ui*dt, vi*dt, wi*dt];
            x = x + advection(1); y = y + advection(2); z = z + advection(3);
            
            if x<=dx || x>=gLen || y<=dx || y>=gLen || z<=dx || z>=gLen
                % if particle moves out system -> reinject
                %% ReinjectAdvection             
                random_index = ceil(Nparticle*rand());
                x = xpt0(random_index);
                y = ypt0(random_index);
                z = zpt0(random_index);
                
            else
                % advection movement within the system
                i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);
                if g(i,j,k) == 1
                    % if particle hits a solid surface -> NO MOVE
                    x = xBA; y = yBA; z = zBA; % restore location
                    advection = [0,0,0]; % update advection displacement
                    
                end
            end
        else 
            advection = [0,0,0]; % update advection displacement          
        end
        
        %% diffusion movement     
        
        % Particle location before diffusion (BD)
        xBD = x; yBD = y; zBD = z;
        
        % Randomize theta & phi
        theta = abs(2*rand()*pi);
        phi = abs(rand()*pi);
        
        % Diffustion displacement
        if inside_solid(p) == true
            zeta = micro_zeta;
        else
            zeta = geo_zeta;
        end
        
        diffusion = [zeta*sin(phi)*cos(theta), zeta*sin(phi)*sin(theta), zeta*cos(phi)];                    
        x = x + diffusion(1); y = y + diffusion(2); z = z + diffusion(3);
        
        if x<=dx || x>=gLen || y<=dx || y>=gLen || z<=dx || z>=gLen
            % if particle jumps out system -> reinject
            %% ReinjectDiffusion                    
            random_index = ceil(Nparticle*rand());
            x = xpt0(random_index);
            y = ypt0(random_index);
            z = zpt0(random_index);
            
        elseif inside_solid(p)
            % if inside_solid (and in system) -> update inside_solid status accordingly
            i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);
            if g(i,j,k)==0
                inside_solid(p) = false;               
            end
            
        else
            %% diffusion movement within the system (as normal)
            i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);
            if g(i,j,k) == 1
                % if hit a surface -> ROLL a chance to STAY INSIDE
                if (rand() <= chance_into_solid)
                    inside_solid(p) = true; % stay inside
                   
                    %% calculate inside location                    
                    i1=ceil(xBD/dx); j1=ceil(yBD/dx); k1=ceil(zBD/dx);
                    
                    if (i ~= i1) %bounce X direction
                        if (diffusion(1) >= 0)
                            insideX = abs(x - (dx*floor(x/dx)));
                            if insideX > inside_limit
                                x = x - insideX + inside_limit;
                            end
                        else
                            insideX = abs(x - (dx*ceil(x/dx)));
                            if insideX > inside_limit
                                x = x + insideX - inside_limit;
                            end
                        end
                        
                    elseif (j ~= j1) %bounce Y direction
                        if (diffusion(2) >= 0)
                            insideY = abs(y - (dx*floor(y/dx)));
                            if insideY > inside_limit
                                y = y - insideY + inside_limit;
                            end
                        else
                            insideY = abs(y - (dx*ceil(y/dx)));
                            if insideY > inside_limit
                                y = y + insideY - inside_limit;
                            end
                        end
                        
                    elseif (k ~= k1) %bounce Z direction
                        if (diffusion(3) >= 0)
                            insideZ = abs(z - (dx*floor(z/dx)));
                            if insideZ > inside_limit
                                z = z - insideZ + inside_limit;
                            end
                        else
                            insideZ = abs(z - (dx*ceil(z/dx)));
                            if insideZ > inside_limit
                                z = z + insideZ - inside_limit;
                            end
                        end
                    end
                    diffusion = [x - xBD, y - zBD, z - zBD]; % update diffusion              
                else
                    %% if hit a surface -> bounce movement (as normal)
                    
                    %% calculate bouncing location                
                    i1=ceil(xBD/dx); j1=ceil(yBD/dx); k1=ceil(zBD/dx);
                    
                    if (i ~= i1) %bounce X direction
                        if (diffusion(1) >= 0)
                            bounceX = - 2*abs(x - (dx*floor(x/dx)));
                        else
                            bounceX =   2*abs(x - (dx*ceil(x/dx)));
                        end
                        x = x + bounceX;
                        
                    elseif (j ~= j1) %bounce Y direction
                        if (diffusion(2) >= 0)
                            bounceY = - 2*abs(y - (dx*floor(y/dx)));
                        else
                            bounceY =   2*abs(y - (dx*ceil(y/dx)));
                        end
                        y = y + bounceY;
                        
                    elseif (k ~= k1) %bounce Z direction
                        if (diffusion(3) >= 0)
                            bounceZ = - 2*abs(z - (dx*floor(z/dx)));
                        else
                            bounceZ =   2*abs(z - (dx*ceil(z/dx)));
                        end
                        z = z + bounceZ;
                        
                    end
                    
                    diffusion = [x - xBD, y - zBD, z - zBD]; % update diffusion
                    
                    if x<=dx || x>=gLen || y<=dx || y>=gLen || z<=dx || z>=gLen
                        % after bouncing: if particle bounces out system -> reinject again
                        %% ReinjectDiffusion
                        random_index = ceil(Nparticle*rand());
                        x = xpt0(random_index);
                        y = ypt0(random_index);
                        z = zpt0(random_index);
                        
                    else
                        %% bounce movement within system
                        i=ceil(x/dx); j=ceil(y/dx); k=ceil(z/dx);
                        if g(i,j,k) == 1
                            % if bounce movement into a solid surface -> NO MOVE
                            x = xBD; y = yBD; z = zBD; % restore location
                            diffusion = [0,0,0]; % update diffusion
                        end
                    end
                end
            end
        end
        
        %% Update particle displacement
        displacement(p,:) = displacement(p,:) + advection + diffusion;
        
        xpt(p) = x; ypt(p) = y; zpt(p) = z; % copy short-name variables back to array(index)
    end % for Nparticle
    variance(t,:) = var(displacement);
end % for Ntimestep

% return value for main flow direction only
run_Npe = Npe(gFlow);
flow_displacement = displacement(:,gFlow);
flow_variance = variance(:,gFlow);

end





