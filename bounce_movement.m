function [x,y,z] = bounce_movement(i,j,k,...
    xBD,yBD,zBD,x,y,z,diffusionX,diffusionY,diffusionZ,dx)

    [iBD,jBD,kBD] = get_ijk (xBD,yBD,zBD,dx);
    
    if (i ~= iBD) %bounce X direction
        if (diffusionX >= 0)
            bounceX = - 2*abs(x - (dx*floor(x/dx)));
        else
            bounceX =   2*abs(x - (dx*ceil(x/dx)));
        end
        x = x + bounceX;
        
    elseif (j ~= jBD) %bounce Y direction
        if (diffusionY >= 0)
            bounceY = - 2*abs(y - (dx*floor(y/dx)));
        else
            bounceY =   2*abs(y - (dx*ceil(y/dx)));
        end
        y = y + bounceY;
        
    elseif (k ~= kBD) %bounce Z direction
        if (diffusionZ >= 0)
            bounceZ = - 2*abs(z - (dx*floor(z/dx)));
        else
            bounceZ =   2*abs(z - (dx*ceil(z/dx)));
        end
        z = z + bounceZ;
        
    end                
end