function out_location = out_system_location (x,y,z,left,right)
    if (x <= left)
        out_location = 1; % OUT_X1;
    elseif (x >= right)
        out_location = 2; % OUT_X2;
    elseif (y <= left)
        out_location = 3; % OUT_Y1;
    elseif (y >= right)
        out_location = 4; % OUT_Y2;
    elseif (z <= left)
        out_location = 5; % OUT_Z1;
    elseif (z >= right)
        out_location = 6; % OUT_Z2;
    else
        out_location = 0;
    end
end