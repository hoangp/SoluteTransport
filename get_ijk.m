function [i,j,k] = get_ijk (x,y,z,dx)
    i = ceil (x / dx);
    j = ceil (y / dx);
    k = ceil (z / dx);
end